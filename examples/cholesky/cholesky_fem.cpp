#include <vector>

#include "a2ddefs.h"
#include "multiphysics/feanalysis.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/integrand_helmholtz.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    using T = double;

    // Number of elements in each dimension
    const int nx = 16, ny = 8, nz = 8;
    const double lx = 1.5, ly = 2.0, lz = 1.0;

    // Set up mesh
    const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
    int ntets = 0, nwedge = 0, npyrmd = 0;
    const int nhex = nx * ny * nz;
    int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
    std::vector<int> hex(8 * nhex);
    std::vector<double> Xloc(3 * nverts);
    MesherBrick3D mesher(nx, ny, nz, lx, ly, lz);
    mesher.set_X_conn<int, double>(Xloc.data(), hex.data());
    MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex.data(), nwedge,
                            wedge, npyrmd, pyrmd);

    // Set up bcs
    auto node_num = [](int i, int j, int k) {
      return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
    };

    const int num_boundary_verts = (ny + 1) * (nz + 1);
    int boundary_verts[num_boundary_verts];

    for (int k = 0, index = 0; k < nz + 1; k++) {
      for (int j = 0; j < ny + 1; j++, index++) {
        boundary_verts[index] = node_num(0, j, k);
      }
    }

    A2D::index_t bc_label =
        conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

    A2D::DirichletBCInfo bcinfo;
    bcinfo.add_boundary_condition(bc_label);

    // Set up the type of element we're using
    const index_t degree = 1;
    const index_t dim = 3;
    const index_t data_size = 1;
    const index_t block_size = dim;
    using Impl_t = typename DirectCholeskyAnalysis<T, block_size>::Impl_t;
    using Vec_t = typename Impl_t::Vec_t;
    constexpr GreenStrainType etype = GreenStrainType::LINEAR;
    using Elem_t = HexTopoElement<Impl_t, etype, degree>;
    using Mat_t = typename DirectCholeskyAnalysis<T, dim>::Mat_t;

    // Create the meshes for the elements
    auto data_mesh = std::make_shared<ElementMesh<Elem_t::DataBasis>>(conn);
    auto geo_mesh = std::make_shared<ElementMesh<Elem_t::GeoBasis>>(conn);
    auto sol_mesh = std::make_shared<ElementMesh<Elem_t::Basis>>(conn);

    // Get the number of different dof for the analysis
    index_t ndata = data_mesh->get_num_dof();
    index_t ngeo = geo_mesh->get_num_dof();
    index_t ndof = sol_mesh->get_num_dof();

    // Create the solution/residual vector for the analysis
    auto data = std::make_shared<Vec_t>(ndata);
    auto geo = std::make_shared<Vec_t>(ngeo);
    auto sol = std::make_shared<Vec_t>(ndof);
    auto res = std::make_shared<Vec_t>(ndof);

    // Create the element assembler
    auto assembler = std::make_shared<ElementAssembler<Impl_t>>();

    // Create the element integrand
    T E = 70.0, nu = 0.3, q = 5.0;
    TopoElasticityIntegrand<T, dim, etype> elem_integrand(E, nu, q);
    assembler->add_element(std::make_shared<Elem_t>(elem_integrand, data_mesh,
                                                    geo_mesh, sol_mesh));

    // Set the geometry
    typename Impl_t::ElementVector<Elem_t::GeoBasis> elem_geo(*geo_mesh, *geo);
    set_geo_from_hex_nodes<Elem_t::GeoBasis>(nhex, hex.data(), Xloc.data(),
                                             elem_geo);

    for (int i = 0; i < ndata; i++) {
      (*data)[i] = T(1.0);
    }

    // Set up the boundary conditions
    auto bcs = std::make_shared<DirichletBCs<T>>();
    bcs->add_bcs(std::make_shared<DirichletBasis<T, Elem_t::Basis>>(
        conn, *sol_mesh, bcinfo, 0.0));

    // Create the assembler object
    auto analysis = std::make_shared<DirectCholeskyAnalysis<T, block_size>>(
        data, geo, sol, res, assembler, bcs);

    analysis->linear_solve();

    // // Create bsr mat and zero bcs rows
    std::shared_ptr<Mat_t> bsr_mat = analysis->get_mat();

    // Convert to csr mat and zero bcs columns
    StopWatch watch;
    CSRMat<T> csr_mat = bsr_to_csr(*bsr_mat);
    double t1 = watch.lap();
    std::printf("bsr->csr time: %12.5e s\n", t1);

    // Convert to csc mat and apply bcs
    CSCMat<T> csc_mat = bsr_to_csc(*bsr_mat);
    const index_t *bc_dofs;
    index_t nbcs = bcs->get_bcs(&bc_dofs);
    csc_mat.zero_columns(nbcs, bc_dofs);
    double t2 = watch.lap();
    std::printf("bsr->csc time: %12.5e s\n", t2 - t1);

    // Create rhs
    std::vector<T> b(csc_mat.nrows);
    for (int i = 0; i < csc_mat.nrows; i++) {
      b[i] = 0.0;
    }
    for (int i = 0; i < csc_mat.nrows; i++) {
      for (int jp = csc_mat.colp[i]; jp < csc_mat.colp[i + 1]; jp++) {
        b[csc_mat.rows[jp]] += csc_mat.vals[jp];
      }
    }

    // Write to mtx
    bsr_mat->write_mtx("bsr_mat.mtx", 1e-12);
    csr_mat.write_mtx("csr_mat.mtx", 1e-12);
    csc_mat.write_mtx("csc_mat.mtx", 1e-12);
    double t3 = watch.lap();
    std::printf("write mtx time: %12.5e s\n", t3 - t2);

    // Perform cholesky factorization
    printf("number of vertices: %d\n", conn.get_num_verts());
    printf("number of edges:    %d\n", conn.get_num_edges());
    printf("number of faces:    %d\n", conn.get_num_bounds());
    printf("number of elements: %d\n", conn.get_num_elements());
    printf("number of dofs:     %d\n", ndof);
    printf("matrix dimension:  (%d, %d)\n", csr_mat.nrows, csr_mat.ncols);
    double t4 = watch.lap();
    SparseCholesky<T> *chol = new SparseCholesky<T>(csc_mat);
    double t5 = watch.lap();
    std::printf("Setup/order/setvalue time: %12.5e s\n", t5 - t4);
    chol->factor();
    double t6 = watch.lap();
    std::printf("Factor time:               %12.5e s\n", t6 - t5);
    chol->solve(b.data());
    double t7 = watch.lap();
    std::printf("Solve time:                %12.5e s\n", t7 - t6);
    T err = 0.0;
    for (int i = 0; i < csc_mat.nrows; i++) {
      err += (1.0 - b[i]) * (1.0 - b[i]);
    }
    std::printf("||x - e||: %25.15e\n", sqrt(err));
  }

  Kokkos::finalize();

  return 0;
}