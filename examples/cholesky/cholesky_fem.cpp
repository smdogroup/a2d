#include <vector>

#include "a2dobjs.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

template <typename T>
class FEProb {
 public:
  static constexpr int spatial_dim = 3;
  static constexpr int data_dim = 1;
  static constexpr int var_dim = 3;
  static constexpr int block_size = var_dim;
  static constexpr int degree = 1;

  using Vec_t = SolutionVector<T>;
  using BSRMat_t = BSRMat<T, block_size, block_size>;

  using Integrand = IntegrandTopoLinearElasticity<T, spatial_dim>;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T, LagrangeL2HexBasis<T, data_dim, degree - 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, var_dim, degree>>;

  using FE =
      FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis>;
  using DataElemVec = ElementVector_Serial<T, DataBasis, Vec_t>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis, Vec_t>;
  using ElemVec = ElementVector_Serial<T, Basis, Vec_t>;

  FEProb(MeshConnectivity3D &conn, int nhex, const int hex[],
         const double Xloc[], DirichletBCInfo &bcinfo, T E = 70e3, T nu = 0.3,
         T q = 5.0)
      : E(E),
        nu(nu),
        q(q),
        datamesh(conn),
        geomesh(conn),
        mesh(conn),
        data(datamesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        sol(mesh.get_num_dof()),
        elem_data(datamesh, data),
        elem_geo(geomesh, geo),
        elem_sol(mesh, sol) {
    // Set geometry
    set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);
  }

  BSRMat_t create_fea_bsr_matrix() {
    // Symbolically create block CSR matrix

    index_t nrows;
    std::vector<index_t> rowp, cols;
    mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    BSRMat_t bsr_mat(nrows, nrows, cols.size(), rowp, cols);

    // Populate the CSR matrix
    FE fe;
    Integrand integrand(E, nu, q);

    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);
    fe.add_jacobian(integrand, elem_data, elem_geo, elem_sol, elem_mat);

    return bsr_mat;
  }

  ElementMesh<Basis> &get_mesh() { return mesh; }

 private:
  T E, nu, q;
  ElementMesh<DataBasis> datamesh;
  ElementMesh<GeoBasis> geomesh;
  ElementMesh<Basis> mesh;
  Vec_t data;
  Vec_t geo;
  Vec_t sol;
  DataElemVec elem_data;
  GeoElemVec elem_geo;
  ElemVec elem_sol;
};

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

    FEProb<T> prob(conn, nhex, hex.data(), Xloc.data(), bcinfo);

    // Create bsr mat and zero bcs rows
    FEProb<T>::BSRMat_t bsr_mat = prob.create_fea_bsr_matrix();
    const index_t *bc_dofs;
    DirichletBCs<FEProb<T>::Basis> bcs(conn, prob.get_mesh(), bcinfo);
    index_t nbcs = bcs.get_bcs(&bc_dofs);
    bsr_mat.zero_rows(nbcs, bc_dofs);

    // Convert to csr mat and zero bcs columns
    StopWatch watch;
    CSRMat<T> csr_mat = bsr_to_csr(bsr_mat);
    double t1 = watch.lap();
    std::printf("bsr->csr time: %12.5e s\n", t1);

    // Convert to csc mat and apply bcs
    CSCMat<T> csc_mat = bsr_to_csc(bsr_mat);
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

    // // Write to mtx
    // bsr_mat.write_mtx("bsr_mat.mtx", 1e-12);
    // csr_mat.write_mtx("csr_mat.mtx", 1e-12);
    // csc_mat.write_mtx("csc_mat.mtx", 1e-12);
    // double t3 = watch.lap();
    // std::printf("write mtx time: %12.5e s\n", t3 - t2);

    // Perform cholesky factorization
    printf("number of vertices: %d\n", conn.get_num_verts());
    printf("number of edges:    %d\n", conn.get_num_edges());
    printf("number of faces:    %d\n", conn.get_num_bounds());
    printf("number of elements: %d\n", conn.get_num_elements());
    printf("number of dofs:     %d\n", prob.get_mesh().get_num_dof());
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