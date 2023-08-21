#include <vector>

#include "a2dobjs.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

/**
 * @tparam T type
 * @tparam degree polynomial degree
 */
template <typename T, index_t degree>
class TopoElasticityAnalysis2D {
 public:
  static constexpr int spatial_dim = 2;
  static constexpr int order = degree + 1;  // Spline order, = degree + 1
  static constexpr int data_degree = degree - 1;
  static constexpr int data_order = data_degree + 1;
  static constexpr int data_dim = 1;
  static constexpr int var_dim = spatial_dim;
  static constexpr int block_size = var_dim;

  using Vec_t = SolutionVector<T>;
  using BSRMat_t = BSRMat<T, block_size, block_size>;

  using Quadrature = QuadGaussQuadrature<order>;
  using DataBasis = FEBasis<T, LagrangeL2QuadBasis<T, data_dim, data_degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, var_dim, degree>>;
  using DataElemVec = ElementVector_Serial<T, DataBasis, Vec_t>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis, Vec_t>;
  using ElemVec = ElementVector_Serial<T, Basis, Vec_t>;

  using TQuadrature = QuadGaussQuadrature<order>;
  using TDataBasis = FEBasis<T>;
  using TGeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using TBasis = FEBasis<T, LagrangeH1QuadBasis<T, var_dim, degree>>;
  using TDataElemVec = ElementVector_Serial<T, TDataBasis, Vec_t>;
  using TGeoElemVec = ElementVector_Serial<T, TGeoBasis, Vec_t>;
  using TElemVec = ElementVector_Serial<T, TBasis, Vec_t>;

  using PDEIntegrand = IntegrandTopoLinearElasticity<T, spatial_dim>;
  using Traction = IntegrandTopoSurfaceTraction<T, spatial_dim>;

  using FE_PDE = FiniteElement<T, PDEIntegrand, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_Traction =
      FiniteElement<T, Traction, TQuadrature, TDataBasis, TGeoBasis, TBasis>;

  TopoElasticityAnalysis2D(MeshConnectivityBase &conn, index_t nquad,
                           const index_t quad[], const double Xloc[],
                           DirichletBCInfo &bcinfo, T E = 70e3, T nu = 0.3,
                           T q = 5.0)
      : E(E),
        nu(nu),
        q(q),
        pde(E, nu, q),
        mesh(conn),
        geomesh(conn),
        datamesh(conn),
        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        data(datamesh.get_num_dof()),
        elem_sol(mesh, sol),
        elem_geo(geomesh, geo),
        elem_data(datamesh, data) {
    // Set geometry
    set_geo_from_quad_nodes<GeoBasis>(nquad, quad, Xloc, elem_geo);
  }

  BSRMat_t create_fea_bsr_matrix() {
    // Symbolically create block CSR matrix
    index_t nrows;
    std::vector<index_t> rowp, cols;
    mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    BSRMat_t bsr_mat(nrows, nrows, cols.size(), rowp, cols);

    // Populate the CSR matrix
    FE_PDE fe;

    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);
    fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);
    bsr_mat.write_mtx("Jacobian.mtx");

    return bsr_mat;
  }

  ElementMesh<Basis> &get_mesh() { return mesh; }

  void tovtk(const std::string filename) {
    A2D::write_quad_to_vtk<3, degree, T, DataBasis, GeoBasis, Basis>(
        pde, elem_data, elem_geo, elem_sol, filename,
        [](index_t k, typename PDEIntegrand::DataSpace &d,
           typename PDEIntegrand::FiniteElementGeometry &g,
           typename PDEIntegrand::FiniteElementSpace &s) {
          if (k == 2) {  // write data
            return (d.template get<0>()).get_value();
          } else {  // write solution components
            auto u = (s.template get<0>()).get_value();
            return u(k);
          }
        });
  }

 private:
  T E, nu, q;
  PDEIntegrand pde;
  ElementMesh<Basis> mesh;
  ElementMesh<GeoBasis> geomesh;
  ElementMesh<DataBasis> datamesh;
  Vec_t sol;
  Vec_t geo;
  Vec_t data;
  ElemVec elem_sol;
  GeoElemVec elem_geo;
  DataElemVec elem_data;
};

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    using T = double;

    // Number of elements in each dimension
    const index_t degree = 2;
    const index_t nx = 32, ny = 32;
    const double lx = 1.5, ly = 2.0;

    // Set up mesh
    const index_t nverts = (nx + 1) * (ny + 1);
    const index_t ntri = 0, nquad = nx * ny;
    index_t *tri = nullptr;
    std::vector<index_t> quad(4 * nquad);
    std::vector<double> Xloc(2 * nverts);
    MesherRect2D mesher(nx, ny, lx, ly);
    bool randomize = true;
    unsigned int seed = 0;
    double fraction = 0.2;
    mesher.set_X_conn<index_t, double>(Xloc.data(), quad.data(), randomize,
                                       seed, fraction);

    MeshConnectivity2D conn(nverts, ntri, tri, nquad, quad.data());
    // Set up bcs
    auto node_num = [](index_t i, index_t j) { return i + j * (nx + 1); };

    const index_t num_boundary_verts = (ny + 1);
    index_t boundary_verts[num_boundary_verts];

    for (index_t j = 0; j < ny + 1; j++) {
      boundary_verts[j] = node_num(0, j);
    }

    A2D::index_t bc_label =
        conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

    A2D::DirichletBCInfo bcinfo;
    bcinfo.add_boundary_condition(bc_label);

    TopoElasticityAnalysis2D<T, degree> prob(conn, nquad, quad.data(),
                                             Xloc.data(), bcinfo);

    // Save mesh
    prob.tovtk("mesh.vtk");

    // Create bsr mat and zero bcs rows
    TopoElasticityAnalysis2D<T, degree>::BSRMat_t bsr_mat =
        prob.create_fea_bsr_matrix();
    const index_t *bc_dofs;
    DirichletBCs<TopoElasticityAnalysis2D<T, degree>::Basis> bcs(
        conn, prob.get_mesh(), bcinfo);
    index_t nbcs = bcs.get_bcs(&bc_dofs);
    bsr_mat.zero_rows(nbcs, bc_dofs);

    // Convert to csr mat and zero bcs columns
    StopWatch watch;
    CSRMat<T> csr_mat = bsr_to_csr(bsr_mat);

    // Convert to csc mat and apply bcs
    CSCMat<T> csc_mat = bsr_to_csc(bsr_mat);
    csc_mat.zero_columns(nbcs, bc_dofs);

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
    // double t2 = watch.lap();
    // bsr_mat.write_mtx("bsr_mat.mtx", 1e-12);
    // csc_mat.write_mtx("csc_mat.mtx", 1e-12);
    // double t3 = watch.lap();
    // printf("write mtx time: %12.5e s\n", t3 - t2);

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
    printf("Setup/order/setvalue time: %12.5e s\n", t5 - t4);
    chol->factor();
    double t6 = watch.lap();
    printf("Factor time:               %12.5e s\n", t6 - t5);
    chol->solve(b.data());
    double t7 = watch.lap();
    printf("Solve time:                %12.5e s\n", t7 - t6);
    T err = 0.0;
    for (int i = 0; i < csc_mat.nrows; i++) {
      err += (1.0 - b[i]) * (1.0 - b[i]);
    }
    printf("||x - e||: %25.15e\n", sqrt(err));
  }

  Kokkos::finalize();

  return 0;
}