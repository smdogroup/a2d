#include <vector>

#include "a2dobjs.h"
#include "multiphysics/elasticity.h"
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

template <typename T>
class TopoElasticityAnalysis2D {
 public:
  static constexpr int spatial_dim = 2;
  static constexpr int degree = 1;          // Polynomial degree
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

  using PDE = TopoLinearElasticity<T, spatial_dim>;
  using Traction = TopoSurfaceTraction<T, spatial_dim>;

  using FE_PDE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_Traction =
      FiniteElement<T, Traction, TQuadrature, TDataBasis, TGeoBasis, TBasis>;

  TopoElasticityAnalysis2D(MeshConnectivityBase &conn, index_t nquad,
                           const index_t quad[], const double Xloc[],
                           DirichletBCInfo &bcinfo, T E = 70e3, T nu = 0.3,
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
    PDE pde(E, nu, q);

    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);
    fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

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
    const index_t nx = 512, ny = 512;
    const double lx = 1.5, ly = 2.0;

    // Set up mesh
    const index_t nverts = (nx + 1) * (ny + 1);
    const index_t ntri = 0, nquad = nx * ny;
    index_t *tri = nullptr;
    std::vector<index_t> quad(4 * nquad);
    std::vector<double> Xloc(2 * nverts);
    MesherRect2D mesher(nx, ny, lx, ly);
    mesher.set_X_conn<index_t, double>(Xloc.data(), quad.data());

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

    TopoElasticityAnalysis2D<T> prob(conn, nquad, quad.data(), Xloc.data(),
                                     bcinfo);

    // Create bsr mat and zero bcs rows
    TopoElasticityAnalysis2D<T>::BSRMat_t bsr_mat =
        prob.create_fea_bsr_matrix();
    const index_t *bc_dofs;
    DirichletBCs<TopoElasticityAnalysis2D<T>::Basis> bcs(conn, prob.get_mesh(),
                                                         bcinfo);
    index_t nbcs = bcs.get_bcs(&bc_dofs);
    bsr_mat.zero_rows(nbcs, bc_dofs);

    // Convert to csr mat and zero bcs columns
    StopWatch watch;
    CSRMat<T> csr_mat = bsr_to_csr(bsr_mat);
    double t1 = watch.lap();
    printf("bsr->csr time: %12.5e s\n", t1);

    // Convert to csc mat and apply bcs
    CSCMat<T> csc_mat = bsr_to_csc(bsr_mat);
    csc_mat.zero_columns(nbcs, bc_dofs);
    double t2 = watch.lap();
    printf("bsr->csc time: %12.5e s\n", t2 - t1);

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
    bsr_mat.write_mtx("bsr_mat.mtx", 1e-12);
    csr_mat.write_mtx("csr_mat.mtx", 1e-12);
    csc_mat.write_mtx("csc_mat.mtx", 1e-12);
    double t3 = watch.lap();
    printf("write mtx time: %12.5e s\n", t3 - t2);

    // Perform cholesky factorization
    printf("number of vertices: %d\n", conn.get_num_verts());
    printf("number of edges:    %d\n", conn.get_num_edges());
    printf("number of faces:    %d\n", conn.get_num_faces());
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