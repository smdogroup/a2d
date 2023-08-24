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
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

template <class... Args>
using ElementVector = ElementVector_Parallel<Args...>;

/**
 * @tparam T type
 * @tparam degree polynomial degree
 */
template <typename T, index_t degree>
class TopoElasticityAnalysis {
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
  using DataElemVec = ElementVector<T, DataBasis, Vec_t>;
  using GeoElemVec = ElementVector<T, GeoBasis, Vec_t>;
  using ElemVec = ElementVector<T, Basis, Vec_t>;

  using TQuadrature = QuadGaussQuadrature<order>;
  using TDataBasis = FEBasis<T>;
  using TGeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using TBasis = FEBasis<T, LagrangeH1QuadBasis<T, var_dim, degree>>;
  using TDataElemVec = ElementVector<T, TDataBasis, Vec_t>;
  using TGeoElemVec = ElementVector<T, TGeoBasis, Vec_t>;
  using TElemVec = ElementVector<T, TBasis, Vec_t>;

  using Integrand_elas = IntegrandTopoLinearElasticity<T, spatial_dim>;
  using Integrand_body = IntegrandTopoBodyForce<T, spatial_dim>;
  using Integrand_traction = IntegrandTopoSurfaceTraction<T, spatial_dim>;

  using FE_elas =
      FiniteElement<T, Integrand_elas, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_body =
      FiniteElement<T, Integrand_body, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_Traction = FiniteElement<T, Integrand_traction, TQuadrature,
                                    TDataBasis, TGeoBasis, TBasis>;

  TopoElasticityAnalysis(MeshConnectivityBase &conn, index_t nquad,
                         const index_t quad[], const double Xloc[],
                         DirichletBCInfo &bcinfo, const T *gravity, T E = 70e3,
                         T nu = 0.3, T q = 5.0)
      : E(E),
        nu(nu),
        q(q),
        integrand_elas(E, nu, q),
        integrand_body(q, gravity),
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

    // Set density
    index_t ndvs = datamesh.get_num_dof();
    for (index_t i = 0; i < ndvs; i++) {
      data[i] = 1.0;
    }
  }

  void add_jacobian() {
    // Symbolically create block CSR matrix
    index_t nrows;
    std::vector<index_t> rowp, cols;
    mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    BSRMat_t bsr_mat(nrows, nrows, cols.size(), rowp, cols);

    // Populate the CSR matrix
    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);
    fe_elas.add_jacobian(integrand_elas, elem_data, elem_geo, elem_sol,
                         elem_mat);
  }

  void add_residual() {
    index_t N = mesh.get_num_dof();
    printf("number of dof: %d\n", N);
    if (N > 20) {
      N = 20;
    }

    Vec_t res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    StopWatch watch;

    double t1 = watch.lap();
    fe_body.add_residual(integrand_body, elem_data, elem_geo, elem_sol,
                         elem_res);
    double t2 = watch.lap();
    printf("add_residual time: %.5f ms\n", (t2 - t1) * 1e3);

    for (index_t i = 0; i < N; i++) {
      printf("res[%2d]: %20.10e\n", i, res[i]);
    }
  }

 private:
  T E, nu, q;
  Integrand_elas integrand_elas;
  Integrand_body integrand_body;
  Integrand_traction integrand_traction;
  ElementMesh<Basis> mesh;
  ElementMesh<GeoBasis> geomesh;
  ElementMesh<DataBasis> datamesh;
  Vec_t sol;
  Vec_t geo;
  Vec_t data;
  ElemVec elem_sol;
  GeoElemVec elem_geo;
  DataElemVec elem_data;
  FE_elas fe_elas;
  FE_body fe_body;
  FE_Traction fe_trac;
};

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    using T = double;

    constexpr bool smoke_test = false;

    // Number of elements in each dimension
    const index_t degree = smoke_test ? 1 : 4;
    const index_t nx = smoke_test ? 2 : 512, ny = smoke_test ? 2 : 512;
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

    index_t bc_label =
        conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

    DirichletBCInfo bcinfo;
    bcinfo.add_boundary_condition(bc_label);

    T gravity[2] = {2.5, -9.8};
    TopoElasticityAnalysis<T, degree> prob(conn, nquad, quad.data(),
                                           Xloc.data(), bcinfo, gravity);

    prob.add_residual();
  }

  Kokkos::finalize();

  return 0;
}