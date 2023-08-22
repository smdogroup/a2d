#ifndef TOPOLOGY_2D_H
#define TOPOLOGY_2D_H

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

template <typename T, int Degree = 1>
class TopoElasticityAnalysis2D {
 public:
  static constexpr int spatial_dim = 2;
  static constexpr int degree = Degree;     // Polynomial degree
  static constexpr int order = degree + 1;  // Spline order, = degree + 1
  static constexpr int data_degree = degree - 1;
  // static constexpr int filter_degree = degree - 1;
  static constexpr int filter_degree = degree;
  static constexpr int data_order = data_degree + 1;
  static constexpr int data_dim = 1;
  static constexpr int var_dim = spatial_dim;
  static constexpr int block_size = var_dim;

  // Data types
  using Vec_t = SolutionVector<T>;
  using BSRMat_t = BSRMat<T, block_size, block_size>;
  using CSCMat_t = CSCMat<T>;

  // Elasticity component
  using Quadrature = QuadGaussQuadrature<order>;
  using DataBasis = FEBasis<T, LagrangeL2QuadBasis<T, data_dim, data_degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, var_dim, degree>>;
  using DataElemVec = ElementVector_Serial<T, DataBasis, Vec_t>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis, Vec_t>;
  using ElemVec = ElementVector_Serial<T, Basis, Vec_t>;

  // Traction component
  using TQuadrature = LineGaussQuadrature<order>;
  using TDataBasis = FEBasis<T>;  // No data related to the traction
  using TGeoBasis = FEBasis<T, LagrangeH1LineBasis<T, spatial_dim, degree>>;
  using TBasis = FEBasis<T, LagrangeH1LineBasis<T, spatial_dim, degree>>;
  using TDataElemVec = ElementVector_Serial<T, TDataBasis, Vec_t>;
  using TGeoElemVec = ElementVector_Serial<T, TGeoBasis, Vec_t>;
  using TElemVec = ElementVector_Serial<T, TBasis, Vec_t>;

  // Filter
  using FSpace = FESpace<T, spatial_dim, H1Space<T, data_dim, spatial_dim>>;
  using FBasis = FEBasis<T, LagrangeH1QuadBasis<T, data_dim, filter_degree>>;
  using FElemVec = ElementVector_Serial<T, FBasis, Vec_t>;

  using Integrand = IntegrandTopoLinearElasticity<T, spatial_dim>;
  using Traction = IntegrandTopoSurfaceTraction<T, spatial_dim>;

  using FE_PDE =
      FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_Traction =
      FiniteElement<T, Traction, TQuadrature, TDataBasis, TGeoBasis, TBasis>;

  TopoElasticityAnalysis2D(MeshConnectivityBase &conn, DirichletBCInfo &bcinfo,
                           index_t nquad, index_t *quad, double *Xloc,
                           const index_t traction_label, const T tx_traction[],
                           T E = 70e3, T nu = 0.3, T q = 5.0)
      : E(E),
        nu(nu),
        q(q),

        mesh(conn),
        geomesh(conn),
        datamesh(conn),
        filter_mesh(conn),
        traction_mesh(traction_label, conn, mesh),
        traction_geomesh(traction_label, conn, geomesh),

        bcs(conn, mesh, bcinfo),

        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        data(datamesh.get_num_dof()),
        filter_data(filter_mesh.get_num_dof()),

        elem_sol(mesh, sol),
        elem_geo(geomesh, geo),
        elem_data(datamesh, data),
        elem_filter_data(filter_mesh, filter_data),

        elem_traction_sol(traction_mesh, sol),
        elem_traction_geo(traction_geomesh, geo),

        integrand(E, nu, q),
        traction_integrand(tx_traction) {
    // Set geometry
    set_geo_from_quad_nodes<GeoBasis>(nquad, quad, Xloc, elem_geo);

    // Symbolically create block CSR matrix
    mesh.template create_block_csr<block_size>(nrows, rowp, cols);
    bsr_mat = BSRMat_t(nrows, nrows, cols.size(), rowp, cols);
    csc_mat = bsr_to_csc(bsr_mat);

    // Set up Cholesky solver but don't set up values yet
    bool set_values = false;
    chol =
        new SparseCholesky(csc_mat, CholOrderingType::ND, nullptr, set_values);
  }

  void set_design_var() {}

  void factor() {
    // Create a new element view of the system matrix
    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);

    // Populate the system matrix
    fe.add_jacobian_new(integrand, elem_data, elem_geo, elem_sol, elem_mat);

    // Apply boundary conditions to each row
    const index_t *bc_dofs;
    index_t nbcs = bcs.get_bcs(&bc_dofs);
    bsr_mat.zero_rows(nbcs, bc_dofs);

    // convert bsr to csc because we want to use Cholesky factorization
    csc_mat = bsr_to_csc(bsr_mat);

    // Apply boundary conditions to each column
    csc_mat.zero_columns(nbcs, bc_dofs);

    // TODO: delete this
    csc_mat.write_mtx("csc_matrix.mtx");

    // Set values to Cholesky solver and factorize
    chol->setValues(csc_mat);
    chol->factor();
  }

  void solve() {
    // Factorize system matrix
    this->factor();

    // Set up right-hand-side
    ElementVector_Empty<ElemVecType::Serial> elem_traction_data;
    A2D::SolutionVector<T> traction_res(mesh.get_num_dof());
    TElemVec elem_traction_res(traction_mesh, traction_res);

    traction.add_residual(traction_integrand, elem_traction_data,
                          elem_traction_geo, elem_traction_sol,
                          elem_traction_res);

    std::vector<T> rhs(sol.get_num_dof());
    for (index_t i = 0; i < sol.get_num_dof(); i++) {
      rhs[i] = -traction_res[i];
      // std::printf("rhs[%4d]: %10.5e\n", i, rhs[i]);  // TODO: delete this
    }

    // Solve
    chol->solve(rhs.data());

    // Record solution
    for (index_t i = 0; i < sol.get_num_dof(); i++) {
      sol[i] = rhs[i];
    }
  }

  void tovtk(const std::string filename) {
    A2D::write_quad_to_vtk<3, degree, T, DataBasis, GeoBasis, Basis>(
        integrand, elem_data, elem_geo, elem_sol, filename,
        [](index_t k, typename Integrand::DataSpace &d,
           typename Integrand::FiniteElementGeometry &g,
           typename Integrand::FiniteElementSpace &s) {
          if (k == 2) {  // write data
            return (d.template get<0>()).get_value();
          } else {  // write solution components
            auto u = (s.template get<0>()).get_value();
            return u(k);
          }
        });
  }

  ElementMesh<Basis> &get_mesh() { return mesh; }
  CSCMat<T> &get_csc_matrix() { return csc_mat; }
  SparseCholesky<T> *get_chol() { return chol; }

 private:
  T E, nu, q;

  ElementMesh<Basis> mesh;
  ElementMesh<GeoBasis> geomesh;
  ElementMesh<DataBasis> datamesh;
  ElementMesh<FBasis> filter_mesh;
  ElementMesh<TBasis> traction_mesh;
  ElementMesh<TGeoBasis> traction_geomesh;

  DirichletBCs<Basis> bcs;

  Vec_t sol;
  Vec_t geo;
  Vec_t data;
  Vec_t filter_data;

  ElemVec elem_sol;
  GeoElemVec elem_geo;
  DataElemVec elem_data;
  FElemVec elem_filter_data;
  TElemVec elem_traction_sol;
  TGeoElemVec elem_traction_geo;

  Integrand integrand;
  Traction traction_integrand;

  FE_PDE fe;
  FE_Traction traction;

  // System matrices
  index_t nrows;
  std::vector<index_t> rowp, cols;
  BSRMat_t bsr_mat;
  CSCMat_t csc_mat;

  // Cholesky solver
  SparseCholesky<T> *chol;
};

#endif  // TOPOLOGY_2D_H