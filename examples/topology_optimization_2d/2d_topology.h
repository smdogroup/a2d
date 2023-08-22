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
  static constexpr int data_order = data_degree + 1;
  static constexpr int data_dim = 1;
  static constexpr int var_dim = spatial_dim;
  static constexpr int block_size = var_dim;

  // Data types
  using Vec_t = SolutionVector<T>;
  using BSRMat_t = BSRMat<T, block_size, block_size>;
  using CSCMat_t = CSCMat<T>;

  // Elasticity component
  using QuadratureElas = QuadGaussQuadrature<order>;
  using DataBasisElas =
      FEBasis<T, LagrangeL2QuadBasis<T, data_dim, data_degree>>;
  using GeoBasisElas = FEBasis<T, LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using BasisElas = FEBasis<T, LagrangeH1QuadBasis<T, var_dim, degree>>;
  using DataElemVecElas = ElementVector_Serial<T, DataBasisElas, Vec_t>;
  using GeoElemVecElas = ElementVector_Serial<T, GeoBasisElas, Vec_t>;
  using ElemVecElas = ElementVector_Serial<T, BasisElas, Vec_t>;

  // Traction component
  using QuadratureTrac = LineGaussQuadrature<order>;
  using DataBasisTrac = FEBasis<T>;  // No data related to the traction
  using GeoBasisTrac = FEBasis<T, LagrangeH1LineBasis<T, spatial_dim, degree>>;
  using BasisTrac = FEBasis<T, LagrangeH1LineBasis<T, spatial_dim, degree>>;
  using DataElemVecTrac = ElementVector_Serial<T, DataBasisTrac, Vec_t>;
  using GeoElemVecTrac = ElementVector_Serial<T, GeoBasisTrac, Vec_t>;
  using ElemVecTrac = ElementVector_Serial<T, BasisTrac, Vec_t>;

  // Integrands
  using IntegrandElas = IntegrandTopoLinearElasticity<T, spatial_dim>;
  using IntegrandTrac = IntegrandTopoSurfaceTraction<T, spatial_dim>;

  // Element operations
  using ElemOpsElas = FiniteElement<T, IntegrandElas, QuadratureElas,
                                    DataBasisElas, GeoBasisElas, BasisElas>;
  using ElemOpsTrac = FiniteElement<T, IntegrandTrac, QuadratureTrac,
                                    DataBasisTrac, GeoBasisTrac, BasisTrac>;

  TopoElasticityAnalysis2D(MeshConnectivityBase &conn, DirichletBCInfo &bcinfo,
                           index_t nquad, index_t *quad, double *Xloc,
                           const index_t traction_label, const T tx_traction[],
                           T E = 70e3, T nu = 0.3, T q = 5.0)
      : E(E),
        nu(nu),
        q(q),

        elem_mesh_elas(conn),
        elem_geomesh_elas(conn),
        elem_datamesh_elas(conn),
        elem_mesh_trac(traction_label, conn, elem_mesh_elas),
        elem_geomesh_trac(traction_label, conn, elem_geomesh_elas),

        bcs(conn, elem_mesh_elas, bcinfo),

        sol(elem_mesh_elas.get_num_dof()),
        geo(elem_geomesh_elas.get_num_dof()),
        data(elem_datamesh_elas.get_num_dof()),

        elem_sol_elas(elem_mesh_elas, sol),
        elem_geo_elas(elem_geomesh_elas, geo),
        elem_data_elas(elem_datamesh_elas, data),

        elem_sol_trac(elem_mesh_trac, sol),
        elem_geo_trac(elem_geomesh_trac, geo),

        integrand_elas(E, nu, q),
        integrand_trac(tx_traction) {
    // Set geometry
    set_geo_from_quad_nodes<GeoBasisElas>(nquad, quad, Xloc, elem_geo_elas);

    // Symbolically create block CSR matrix
    elem_mesh_elas.template create_block_csr<block_size>(nrows, rowp, cols);
    bsr_mat = BSRMat_t(nrows, nrows, cols.size(), rowp, cols);
    csc_mat = bsr_to_csc(bsr_mat);

    // Set up Cholesky solver but don't set up values yet
    bool set_values = false;
    chol =
        new SparseCholesky(csc_mat, CholOrderingType::ND, nullptr, set_values);
  }

  template <class VecType>
  void set_design_vars(const VecType &xvec) {
    for (index_t i = 0; i < elem_datamesh_elas.get_num_dof(); i++) {
      data[i] = xvec[i];
    }
  }
  void set_design_vars(T xval) {
    for (index_t i = 0; i < elem_datamesh_elas.get_num_dof(); i++) {
      data[i] = xval;
    }
  }

  void factor() {
    // Create a new element view of the system matrix
    ElementMat_Serial<T, BasisElas, BSRMat_t> elem_mat(elem_mesh_elas, bsr_mat);

    // Populate the system matrix
    elem_ops_elas.add_jacobian_new(integrand_elas, elem_data_elas,
                                   elem_geo_elas, elem_sol_elas, elem_mat);

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
    ElementVector_Empty<ElemVecType::Serial> elem_data_trac;
    SolutionVector<T> traction_res(elem_mesh_elas.get_num_dof());
    ElemVecTrac elem_res_trac(elem_mesh_trac, traction_res);

    elem_ops_trac.add_residual(integrand_trac, elem_data_trac, elem_geo_trac,
                               elem_sol_trac, elem_res_trac);

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
    write_quad_to_vtk<3, degree, T, DataBasisElas, GeoBasisElas, BasisElas,
                      IntegrandElas>(
        elem_data_elas, elem_geo_elas, elem_sol_elas, filename,
        [](index_t k, typename IntegrandElas::DataSpace &d,
           typename IntegrandElas::FiniteElementGeometry &g,
           typename IntegrandElas::FiniteElementSpace &s) {
          if (k == 2) {  // write data
            return (d.template get<0>()).get_value();
          } else {  // write solution components
            auto u = (s.template get<0>()).get_value();
            return u(k);
          }
        });
  }

  ElementMesh<BasisElas> &get_mesh() { return elem_mesh_elas; }
  CSCMat_t &get_csc_matrix() { return csc_mat; }
  SparseCholesky<T> *get_chol() { return chol; }

 private:
  T E, nu, q;

  ElementMesh<BasisElas> elem_mesh_elas;
  ElementMesh<GeoBasisElas> elem_geomesh_elas;
  ElementMesh<DataBasisElas> elem_datamesh_elas;

  ElementMesh<BasisTrac> elem_mesh_trac;
  ElementMesh<GeoBasisTrac> elem_geomesh_trac;

  DirichletBCs<BasisElas> bcs;

  Vec_t sol;
  Vec_t geo;
  Vec_t data;

  ElemVecElas elem_sol_elas;
  GeoElemVecElas elem_geo_elas;
  DataElemVecElas elem_data_elas;

  ElemVecTrac elem_sol_trac;
  GeoElemVecTrac elem_geo_trac;

  IntegrandElas integrand_elas;
  IntegrandTrac integrand_trac;

  ElemOpsElas elem_ops_elas;
  ElemOpsTrac elem_ops_trac;

  // System matrices
  index_t nrows;
  std::vector<index_t> rowp, cols;
  BSRMat_t bsr_mat;
  CSCMat_t csc_mat;

  // Cholesky solver
  SparseCholesky<T> *chol;
};

#endif  // TOPOLOGY_2D_H