#ifndef TOPOLOGY_HEAT_H
#define TOPOLOGY_HEAT_H

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dprofiler.h"

/**
 * @brief Performs heat conduction analysis for topology optimization.
 */
template <typename T, A2D::index_t degree, A2D::index_t filter_degree>
class TopoHeatAnalysis {
 public:
  // Alias templates
  template <class... Args>
  using ElementVector = A2D::ElementVector_Serial<Args...>;

  // Basic types
  using I = A2D::index_t;

  // Magic integers
  static constexpr I spatial_dim = 3;  // spatial dimension
  static constexpr I var_dim = 1;      // dimension of the PDE solution variable
  static constexpr I data_dim = 1;     // dimension of material data
  static constexpr I low_degree = 1;   // low order preconditinoer mesh degree
  static constexpr I block_size = var_dim;  // block size for BSR matrix

  // The type of solution vector to use
  using BasisVecType = A2D::SolutionVector<T>;

  // Quadrature, basis and element views for original mesh
  using Quadrature = A2D::HexGaussQuadrature<degree + 1>;
  using DataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, degree - 1>>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using DataElemVec = ElementVector<T, DataBasis, BasisVecType>;
  using GeoElemVec = ElementVector<T, GeoBasis, BasisVecType>;
  using ElemVec = ElementVector<T, Basis, BasisVecType>;

  // Quadrature, basis and element views for low order preconditioner mesh
  using LOrderQuadrature = A2D::HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, low_degree - 1>>;
  using LOrderGeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, low_degree>>;
  using LOrderBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, low_degree>>;
  using LOrderDataElemVec = ElementVector<T, LOrderDataBasis, BasisVecType>;
  using LOrderGeoElemVec = ElementVector<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = ElementVector<T, LOrderBasis, BasisVecType>;

  // Block compressed row sparse matrix
  using BSRMatType = A2D::BSRMat<T, block_size, block_size>;

  /* Problem specific types */

  // Problem PDE
  using PDE = A2D::HeatConduction<T, spatial_dim>;
  using AdjRHS = A2D::AdjRHS<T, spatial_dim>;

  using FE_PDE =
      A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_AdjRHS =
      A2D::FiniteElement<T, AdjRHS, Quadrature, DataBasis, GeoBasis, Basis>;

  // Finite element functional for low order preconditioner mesh
  using LOrderFE = A2D::FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                      LOrderGeoBasis, LOrderBasis>;

  // Matrix-free operator
  using MatFree =
      A2D::MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Algebraic multigrid solver
  static constexpr I null_size = 1;
  using BSRMatAmgType = A2D::BSRMatAmg<T, block_size, null_size>;

  // Filter information
  // Use the Gauss quadrature points here so that the filter can be evaluated at
  // the Gauss points for the DataBasis
  using FilterQuadrature = A2D::HexGaussQuadrature<degree>;

  // Use a continuous H1 space for the filter
  using FilterSpace =
      A2D::FESpace<T, spatial_dim, A2D::H1Space<T, data_dim, spatial_dim>>;

  // Use Bernstein points for the design vector on a lower-order mesh
  using FilterBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, data_dim, filter_degree,
                                              A2D::BERNSTEIN_INTERPOLATION>>;

  // The design variable element mesh
  using FilterElemVec = ElementVector<T, FilterBasis, BasisVecType>;

  // Functional definitions
  using VolumePDE = A2D::TopoVolume<T, var_dim, spatial_dim, PDE>;

  using VolumeFunctional =
      A2D::FiniteElement<T, VolumePDE, Quadrature, DataBasis, GeoBasis, Basis>;

  /**
   * @brief The heat conduction topology analysis class.
   *
   * @param conn mesh connectivity
   * @param bcinfo boundary condition information
   * @param kappa heat conduction coefficient
   * @param heat_source source term value
   * @param bc_temp boundary temperature
   * @param q RAMP parameter
   * @param verbose if execute verbosely
   * @param amg_nlevels number of algebraic multi-grid levels
   * @param cg_it max number of iterations for conjugate gradient linear solver
   * @param cg_rtol relative error tolerance for conjugate gradient solver
   * @param cg_atol absolute error tolerance for conjugate gradient solver
   */
  TopoHeatAnalysis(A2D::MeshConnectivity3D &conn, A2D::DirichletBCInfo &bcinfo,
                   T kappa, T heat_source, T bc_temp, T q, bool verbose,
                   int amg_nlevels, int cg_it, double cg_rtol,
                   double cg_atol)
      :  // Material parameters and penalization
        kappa(kappa),
        heat_source(heat_source),
        bc_temp(bc_temp),
        q(q),

        // Meshes for the solution, geometry and data
        mesh(conn),
        geomesh(conn),
        datamesh(conn),
        filtermesh(conn),

        bcs(conn, mesh, bcinfo),

        // Project the meshes onto the low-order meshes
        lorder_mesh(mesh),
        lorder_geomesh(geomesh),
        lorder_datamesh(datamesh),

        // Solution, adjoint, geometry and data vectors
        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        data(datamesh.get_num_dof()),
        filter_data(filtermesh.get_num_dof()),

        // Element-level views of the solution geometry and data
        elem_data(datamesh, data),
        elem_geo(geomesh, geo),
        elem_sol(mesh, sol),
        elem_filter_data(filtermesh, filter_data),

        // Low-order views of the solution geometry and data
        lorder_elem_data(lorder_datamesh, data),
        lorder_elem_geo(lorder_geomesh, geo),
        lorder_elem_sol(lorder_mesh, sol),

        pde(kappa, q, heat_source),
        adjrhs(heat_source),

        B("B", sol.get_num_dof() / block_size),

        verbose(verbose),
        amg_nlevels(amg_nlevels),
        cg_it(cg_it),
        cg_rtol(cg_rtol),
        cg_atol(cg_atol) {
    // Initialize the data
    for (I i = 0; i < data.get_num_dof(); i++) {
      data[i] = 1.0;
    }

    // Create near null sapce
    A2D::BLAS::fill(B, T(1.0));

    // Create the matrix for the low-order mesh
    I nrows;
    std::vector<I> rowp, cols;
    lorder_mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    // Create the shared pointer
    mat = std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);
  }

  GeoElemVec &get_geometry() { return elem_geo; }

  /**
   * @brief This is not used for heat conduction problem hence near null space B
   * is always the vector of ones.
   */
  void reset_geometry() {}

  /**
   * @brief Solve the governing equations and set the new solution vector
   *
   */
  void solve() {
    A2D::Timer timer("TopoHeatAnalysis::solve()");

    // Create a view of the low-order element matrix
    A2D::ElementMat_Serial<T, LOrderBasis, BSRMatType> elem_mat(lorder_mesh,
                                                                *mat);

    // Initialie the Jacobian matrix
    lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                           lorder_elem_sol, elem_mat);

    // Apply the boundary conditions
    const I *bc_dofs;
    I nbcs = bcs.get_bcs(&bc_dofs);
    mat->zero_rows(nbcs, bc_dofs);

    // Initialize the matrix-free data
    matfree.initialize(pde, elem_data, elem_geo, elem_sol);

    // Allocate space for temporary variables with the matrix-vector code
    A2D::SolutionVector<T> xvec(mesh.get_num_dof());
    A2D::SolutionVector<T> yvec(mesh.get_num_dof());
    ElemVec elem_xvec(mesh, xvec);
    ElemVec elem_yvec(mesh, yvec);

    auto mat_vec = [&](A2D::MultiArrayNew<T *[block_size]> &in,
                       A2D::MultiArrayNew<T *[block_size]> &out) -> void {
      xvec.zero();
      yvec.zero();
      for (I i = 0; i < xvec.get_num_dof(); i++) {
        xvec[i] = in(i / block_size, i % block_size);
      }
      matfree.add_jacobian_vector_product(elem_xvec, elem_yvec);

      for (I i = 0; i < yvec.get_num_dof(); i++) {
        out(i / block_size, i % block_size) = yvec[i];
      }

      // Set the boundary conditions as equal to the inputs
      const I *bc_dofs;
      I nbcs = bcs.get_bcs(&bc_dofs);
      for (I i = 0; i < nbcs; i++) {
        I dof = bc_dofs[i];

        out(dof / block_size, dof % block_size) =
            in(dof / block_size, dof % block_size);
      }
    };

    // Allocate the solver - we should add some of these as solver options
    double omega = 4.0 / 3.0;
    double epsilon = 0.0;
    bool print_info = false;
    if (verbose) {
      print_info = true;
    }
    BSRMatAmgType amg(amg_nlevels, omega, epsilon, mat, B, print_info);

    // Create the solution and right-hand-side vectors
    I size = sol.get_num_dof() / block_size;
    A2D::MultiArrayNew<T *[block_size]> sol_vec("sol_vec", size);
    A2D::MultiArrayNew<T *[block_size]> rhs_vec("rhs_vec", size);

    A2D::SolutionVector<T> res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    sol.zero();
    fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

    for (I i = 0; i < sol.get_num_dof(); i++) {
      rhs_vec(i / block_size, i % block_size) = -res[i];
    }

    // Zero out the boundary conditions
    for (I i = 0; i < nbcs; i++) {
      I dof = bc_dofs[i];
      rhs_vec(dof / block_size, dof % block_size) = bc_temp;
    }

    // Solve the problem
    I monitor = 0;
    if (verbose) {
      monitor = 5;
    }

    bool succ =
        amg.cg(mat_vec, rhs_vec, sol_vec, monitor, cg_it, cg_rtol, cg_atol);
    if (!succ) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "%s:%d: CG failed to converge after %d iterations given "
                    "rtol=%.1e, atol=%.1e",
                    __FILE__, __LINE__, cg_it, cg_rtol, cg_atol);
      throw std::runtime_error(msg);
    }

    // Record the solution
    for (I i = 0; i < sol.get_num_dof(); i++) {
      sol[i] = sol_vec(i / block_size, i % block_size);
    }
  }

  /**
   * @brief Interpolate to the data space from the design variables
   */
  template <class VecType>
  void filter_interp(const VecType &xvec) {
    // Copy the design variables to the filter data
    for (A2D::index_t i = 0; i < filtermesh.get_num_dof(); i++) {
      filter_data[i] = xvec[i];
    }

    // Loop over the elements and interpolate the values to the refined data
    // mesh
    const A2D::index_t num_elements = elem_filter_data.get_num_elements();
    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Interpolate the data from the Bernstein filter
      typename FilterElemVec::FEDof filter_dof(i, elem_filter_data);
      elem_filter_data.get_element_values(i, filter_dof);

      // Interpolate from the filter degrees of freedom to the high-order
      // quadrature points
      A2D::QptSpace<FilterQuadrature, FilterSpace> qdata;
      FilterBasis::template interp(filter_dof, qdata);

      typename DataElemVec::FEDof data_dof(i, elem_data);

      // Set the data values into the data space for the high-order GLL mesh
      for (A2D::index_t j = 0; j < DataBasis::ndof; j++) {
        data_dof[j] = qdata.get(j)[0];
      }

      elem_data.set_element_values(i, data_dof);
    }
  }

  /**
   * @brief Add the values back to the design variable vector
   */
  template <class VecType, class DataDerivElemVec>
  void filter_add(DataDerivElemVec &elem_dfdx, VecType &dfdx) {
    ElementVector<T, FilterBasis, VecType> filter_dfdx(filtermesh, dfdx);

    // Loop over the elements and interpolate the values to the refined data
    // mesh
    const A2D::index_t num_elements = elem_filter_data.get_num_elements();
    for (A2D::index_t i = 0; i < num_elements; i++) {
      typename DataDerivElemVec::FEDof data_dof(i, elem_dfdx);
      elem_dfdx.get_element_values(i, data_dof);

      // Set the data values into the data space for the high-order GLL mesh
      A2D::QptSpace<FilterQuadrature, FilterSpace> qdata;
      for (A2D::index_t j = 0; j < DataBasis::ndof; j++) {
        qdata.get(j)[0] = data_dof[j];
      }

      // Add the derivative from the Bernstein filter
      typename ElementVector<T, FilterBasis, VecType>::FEDof filter_dof(
          i, filter_dfdx);
      FilterBasis::template add(qdata, filter_dof);
      filter_dfdx.add_element_values(i, filter_dof);
    }
  }

  /**
   * @brief Get the number of design variables
   */
  I get_num_design_vars() { return filtermesh.get_num_dof(); }

  /**
   * @brief Get the number of degrees of freedom
   */
  I get_num_dofs() { return mesh.get_num_dof(); }

  /**
   * @brief Get the number of finite elements
   */
  I get_num_elements() { return mesh.get_num_elements(); }

  /**
   * @brief Set the design variable values into the topology
   *
   * @tparam VecType The type of design vector
   * @param xvec The vector of design variable values
   */
  template <class VecType>
  void set_design_vars(const VecType &xvec) {
    filter_interp(xvec);
  }

  /**
   * @brief Evaluate the compliance
   */
  T eval_compliance() {
    A2D::Timer timer("TopoHeatAnalysis::eval_compliance()");
    return fe.integrate(pde, elem_data, elem_geo, elem_sol);
  }

  /**
   * @brief Evaluate the compliance gradient
   *
   */
  template <class VecType>
  void add_compliance_gradient(VecType &dfdx) {
    A2D::Timer timer("TopoHeatAnalysis::add_compliance_gradient()");
    BasisVecType dfdrho(datamesh.get_num_dof());
    DataElemVec elem_dfdrho(datamesh, dfdrho);

    BasisVecType adjoint(mesh.get_num_dof());
    for (I i = 0; i < adjoint.get_num_dof(); i++) {
      adjoint[i] = -sol[i];
    }
    ElemVec elem_adjoint(mesh, adjoint);

    fe.add_adjoint_residual_data_derivative(pde, elem_data, elem_geo, elem_sol,
                                            elem_adjoint, elem_dfdrho);
    filter_add(elem_dfdrho, dfdx);
  }

  /**
   * @brief Evaluate the volume
   *
   * @return The volume of the parametrized topology
   */
  T eval_volume() {
    A2D::Timer timer("TopoHeatAnalysis::eval_volume()");
    VolumeFunctional functional;
    VolumePDE volume;
    T vol = functional.integrate(volume, elem_data, elem_geo, elem_sol);

    return vol;
  }

  /**
   * @brief Add the volume gradient to the specified functional
   *
   * @tparam VecType The derivative vector type
   * @param dfdx The derivative value
   */
  template <class VecType>
  void add_volume_gradient(VecType &dfdx) {
    A2D::Timer timer("TopoHeatAnalysis::add_volume_gradient()");
    VolumeFunctional functional;
    VolumePDE volume;

    // Create the element-view for the derivative
    BasisVecType dfdrho(datamesh.get_num_dof());
    DataElemVec elem_dfdrho(datamesh, dfdrho);
    functional.add_data_derivative(volume, elem_data, elem_geo, elem_sol,
                                   elem_dfdrho);

    filter_add(elem_dfdrho, dfdx);
  }

  void tovtk(const std::string filename) {
    A2D::write_hex_to_vtk<2, degree, T, DataBasis, GeoBasis, Basis>(
        pde, elem_data, elem_geo, elem_sol, filename,
        [](I k, typename PDE::DataSpace &d,
           typename PDE::FiniteElementGeometry &g,
           typename PDE::FiniteElementSpace &s) {
          if (k == 0) {
            return (s.template get<0>()).get_value();  // state
          } else {
            return (d.template get<0>()).get_value();  // design
          }
        });
  }

 private:
  T kappa, heat_source, bc_temp, q;

  A2D::ElementMesh<Basis> mesh;
  A2D::ElementMesh<GeoBasis> geomesh;
  A2D::ElementMesh<DataBasis> datamesh;
  A2D::ElementMesh<FilterBasis> filtermesh;

  A2D::DirichletBCs<Basis> bcs;

  A2D::ElementMesh<LOrderBasis> lorder_mesh;
  A2D::ElementMesh<LOrderGeoBasis> lorder_geomesh;
  A2D::ElementMesh<LOrderDataBasis> lorder_datamesh;

  A2D::SolutionVector<T> sol;
  A2D::SolutionVector<T> geo;
  A2D::SolutionVector<T> data;
  A2D::SolutionVector<T> filter_data;

  DataElemVec elem_data;
  GeoElemVec elem_geo;
  ElemVec elem_sol;
  FilterElemVec elem_filter_data;

  LOrderDataElemVec lorder_elem_data;
  LOrderGeoElemVec lorder_elem_geo;
  LOrderElemVec lorder_elem_sol;

  PDE pde;
  AdjRHS adjrhs;

  FE_PDE fe;
  FE_AdjRHS fea;

  LOrderFE lorder_fe;
  MatFree matfree;

  // The near null-space to an appropriate vector
  A2D::MultiArrayNew<T *[block_size][null_size]> B;

  // System matrix
  std::shared_ptr<BSRMatType> mat;

  // If we print detailed info to stdout
  bool verbose;

  // AMG settings
  int amg_nlevels;

  // CG settings
  int cg_it;
  double cg_rtol, cg_atol;
};

void test_heat_analysis(int argc, char *argv[]) {
  // Magic numbers
  const int degree = 1;            // polynomial degree, 1 = linear
  const int order = 2;             // order = degree + 1
  const int spatial_dim = 3;       // spatial dimension
  const int var_dim = 1;           // solution variable dimension
  const int data_dim = 1;          // design variable dimension
  const int block_size = var_dim;  // block size for BSR matrix

  /* Types */

  // Basic type
  using T = double;
  using I = A2D::index_t;

  // Quadrature and basis
  using Quadrature = A2D::HexGaussQuadrature<order>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using DataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, degree - 1>>;

  // Block sparse compressed row matrix
  using BSRMatType = A2D::BSRMat<T, block_size, block_size>;

  // Element-centric views
  using GlobalVecType = A2D::SolutionVector<T>;
  using ElemVec = A2D::ElementVector_Serial<T, Basis, GlobalVecType>;
  using GeoElemVec = A2D::ElementVector_Serial<T, GeoBasis, GlobalVecType>;
  using DataElemVec = A2D::ElementVector_Serial<T, DataBasis, GlobalVecType>;
  using ElemMat = A2D::ElementMat_Serial<T, Basis, BSRMatType>;

  // Physics and functional
  using PDE = A2D::HeatConduction<T, spatial_dim>;
  using FE = A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  /* Load mesh and boundary vertices from vtk */

  // Load vtk
  std::string vtk_name = "3d_hex_1.vtk";
  if (argc > 1) {
    vtk_name = argv[1];
  }
  A2D::ReadVTK3D<I, T> readvtk(vtk_name);

  // Get connectivity for each element type
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;

  // Construct connectivity
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Extract boundary vertices
  std::vector<int> ids{100};
  std::vector<I> bc_verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Construct A2D mesh objects
  A2D::ElementMesh<Basis> mesh(conn);
  A2D::ElementMesh<GeoBasis> geomesh(conn);
  A2D::ElementMesh<DataBasis> datamesh(conn);

  /* Construct solution, data and X vectors and their element views */

  // Allocate vectors
  GlobalVecType sol(mesh.get_num_dof());
  GlobalVecType res(mesh.get_num_dof());
  GlobalVecType data(datamesh.get_num_dof());
  GlobalVecType geo(geomesh.get_num_dof());

  // Allocate elem views
  ElemVec elem_sol(mesh, sol);
  ElemVec elem_res(mesh, res);
  DataElemVec elem_data(datamesh, data);
  GeoElemVec elem_geo(geomesh, geo);

  // Populate data and geo
  data.fill(1.0);
  A2D::set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  /* Construct the Jacobian matrix */

  // Construct the nonzero pattern
  I nrows;
  std::vector<I> rowp, cols;
  mesh.create_block_csr<block_size>(nrows, rowp, cols);
  std::shared_ptr<BSRMatType> mat =
      std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);

  // Populate the stiffness matrix via element view
  ElemMat elem_mat(mesh, *mat);
  FE fe;
  T kappa = 1.0, q = 0.0, heat_source = 1.0, bc_temp = 0.0;
  PDE pde(kappa, q, heat_source);
  fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

  mat->write_mtx("heat_jacobian_nobc.mtx");

  // Apply boundary condition
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data()));
  A2D::DirichletBCs<Basis> bcs(conn, mesh, bcinfo);
  const I *bc_dofs;
  I nbcs = bcs.get_bcs(&bc_dofs);
  mat->zero_rows(nbcs, bc_dofs);

  // Write the Jacobian matrix to mtx for visualization
  mat->write_mtx("heat_jacobian_withbc.mtx");

  /* Construct the residual and rhs */

  fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

  I size = sol.get_num_dof() / block_size;
  A2D::MultiArrayNew<T *[block_size]> sol_vec("sol_vec", size);
  A2D::MultiArrayNew<T *[block_size]> rhs_vec("rhs_vec", size);

  for (I i = 0; i < sol.get_num_dof(); i++) {
    rhs_vec(i / block_size, i % block_size) = -res[i];
  }

  // Set the temperature boundary conditions
  for (I i = 0; i < nbcs; i++) {
    I dof = bc_dofs[i];
    rhs_vec(dof / block_size, dof % block_size) = bc_temp;
  }

  /* Construct the AMG solver and solve the problem */

  // Create amg
  A2D::index_t num_levels = 3;
  double omega = 4.0 / 3.0;
  double epsilon = 0.0;
  bool print_info = true;
  const int null_size = 1;
  A2D::MultiArrayNew<T *[block_size][null_size]> B(
      "B", sol.get_num_dof() / block_size);
  A2D::BLAS::fill(B, 1.0);
  A2D::BSRMatAmg<T, block_size, null_size> amg(num_levels, omega, epsilon, mat,
                                               B, print_info);

  // Solve
  auto mat_vec = [&](A2D::MultiArrayNew<T *[block_size]> &in,
                     A2D::MultiArrayNew<T *[block_size]> &out) -> void {
    A2D::BSRMatVecMult<T, block_size, block_size>(*mat, in, out);
  };

  amg.cg(mat_vec, rhs_vec, sol_vec, 5, 100);

  // Record the solution
  for (I i = 0; i < sol.get_num_dof(); i++) {
    sol[i] = sol_vec(i / block_size, i % block_size);
  }

  // Write result to vtk
  A2D::write_hex_to_vtk<2, degree, T, DataBasis, GeoBasis, Basis>(
      pde, elem_data, elem_geo, elem_sol, "heat_analysis.vtk",
      [](A2D::index_t k, typename PDE::DataSpace &d,
         typename PDE::FiniteElementGeometry &g,
         typename PDE::FiniteElementSpace &s) {
        if (k == 0) {
          return (s.template get<0>()).get_value();  // state
        } else {
          return (d.template get<0>()).get_value();  // design
        }
      });
}

#endif