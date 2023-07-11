#ifndef TOPOLOGY_ELASTICITY_H
#define TOPOLOGY_ELASTICITY_H

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/lagrange_quad_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dprofiler.h"

/**
 * @brief Performs elasticity analysis for topology optimization.
 */
template <typename T, A2D::index_t degree, A2D::index_t filter_degree>
class TopoElasticityAnalysis {
 public:
  // Alias templates
  template <class... Args>
  using ElementVector = A2D::ElementVector_Parallel<Args...>;
  using ElementVectorEmpty =
      A2D::ElemenetVector_Empty<A2D::ElemVecType::Parallel>;

  // Basic types
  using I = A2D::index_t;

  // Magic integers
  static constexpr I spatial_dim = 3;  // spatial dimension
  static constexpr I var_dim = 3;      // dimension of the PDE solution variable
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

  // Quadrature, basis and element views for high-order surface traction
  using TractionQuadrature = A2D::QuadGaussQuadrature<degree + 1>;
  using TractionDataBasis = A2D::FEBasis<T>;
  using TractionGeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using TractionBasis =
      A2D::FEBasis<T, A2D::LagrangeH1QuadBasis<T, spatial_dim, degree>>;
  using TractionDataElemVec = ElementVector<T, TractionDataBasis, BasisVecType>;
  using TractionGeoElemVec = ElementVector<T, TractionGeoBasis, BasisVecType>;
  using TractionElemVec = ElementVector<T, TractionBasis, BasisVecType>;

  // Block compressed row sparse matrix
  using BSRMatType = A2D::BSRMat<T, block_size, block_size>;

  /* Problem specific types */

  // Problem PDE
  using PDE = A2D::TopoLinearElasticity<T, spatial_dim>;
  using BodyForce = A2D::TopoBodyForce<T, spatial_dim>;
  using TractionPDE = A2D::TopoSurfaceTraction<T, spatial_dim>;

  // Finite element functional
  using FE_PDE =
      A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_BodyForce =
      A2D::FiniteElement<T, BodyForce, Quadrature, DataBasis, GeoBasis, Basis>;
  using FE_Traction =
      A2D::FiniteElement<T, TractionPDE, TractionQuadrature, TractionDataBasis,
                         TractionGeoBasis, TractionBasis>;

  // Finite element functional for low order preconditioner mesh
  using LOrderFE = A2D::FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                      LOrderGeoBasis, LOrderBasis>;

  // Matrix-free operator
  using MatFree =
      A2D::MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Algebraic multigrid solver
  static constexpr I null_size = 6;
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
  using AggregationFunctional =
      A2D::FiniteElement<T, A2D::TopoVonMisesAggregation<T, spatial_dim>,
                         Quadrature, DataBasis, GeoBasis, Basis>;

  /**
   * @brief The elasticity topology analysis class.
   *
   * @param conn mesh connectivity
   * @param bcinfo boundary condition information
   * @param E Young's modulus
   * @param nu Poisson's ratio
   * @param q RAMP penalization parameter
   * @param tx_bodyforce the body force components
   * @param traction_label label of nodes to which the traction force is applied
   * @param tx_traction traction force components
   * @param x0_torque origin the torque vector is about
   * @param tx_torque the torque vector
   * @param verbose if executed verbosely
   * @param amg_nlevels number of algebraic multi-grid levels
   * @param cg_it max number of iterations for conjugate gradient linear solver
   * @param cg_rtol relative error tolerance for conjugate gradient solver
   * @param cg_atol absolute error tolerance for conjugate gradient solver
   */
  TopoElasticityAnalysis(A2D::MeshConnectivity3D &conn,
                         A2D::DirichletBCInfo &bcinfo, T E, T nu, T q,
                         const T tx_bodyforce[], A2D::index_t traction_label,
                         const T tx_traction[], const T x0_torque[],
                         const T tx_torque[], bool verbose, int amg_nlevels,
                         int cg_it, double cg_rtol,
                         double cg_atol)
      :  // Material parameters and penalization
        E(E),
        nu(nu),
        q(q),

        // Meshes for the solution, geometry and data
        mesh(conn),
        geomesh(conn),
        datamesh(conn),
        filtermesh(conn),

        traction_mesh(traction_label, conn, mesh),
        traction_geomesh(traction_label, conn, geomesh),

        // Boundary conditions
        bcs(conn, mesh, bcinfo),

        // Project the meshes onto the low-order meshes
        lorder_mesh(mesh),
        lorder_geomesh(geomesh),
        lorder_datamesh(datamesh),

        // Solution, geometry and data vectors
        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        data(datamesh.get_num_dof()),
        filter_data(filtermesh.get_num_dof()),

        // Element-level views of the solution geometry and data
        elem_sol(mesh, sol),
        elem_geo(geomesh, geo),
        elem_data(datamesh, data),
        elem_filter_data(filtermesh, filter_data),

        // Element view of the traction class
        elem_traction_sol(traction_mesh, sol),
        elem_traction_geo(traction_geomesh, geo),

        // Low-order views of the solution geometry and data
        lorder_elem_sol(lorder_mesh, sol),
        lorder_elem_geo(lorder_geomesh, geo),
        lorder_elem_data(lorder_datamesh, data),

        pde(E, nu, q),
        bodyforce(q, tx_bodyforce),
        traction_pde(tx_traction, tx_torque, x0_torque),

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

    // Create the matrix for the low-order mesh
    I nrows;
    std::vector<I> rowp, cols;
    lorder_mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    // Create the shared pointer
    mat = std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);
  }

  GeoElemVec &get_geometry() { return elem_geo; }

  void reset_geometry() {
    A2D::SolutionVector<T> x(mesh.get_num_dof());
    A2D::SolutionVector<T> y(mesh.get_num_dof());
    A2D::SolutionVector<T> z(mesh.get_num_dof());

    ElemVec elem_x(mesh, x);
    ElemVec elem_y(mesh, y);
    ElemVec elem_z(mesh, z);

    A2D::DOFCoordinates<T, PDE, GeoBasis, Basis> coords;
    coords.get_dof_coordinates(elem_geo, elem_x, elem_y, elem_z);

    // Initialize the near null-space to an appropriate vector
    for (I i = 0; i < B.extent(0); i++) {
      B(i, 0, 0) = 1.0;
      B(i, 1, 1) = 1.0;
      B(i, 2, 2) = 1.0;

      // Rotation about the x-axis
      B(i, 1, 3) = z[3 * i + 2];
      B(i, 2, 3) = -y[3 * i + 1];

      // Rotation about the y-axis
      B(i, 0, 4) = z[3 * i + 2];
      B(i, 2, 4) = -x[3 * i];

      // Rotation about the z-axis
      B(i, 0, 5) = y[3 * i + 1];
      B(i, 1, 5) = -x[3 * i];
    }

    // Zero out the boundary conditions
    const I *bc_dofs;
    I nbcs = bcs.get_bcs(&bc_dofs);
    for (I i = 0; i < nbcs; i++) {
      I dof = bc_dofs[i];
      for (I j = 0; j < null_size; j++) {
        B(dof / block_size, dof % block_size, j) = 0.0;
      }
    }
  }

  /**
   * @brief Solve the governing equations and set the new solution vector
   *
   */
  void solve() {
    A2D::Timer timer("TopoElasticityAnalysis::solve()");
    // Create a view of the low-order element matrix
    A2D::ElementMat_Serial<T, LOrderBasis, BSRMatType> elem_mat(lorder_mesh,
                                                                *mat);

    // Initialie the Jacobian matrix
    lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                           lorder_elem_sol, elem_mat);
    mat->write_mtx("low_order_Jacobian.mtx");  // TODO: delete this

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

    // Allocate the solver
    // TODO: we should add some of these as solver options
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

    // Zero the solution
    sol.zero();

    // No data associated with the traction class
    ElementVectorEmpty elem_traction_data;

    // Assemble the force contribution
    A2D::SolutionVector<T> traction_res(mesh.get_num_dof());
    TractionElemVec elem_traction_res(traction_mesh, traction_res);
    traction.add_residual(traction_pde, elem_traction_data, elem_traction_geo,
                          elem_traction_sol, elem_traction_res);

    // Assemble the body force contribution
    A2D::SolutionVector<T> res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    feb.add_residual(bodyforce, elem_data, elem_geo, elem_sol, elem_res);

    // Set the right-hand-side
    for (I i = 0; i < sol.get_num_dof(); i++) {
      rhs_vec(i / block_size, i % block_size) = -traction_res[i] - res[i];
    }

    // Zero out the boundary conditions
    for (A2D::index_t i = 0; i < nbcs; i++) {
      A2D::index_t dof = bc_dofs[i];
      rhs_vec(dof / block_size, dof % block_size) = 0.0;
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

    if constexpr (decltype(elem_filter_data)::evtype ==
                  A2D::ElemVecType::Parallel) {
      elem_filter_data.get_values();
    }

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Interpolate the data from the Bernstein filter
      typename FilterElemVec::FEDof filter_dof(i, elem_filter_data);
      if constexpr (decltype(elem_filter_data)::evtype ==
                    A2D::ElemVecType::Serial) {
        elem_filter_data.get_element_values(i, filter_dof);
      }

      // Interpolate from the filter degrees of freedom to the high-order
      // quadrature points
      A2D::QptSpace<FilterQuadrature, FilterSpace> qdata;
      FilterBasis::template interp(filter_dof, qdata);

      typename DataElemVec::FEDof data_dof(i, elem_data);

      // Set the data values into the data space for the high-order GLL mesh
      for (A2D::index_t j = 0; j < DataBasis::ndof; j++) {
        data_dof[j] = qdata.get(j)[0];
      }

      if constexpr (decltype(elem_data)::evtype == A2D::ElemVecType::Serial) {
        elem_data.set_element_values(i, data_dof);
      }
    }
    if constexpr (decltype(elem_data)::evtype == A2D::ElemVecType::Parallel) {
      elem_data.set_values();
    }
  }

  /**
   * @brief Add the values back to the design variable vector
   */
  template <A2D::ElemVecType evtype, class VecType, class DataDerivElemVec>
  void filter_add(A2D::ElementVectorBase<evtype, DataDerivElemVec> &elem_dfdx,
                  VecType &dfdx) {
    ElementVector<T, FilterBasis, VecType> filter_dfdx(filtermesh, dfdx);

    // Loop over the elements and interpolate the values to the refined data
    // mesh
    const A2D::index_t num_elements = elem_filter_data.get_num_elements();

    if constexpr (evtype == A2D::ElemVecType::Parallel) {
      elem_dfdx.get_values();
      filter_dfdx.get_zero_values();
    }
    for (A2D::index_t i = 0; i < num_elements; i++) {
      typename DataDerivElemVec::FEDof data_dof(i, elem_dfdx);
      if constexpr (evtype == A2D::ElemVecType::Serial) {
        elem_dfdx.get_element_values(i, data_dof);
      }

      // Set the data values into the data space for the high-order GLL mesh
      A2D::QptSpace<FilterQuadrature, FilterSpace> qdata;
      for (A2D::index_t j = 0; j < DataBasis::ndof; j++) {
        qdata.get(j)[0] = data_dof[j];
      }

      // Add the derivative from the Bernstein filter
      typename ElementVector<T, FilterBasis, VecType>::FEDof filter_dof(
          i, filter_dfdx);
      FilterBasis::template add(qdata, filter_dof);
      if constexpr (evtype == A2D::ElemVecType::Serial) {
        filter_dfdx.add_element_values(i, filter_dof);
      }
    }
    if constexpr (evtype == A2D::ElemVecType::Parallel) {
      filter_dfdx.add_values();
    }
  }

  /**
   * @brief Get the number of design variables
   *
   * @return The number of design variables
   */
  A2D::index_t get_num_design_vars() { return filtermesh.get_num_dof(); }

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
    A2D::Timer timer("TopoElasticityAnalysis::eval_compliance()");
    return fe.integrate(pde, elem_data, elem_geo, elem_sol);
  }

  /**
   * @brief Evaluate the compliance gradient
   *
   */
  template <class VecType>
  void add_compliance_gradient(VecType &dfdx) {
    A2D::Timer timer("TopoElasticityAnalysis::add_compliance_gradient()");

    BasisVecType dfdrho(datamesh.get_num_dof());
    DataElemVec elem_dfdrho(datamesh, dfdrho);

    BasisVecType adjoint(mesh.get_num_dof());
    for (I i = 0; i < adjoint.get_num_dof(); i++) {
      adjoint[i] = -0.5 * sol[i];
    }
    ElemVec elem_adjoint(mesh, adjoint);

    fe.add_adjoint_residual_data_derivative(pde, elem_data, elem_geo, elem_sol,
                                            elem_adjoint, elem_dfdrho);

    for (I i = 0; i < adjoint.get_num_dof(); i++) {
      adjoint[i] = -sol[i];
    }
    feb.add_adjoint_residual_data_derivative(
        bodyforce, elem_data, elem_geo, elem_sol, elem_adjoint, elem_dfdrho);
    filter_add(elem_dfdrho, dfdx);
  }

  /**
   * @brief Evaluate the volume
   *
   * @return The volume of the parametrized topology
   */
  T eval_volume() {
    A2D::Timer timer("TopoElasticityAnalysis::eval_volume()");
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
    A2D::Timer timer("TopoElasticityAnalysis::add_volume_gradient()");
    VolumeFunctional functional;
    VolumePDE volume;

    // Create the element-view for the derivative
    BasisVecType dfdrho(datamesh.get_num_dof());
    DataElemVec elem_dfdrho(datamesh, dfdrho);
    functional.add_data_derivative(volume, elem_data, elem_geo, elem_sol,
                                   elem_dfdrho);

    filter_add(elem_dfdrho, dfdx);
  }

  /**
   * @brief Evaluate the aggregated stress over the domain
   *
   * @param design_stress The design stress value
   * @param ks_penalty The KS penalty parameter
   * @return The aggregation functional value
   */
  T eval_aggregation(T design_stress, T ks_penalty) {
    A2D::Timer timer("TopoElasticityAnalysis::eval_aggregation()");
    AggregationFunctional functional;
    A2D::TopoVonMisesAggregation<T, spatial_dim> aggregation(
        E, nu, q, design_stress, ks_penalty);

    T max_value = functional.max(aggregation, elem_data, elem_geo, elem_sol);
    aggregation.set_max_failure_index(max_value);

    T integral =
        functional.integrate(aggregation, elem_data, elem_geo, elem_sol);
    T value = aggregation.evaluate_functional(integral);

    return value;
  }

  /**
   * @brief Evaluate the aggregated stress over the domain
   *
   * @param design_stress The design stress value
   * @param ks_penalty The KS penalty parameter
   * @return The aggregation functional value
   */
  template <class VecType>
  void add_aggregation_gradient(T design_stress, T ks_penalty, VecType &dfdx) {
    A2D::Timer timer("TopoElasticityAnalysis::add_aggregation_gradient()");
    AggregationFunctional functional;
    A2D::TopoVonMisesAggregation<T, spatial_dim> aggregation(
        E, nu, q, design_stress, ks_penalty);

    T max_value = functional.max(aggregation, elem_data, elem_geo, elem_sol);
    aggregation.set_max_failure_index(max_value);

    T integral =
        functional.integrate(aggregation, elem_data, elem_geo, elem_sol);
    T value = aggregation.evaluate_functional(integral);

    // Create the element-view for the derivative
    BasisVecType dfdrho(datamesh.get_num_dof());
    DataElemVec elem_dfdrho(datamesh, dfdrho);

    functional.add_data_derivative(aggregation, elem_data, elem_geo, elem_sol,
                                   elem_dfdrho);

    A2D::SolutionVector<T> dfdu(mesh.get_num_dof());
    ElementVector<T, Basis, A2D::SolutionVector<T>> elem_dfdu(mesh, dfdu);

    // Set up and solve the adjoint equations...
    functional.add_residual(aggregation, elem_data, elem_geo, elem_sol,
                            elem_dfdu);

    //  Create a view of the low-order element matrix
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

    // Set the right-hand-side
    for (I i = 0; i < dfdu.get_num_dof(); i++) {
      rhs_vec(i / block_size, i % block_size) = -dfdu[i];
    }

    // Zero out the boundary conditions
    for (I i = 0; i < nbcs; i++) {
      I dof = bc_dofs[i];
      rhs_vec(dof / block_size, dof % block_size) = 0.0;
    }

    // Solve the problem
    I monitor = 0;
    if (verbose) {
      monitor = 5;
    }
    bool succ =
        amg.cg(mat_vec, rhs_vec, sol_vec, monitor, cg_it, cg_rtol, cg_atol);
    if (!succ) {
      throw std::runtime_error("CG failed to converge!");
    }

    // Record the solution back into the right-hand-side vector
    for (I i = 0; i < sol.get_num_dof(); i++) {
      dfdu[i] = sol_vec(i / block_size, i % block_size);
    }

    fe.add_adjoint_residual_data_derivative(pde, elem_data, elem_geo, elem_sol,
                                            elem_dfdu, elem_dfdrho);
    feb.add_adjoint_residual_data_derivative(bodyforce, elem_data, elem_geo,
                                             elem_sol, elem_dfdu, elem_dfdrho);
    filter_add(elem_dfdrho, dfdx);
  }

  void tovtk(const std::string filename) {
    A2D::write_hex_to_vtk<4, degree, T, DataBasis, GeoBasis, Basis>(
        pde, elem_data, elem_geo, elem_sol, filename,
        [](I k, typename PDE::DataSpace &d,
           typename PDE::FiniteElementGeometry &g,
           typename PDE::FiniteElementSpace &s) {
          if (k == 3) {
            return (d.template get<0>()).get_value();
          } else {
            auto u = (s.template get<0>()).get_value();
            return u(k);
          }
        });
  }

 private:
  T E, nu, q;

  A2D::ElementMesh<Basis> mesh;
  A2D::ElementMesh<GeoBasis> geomesh;
  A2D::ElementMesh<DataBasis> datamesh;
  A2D::ElementMesh<FilterBasis> filtermesh;
  A2D::ElementMesh<TractionBasis> traction_mesh;
  A2D::ElementMesh<TractionGeoBasis> traction_geomesh;

  A2D::DirichletBCs<Basis> bcs;

  A2D::ElementMesh<LOrderBasis> lorder_mesh;
  A2D::ElementMesh<LOrderGeoBasis> lorder_geomesh;
  A2D::ElementMesh<LOrderDataBasis> lorder_datamesh;

  A2D::SolutionVector<T> sol;
  A2D::SolutionVector<T> geo;
  A2D::SolutionVector<T> data;
  A2D::SolutionVector<T> filter_data;

  ElemVec elem_sol;
  GeoElemVec elem_geo;
  DataElemVec elem_data;
  FilterElemVec elem_filter_data;
  TractionElemVec elem_traction_sol;
  TractionGeoElemVec elem_traction_geo;

  LOrderElemVec lorder_elem_sol;
  LOrderGeoElemVec lorder_elem_geo;
  LOrderDataElemVec lorder_elem_data;

  PDE pde;
  BodyForce bodyforce;
  TractionPDE traction_pde;

  FE_PDE fe;
  FE_BodyForce feb;
  FE_Traction traction;

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

#endif