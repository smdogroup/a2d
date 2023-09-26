#ifndef A2D_ANALSYS_H
#define A2D_ANALSYS_H

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/febase.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"

namespace A2D {

template <typename T, index_t block_size>
class DirectCholeskyAnalysis {
 public:
  // Set the implementation type
  class Impl_t {
   public:
    using type = T;
    using Vec_t = SolutionVector<T>;
    using Mat_t = BSRMat<T, block_size, block_size>;

    template <class Basis>
    using ElementVector = ElementVector_Serial<type, Basis, Vec_t>;

    template <class Basis>
    using ElementMatrix = ElementMat_Serial<type, Basis, Mat_t>;
  };

  // Data types
  using Vec_t = typename Impl_t::Vec_t;
  using Mat_t = typename Impl_t::Mat_t;
  using CSCMat_t = CSCMat<T>;

  DirectCholeskyAnalysis(index_t ndata, index_t ngeo, index_t ndof,
                         std::shared_ptr<ElementAssembler<Impl_t>> assembler,
                         std::shared_ptr<DirichletBCs<T>> bcs)
      : ndata(ndata),
        ngeo(ngeo),
        ndof(ndof),
        assembler(assembler),
        bcs(bcs),
        data(ndata),
        geo(ngeo),
        sol(ndof),
        res(ndof) {
    // Create the system matrix
    index_t nrows;
    std::vector<index_t> rowp, cols;
    assembler->get_bsr_data(block_size, nrows, rowp, cols);
    bsr_mat = std::make_shared<Mat_t>(nrows, nrows, cols.size(), rowp, cols);

    // Create the sparse matrix
    csc_mat = bsr_to_csc(*bsr_mat);

    // Set up Cholesky solver but don't set up values yet
    bool set_values = false;
    chol =
        new SparseCholesky(csc_mat, CholOrderingType::ND, nullptr, set_values);

    // Set the boundary conditions
    bcs->set_bcs(sol);
  }

  /**
   * @brief Get the data vector object
   *
   * @return Vec_t& Reference to the data
   */
  Vec_t &get_data() { return data; }

  /**
   * @brief Get the geometry vector object
   *
   * @return Vec_t& Reference to the data
   */
  Vec_t &get_geo() { return geo; }

  /**
   * @brief Get the solution vector object
   *
   * @return Vec_t& Reference to the data
   */
  Vec_t &get_sol() { return sol; }

  /**
   * @brief Compute the residual and store it
   */
  void residual() {
    res.zero();
    assembler->add_residual(T(1.0), data, geo, sol, res);

    // Apply boundary conditions
    bcs->zero_bcs(res);
  }

  /**
   * @brief Compute the Jacobian and store it
   */
  void jacobian() {
    // Compute the BSR matrix
    bsr_mat->zero();
    assembler->add_jacobian(T(1.0), data, geo, sol, *bsr_mat);

    // Zero the rows corresponding to the boundary conditions
    const index_t *bc_dofs;
    index_t nbcs = bcs->get_bcs(&bc_dofs);
    bsr_mat->zero_rows(nbcs, bc_dofs);
  }

  /**
   * @brief Factor the Jacobian
   */
  void factor() {
    // This is a problem here -- matrix is passed by value, but should be passed
    // by reference to avoid copying. Convert bsr to csc because we want to use
    // Cholesky factorization
    csc_mat = bsr_to_csc(*bsr_mat);

    // Apply boundary conditions to each column
    const index_t *bc_dofs;
    index_t nbcs = bcs->get_bcs(&bc_dofs);
    csc_mat.zero_columns(nbcs, bc_dofs);

    // Set values to Cholesky solver and factorize
    chol->setValues(csc_mat);
    chol->factor();
  }

  /**
   * @brief Perform a linear solve to obtain the solution
   */
  void linear_solve() {
    residual();
    jacobian();
    factor();

    // Find the solution
    chol->solve(res.data());

    // Apply the boundary condition values
    bcs->set_bcs(sol);

    // Update the solution
    for (index_t i = 0; i < sol.get_num_dof(); i++) {
      sol[i] -= res[i];
    }
  }

  void nonlinear_solve() {}

  T evaluate(FunctionalBase<Impl_t> &func) {
    return func.evaluate(data, geo, sol);
  }

  void add_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt, T alpha,
                      Vec_t &dfdx) {
    func.add_derivative(wrt, alpha, data, geo, sol, dfdx);
  }

  void eval_adjoint_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt,
                               Vec_t &dfdx) {
    // This doesn't make sense for the adjoint method
    if (wrt == FEVarType::STATE) {
      return;
    }

    // Add the contributions to the derivative from the partial
    dfdx.zero();
    func.add_derivative(wrt, T(1.0), data, geo, sol, dfdx);

    // Compute the derivative of the function wrt state and solve the adjoint
    // equations
    res.zero();
    func.add_derivative(FEVarType::STATE, T(-1.0), data, geo, sol, res);
    chol->solve(res.data());

    // Add the terms from the total derivative
    assembler->add_adjoint_res_product(wrt, T(1.0), data, geo, sol, res, dfdx);
  }

 private:
  // The solution information
  index_t ndata, ngeo, ndof;

  // Element matrix assembler object
  std::shared_ptr<ElementAssembler<Impl_t>> assembler;

  // Boundary conditions
  std::shared_ptr<DirichletBCs<T>> bcs;

  // Vectors of data, geometry and solution degrees of freedom
  Vec_t data;
  Vec_t geo;
  Vec_t sol;
  Vec_t res;

  // System matrices
  std::shared_ptr<Mat_t> bsr_mat;
  CSCMat_t csc_mat;

  // Cholesky solver
  SparseCholesky<T> *chol;
};

}  // namespace A2D

#endif  // A2D_ANALSYS_H