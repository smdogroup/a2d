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

template <class Impl>
class Analysis {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;

  virtual ~Analysis() {}

  virtual std::shared_ptr<Vec_t> get_data() = 0;
  virtual std::shared_ptr<Vec_t> get_geo() = 0;
  virtual std::shared_ptr<Vec_t> get_sol() = 0;
  virtual std::shared_ptr<Vec_t> get_res() = 0;

  virtual void linear_solve() = 0;
  virtual void nonlinear_solve() = 0;

  virtual T evaluate(FunctionalBase<Impl> &func) = 0;
  virtual void add_derivative(FunctionalBase<Impl> &func, FEVarType wrt,
                              T alpha, Vec_t &dfdx) = 0;
  virtual void eval_adjoint_derivative(FunctionalBase<Impl> &func,
                                       FEVarType wrt, Vec_t &dfdx) = 0;
  virtual void eval_adjoint_derivative(FEVarType wrt, Vec_t &dfdx) = 0;
  virtual void to_vtk(const std::string prefix) {}
};

template <typename T, index_t block_size>
class DirectCholeskyAnalysisImpl {
 public:
  using type = T;
  using Vec_t = SolutionVector<T>;
  using Mat_t = BSRMat<T, block_size, block_size>;

  template <class Basis>
  using ElementVector = ElementVector_Serial<type, Basis, Vec_t>;

  template <class Basis>
  using ElementMatrix = ElementMat_Serial<type, Basis, Mat_t>;
};

template <typename T, index_t block_size>
class DirectCholeskyAnalysis
    : public Analysis<DirectCholeskyAnalysisImpl<T, block_size>> {
 public:
  // Set the implementation type
  using Impl_t = DirectCholeskyAnalysisImpl<T, block_size>;

  // Data types
  using Vec_t = typename Impl_t::Vec_t;
  using Mat_t = typename Impl_t::Mat_t;
  using CSCMat_t = CSCMat<T>;

  DirectCholeskyAnalysis(std::shared_ptr<Vec_t> data,
                         std::shared_ptr<Vec_t> geo, std::shared_ptr<Vec_t> sol,
                         std::shared_ptr<Vec_t> res,
                         std::shared_ptr<ElementAssembler<Impl_t>> assembler,
                         std::shared_ptr<DirichletBCs<T>> bcs = nullptr)
      : data(data),
        geo(geo),
        sol(sol),
        res(res),
        assembler(assembler),
        bcs(bcs) {
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
    if (bcs != nullptr) {
      bcs->set_bcs(*sol);
    }
  }
  ~DirectCholeskyAnalysis() {}

  /**
   * @brief Get the data vector object
   *
   * @return std::shared_ptr<Vec_t>
   */
  std::shared_ptr<Vec_t> get_data() { return data; }

  /**
   * @brief Get the geometry vector object
   *
   * @return std::shared_ptr<Vec_t>
   */
  std::shared_ptr<Vec_t> get_geo() { return geo; }

  /**
   * @brief Get the solution vector object
   *
   * @return std::shared_ptr<Vec_t>
   */
  std::shared_ptr<Vec_t> get_sol() { return sol; }

  /**
   * @brief Get the residual vector - stores residuals/dfdu
   *
   * @return std::shared_ptr<Vec_t>
   */
  std::shared_ptr<Vec_t> get_res() { return res; }

  /**
   * @brief Get the underlying matrix
   *
   */
  std::shared_ptr<Mat_t> get_mat() { return bsr_mat; }

  /**
   * @brief Compute the residual and store it
   */
  void residual() {
    res->zero();
    assembler->add_residual(T(1.0), *data, *geo, *sol, *res);

    // Apply boundary conditions
    if (bcs != nullptr) {
      bcs->zero_bcs(*res);
    }
  }

  /**
   * @brief Compute the Jacobian and store it
   */
  void jacobian() {
    // Compute the BSR matrix
    bsr_mat->zero();
    assembler->add_jacobian(T(1.0), *data, *geo, *sol, *bsr_mat);

    // Zero the rows corresponding to the boundary conditions
    if (bcs != nullptr) {
      const index_t *bc_dofs;
      index_t nbcs = bcs->get_bcs(&bc_dofs);
      bsr_mat->zero_rows(nbcs, bc_dofs);
    }
  }

  /**
   * @brief Factor the Jacobian
   */
  void factor() {
    // This is a problem here -- matrix is passed by value, but should be
    // passed by reference to avoid copying. Convert bsr to csc because we
    // want to use Cholesky factorization
    csc_mat = bsr_to_csc(*bsr_mat);

    // Apply boundary conditions to each column
    if (bcs != nullptr) {
      const index_t *bc_dofs;
      index_t nbcs = bcs->get_bcs(&bc_dofs);
      csc_mat.zero_columns(nbcs, bc_dofs);
    }

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
    chol->solve(res->data());

    // Apply the boundary condition values
    if (bcs != nullptr) {
      bcs->set_bcs(*sol);
    }

    // The residual
    sol->axpy(-1.0, *res);
  }

  void nonlinear_solve() {}

  T evaluate(FunctionalBase<Impl_t> &func) {
    return func.evaluate(*data, *geo, *sol);
  }

  void add_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt, T alpha,
                      Vec_t &dfdx) {
    func.add_derivative(wrt, alpha, *data, *geo, *sol, dfdx);
  }

  void eval_adjoint_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt,
                               Vec_t &dfdx) {
    // This doesn't make sense for the adjoint method
    if (wrt == FEVarType::STATE) {
      return;
    }

    // Add the contributions to the derivative from the partial
    dfdx.zero();
    func.add_derivative(wrt, T(1.0), *data, *geo, *sol, dfdx);

    // Compute the derivative of the function wrt state and solve the adjoint
    // equations
    res->zero();
    func.add_derivative(FEVarType::STATE, T(1.0), *data, *geo, *sol, *res);

    // Apply boundary conditions
    if (bcs != nullptr) {
      bcs->zero_bcs(*res);
    }

    // Solve the adjoint system
    chol->solve(res->data());

    // Add the terms from the total derivative
    assembler->add_adjoint_res_product(wrt, T(-1.0), *data, *geo, *sol, *res,
                                       dfdx);
  }

  void eval_adjoint_derivative(FEVarType wrt, Vec_t &dfdx) {
    // This doesn't make sense for the adjoint method
    if (wrt == FEVarType::STATE) {
      return;
    }

    // Apply boundary conditions
    if (bcs != nullptr) {
      bcs->zero_bcs(*res);
    }

    // Solve the adjoint system
    chol->solve(res->data());

    // Add the terms from the total derivative
    assembler->add_adjoint_res_product(wrt, T(-1.0), *data, *geo, *sol, *res,
                                       dfdx);
  }

  void to_vtk(const std::string prefix) {
    assembler->to_vtk(*data, *geo, *sol, prefix);
  }

 private:
  // Vectors of data, geometry and solution degrees of freedom
  std::shared_ptr<Vec_t> data;
  std::shared_ptr<Vec_t> geo;
  std::shared_ptr<Vec_t> sol;
  std::shared_ptr<Vec_t> res;

  // Element matrix assembler object
  std::shared_ptr<ElementAssembler<Impl_t>> assembler;

  // Boundary conditions
  std::shared_ptr<DirichletBCs<T>> bcs;

  // System matrices
  std::shared_ptr<Mat_t> bsr_mat;
  CSCMat_t csc_mat;

  // Cholesky solver
  SparseCholesky<T> *chol;
};

// template <typename T, index_t block_size>
// class MatFreeAmgAnalysisImpl {
//  public:
//   using type = T;
//   using Vec_t = SolutionVector<T>;
//   using Mat_t = BSRMat<T, block_size, block_size>;

//   template <class Basis>
//   using ElementVector = ElementVector_Serial<type, Basis, Vec_t>;

//   template <class Basis>
//   using ElementMatrix = ElementMat_Serial<type, Basis, Mat_t>;
// };

// template <typename T, index_t block_size, index_t null_size>
// class MatFreeAmgAnalysis
//     : public Analysis<MatFreeAmgAnalysisImpl<T, block_size>> {
//  public:
//   // Set the implementation type
//   using Impl_t = MatFreeAmgAnalysisImpl<T, block_size>;

//   // Data types
//   using Vec_t = typename Impl_t::Vec_t;
//   using Mat_t = typename Impl_t::Mat_t;
//   using Amg_t = BSRMatAmg<T, block_size, null_size>;

//   MatFreeAmgAnalysis(std::shared_ptr<Vec_t> data, std::shared_ptr<Vec_t> geo,
//                      std::shared_ptr<Vec_t> sol, std::shared_ptr<Vec_t> res,
//                      std::shared_ptr<ElementAssembler<Impl_t>> assembler,
//                      std::shared_ptr<ElementAssembler<Impl_t>> low_assembler,
//                      std::shared_ptr<DirichletBCs<T>> bcs = nullptr,
//                      int amg_nlevels = 2, int cg_it = 100,
//                      double cg_rtol = 1e-8, double cg_atol = 1e-30)
//       : data(data),
//         geo(geo),
//         sol(sol),
//         res(res),
//         assembler(assembler),
//         low_assembler(low_assembler),
//         bcs(bcs),
//         amg_nlevels(amg_nlevels),
//         cg_it(cg_it),
//         cg_rtol(cg_rtol),
//         cg_atol(cg_atol) {
//     // Create the system matrix
//     index_t nrows;
//     std::vector<index_t> rowp, cols;
//     low_assembler->get_bsr_data(block_size, nrows, rowp, cols);
//     bsr_mat = std::make_shared<Mat_t>(nrows, nrows, cols.size(), rowp, cols);

//     // Create the sparse matrix
//     csc_mat = bsr_to_csc(*bsr_mat);

//     // Set up Cholesky solver but don't set up values yet
//     bool set_values = false;
//     chol =
//         new SparseCholesky(csc_mat, CholOrderingType::ND, nullptr,
//         set_values);

//     // Set the boundary conditions
//     if (bcs != nullptr) {
//       bcs->set_bcs(*sol);
//     }
//   }
//   ~MatFreeAmgAnalysis() {}

//   std::shared_ptr<Vec_t> get_data() { return data; }
//   std::shared_ptr<Vec_t> get_geo() { return geo; }
//   std::shared_ptr<Vec_t> get_sol() { return sol; }
//   std::shared_ptr<Vec_t> get_res() { return res; }
//   std::shared_ptr<Mat_t> get_mat() { return bsr_mat; }

//   /**
//    * @brief Compute the residual and store it
//    */
//   void residual() {
//     res->zero();
//     assembler->add_residual(T(1.0), *data, *geo, *sol, *res);

//     // Apply boundary conditions
//     if (bcs != nullptr) {
//       bcs->zero_bcs(*res);
//     }
//   }

//   /**
//    * @brief Compute the Jacobian and store it
//    */
//   void jacobian() {
//     // Compute the BSR matrix
//     bsr_mat->zero();
//     low_assembler->add_jacobian(T(1.0), *data, *geo, *sol, *bsr_mat);

//     // Zero the rows corresponding to the boundary conditions
//     if (bcs != nullptr) {
//       const index_t *bc_dofs;
//       index_t nbcs = bcs->get_bcs(&bc_dofs);
//       bsr_mat->zero_rows(nbcs, bc_dofs);
//     }
//   }

//   /**
//    * @brief Factor the Jacobian
//    */
//   void factor() {
//     // This is a problem here -- matrix is passed by value, but should be
//     // passed by reference to avoid copying. Convert bsr to csc because we
//     // want to use Cholesky factorization
//     csc_mat = bsr_to_csc(*bsr_mat);

//     // Apply boundary conditions to each column
//     if (bcs != nullptr) {
//       const index_t *bc_dofs;
//       index_t nbcs = bcs->get_bcs(&bc_dofs);
//       csc_mat.zero_columns(nbcs, bc_dofs);
//     }

//     // Set values to Cholesky solver and factorize
//     chol->setValues(csc_mat);
//     chol->factor();
//   }

//   /**
//    * @brief Perform a linear solve to obtain the solution
//    */
//   void linear_solve() {
//     residual();
//     jacobian();
//     factor();

//     // Find the solution
//     chol->solve(res->data());

//     // Apply the boundary condition values
//     if (bcs != nullptr) {
//       bcs->set_bcs(*sol);
//     }

//     // The residual
//     sol->axpy(-1.0, *res);
//   }

//   void nonlinear_solve() {}

//   T evaluate(FunctionalBase<Impl_t> &func) {
//     return func.evaluate(*data, *geo, *sol);
//   }

//   void add_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt, T alpha,
//                       Vec_t &dfdx) {
//     func.add_derivative(wrt, alpha, *data, *geo, *sol, dfdx);
//   }

//   void eval_adjoint_derivative(FunctionalBase<Impl_t> &func, FEVarType wrt,
//                                Vec_t &dfdx) {
//     // This doesn't make sense for the adjoint method
//     if (wrt == FEVarType::STATE) {
//       return;
//     }

//     // Add the contributions to the derivative from the partial
//     dfdx.zero();
//     func.add_derivative(wrt, T(1.0), *data, *geo, *sol, dfdx);

//     // Compute the derivative of the function wrt state and solve the adjoint
//     // equations
//     res->zero();
//     func.add_derivative(FEVarType::STATE, T(1.0), *data, *geo, *sol, *res);

//     // Apply boundary conditions
//     if (bcs != nullptr) {
//       bcs->zero_bcs(*res);
//     }

//     // Solve the adjoint system
//     chol->solve(res->data());

//     // Add the terms from the total derivative
//     assembler->add_adjoint_res_product(wrt, T(-1.0), *data, *geo, *sol, *res,
//                                        dfdx);
//   }

//   void eval_adjoint_derivative(FEVarType wrt, Vec_t &dfdx) {
//     // This doesn't make sense for the adjoint method
//     if (wrt == FEVarType::STATE) {
//       return;
//     }

//     // Apply boundary conditions
//     if (bcs != nullptr) {
//       bcs->zero_bcs(*res);
//     }

//     // Solve the adjoint system
//     chol->solve(res->data());

//     // Add the terms from the total derivative
//     assembler->add_adjoint_res_product(wrt, T(-1.0), *data, *geo, *sol, *res,
//                                        dfdx);
//   }

//   void to_vtk(const std::string prefix) {
//     assembler->to_vtk(*data, *geo, *sol, prefix);
//   }

//  private:
//   // Vectors of data, geometry and solution degrees of freedom
//   std::shared_ptr<Vec_t> data;
//   std::shared_ptr<Vec_t> geo;
//   std::shared_ptr<Vec_t> sol;
//   std::shared_ptr<Vec_t> res;

//   // Element matrix assembler object
//   std::shared_ptr<ElementAssembler<Impl_t>> assembler;
//   std::shared_ptr<ElementAssembler<Impl_t>> low_assembler;

//   // Boundary conditions
//   std::shared_ptr<DirichletBCs<T>> bcs;

//   // System matrices
//   std::shared_ptr<Mat_t> bsr_mat;

//   // Preconditioner
//   std::shared_ptr<Amg_t> amg;
// };

}  // namespace A2D

#endif  // A2D_ANALSYS_H