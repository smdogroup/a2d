#ifndef A2D_FE_BASE_H
#define A2D_FE_BASE_H

#include "multiphysics/fesolution.h"

namespace A2D {

/**
 * @brief Element base class.
 *
 * This defines an element that is compatible with the given PDE. You can
 * set the node locations into the element, add the non-zero-pattern of the
 * Jacobian matrix via the add_node_set as well as adding the residual and
 * Jacobian contributions
 *
 * @tparam T data type
 */
template <typename T>
class ElementBase {
 public:
  virtual ~ElementBase() {}
  virtual void set_geo(SolutionVector<T>& X) = 0;
  virtual void set_solution(SolutionVector<T>& U) = 0;
  virtual void add_residual(SolutionVector<T>& res) = 0;
  // virtual void add_jacobian_vector_product(SolutionVector<T>& x,
  //                                          SolutionVector<T>& y) = 0;

  // virtual void add_dof_set(Kokkos::UnorderedMap<COO<I>, void>& node_set) = 0;
  // virtual void add_jacobian(typename PDEInfo::SparseMat& jac) = 0;
  virtual void add_adjoint_dfdnodes(SolutionVector<T>& psi,
                                    SolutionVector<T>& dfdx) {}
};

/**
 * @brief Base class for a constitutive object.
 *
 * These objects provide a way to set the design variables into the elements.
 * The mapping from design variables to element data is not defined and depends
 * on the nature of the design variable.
 *
 * @tparam T data type
 */
template <typename T>
class ConstitutiveBase {
 public:
  virtual void set_design_vars(SolutionVector<T>& x) {}
  virtual void add_adjoint_dfdx(SolutionVector<T>& psi,
                                SolutionVector<T>& dfdx) {}
};

/*
  ElementFunction base class.

  This class enables the computation of a function that sums up over all
  the elements within the element.
*/
template <typename T>
class ElementFunctional {
 public:
  virtual ~ElementFunctional() {}
  virtual T eval_functional() = 0;
  virtual void add_dfdu(SolutionVector<T>& dfdu) {}
  virtual void add_dfdx(SolutionVector<T>& dfdx) {}
  virtual void add_dfdnodes(SolutionVector<T>& dfdx) {}
};

/*
  Container class for all functionals.

  The functional must be the result of a sum over the ElementFunctions in the
  container.
*/
// template <typename I, typename T, class PDE>
// class Functional {
//  public:
//   Functional() {}
//   virtual ~Functional() {}

//   void add_functional(
//       std::shared_ptr<ElementFunctional<I, T, PDE>> functional) {
//     functionals.push_back(functional);
//   }

//   /*
//     Evaluate the functional by summing the values over all components
//   */
//   T eval_functional() {
//     T value = 0.0;
//     for (auto it = functionals.begin(); it != functionals.end(); it++) {
//       value += (*it)->eval_functional();
//     }
//     return value;
//   }

//   /*
//     Compute the derivative of the functional w.r.t. state variables
//   */
//   void eval_dfdu(std::shared_ptr<typename PDE::SolutionArray> dfdu) {
//     A2D::BLAS::zero(*dfdu);
//     for (auto it = functionals.begin(); it != functionals.end(); it++) {
//       (*it)->add_dfdu(*dfdu);
//     }
//   }

//   /*
//     Compute the derivative of the functional w.r.t. design variables
//   */
//   void eval_dfdx(std::shared_ptr<typename PDE::DesignArray> dfdx) {
//     A2D::BLAS::zero(*dfdx);
//     for (auto it = functionals.begin(); it != functionals.end(); it++) {
//       (*it)->add_dfdx(*dfdx);
//     }
//   }

//   /*
//     Compute the derivative of the functional w.r.t. nodes
//   */
//   void eval_dfdnodes(std::shared_ptr<typename PDE::NodeArray> dfdx) {
//     A2D::BLAS::zero(*dfdx);
//     for (auto it = functionals.begin(); it != functionals.end(); it++) {
//       (*it)->add_dfdnodes(*dfdx);
//     }
//   }

//  private:
//   std::list<std::shared_ptr<ElementFunctional<I, T, PDE>>> functionals;
// };

}  // namespace A2D

#endif  //  A2D_FE_BASE_H