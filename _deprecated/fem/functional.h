#ifndef A2D_FUNCTIONAL_H
#define A2D_FUNCTIONAL_H

#include <list>
#include <memory>

#include "element.h"

namespace A2D {

/*
  ElementFunction base class.

  This class enables the computation of a function that sums up over all
  the elements within the element.
*/
template <typename I, typename T, class PDE>
class ElementFunctional {
 public:
  virtual ~ElementFunctional() {}
  virtual T eval_functional() = 0;
  virtual void add_dfdu(typename PDE::SolutionArray& dfdu) {}
  virtual void add_dfdx(typename PDE::DesignArray& dfdx) {}
  virtual void add_dfdnodes(typename PDE::NodeArray& dfdx) {}
};

/*
  Container class for all functionals.

  The functional must be the result of a sum over the ElementFunctions in the
  container.
*/
template <typename I, typename T, class PDE>
class Functional {
 public:
  Functional() {}
  virtual ~Functional() {}

  void add_functional(
      std::shared_ptr<ElementFunctional<I, T, PDE>> functional) {
    functionals.push_back(functional);
  }

  /*
    Evaluate the functional by summing the values over all components
  */
  T eval_functional() {
    T value = 0.0;
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      value += (*it)->eval_functional();
    }
    return value;
  }

  /*
    Compute the derivative of the functional w.r.t. state variables
  */
  void eval_dfdu(std::shared_ptr<typename PDE::SolutionArray> dfdu) {
    A2D::BLAS::zero(*dfdu);
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      (*it)->add_dfdu(*dfdu);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. design variables
  */
  void eval_dfdx(std::shared_ptr<typename PDE::DesignArray> dfdx) {
    A2D::BLAS::zero(*dfdx);
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      (*it)->add_dfdx(*dfdx);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. nodes
  */
  void eval_dfdnodes(std::shared_ptr<typename PDE::NodeArray> dfdx) {
    A2D::BLAS::zero(*dfdx);
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      (*it)->add_dfdnodes(*dfdx);
    }
  }

 private:
  std::list<std::shared_ptr<ElementFunctional<I, T, PDE>>> functionals;
};

}  // namespace A2D

#endif  // A2D_FUNCTIONAL_H