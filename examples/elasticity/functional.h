#ifndef A2D_FUNCTIONAL_H
#define A2D_FUNCTIONAL_H

#include <list>

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

  void add_functional(ElementFunctional<I, T, PDE>* functional) {
    functionals.push_back(functional);
  }

  /*
    Evaluate the functional by summing the values over all components
  */
  T eval_functional() {
    T value = 0.0;
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      ElementFunctional<I, T, PDE>* functional = *it;
      value += functional->eval_functional();
    }
    return value;
  }

  /*
    Compute the derivative of the functional w.r.t. state variables
  */
  void eval_dfdu(typename PDE::SolutionArray& dfdu) {
    dfdu.zero();
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      ElementFunctional<I, T, PDE>* functional = *it;
      functional->add_dfdu(dfdu);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. design variables
  */
  void eval_dfdx(typename PDE::DesignArray& dfdx) {
    dfdx.zero();
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      ElementFunctional<I, T, PDE>* functional = *it;
      functional->add_dfdx(dfdx);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. nodes
  */
  void eval_dfdnodes(typename PDE::NodeArray& dfdx) {
    dfdx.zero();
    for (auto it = functionals.begin(); it != functionals.end(); it++) {
      ElementFunctional<I, T, PDE>* functional = *it;
      functional->eval_dfdnodes(dfdx);
    }
  }

 private:
  std::list<ElementFunctional<I, T, PDE>*> functionals;
};

}  // namespace A2D

#endif  // A2D_FUNCTIONAL_H