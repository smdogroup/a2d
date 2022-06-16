#ifndef A2D_FUNCTIONAL_H
#define A2D_FUNCTIONAL_H

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
  virtual void add_dfdu(typename PDE::SolutionArray& dfdu) = 0;
  virtual void add_dfdx(typename PDE::NodeArray& dfdx) = 0;
};

/*
  Container class for all functionals.

  The functional must be the result of a sum over the ElementFunctions in the
  container.
*/
template <typename I, typename T, class PDE>
class Functional {
 public:
  Functional() { num_functionals = 0; }
  virtual ~Functional() {}

  void add_functional(ElementFunctional<I, T, PDE>* functional) {
    functionals[num_functionals] = functional;
    num_functionals++;
  }

  index_t get_num_functionals() { return num_functionals; }
  ElementFunctional<I, T, PDE>* get_functional(index_t index) {
    return functionals[index];
  }

  /*
    Evaluate the functional by summing the values over all components
  */
  T eval_functional(Functional& functional) {
    T value = 0.0;
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      value += functional.eval_functional();
    }
  }

  /*
    Compute the derivative of the functional w.r.t. state variables
  */
  void eval_dfdu(Functional& functional, typename PDE::SolutionArray& dfdu) {
    dfdu.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdu(dfdu);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. design variables
  */
  void eval_dfdx(Functional& functional, typename PDE::DesignArray& dfdx) {
    dfdx.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdx(dfdx);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. nodes
  */
  void eval_dfdnodes(Functional& functional, typename PDE::NodeArray& dfdx) {
    dfdx.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdx(dfdx);
    }
  }

 private:
  index_t num_functionals;
  ElementFunctional<I, T, PDE>* functionals[10];
};

// /*
//   A minimal set of data needed for the ElementFunctional class
// */
// template <typename I, typename T, class PDE, class Basis>
// class FunctionalBasis : public ElementFunctional<I, T, PDE> {
//  public:
//   FunctionalBasis(ElementBasis<I, T, PDE, Basis>& element) : element(element)
//   {}

//   ElementBasis<I, T, PDE, Basis>& element;
// };

}  // namespace A2D

#endif  // A2D_FUNCTIONAL_H