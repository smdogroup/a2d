#ifndef A2D_INTEGRAND_TEST_H
#define A2D_INTEGRAND_TEST_H

#include "ad/a2dtest.h"

namespace A2D {

namespace Test {

template <template <typename> class Integrand, typename T>
class A2DIntegrandAnalysisTest
    : public A2DTest<T, T, typename Integrand<T>::FiniteElementSpace> {
 public:
  // Set the solution space data
  using DataSpace = typename Integrand<T>::DataSpace;
  using FiniteElementSpace = typename Integrand<T>::FiniteElementSpace;
  using FiniteElementGeometry = typename Integrand<T>::FiniteElementGeometry;
  using SolutionMapping = typename Integrand<T>::SolutionMapping;
  using QMatType = typename Integrand<T>::QMatType;

  // Set the inputs and outputs
  using Inputs = FiniteElementSpace;
  using Output = T;

  A2DIntegrandTest(Integrand<T>& integrand)
      : A2DTest<T, T, typename Integrand<T>::FiniteElementSpace>(
            TestType::SECOND_ORDER),
        integrand(integrand) {
    // Set values for the data and the geometry
    this->set_rand(DataSpace::ncomp, data);
    this->set_rand(FiniteElementGeometry::ncomp, gref);
  }

  // Finite element geometry and data space
  DataSpace data;
  FiniteElementGeometry gref;

  // The integrand that we will test
  Integrand<T> integrand;

  /**
   * @brief Get the name of the test
   *
   * @returns The name of the test
   */
  std::string name() { return std::string("A2DIntegrandTest"); }

  /**
   * @brief Evaluate the output as a function of the inputs
   *
   * @param x Variable tuple for the input
   * @return The output value
   */
  VarTuple<T, Output> eval(const VarTuple<T, Inputs>& x) {
    // Get the solution/geometry in the reference domain
    FiniteElementSpace sref;

    // Get the input values
    x.get_values(sref);

    // Initialize the transform object
    T detJ;
    SolutionMapping transform(gref, detJ);

    // Transform from the reference element to the physical space
    FiniteElementSpace s;
    transform.transform(sref, s);

    // Compute the coefficients for the weak form of the Integrand
    FiniteElementSpace coef;
    T value = integrand.integrand(detJ, data, gref, s);

    return MakeVarTuple<T>(value);
  }

  /**
   * @brief Compute the derivative of the output as a function of the inputs
   *
   * @param seed The seed value for the input
   * @param x Variable tuple for the input
   * @param g Derivative of the output w.r.t. the input
   */
  void deriv(const VarTuple<T, Output>& seed, const VarTuple<T, Inputs>& x,
             VarTuple<T, Inputs>& g) {
    // Get the solution/geometry in the reference domain
    FiniteElementSpace sref;

    // Get the input values
    x.get_values(sref);

    // Initialize the transform object
    T detJ;
    SolutionMapping transform(gref, detJ);

    // Transform from the reference element to the physical space
    FiniteElementSpace s;
    transform.transform(sref, s);

    // Compute the coefficients for the weak form of the Integrand
    FiniteElementSpace coef;
    integrand.residual(detJ, data, gref, s, coef);

    // Transform the coefficients back to the reference element
    FiniteElementSpace cref;
    transform.rtransform(coef, cref);

    // Set the coefficients as the gradient
    for (int i = 0; i < FiniteElementSpace::ncomp; i++) {
      g[i] = cref[i];
    }
  }

  /**
   * @brief Compute the derivative of the output as a function of the inputs
   *
   * @param x Variable tuple for the input
   * @param g Derivative of the output w.r.t. the input
   */
  void hprod(const VarTuple<T, Output>& seed, const VarTuple<T, Output>& hval,
             const VarTuple<T, Inputs>& x, const VarTuple<T, Inputs>& p,
             VarTuple<T, Inputs>& h) {
    // Get the solution/geometry in the reference domain
    FiniteElementSpace sref;

    // Get the input values
    x.get_values(sref);

    // Initialize the transform object
    T detJ;
    SolutionMapping transform(gref, detJ);

    // Transform from the reference element to the physical space
    FiniteElementSpace s;
    transform.transform(sref, s);

    // Compute the Jacobian for the weak form at the quadrature point
    QMatType jac_ref, jac;
    integrand.jacobian(detJ, data, gref, s, jac);

    // Transform second derivatives from w.r.t. x to w.r.t. xi
    transform.template jtransform<FiniteElementSpace>(jac, jac_ref);

    for (int i = 0; i < FiniteElementSpace::ncomp; i++) {
      h[i] = T(0.0);
      for (int j = 0; j < FiniteElementSpace::ncomp; j++) {
        h[i] += jac_ref(i, j) * p[j];
      }
    }
  }
};

}  // namespace Test

}  // namespace A2D

#endif  // A2D_INTEGRAND_TEST_H