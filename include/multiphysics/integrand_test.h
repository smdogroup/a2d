#ifndef A2D_INTEGRAND_TEST_H
#define A2D_INTEGRAND_TEST_H

#include "ad/a2dtest.h"

namespace A2D {

namespace Test {

template <template <typename> class Integrand, typename T>
class A2DIntegrandAnalysisTest
    : public A2DTest<T, T, typename Integrand<T>::DataSpace,
                     typename Integrand<T>::FiniteElementGeometry,
                     typename Integrand<T>::FiniteElementSpace> {
 public:
  // Set the solution space data
  using DataSpace = typename Integrand<T>::DataSpace;
  using Geometry = typename Integrand<T>::FiniteElementGeometry;
  using Space = typename Integrand<T>::FiniteElementSpace;

  static constexpr FEVarType DATA = FEVarType::DATA;
  static constexpr FEVarType GEOMETRY = FEVarType::GEOMETRY;
  static constexpr FEVarType STATE = FEVarType::STATE;

  // Set the inputs and outputs
  using Input = VarTuple<T, DataSpace, Geometry, Space>;
  using Output = VarTuple<T, T>;

  A2DIntegrandAnalysisTest(Integrand<T>& integrand)
      : A2DTest<T, T, DataSpace, Geometry, Space>(TestType::SECOND_ORDER),
        integrand(integrand) {}

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
  Output eval(const Input& x) {
    DataSpace data;
    Geometry geo;
    Space sref;
    x.get_values(data, geo, sref);

    // Compute the coefficients for the weak form of the Integrand
    double weight(0.1378);
    T value = integrand.integrand(weight, data, geo, sref);

    return MakeVarTuple<T>(value);
  }

  /**
   * @brief Compute the derivative of the output as a function of the inputs
   *
   * @param seed The seed value for the input
   * @param x Variable tuple for the input
   * @param g Derivative of the output w.r.t. the input
   */
  void deriv(const Output& seed, const Input& x, Input& g) {
    DataSpace data, bdata;
    Geometry geo, bgeo;
    Space sref, bsref;
    x.get_values(data, geo, sref);

    double weight(0.1378);
    integrand.template residual<DATA>(weight, data, geo, sref, bdata);
    integrand.template residual<GEOMETRY>(weight, data, geo, sref, bgeo);
    integrand.template residual<STATE>(weight, data, geo, sref, bsref);

    g.set_values(bdata, bgeo, bsref);
  }

  /**
   * @brief Compute the derivative of the output as a function of the inputs
   *
   * @param x Variable tuple for the input
   * @param g Derivative of the output w.r.t. the input
   */
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    DataSpace data, pdata, bdata, hdata;
    Geometry geo, pgeo, bgeo, hgeo;
    Space sref, psref, bsref, hsref;
    x.get_values(data, geo, sref);
    p.get_values(pdata, pgeo, psref);

    double weight(0.1378);
    if constexpr (DataSpace::ncomp > 0) {
      integrand.template jacobian_product<DATA, DATA>(weight, data, geo, sref,
                                                      pdata, hdata);
      integrand.template jacobian_product<DATA, GEOMETRY>(weight, data, geo,
                                                          sref, pgeo, bdata);
      for (index_t i = 0; i < DataSpace::ncomp; i++) {
        hdata[i] += bdata[i];
      }
      integrand.template jacobian_product<DATA, STATE>(weight, data, geo, sref,
                                                       psref, bdata);
      for (index_t i = 0; i < DataSpace::ncomp; i++) {
        hdata[i] += bdata[i];
      }
    }

    if constexpr (Geometry::ncomp > 0) {
      integrand.template jacobian_product<GEOMETRY, DATA>(weight, data, geo,
                                                          sref, pdata, hgeo);
      integrand.template jacobian_product<GEOMETRY, GEOMETRY>(weight, data, geo,
                                                              sref, pgeo, bgeo);
      for (index_t i = 0; i < Geometry::ncomp; i++) {
        hgeo[i] += bgeo[i];
      }
      integrand.template jacobian_product<GEOMETRY, STATE>(weight, data, geo,
                                                           sref, psref, bgeo);
      for (index_t i = 0; i < Geometry::ncomp; i++) {
        hgeo[i] += bgeo[i];
      }
    }

    if constexpr (Space::ncomp > 0) {
      integrand.template jacobian_product<STATE, DATA>(weight, data, geo, sref,
                                                       pdata, hsref);
      integrand.template jacobian_product<STATE, GEOMETRY>(weight, data, geo,
                                                           sref, pgeo, bsref);
      for (index_t i = 0; i < Space::ncomp; i++) {
        hsref[i] += bsref[i];
      }
      integrand.template jacobian_product<STATE, STATE>(weight, data, geo, sref,
                                                        psref, bsref);
      for (index_t i = 0; i < Space::ncomp; i++) {
        hsref[i] += bsref[i];
      }
    }

    h.set_values(hdata, hgeo, hsref);
  }
};

}  // namespace Test

}  // namespace A2D

#endif  // A2D_INTEGRAND_TEST_H