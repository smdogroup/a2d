#ifndef A2D_GREEN_STRAIN_H
#define A2D_GREEN_STRAIN_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "ad/core/a2dgreenstraincore.h"

namespace A2D {

enum class GreenStrain { LINEAR, NONLINEAR };

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrain(Mat<T, N, N>& Ux, SymMat<T, N>& E) {
  if constexpr (etype == GreenStrain::LINEAR) {
    LinearGreenStrainCore<T, N>(get_data(Ux), get_data(E));
  } else {
    NonlinearGreenStrainCore<T, N>(get_data(Ux), get_data(E));
  }
}

template <GreenStrain etype, typename T, int N, ADorder order>
class MatGreenStrainExpr {
 public:
  using Utype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, N>>;
  using Etype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, N>>;

  KOKKOS_FUNCTION MatGreenStrainExpr(Utype& Ux, Etype& E) : Ux(Ux), E(E) {
    if constexpr (etype == GreenStrain::LINEAR) {
      LinearGreenStrainCore<T, N>(get_data(Ux), get_data(E));
    } else {
      NonlinearGreenStrainCore<T, N>(get_data(Ux), get_data(E));
    }
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (etype == GreenStrain::LINEAR) {
      LinearGreenStrainForwardCore<T, N>(GetSeed<seed>::get_data(Ux),
                                         GetSeed<seed>::get_data(E));
    } else {
      NonlinearGreenStrainForwardCore<T, N>(get_data(Ux),
                                            GetSeed<seed>::get_data(Ux),
                                            GetSeed<seed>::get_data(E));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    if constexpr (etype == GreenStrain::LINEAR) {
      LinearGreenStrainReverseCore<T, N>(GetSeed<ADseed::b>::get_data(E),
                                         GetSeed<ADseed::b>::get_data(Ux));

    } else {
      NonlinearGreenStrainReverseCore<T, N>(get_data(Ux),
                                            GetSeed<ADseed::b>::get_data(E),
                                            GetSeed<ADseed::b>::get_data(Ux));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (etype == GreenStrain::LINEAR) {
      LinearGreenStrainHReverseCore<T, N>(GetSeed<ADseed::h>::get_data(E),
                                          GetSeed<ADseed::h>::get_data(Ux));

    } else {
      NonlinearGreenStrainHReverseCore<T, N>(
          get_data(Ux), GetSeed<ADseed::p>::get_data(Ux),
          GetSeed<ADseed::b>::get_data(E), GetSeed<ADseed::h>::get_data(E),
          GetSeed<ADseed::h>::get_data(Ux));
    }
  }

  Utype& Ux;
  Etype& E;
};

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION auto MatGreenStrain(ADMat<Mat<T, N, N>>& Ux,
                                    ADMat<SymMat<T, N>>& E) {
  return MatGreenStrainExpr<etype, T, N, ADorder::FIRST>(Ux, E);
}

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION auto MatGreenStrain(A2DMat<Mat<T, N, N>>& Ux,
                                    A2DMat<SymMat<T, N>>& E) {
  return MatGreenStrainExpr<etype, T, N, ADorder::SECOND>(Ux, E);
}

namespace Test {

template <GreenStrain etype, typename T, int N>
class MatGreenStrainTest : public A2DTest<T, SymMat<T, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatGreenStrain<";
    if (etype == GreenStrain::LINEAR) {
      s << "LINEAR," << N << ">";
    } else {
      s << "NONLINEAR," << N << ">";
    }

    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> Ux;
    SymMat<T, N> E;
    x.get_values(Ux);
    MatGreenStrain<etype>(Ux, E);
    return MakeVarTuple<T>(E);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    Mat<T, N, N> Ux0, Uxb;
    SymMat<T, N> E0, Eb;
    ADMat<Mat<T, N, N>> Ux(Ux0, Uxb);
    ADMat<SymMat<T, N>> E(E0, Eb);

    x.get_values(Ux0);
    auto op = MatGreenStrain<etype>(Ux, E);
    auto stack = MakeStack(op);
    seed.get_values(Eb);
    stack.reverse();
    g.set_values(Uxb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DMat<Mat<T, N, N>> Ux;
    A2DMat<SymMat<T, N>> E;

    x.get_values(Ux.value());
    p.get_values(Ux.pvalue());

    auto op = MatGreenStrain<etype>(Ux, E);
    auto stack = MakeStack(op);

    seed.get_values(E.bvalue());
    hval.get_values(E.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(Ux.hvalue());
  }
};

bool MatGreenStrainTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatGreenStrainTest<GreenStrain::LINEAR, Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatGreenStrainTest<GreenStrain::NONLINEAR, Tc, 2> test2;
  passed = passed && Run(test2, component, write_output);

  MatGreenStrainTest<GreenStrain::LINEAR, Tc, 3> test3;
  passed = passed && Run(test3, component, write_output);
  MatGreenStrainTest<GreenStrain::NONLINEAR, Tc, 3> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_GREEN_STRAIN_H