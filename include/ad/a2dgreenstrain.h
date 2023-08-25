#ifndef A2D_GREEN_STRAIN_H
#define A2D_GREEN_STRAIN_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "ad/core/a2dgemmcore.h"

namespace A2D {

enum class GreenStrain { LINEAR, NONLINEAR };

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrainCore(const T Ux[], T E[]) {
  static_assert(N == 2 || N == 3,
                "MatGreenStrainCore must use N == 2 or N == 3");

  if constexpr (etype == GreenStrain::LINEAR) {
    // E = 0.5 * (Ux + Ux^{T})
    if constexpr (N == 2) {
      E[0] = Ux[0];

      E[1] = 0.5 * (Ux[1] + Ux[2]);
      E[2] = Ux[3];
    } else {
      E[0] = Ux[0];

      E[1] = 0.5 * (Ux[1] + Ux[3]);
      E[2] = Ux[4];

      E[3] = 0.5 * (Ux[2] + Ux[6]);
      E[4] = 0.5 * (Ux[5] + Ux[7]);
      E[5] = Ux[8];
    }
  } else {
    // E = 0.5 * (Ux + Ux^{T} + Ux^{T} * Ux)
    if constexpr (N == 2) {
      E[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);

      E[1] = 0.5 * (Ux[1] + Ux[2] + Ux[0] * Ux[1] + Ux[2] * Ux[3]);
      E[2] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
    } else {
      E[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);

      E[1] =
          0.5 * (Ux[1] + Ux[3] + Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
      E[2] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);

      E[3] =
          0.5 * (Ux[2] + Ux[6] + Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
      E[4] =
          0.5 * (Ux[5] + Ux[7] + Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
      E[5] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);
    }
  }
}

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrainForwardCore(const T Ux[], const T Ud[],
                                               T E[]) {
  static_assert(N == 2 || N == 3,
                "MatGreenStrainForwardCore must use N == 2 or N == 3");

  if constexpr (etype == GreenStrain::LINEAR) {
    // E = 0.5 * (Ux + Ux^{T})
    if constexpr (N == 2) {
      E[0] = Ud[0];

      E[1] = 0.5 * (Ud[1] + Ud[2]);
      E[2] = Ud[3];
    } else {
      E[0] = Ud[0];

      E[1] = 0.5 * (Ud[1] + Ud[3]);
      E[2] = Ud[4];

      E[3] = 0.5 * (Ud[2] + Ud[6]);
      E[4] = 0.5 * (Ud[5] + Ud[7]);
      E[5] = Ud[8];
    }
  } else {
    // E = 0.5 * (Ux + Ux^{T} + Ux^{T} * Ux)
    if constexpr (N == 2) {
      E[0] = Ud[0] + Ux[0] * Ud[0] + Ux[2] * Ud[2];

      E[1] = 0.5 * (Ud[1] + Ud[2] + Ux[0] * Ud[1] + Ux[2] * Ud[3] +
                    Ud[0] * Ux[1] + Ud[2] * Ux[3]);
      E[2] = Ud[2] + Ux[1] * Ud[1] + Ux[3] * Ud[3];

    } else {
      E[0] = Ud[0] + Ux[0] * Ud[0] + Ux[3] * Ud[3] + Ux[6] * Ud[6];

      E[1] =
          0.5 * (Ud[1] + Ud[3] + Ux[0] * Ud[1] + Ux[3] * Ud[4] + Ux[6] * Ud[7] +
                 Ud[0] * Ux[1] + Ud[3] * Ux[4] + Ud[6] * Ux[7]);
      E[2] = Ud[4] + Ux[1] * Ud[1] + Ux[4] * Ud[4] + Ux[7] * Ud[7];

      E[3] =
          0.5 * (Ud[2] + Ud[6] + Ux[0] * Ud[2] + Ux[3] * Ud[5] + Ux[6] * Ud[8] +
                 Ud[0] * Ux[2] + Ud[3] * Ux[5] + Ud[6] * Ux[8]);
      E[4] =
          0.5 * (Ud[5] + Ud[7] + Ux[1] * Ud[2] + Ux[4] * Ud[5] + Ux[7] * Ud[8] +
                 Ud[1] * Ux[2] + Ud[4] * Ux[5] + Ud[7] * Ux[8]);
      E[5] = Ud[8] + Ux[2] * Ud[2] + Ux[5] * Ud[5] + Ux[8] * Ud[8];
    }
  }
}

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrainReverseCore(const T Ux[], const T Eb[],
                                               T Ub[]) {
  static_assert(N == 2 || N == 3,
                "MatGreenStrainReverseCore must use N == 2 or N == 3");

  if constexpr (etype == GreenStrain::LINEAR) {
    // E = 0.5 * (Ux + Ux^{T})
    if constexpr (N == 2) {
      Ub[0] += Eb[0];
      Ub[1] += 0.5 * Eb[1];

      Ub[2] += 0.5 * Eb[1];
      Ub[3] += Eb[2];
    } else {
      // Uxb = Eb
      Ub[0] += Eb[0];
      Ub[1] += 0.5 * Eb[1];
      Ub[2] += 0.5 * Eb[3];

      Ub[3] += 0.5 * Eb[1];
      Ub[4] += Eb[2];
      Ub[5] += 0.5 * Eb[4];

      Ub[6] += 0.5 * Eb[3];
      Ub[7] += 0.5 * Eb[4];
      Ub[8] += Eb[5];
    }
  } else {
    // Uxb = (I + Ux) * Eb
    if constexpr (N == 2) {
      // Uxb = (I + Ux) * Eb
      Ub[0] += (Ux[0] + 1.0) * Eb[0] + 0.5 * Ux[1] * Eb[1];
      Ub[1] += 0.5 * (Ux[0] + 1.0) * Eb[1] + Ux[1] * Eb[2];

      Ub[2] += Ux[2] * Eb[0] + 0.5 * (Ux[3] + 1.0) * Eb[1];
      Ub[3] += 0.5 * Ux[2] * Eb[1] + (Ux[3] + 1.0) * Eb[2];
    } else {
      Ub[0] +=
          (Ux[0] + 1.0) * Eb[0] + 0.5 * Ux[1] * Eb[1] + 0.5 * Ux[2] * Eb[3];
      Ub[1] +=
          0.5 * (Ux[0] + 1.0) * Eb[1] + Ux[1] * Eb[2] + 0.5 * Ux[2] * Eb[4];
      Ub[2] +=
          0.5 * (Ux[0] + 1.0) * Eb[3] + 0.5 * Ux[1] * Eb[4] + Ux[2] * Eb[5];

      Ub[3] +=
          Ux[3] * Eb[0] + 0.5 * (Ux[4] + 1.0) * Eb[1] + 0.5 * Ux[5] * Eb[3];
      Ub[4] +=
          0.5 * Ux[3] * Eb[1] + (Ux[4] + 1.0) * Eb[2] + 0.5 * Ux[5] * Eb[4];
      Ub[5] +=
          0.5 * Ux[3] * Eb[3] + 0.5 * (Ux[4] + 1.0) * Eb[4] + Ux[5] * Eb[5];

      Ub[6] +=
          Ux[6] * Eb[0] + 0.5 * Ux[7] * Eb[1] + 0.5 * (Ux[8] + 1.0) * Eb[3];
      Ub[7] +=
          0.5 * Ux[6] * Eb[1] + Ux[7] * Eb[2] + 0.5 * (Ux[8] + 1.0) * Eb[4];
      Ub[8] +=
          0.5 * Ux[6] * Eb[3] + 0.5 * Ux[7] * Eb[4] + (Ux[8] + 1.0) * Eb[5];
    }
  }
}

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrainHReverseCore(const T Ux[], const T Up[],
                                                const T Eb[], const T Eh[],
                                                T Uh[]) {
  static_assert(N == 2 || N == 3,
                "MatGreenStrainReverseCore must use N == 2 or N == 3");

  if constexpr (etype == GreenStrain::LINEAR) {
    if constexpr (N == 2) {
      Uh[0] += Eh[0];
      Uh[1] += 0.5 * Eh[1];

      Uh[2] += 0.5 * Eh[1];
      Uh[3] += Eh[2];
    } else {
      Uh[0] += Eh[0];
      Uh[1] += 0.5 * Eh[1];
      Uh[2] += 0.5 * Eh[3];

      Uh[3] += 0.5 * Eh[1];
      Uh[4] += Eh[2];
      Uh[5] += 0.5 * Eh[4];

      Uh[6] += 0.5 * Eh[3];
      Uh[7] += 0.5 * Eh[4];
      Uh[8] += Eh[5];
    }
  } else {
    // if constexpr (N == 2){
    // Uxh(0, 0) += Uxp(0, 0) * Eb(0, 0) + 0.5 * Uxp(0, 1) * Eb(0, 1);
    // Uxh(0, 1) += 0.5 * Uxp(0, 0) * Eb(0, 1) + Uxp(0, 1) * Eb(1, 1);

    // Uxh(1, 0) += Uxp(1, 0) * Eb(0, 0) + 0.5 * Uxp(1, 1) * Eb(0, 1);
    // Uxh(1, 1) += 0.5 * Uxp(1, 0) * Eb(0, 1) + Uxp(1, 1) * Eb(1, 1);

    // Uxh(0, 0) += (Ux(0, 0) + 1.0) * Eh(0, 0) + 0.5 * Ux(0, 1) * Eh(0, 1);
    // Uxh(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eh(0, 1) + Ux(0, 1) * Eh(1, 1);

    // Uxh(1, 0) += Ux(1, 0) * Eh(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eh(0, 1);
    // Uxh(1, 1) += 0.5 * Ux(1, 0) * Eh(0, 1) + (Ux(1, 1) + 1.0) * Eh(1, 1);
  }
}

template <GreenStrain etype, typename T, int N>
KOKKOS_FUNCTION void MatGreenStrain(Mat<T, N, N>& Ux, SymMat<T, N>& E) {
  MatGreenStrainCore<etype, T, N>(get_data(Ux), get_data(E));
}

template <GreenStrain etype, typename T, int N, ADorder order>
class MatGreenStrainExpr {
 public:
  using Utype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, N>>;
  using Etype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, N>>;

  KOKKOS_FUNCTION MatGreenStrainExpr(Utype& Ux, Etype& E) : Ux(Ux), E(E) {
    MatGreenStrainCore<etype, T, N>(get_data(Ux), get_data(E));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    MatGreenStrainForwardCore<etype, T, N>(
        get_data(Ux), GetSeed<seed>::get_data(Ux), GetSeed<seed>::get_data(E));
  }

  KOKKOS_FUNCTION void reverse() {
    MatGreenStrainReverseCore<etype, T, N>(get_data(Ux),
                                           GetSeed<ADseed::b>::get_data(E),
                                           GetSeed<ADseed::b>::get_data(Ux));
  }

  KOKKOS_FUNCTION void hreverse() {
    MatGreenStrainHReverseCore<etype, T, N>(
        get_data(Ux), GetSeed<ADseed::p>::get_data(Ux),
        GetSeed<ADseed::b>::get_data(E), GetSeed<ADseed::h>::get_data(E),
        GetSeed<ADseed::h>::get_data(Ux));
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

void MatGreenStrainTestAll() {
  using Tc = std::complex<double>;

  MatGreenStrainTest<GreenStrain::LINEAR, Tc, 2> test1;
  Run(test1);
  MatGreenStrainTest<GreenStrain::NONLINEAR, Tc, 2> test2;
  Run(test2);

  MatGreenStrainTest<GreenStrain::LINEAR, Tc, 3> test3;
  Run(test3);
  MatGreenStrainTest<GreenStrain::NONLINEAR, Tc, 3> test4;
  Run(test4);
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_GREEN_STRAIN_H