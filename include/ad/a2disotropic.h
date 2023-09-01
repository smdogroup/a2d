#ifndef A2D_ISOTROPIC_H
#define A2D_ISOTROPIC_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dscalar.h"
#include "a2dstack.h"
#include "a2dtest.h"
// #include "ad/core/a2disotropiccore.h"

namespace A2D {

template <typename T, int N>
KOKKOS_FUNCTION void SymIsotropicCore(const T mu, const T lambda, const T E[],
                                      T S[]) {
  static_assert(N == 2 || N == 3, "SymIsotropicCore must use N == 2 or N == 3");

  if constexpr (N == 2) {
    T tr = lambda * (E[0] + E[2]);
    T mu2 = 2.0 * mu;
    S[0] = mu2 * E[0] + tr;
    S[1] = mu2 * E[1];
    S[2] = mu2 * E[2] + tr;
  } else {
    T tr = lambda * (E[0] + E[2] + E[5]);
    T mu2 = 2.0 * mu;
    S[0] = mu2 * E[0] + tr;
    S[1] = mu2 * E[1];
    S[2] = mu2 * E[2] + tr;
    S[3] = mu2 * E[3];
    S[4] = mu2 * E[4];
    S[5] = mu2 * E[5] + tr;
  }
}

template <typename T, int N>
KOKKOS_FUNCTION void SymIsotropicAddCore(const T mu, const T lambda,
                                         const T E[], T S[]) {
  static_assert(N == 2 || N == 3, "SymIsotropicCore must use N == 2 or N == 3");

  if constexpr (N == 2) {
    T tr = lambda * (E[0] + E[2]);
    T mu2 = 2.0 * mu;
    S[0] += mu2 * E[0] + tr;
    S[1] += mu2 * E[1];
    S[2] += mu2 * E[2] + tr;
  } else {
    T tr = lambda * (E[0] + E[2] + E[5]);
    T mu2 = 2.0 * mu;
    S[0] += mu2 * E[0] + tr;
    S[1] += mu2 * E[1];
    S[2] += mu2 * E[2] + tr;
    S[3] += mu2 * E[3];
    S[4] += mu2 * E[4];
    S[5] += mu2 * E[5] + tr;
  }
}

template <typename T, int N>
KOKKOS_FUNCTION void SymIsotropicReverseCoefCore(const T E[], const T Sb[],
                                                 T& mu, T& lambda) {
  static_assert(N == 2 || N == 3, "SymIsotropicCore must use N == 2 or N == 3");

  if constexpr (N == 2) {
    lambda += (Sb[0] + Sb[2]) * (E[0] + E[2]);
    mu += 2.0 * (Sb[0] * E[0] + Sb[1] * E[1] + Sb[2] * E[2]);
  } else {
    lambda += (Sb[0] + Sb[2] + Sb[5]) * (E[0] + E[2] + E[5]);
    mu += 2.0 * (Sb[0] * E[0] + Sb[1] * E[1] + Sb[2] * E[2] + Sb[3] * E[3] +
                 Sb[4] * E[4] + Sb[5] * E[5]);
  }
}

template <typename T, int N>
KOKKOS_FUNCTION void SymIsotropic(const T mu, const T lambda,
                                  const SymMat<T, N>& E, SymMat<T, N>& S) {
  SymIsotropicCore<T, N>(mu, lambda, get_data(E), get_data(S));
}

template <typename T, int N, ADorder order,
          ADiffType mudiff = ADiffType::ACTIVE,
          ADiffType Ediff = ADiffType::ACTIVE>
class SymIsotropicExpr {
 public:
  using mutype = ADScalarInputType<mudiff, order, T>;
  using Etype = ADMatType<Ediff, order, SymMat<T, N>>;
  using Stype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, N>>;

  KOKKOS_FUNCTION SymIsotropicExpr(mutype mu, mutype lambda, Etype& E, Stype& S)
      : mu(mu), lambda(lambda), E(E), S(S) {}

  KOKKOS_FUNCTION void eval() {
    SymIsotropicCore<T, N>(get_data(mu), get_data(lambda), get_data(E),
                           get_data(S));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (mudiff == ADiffType::ACTIVE && Ediff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(GetSeed<seed>::get_data(mu),
                             GetSeed<seed>::get_data(lambda), get_data(E),
                             GetSeed<seed>::get_data(S));
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<seed>::get_data(E),
                                GetSeed<seed>::get_data(S));
    } else if constexpr (mudiff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(GetSeed<seed>::get_data(mu),
                             GetSeed<seed>::get_data(lambda), get_data(E),
                             GetSeed<seed>::get_data(S));
    } else if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(get_data(mu), get_data(lambda),
                             GetSeed<seed>::get_data(E),
                             GetSeed<seed>::get_data(S));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<ADseed::b>::get_data(S),
                                GetSeed<ADseed::b>::get_data(E));
    }
    if constexpr (mudiff == ADiffType::ACTIVE) {
      SymIsotropicReverseCoefCore<T, N>(get_data(E),
                                        GetSeed<ADseed::b>::get_data(S),
                                        GetSeed<ADseed::b>::get_data(mu),
                                        GetSeed<ADseed::b>::get_data(lambda));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<ADseed::h>::get_data(S),
                                GetSeed<ADseed::h>::get_data(E));
    }
    if constexpr (mudiff == ADiffType::ACTIVE) {
      SymIsotropicReverseCoefCore<T, N>(get_data(E),
                                        GetSeed<ADseed::h>::get_data(S),
                                        GetSeed<ADseed::h>::get_data(mu),
                                        GetSeed<ADseed::h>::get_data(lambda));
    }
    if constexpr (Ediff == ADiffType::ACTIVE && mudiff == ADiffType::ACTIVE) {
      SymIsotropicAddCore<T, N>(GetSeed<ADseed::p>::get_data(mu),
                                GetSeed<ADseed::p>::get_data(lambda),
                                GetSeed<ADseed::b>::get_data(S),
                                GetSeed<ADseed::h>::get_data(E));
      SymIsotropicReverseCoefCore<T, N>(GetSeed<ADseed::p>::get_data(E),
                                        GetSeed<ADseed::b>::get_data(S),
                                        GetSeed<ADseed::h>::get_data(mu),
                                        GetSeed<ADseed::h>::get_data(lambda));
    }
  }

  mutype mu;
  mutype lambda;
  Etype& E;
  Stype& S;
};

template <typename T, int N>
KOKKOS_FUNCTION auto SymIsotropic(ADScalar<T>& mu, ADScalar<T>& lambda,
                                  ADMat<SymMat<T, N>>& E,
                                  ADMat<SymMat<T, N>>& S) {
  return SymIsotropicExpr<T, N, ADorder::FIRST>(mu, lambda, E, S);
}

template <typename T, int N>
KOKKOS_FUNCTION auto SymIsotropic(const T mu, const T lambda,
                                  ADMat<SymMat<T, N>>& E,
                                  ADMat<SymMat<T, N>>& S) {
  return SymIsotropicExpr<T, N, ADorder::FIRST, ADiffType::PASSIVE>(mu, lambda,
                                                                    E, S);
}

template <typename T, int N>
KOKKOS_FUNCTION auto SymIsotropic(A2DScalar<T>& mu, A2DScalar<T>& lambda,
                                  A2DMat<SymMat<T, N>>& E,
                                  A2DMat<SymMat<T, N>>& S) {
  return SymIsotropicExpr<T, N, ADorder::SECOND>(mu, lambda, E, S);
}

template <typename T, int N>
KOKKOS_FUNCTION auto SymIsotropic(const T mu, const T lambda,
                                  A2DMat<SymMat<T, N>>& E,
                                  A2DMat<SymMat<T, N>>& S) {
  return SymIsotropicExpr<T, N, ADorder::SECOND, ADiffType::PASSIVE>(mu, lambda,
                                                                     E, S);
}

namespace Test {

template <typename T, int N>
class SymIsotropicConstTest : public A2DTest<T, SymMat<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  std::string name() {
    std::stringstream s;
    s << "SymIsotropic<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    SymMat<T, N> E, S;
    x.get_values(E);
    SymIsotropic(T(0.314), T(0.731), E, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    SymMat<T, N> E0, Eb, S0, Sb;
    ADMat<SymMat<T, N>> E(E0, Eb), S(S0, Sb);

    x.get_values(E.value());
    auto op = SymIsotropic(T(0.314), T(0.731), E, S);
    auto stack = MakeStack(op);
    seed.get_values(Sb);
    stack.reverse();
    g.set_values(E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DMat<SymMat<T, N>> E, S;

    x.get_values(E.value());
    p.get_values(E.pvalue());
    auto op = SymIsotropic(T(0.314), T(0.731), E, S);
    auto stack = MakeStack(op);
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(E.hvalue());
  }
};

template <typename T, int N>
class SymIsotropicTest : public A2DTest<T, SymMat<T, N>, T, T, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, T, T, SymMat<T, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  std::string name() {
    std::stringstream s;
    s << "SymIsotropic<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T mu, lambda;
    SymMat<T, N> E, S;
    x.get_values(mu, lambda, E);
    SymIsotropic(mu, lambda, E, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADScalar<T> mu, lambda;
    SymMat<T, N> E0, Eb, S0, Sb;
    ADMat<SymMat<T, N>> E(E0, Eb), S(S0, Sb);

    x.get_values(mu.value, lambda.value, E.value());
    auto op = SymIsotropic(mu, lambda, E, S);
    auto stack = MakeStack(op);
    seed.get_values(Sb);
    stack.reverse();
    g.set_values(mu.bvalue, lambda.bvalue, E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DScalar<T> mu, lambda;
    A2DMat<SymMat<T, N>> E, S;

    x.get_values(mu.value, lambda.value, E.value());
    p.get_values(mu.pvalue, lambda.pvalue, E.pvalue());

    auto op = SymIsotropic(mu, lambda, E, S);
    auto stack = MakeStack(op);
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(mu.hvalue, lambda.hvalue, E.hvalue());
  }
};

bool SymIsotropicTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  SymIsotropicConstTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymIsotropicConstTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  SymIsotropicTest<Tc, 2> test3;
  passed = passed && Run(test3, component, write_output);
  SymIsotropicTest<Tc, 3> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_ISOTROPIC_H