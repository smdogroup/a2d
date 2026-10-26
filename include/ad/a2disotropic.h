#ifndef A2D_ISOTROPIC_H
#define A2D_ISOTROPIC_H

#include "../a2ddefs.h"
#include "../a2dmat.h"

namespace A2D {

template <typename T, int N>
A2D_FUNCTION void SymIsotropicCore(const T mu, const T lambda, const T E[],
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
A2D_FUNCTION void SymIsotropicAddCore(const T mu, const T lambda, const T E[],
                                      T S[]) {
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
A2D_FUNCTION void SymIsotropicReverseCoefCore(const T E[], const T Sb[], T& mu,
                                              T& lambda) {
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
A2D_FUNCTION void SymIsotropic(const T mu, const T lambda,
                               const SymMat<T, N>& E, SymMat<T, N>& S) {
  SymIsotropicCore<T, N>(mu, lambda, get_data(E), get_data(S));
}

template <class mutype, class lamtype, class Etype, class Stype>
class SymIsotropicExpr {
 public:
  // Extract the numeric type to use
  using T = object_scalar_t<Stype>;

  // Extract the dimensions of the matrices
  static constexpr int N = get_symmatrix_size<Etype>::size;
  static constexpr int M = get_symmatrix_size<Stype>::size;

  // Assert that the matrices are the right size
  static_assert(N == M, "Matrices must be the same size");

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<Stype>::order;

  // Get the types of the different objects
  static constexpr ADiffType mudiff = get_diff_type<mutype>::diff_type;
  static constexpr ADiffType lamdiff = get_diff_type<lamtype>::diff_type;
  static constexpr ADiffType Ediff = get_diff_type<Etype>::diff_type;

  A2D_FUNCTION SymIsotropicExpr(mutype mu, lamtype lambda, Etype& E, Stype& S)
      : mu(mu), lambda(lambda), E(E), S(S) {}

  A2D_FUNCTION void eval() {
    SymIsotropicCore<T, N>(get_data(mu), get_data(lambda), get_data(E),
                           get_data(S));
  }

  A2D_FUNCTION void bzero() { S.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (mudiff == ADiffType::ACTIVE && lamdiff == ADiffType::ACTIVE &&
                  Ediff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(GetSeed<seed>::get_data(mu),
                             GetSeed<seed>::get_data(lambda), get_data(E),
                             GetSeed<seed>::get_data(S));
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<seed>::get_data(E),
                                GetSeed<seed>::get_data(S));
    } else if constexpr (mudiff == ADiffType::ACTIVE &&
                         lamdiff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(GetSeed<seed>::get_data(mu),
                             GetSeed<seed>::get_data(lambda), get_data(E),
                             GetSeed<seed>::get_data(S));
    } else if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicCore<T, N>(get_data(mu), get_data(lambda),
                             GetSeed<seed>::get_data(E),
                             GetSeed<seed>::get_data(S));
    }
  }

  A2D_FUNCTION void reverse() {
    if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<ADseed::b>::get_data(S),
                                GetSeed<ADseed::b>::get_data(E));
    }
    if constexpr (mudiff == ADiffType::ACTIVE || lamdiff == ADiffType::ACTIVE) {
      T mu0(0.0), lam0(0.0);
      SymIsotropicReverseCoefCore<T, N>(
          get_data(E), GetSeed<ADseed::b>::get_data(S), mu0, lam0);
      if constexpr (mudiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::b>::get_data(mu) += mu0;
      }
      if constexpr (lamdiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::b>::get_data(lambda) += lam0;
      }
    }
  }

  A2D_FUNCTION void hzero() { S.hzero(); }

  A2D_FUNCTION void hreverse() {
    if constexpr (Ediff == ADiffType::ACTIVE) {
      SymIsotropicAddCore<T, N>(get_data(mu), get_data(lambda),
                                GetSeed<ADseed::h>::get_data(S),
                                GetSeed<ADseed::h>::get_data(E));
    }
    if constexpr (mudiff == ADiffType::ACTIVE || lamdiff == ADiffType::ACTIVE) {
      T mu0(0.0), lam0(0.0);
      SymIsotropicReverseCoefCore<T, N>(
          get_data(E), GetSeed<ADseed::h>::get_data(S), mu0, lam0);
      if constexpr (mudiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::h>::get_data(mu) += mu0;
      }
      if constexpr (lamdiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::h>::get_data(lambda) += lam0;
      }
    }
    if constexpr (Ediff == ADiffType::ACTIVE &&
                  (mudiff == ADiffType::ACTIVE ||
                   lamdiff == ADiffType::ACTIVE)) {
      T pmu(0.0), plam(0.0), hmu(0.0), hlam(0.0);
      if constexpr (mudiff == ADiffType::ACTIVE) {
        pmu = GetSeed<ADseed::p>::get_data(mu);
      }
      if constexpr (lamdiff == ADiffType::ACTIVE) {
        plam = GetSeed<ADseed::p>::get_data(lambda);
      }
      SymIsotropicAddCore<T, N>(pmu, plam, GetSeed<ADseed::b>::get_data(S),
                                GetSeed<ADseed::h>::get_data(E));
      SymIsotropicReverseCoefCore<T, N>(GetSeed<ADseed::p>::get_data(E),
                                        GetSeed<ADseed::b>::get_data(S), hmu,
                                        hlam);
      if constexpr (mudiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::h>::get_data(mu) += hmu;
      }
      if constexpr (lamdiff == ADiffType::ACTIVE) {
        GetSeed<ADseed::h>::get_data(lambda) += hlam;
      }
    }
  }

  mutype mu;
  mutype lambda;
  Etype& E;
  Stype& S;
};

template <class mutype, class lamtype, class Etype, class Stype>
A2D_FUNCTION auto SymIsotropic(ADObj<mutype>& mu, ADObj<lamtype>& lambda,
                               ADObj<Etype>& E, ADObj<Stype>& S) {
  return SymIsotropicExpr<ADObj<mutype>&, ADObj<lamtype>&, ADObj<Etype>,
                          ADObj<Stype>>(mu, lambda, E, S);
}

template <class mutype, class lamtype, class Etype, class Stype>
A2D_FUNCTION auto SymIsotropic(mutype mu, lamtype lambda, ADObj<Etype>& E,
                               ADObj<Stype>& S) {
  return SymIsotropicExpr<mutype, lamtype, ADObj<Etype>, ADObj<Stype>>(
      mu, lambda, E, S);
}

template <class mutype, class lamtype, class Etype, class Stype>
A2D_FUNCTION auto SymIsotropic(A2DObj<mutype>& mu, A2DObj<lamtype>& lambda,
                               A2DObj<Etype>& E, A2DObj<Stype>& S) {
  return SymIsotropicExpr<A2DObj<mutype>&, A2DObj<lamtype>&, A2DObj<Etype>,
                          A2DObj<Stype>>(mu, lambda, E, S);
}

template <class mutype, class lamtype, class Etype, class Stype>
A2D_FUNCTION auto SymIsotropic(mutype mu, lamtype lambda, A2DObj<Etype>& E,
                               A2DObj<Stype>& S) {
  return SymIsotropicExpr<mutype, lamtype, A2DObj<Etype>, A2DObj<Stype>>(
      mu, lambda, E, S);
}

}  // namespace A2D

#endif  // A2D_ISOTROPIC_H
