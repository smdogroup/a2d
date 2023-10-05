#ifndef A2D_SYM_MAT_EIGS_H
#define A2D_SYM_MAT_EIGS_H

#include "a2ddefs.h"

namespace A2D {

template <typename T>
KOKKOS_FUNCTION void SymEigs2x2(const T* S, T* eigs, T* Q) {
  T tr = S[0] + S[2];
  T diff = S[0] - S[2];
  T discrm = sqrt(diff * diff + 4.0 * S[1] * S[1]);
  T det = S[0] * S[2] - S[1] * S[1];

  // Compute the eigenvalues such that eigs[0] <= eigs[1]
  if (std::real(tr) > 0.0) {
    eigs[1] = 0.5 * (tr + discrm);
    eigs[0] = det / eigs[1];
  } else if (std::real(tr) < 0.0) {
    eigs[0] = 0.5 * (tr - discrm);
    eigs[1] = det / eigs[0];
  } else {
    eigs[0] = -0.5 * discrm;
    eigs[1] = 0.5 * discrm;
  }

  // Compute the eigenvector components
  T u = 1.0, v = 0.0;
  if (std::real(S[1]) != 0.0) {
    if (std::real(diff) > 0.0) {
      T a = 0.5 * (diff + discrm);
      T inv = 1.0 / sqrt(a * a + S[1] * S[1]);
      u = inv * S[1];
      v = -inv * a;
    } else {
      T a = 0.5 * (diff - discrm);
      T inv = 1.0 / sqrt(a * a + S[1] * S[1]);
      u = inv * a;
      v = inv * S[1];
    }
  }

  // Set the eigenvector components
  Q[0] = u;
  Q[1] = -v;
  Q[2] = v;
  Q[3] = u;
}

/**
 * Compute the eigenvalues and eigenvectors of a symmetric eigenvalue problem
 */
template <typename T, int N>
KOKKOS_FUNCTION void SymEigs(const SymMat<T, N>& S, Vec<T, N>& eigs,
                             Mat<T, N, N>& Q) {
  if constexpr (N == 2) {
    SymEigs2x2(get_data(S), get_data(eigs), get_data(Q));
  }
}

template <class Stype, class etype, class Qtype>
class SymEigsExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<etype>::type T;

  // Extract the dimensions of the underlying matrices
  static constexpr int N = get_symmatrix_size<Stype>::size;

  // Make sure we have the correct size
  static constexpr int K = get_vec_size<etype>::size;
  static constexpr int P = get_matrix_rows<Qtype>::size;
  static constexpr int R = get_matrix_columns<Qtype>::size;
  static_assert(K == N, "Vector of eigenvalues must be correct size");
  static_assert(P == N && R == N, "Eigenvector matrix must be correct size");

  SymEigsExpr(Stype& S, etype& eigs, Qtype& Q) : S(S), eigs(eigs), Q(Q) {}

  KOKKOS_FUNCTION void eval() {
    if constexpr (N == 2) {
      SymEigs2x2(get_data(S), get_data(eigs), get_data(Q));
    }
  }

  KOKKOS_FUNCTION void bzero() {
    Q.bzero();
    eigs.bzero();
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
  }

  KOKKOS_FUNCTION void reverse() {}

  KOKKOS_FUNCTION void hzero() {
    eigs.hzero();
    Q.hzero();
  }

  KOKKOS_FUNCTION void hreverse() {}

 private:
  Stype& S;
  etype& eigs;
  Qtype& Q;
};

template <class Stype, class etype, class Qtype>
KOKKOS_FUNCTION auto SymEigs(ADObj<Stype>& S, ADObj<etype>& eigs,
                             ADObj<Qtype>& Q) {
  return SymEigsExpr<ADObj<Stype>, ADObj<etype>, ADObj<Qtype>>(S, eigs, Q);
}

template <class Stype, class etype, class Qtype>
KOKKOS_FUNCTION auto SymEigs(A2DObj<Stype>& S, A2DObj<etype>& eigs,
                             A2DObj<Qtype>& Q) {
  return SymEigsExpr<A2DObj<Stype>, A2DObj<etype>, A2DObj<Qtype>>(S, eigs, Q);
}

namespace Test {

template <typename T, int N>
class SymEigsTest
    : public A2DTest<T, VarTuple<T, Vec<T, N>, Mat<T, N, N>>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, VarTuple<T, Vec<T, N>, Mat<T, N, N>>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymEigsTest<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    SymMat<T, N> S;
    Mat<T, N, N> Q;
    Vec<T, N> eigs;
    x.get_values(S);
    SymEigs(S, eigs, Q);
    return MakeVarTuple<T>(eigs, Q);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<SymMat<T, N>> S;
    ADObj<Mat<T, N, N>> Q;
    ADObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    auto stack = MakeStack(SymEigs(S, eigs, Q));
    auto sv = MakeVarTuple<T>(eigs.bvalue(), Q.bvalue());
    seed.get_values(sv);
    stack.reverse();
    g.set_values(S.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> S;
    A2DObj<Mat<T, N, N>> Q;
    A2DObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    p.get_values(S.pvalue());
    auto stack = MakeStack(SymEigs(S, eigs, Q));
    auto sv = MakeVarTuple<T>(eigs.bvalue(), Q.bvalue());
    auto hv = MakeVarTuple<T>(eigs.hvalue(), Q.hvalue());
    seed.get_values(sv);
    hval.get_values(hv);
    stack.hproduct();
    h.set_values(S.hvalue());
  }
};

bool SymEigsTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 2> test1;
    passed = passed && Run(test1, component, write_output);
  }

  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 3> test1;
    passed = passed && Run(test1, component, write_output);
  }

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYM_MAT_EIGS_H
