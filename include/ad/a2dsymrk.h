#ifndef A2D_SYMMAT_RK_H
#define A2D_SYMMAT_RK_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dtest.h"
#include "ad/core/a2dgemmcore.h"

namespace A2D {

/*
  Compute S = A * A^{T}  or S = A^{T} * A
*/
template <typename T, int N, int K, MatOp op = MatOp::NORMAL,
          bool additive = false>
KOKKOS_FUNCTION void SymMatRKCore(const T A[], T S[]) {
  if constexpr (op == MatOp::NORMAL) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        const T* a = &A[K * i];
        const T* b = &A[K * j];

        T val = 0.0;
        for (int k = 0; k < K; k++) {
          val += a[0] * b[0];
          a++, b++;
        }
        if constexpr (additive) {
          S[0] += val;
        } else {
          S[0] = val;
        }
        S++;
      }
    }
  } else {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j <= i; j++) {
        const T* a = &A[i];
        const T* b = &A[j];

        T val = 0.0;
        for (int k = 0; k < N; k++) {
          val += a[0] * b[0];
          a += K, b += K;
        }
        if constexpr (additive) {
          S[0] += val;
        } else {
          S[0] = val;
        }
        S++;
      }
    }
  }
}

/*
  Compute S = alpha * A * A^{T}  or S = alpha * A^{T} * A
*/
template <typename T, int N, int K, MatOp op = MatOp::NORMAL,
          bool additive = false>
KOKKOS_FUNCTION void SymMatRKCoreScale(const T alpha, const T A[], T S[]) {
  if constexpr (op == MatOp::NORMAL) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        const T* a = &A[K * i];
        const T* b = &A[K * j];

        T val = 0.0;
        for (int k = 0; k < K; k++) {
          val += a[0] * b[0];
          a++, b++;
        }
        if constexpr (additive) {
          S[0] += alpha * val;
        } else {
          S[0] = alpha * val;
        }
        S++;
      }
    }
  } else {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j <= i; j++) {
        const T* a = &A[i];
        const T* b = &A[j];

        T val = 0.0;
        for (int k = 0; k < N; k++) {
          val += a[0] * b[0];
          a += K, b += K;
        }
        if constexpr (additive) {
          S[0] += alpha * val;
        } else {
          S[0] = alpha * val;
        }
        S++;
      }
    }
  }
}

/*
  Compute S = A * B^{T} + B * A^{T}  or  S = A^{T} * B + B^{T} * A^{T}
*/
template <typename T, int N, int K, MatOp op = MatOp::NORMAL,
          bool additive = false>
KOKKOS_FUNCTION void SymMatR2KCore(const T A[], const T B[], T S[]) {
  if constexpr (op == MatOp::NORMAL) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        T val = 0.0;

        const T* a = &A[K * i];
        const T* b = &B[K * j];
        for (int k = 0; k < K; k++) {
          val += a[0] * b[0];
          a++, b++;
        }

        a = &A[K * j];
        b = &B[K * i];
        for (int k = 0; k < K; k++) {
          val += a[0] * b[0];
          a++, b++;
        }

        if constexpr (additive) {
          S[0] += val;
        } else {
          S[0] = val;
        }
        S++;
      }
    }
  } else {
    for (int i = 0; i < K; i++) {
      for (int j = 0; j <= i; j++) {
        T val = 0.0;

        const T* a = &A[i];
        const T* b = &B[j];
        for (int k = 0; k < N; k++) {
          val += a[0] * b[0];
          a += K, b += K;
        }

        a = &A[j];
        b = &B[i];
        for (int k = 0; k < N; k++) {
          val += a[0] * b[0];
          a += K, b += K;
        }

        if constexpr (additive) {
          S[0] += val;
        } else {
          S[0] = val;
        }
        S++;
      }
    }
  }
}

template <typename T, int N, int K, MatOp op = MatOp::NORMAL>
KOKKOS_FUNCTION void SymMatRKCoreReverse(const T A[], const T Sb[], T Ab[]) {
  if constexpr (op == MatOp::NORMAL) {
    // Ab = Sb * A
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < K; j++) {
        int k = 0;
        const T* s = &Sb[i * (i + 1) / 2];
        const T* a = &A[j];

        T val = 0.0;
        for (; k < i; k++) {
          val += s[0] * a[0];
          a += K, s++;
        }

        for (; k < N; k++) {
          val += s[0] * a[0];
          a += K, s += k + 1;
        }

        val += A[K * i + j] * Sb[i + i * (i + 1) / 2];

        Ab[0] += val;
        Ab++;
      }
    }
  } else {  // op == MatOp::TRANSPOSE
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < K; j++) {
        int k = 0;
        const T* a = &A[K * i];
        const T* s = &Sb[j * (j + 1) / 2];

        T val = 0.0;
        for (; k < j; k++) {
          val += s[0] * a[0];
          a++, s++;
        }

        for (; k < K; k++) {
          val += s[0] * a[0];
          a++, s += k + 1;
        }

        val += A[K * i + j] * Sb[j + j * (j + 1) / 2];

        Ab[0] += val;
        Ab++;
      }
    }
  }
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION void SymMatRK(const Mat<T, N, K>& A, SymMat<T, P>& S) {
  static_assert(
      (op == MatOp::NORMAL && P == N) || (op == MatOp::TRANSPOSE && K == P),
      "SymMatRK matrix dimensions must agree");
  SymMatRKCore<T, N, K, op>(get_data(A), get_data(S));
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION void SymMatRK(const T alpha, const Mat<T, N, K>& A,
                              SymMat<T, P>& S) {
  static_assert(
      (op == MatOp::NORMAL && P == N) || (op == MatOp::TRANSPOSE && K == P),
      "SymMatRK matrix dimensions must agree");
  SymMatRKCoreScale<T, N, K, op>(get_data(alpha), get_data(A), get_data(S));
}

template <typename T, int N, int K, int P, ADorder order, MatOp op>
class SymMatRKExpr {
 private:
  using Atype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, K>>;
  using Stype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, P>>;

 public:
  KOKKOS_FUNCTION SymMatRKExpr(Atype& A, Stype& S) : A(A), S(S) {
    static_assert(
        (op == MatOp::NORMAL && P == N) || (op == MatOp::TRANSPOSE && K == P),
        "SymMatRK matrix dimensions must agree");
    SymMatRKCore<T, N, K, op>(get_data(A), get_data(S));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    SymMatR2KCore<T, N, K, op>(get_data(A), GetSeed<seed>::get_data(A),
                               GetSeed<seed>::get_data(S));
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    SymMatRKCoreReverse<T, N, K, op>(get_data(A), GetSeed<seed>::get_data(S),
                                     GetSeed<seed>::get_data(A));
  }

  KOKKOS_FUNCTION void hreverse() {
    SymMatRKCoreReverse<T, N, K, op>(get_data(A),
                                     GetSeed<ADseed::h>::get_data(S),
                                     GetSeed<ADseed::h>::get_data(A));

    SymMatRKCoreReverse<T, N, K, op>(GetSeed<ADseed::p>::get_data(A),
                                     GetSeed<ADseed::b>::get_data(S),
                                     GetSeed<ADseed::h>::get_data(A));
  }

  Atype& A;
  Stype& S;
};

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(ADMat<Mat<T, N, K>>& A, ADMat<SymMat<T, P>>& S) {
  return SymMatRKExpr<T, N, K, P, ADorder::FIRST, op>(A, S);
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(A2DMat<Mat<T, N, K>>& A,
                              A2DMat<SymMat<T, P>>& S) {
  return SymMatRKExpr<T, N, K, P, ADorder::SECOND, op>(A, S);
}

template <typename T, int N, int K, int P, ADorder order, ADiffType ada,
          ADiffType adA, MatOp op>
class SymMatRKScaleExpr {
 private:
  using atype = ADScalarInputType<ada, order, T>;
  using Atype = ADMatType<adA, order, Mat<T, N, K>>;
  using Stype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, P>>;

 public:
  KOKKOS_FUNCTION SymMatRKScaleExpr(atype& alpha, Atype& A, Stype& S)
      : alpha(alpha), A(A), S(S) {
    static_assert(
        (op == MatOp::NORMAL && P == N) || (op == MatOp::TRANSPOSE && K == P),
        "SymMatRK matrix dimensions must agree");
    SymMatRKCore<T, N, K, op>(get_data(A), get_data(S));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    SymMatR2KCore<T, N, K, op>(get_data(A), GetSeed<seed>::get_data(A),
                               GetSeed<seed>::get_data(S));
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    SymMatRKCoreReverse<T, N, K, op>(get_data(A), GetSeed<seed>::get_data(S),
                                     GetSeed<seed>::get_data(A));
  }

  KOKKOS_FUNCTION void hreverse() {
    SymMatRKCoreReverse<T, N, K, op>(get_data(A),
                                     GetSeed<ADseed::h>::get_data(S),
                                     GetSeed<ADseed::h>::get_data(A));

    SymMatRKCoreReverse<T, N, K, op>(GetSeed<ADseed::p>::get_data(A),
                                     GetSeed<ADseed::b>::get_data(S),
                                     GetSeed<ADseed::h>::get_data(A));
  }
};

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(ADScalar<T>& alpha, ADMat<Mat<T, N, K>>& A,
                              ADMat<SymMat<T, P>>& S) {
  return SymMatRKScaleExpr<T, N, K, P, ADorder::FIRST, op>(alpha, A, S);
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(const T alpha, ADMat<Mat<T, N, K>>& A,
                              ADMat<SymMat<T, P>>& S) {
  return SymMatRKScaleExpr<T, N, K, P, ADorder::FIRST, op>(alpha, A, S);
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(ADScalar<T>& alpha, const Mat<T, N, K>& A,
                              ADMat<SymMat<T, P>>& S) {
  return SymMatRKScaleExpr<T, N, K, P, ADorder::FIRST, op>(alpha, A, S);
}

template <MatOp op, typename T, int N, int K, int P>
KOKKOS_FUNCTION auto SymMatRK(A2DMat<Mat<T, N, K>>& A,
                              A2DMat<SymMat<T, P>>& S) {
  return SymMatRKExpr<T, N, K, P, ADorder::SECOND, op>(A, S);
}

namespace Test {

template <MatOp op, typename T, int N, int M, int P>
class SymMatRKTest : public A2DTest<T, SymMat<T, P>, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>>;
  using Output = VarTuple<T, SymMat<T, P>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatRKTest<";
    if (op == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    s << N << "," << M << "," << P << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, M> A;
    SymMat<T, P> S;
    x.get_values(A);
    SymMatRK<op>(A, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    Mat<T, N, M> A0, Ab;
    SymMat<T, P> S0, Sb;
    ADMat<Mat<T, N, M>> A(A0, Ab);
    ADMat<SymMat<T, P>> S(S0, Sb);

    x.get_values(A0);
    auto mult = SymMatRK<op>(A, S);
    auto stack = MakeStack(mult);
    seed.get_values(Sb);
    stack.reverse();
    g.set_values(Ab);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DMat<Mat<T, N, M>> A;
    A2DMat<SymMat<T, P>> S;

    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto mult = SymMatRK<op>(A, S);
    auto stack = MakeStack(mult);

    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N, int K>
bool SymMatRKTestHelper(bool component = false, bool write_output = true) {
  const MatOp NORMAL = MatOp::NORMAL;
  const MatOp TRANSPOSE = MatOp::TRANSPOSE;
  using Tc = std::complex<T>;

  bool passed = true;
  SymMatRKTest<NORMAL, Tc, N, K, N> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatRKTest<TRANSPOSE, Tc, N, K, K> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

bool SymMatRKTestAll(bool component = false, bool write_output = true) {
  bool passed = true;
  passed = passed && SymMatRKTestHelper<double, 3, 3>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 2, 4>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 5, 2>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 7, 5>(component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYMMAT_RK_H