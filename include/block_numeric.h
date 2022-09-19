#ifndef A2D_BLOCK_NUMERIC_H
#define A2D_BLOCK_NUMERIC_H

#include <complex>

#include "a2dobjs.h"
#include "multiarray.h"

namespace A2D {

double absfunc(std::complex<double> a) {
  if (a.real() >= 0.0) {
    return a.real();
  } else {
    return -a.real();
  }
}

double absfunc(double a) {
  if (a >= 0.0) {
    return a;
  } else {
    return -a;
  }
}

double RealPart(double a) { return a; }

double RealPart(std::complex<double> a) { return a.real(); }

/*
  Compute y = A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemv(const AType& A, const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) = prod;
  }
}

/*
  Compute y += A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvAdd(const AType& A, const xType& x,
                                      yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) += prod;
  }
}

/*
  Compute y -= A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvSub(const AType& A, const xType& x,
                                      yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) -= prod;
  }
}

/*
  Compute y = scale * A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvScale(T scale, const AType& A, const xType& x,
                                        yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) = scale * prod;
  }
}

/*
  Compute y += scale * A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvAddScale(T scale, const AType& A,
                                           const xType& x, yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(i, j) * x(j);
    }
    y(i) += scale * prod;
  }
}

/*
  Compute: C = A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemm(const AType& A, const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) = prod;
    }
  }
}

/*
  Compute: C += A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmAdd(const AType& A, const BType& B,
                                      CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) += prod;
    }
  }
}

/*
  Compute: C -= A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmSub(const AType& A, const BType& B,
                                      CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) -= prod;
    }
  }
}

/*
  Compute: C = scale * A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmScale(T scale, const AType& A, const BType& B,
                                        CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) = scale * prod;
    }
  }
}

/*
  Compute: C += scale * A * B

  A in M x N
  B in N x P
  C in M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmAddScale(T scale, const AType& A,
                                           const BType& B, CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(k, j);
      }
      C(i, j) += scale * prod;
    }
  }
}

extern "C" {
extern void dgelss_(int* m, int* n, int* nrhs, double* a, int* lda, double* b,
                    int* ldb, double* s, double* rcond, int* rank, double* work,
                    int* lwork, int* info);

extern void zgelss_(int* m, int* n, int* nrhs, std::complex<double>* a,
                    int* lda, std::complex<double>* b, int* ldb, double* s,
                    double* rcond, int* rank, std::complex<double>* work,
                    int* lwork, double* rwork, int* info);
}

/*
  Compute the pseudo-inverse when A is singular: Ainv = A^{-1}
*/
template <index_t N>
int blockPseudoInverse(A2D::MultiArraySlice<double, N, N>& A,
                       A2D::Mat<double, N, N>& Ainv) {
  // Populate the diaginal matrix
  for (int ii = 0; ii < N; ii++) {
    Ainv(ii, ii) = 1.0;
  }

  int fail = 1;

  int m = N;
  int n = N;
  int nrhs = N;
  double* a = A.get_pointer();
  int lda = N;
  double* b = Ainv.A;
  int ldb = N;
  double s[N];
  double rcond = 1e-12;
  int rank;
  int lwork = 5 * N;
  double work[5 * N];
  dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
          &fail);

  return fail;
}

template <index_t N>
int blockPseudoInverse(A2D::MultiArraySlice<std::complex<double>, N, N>& A,
                       A2D::Mat<std::complex<double>, N, N>& Ainv) {
  // Populate the diaginal matrix
  for (int ii = 0; ii < N; ii++) {
    Ainv(ii, ii) = 1.0;
  }

  int fail = 1;

  int m = N;
  int n = N;
  int nrhs = N;
  std::complex<double>* a = A.get_pointer();
  int lda = N;
  std::complex<double>* b = Ainv.A;
  int ldb = N;
  double s[N];
  double rcond = 1e-12;
  int rank;
  int lwork = 5 * N;
  std::complex<double> work[5 * N];
  double rwork[5 * N];
  zgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
          rwork, &fail);

  return fail;
}

/*
  Compute: Ainv = A^{-1} with pivoting
*/
template <typename T, int N, class AType, class AinvType, class IType>
int blockInverse(AType& A, AinvType& Ainv, IType& ipiv) {
  int fail = 0;

  for (int k = 0; k < N - 1; k++) {
    // Find the maximum value and use it as the pivot
    int r = k;
    T maxv = A(k, k);
    for (int j = k + 1; j < N; j++) {
      T t = A(j, k);
      if (absfunc(t) > absfunc(maxv)) {
        maxv = t;
        r = j;
      }
    }

    ipiv(k) = r;

    // If a swap is required, swap the rows
    if (r != k) {
      for (int j = 0; j < N; j++) {
        T t = A(k, j);
        A(k, j) = A(r, j);
        A(r, j) = t;
      }
    }

    if (absfunc(A(k, k)) == 0.0) {
      fail = k + 1;
      return fail;
    }

    for (int i = k + 1; i < N; i++) {
      A(i, k) = A(i, k) / A(k, k);
    }

    for (int i = k + 1; i < N; i++) {
      for (int j = k + 1; j < N; j++) {
        A(i, j) -= A(i, k) * A(k, j);
      }
    }
  }

  // Now, compute the matrix-inverse
  for (int k = 0; k < N; k++) {
    int ip = k;
    for (int i = 0; i < N - 1; i++) {
      if (ip == ipiv(i)) {
        ip = i;
      } else if (ip == i) {
        ip = ipiv(i);
      }
    }

    for (int i = 0; i < ip; i++) {
      Ainv(i, k) = 0.0;
    }

    Ainv(ip, k) = 1.0;

    for (int i = ip + 1; i < N; i++) {
      Ainv(i, k) = 0.0;
      for (int j = ip; j < i; j++) {
        Ainv(i, k) -= A(i, j) * Ainv(j, k);
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int j = i + 1; j < N; j++) {
        Ainv(i, k) -= A(i, j) * Ainv(j, k);
      }
      Ainv(i, k) = Ainv(i, k) / A(i, i);
    }
  }

  return fail;
}

}  // namespace A2D

#endif  // A2D_BLOCK_NUMERIC_H
