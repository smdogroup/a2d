#ifndef A2D_BLOCK_NUMERIC_H
#define A2D_BLOCK_NUMERIC_H

#include "a2dobjs.h"
#include "array.h"

namespace A2D {

/*
  Compute y = A * x
*/
template <typename T, int M, int N, class AType, class xType, class yType>
KOKKOS_FUNCTION void blockGemv(const AType& A, const xType& x, yType& y) {
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
KOKKOS_FUNCTION void blockGemvAdd(const AType& A, const xType& x,
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
KOKKOS_FUNCTION void blockGemvSub(const AType& A, const xType& x,
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
KOKKOS_FUNCTION void blockGemvScale(T scale, const AType& A, const xType& x,
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
KOKKOS_FUNCTION void blockGemvAddScale(T scale, const AType& A,
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
KOKKOS_FUNCTION void blockGemm(const AType& A, const BType& B, CType& C) {
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
KOKKOS_FUNCTION void blockGemmAdd(const AType& A, const BType& B,
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
KOKKOS_FUNCTION void blockGemmSub(const AType& A, const BType& B,
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
KOKKOS_FUNCTION void blockGemmScale(T scale, const AType& A, const BType& B,
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
KOKKOS_FUNCTION void blockGemmAddScale(T scale, const AType& A,
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
extern void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
extern void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work,
                    int* lwork, int* info);
extern void zgelss_(int* m, int* n, int* nrhs, void* a, int* lda, void* b,
                    int* ldb, double* s, double* rcond, int* rank, void* work,
                    int* lwork, double* rwork, int* info);
extern void zgetrf_(int* m, int* n, void* a, int* lda, int* ipiv, int* info);
extern void zgetri_(int* n, void* a, int* lda, int* ipiv, void* work,
                    int* lwork, int* info);
}

/**
 * @brief Compute the pseudo-inverse when A is singular: Ainv = A^{-1}
 *
 * Note1: different routines get invoked depending on the input type T.
 *
 * If T is complex<double>, then call ZGETRF and ZGETRI to compute real inverse;
 * If T is real, then call DGELSS to compute pseudo-inverse. (ZGELSS doesn't
 * converge to the precision required by the small complex step, e.g. h = 1e-30)
 *
 * Note2: LAPACK routines require input matrices to store in continuous memory
 * chunk. If FLayout is used, need to prepare the matrix first.
 */
template <typename T, int N, class AType>
int blockPseudoInverse(AType& A, Mat<T, N, N>& Ainv) {
  // Populate the diaginal matrix
  for (int ii = 0; ii < N; ii++) {
    Ainv(ii, ii) = 1.0;
  }

  // Decide which branch to go
  const bool is_layout_right =
      std::is_same<typename AType::array_layout, Kokkos::LayoutRight>::value;

  // Set parameters
  int m = N;
  int n = N;
  int nrhs = N;
  int lda = N;
  T* b = Ainv.data();
  int ldb = N;
  double rcond = -1;
  int rank;
  int lwork = 5 * N;
  T work[5 * N];
  int fail = -1;

  // Prepare A
  T* a;
  if constexpr (!is_layout_right) {
    a = new T[N * N];
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        a[i * N + j] = A(i, j);
      }
    }
  } else {
    a = A.data();
  }

  if constexpr (is_complex<T>::value) {
    int ipiv[N];
    zgetrf_(&m, &n, a, &lda, ipiv, &fail);
    zgetri_(&n, a, &lda, ipiv, work, &lwork, &fail);
    for (int ii = 0; ii != N * N; ii++) {
      b[ii] = a[ii];
    }
  } else {
    double s[N];
    dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
            &fail);
  }

  if constexpr (!is_layout_right) {
    delete[] a;
  }

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
