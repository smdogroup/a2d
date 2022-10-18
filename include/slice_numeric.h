/**
 * @file slice_numeric.h
 * @author Yicong Fu (yfu97@gatech.edu)
 * @brief Similar to block_numeric.h but perform operations slices of bigger
 * matrices.
 * @date 2022-10-18
 *
 * Note: This is supposed to be a set of temporary implementation and should be
 * replaced by using kokkos-kernels in a more general (and performant) way.
 */

#ifndef A2D_SLICE_NUMERIC_H
#define A2D_SLICE_NUMERIC_H

#include "a2dobjs.h"

namespace A2D {

/*
  Compute y[k, :] = A[i, :, :] * x[j, :]

  A is I x M x N
  x is J x N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvSlice(const AType& A, const int Ai,
                                        const xType& x, const int xj, yType& y,
                                        const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(xj, j);
    }
    y(yk, i) = prod;
  }
}

/*
  Compute y[k, :] = A[i, :, :] * x[:]

  A is I x M x N
  x is N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvSlice(const AType& A, const int Ai,
                                        const xType& x, yType& y,
                                        const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(j);
    }
    y(yk, i) = prod;
  }
}

/*
  Compute y[k, :] += A[i, :, :] * x[j, :]

  A is I x M x N
  x is J x N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvAddSlice(const AType& A, const int Ai,
                                           const xType& x, const int xj,
                                           yType& y, const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(xj, j);
    }
    y(yk, i) += prod;
  }
}

/*
  Compute y[k, :] += scale * A[i, :, :] * x[j, :]

  A is I x M x N
  x is J x N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvAddScaleSlice(const T scale, const AType& A,
                                                const int Ai, const xType& x,
                                                const int xj, yType& y,
                                                const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(xj, j);
    }
    y(yk, i) += scale * prod;
  }
}

/*
  Compute y[k, :] += scale * A[i, :, :] * x[:]

  A is I x M x N
  x is N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvAddScaleSlice(const T scale, const AType& A,
                                                const int Ai, const xType& x,
                                                yType& y, const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(j);
    }
    y(yk, i) += scale * prod;
  }
}

/*
  Compute y[k, :] -= A[i, :, :] * x[j, :]

  A is I x M x N
  x is J x N
  y is K x M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvSubSlice(const AType& A, const int Ai,
                                           const xType& x, const int xj,
                                           yType& y, const int yk) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(xj, j);
    }
    y(yk, i) -= prod;
  }
}

/*
  Compute y[:] -= A[i, :, :] * x[j, :]

  A is I x M x N
  x is J x N
  y is M
*/
template <typename T, int M, int N, class AType, class xType, class yType>
A2D_INLINE_FUNCTION void blockGemvSubSlice(const AType& A, const int Ai,
                                           const xType& x, const int xj,
                                           yType& y) {
  for (int i = 0; i < M; i++) {
    T prod = 0.0;
    for (int j = 0; j < N; j++) {
      prod += A(Ai, i, j) * x(xj, j);
    }
    y(i) -= prod;
  }
}

/*
  Compute: C[k, :, :] = A[i, :, :] * B[j, :, :]

  A is I x M x N
  B is J x N x P
  C is K x M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmSlice(const AType& A, const int Ai,
                                        const BType& B, const int Bj, CType& C,
                                        const int Ck) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(Ai, i, k) * B(Bj, k, j);
      }
      C(Ck, i, j) = prod;
    }
  }
}

/*
  Compute: C[:, :] = A[i, :, :] * B[j, :, :]

  A is I x M x N
  B is J x N x P
  C is M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmSlice(const AType& A, const int Ai,
                                        const BType& B, const int Bj,
                                        CType& C) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(Ai, i, k) * B(Bj, k, j);
      }
      C(i, j) = prod;
    }
  }
}

/*
  Compute: C[k, :, :] += A[i, :, :] * B[j, :, :]

  A is I x M x N
  B is J x N x P
  C is K x M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmAddSlice(const AType& A, const int Ai,
                                           const BType& B, const int Bj,
                                           CType& C, const int Ck) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(Ai, i, k) * B(Bj, k, j);
      }
      C(Ck, i, j) += prod;
    }
  }
}

/*
  Compute: C[k, :, :] += scale A[i, :, :] * B[j, :, :]

  A is I x M x N
  B is J x N x P
  C is K x M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmAddScaleSlice(const T scale, const AType& A,
                                                const int Ai, const BType& B,
                                                const int Bj, CType& C,
                                                const int Ck) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(Ai, i, k) * B(Bj, k, j);
      }
      C(Ck, i, j) += scale * prod;
    }
  }
}

/*
  Compute: C[k, :, :] -= A[i, :, :] * B[j, :, :]

  A is I x M x N
  B is J x N x P
  C is K x M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmSubSlice(const AType& A, const int Ai,
                                           const BType& B, const int Bj,
                                           CType& C, const int Ck) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(Ai, i, k) * B(Bj, k, j);
      }
      C(Ck, i, j) -= prod;
    }
  }
}

/*
  Compute: C[k, :, :] -= A[:, :] * B[j, :, :]

  A is M x N
  B is J x N x P
  C is K x M x P
*/
template <typename T, int M, int N, int P, class AType, class BType,
          class CType>
A2D_INLINE_FUNCTION void blockGemmSubSlice(const AType& A, const BType& B,
                                           const int Bj, CType& C,
                                           const int Ck) {
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < P; j++) {
      T prod = 0.0;
      for (int k = 0; k < N; k++) {
        prod += A(i, k) * B(Bj, k, j);
      }
      C(Ck, i, j) -= prod;
    }
  }
}

/*
  Compute: Ainv = A^{-1} with pivoting
*/
template <typename T, int N, class AType, class AinvType, class IType>
int blockInverseSlice(AType& A, const int Ai, AinvType& Ainv, IType& ipiv) {
  int fail = 0;

  for (int k = 0; k < N - 1; k++) {
    // Find the maximum value and use it as the pivot
    int r = k;
    T maxv = A(Ai, k, k);
    for (int j = k + 1; j < N; j++) {
      T t = A(Ai, j, k);
      if (absfunc(t) > absfunc(maxv)) {
        maxv = t;
        r = j;
      }
    }

    ipiv(k) = r;

    // If a swap is required, swap the rows
    if (r != k) {
      for (int j = 0; j < N; j++) {
        T t = A(Ai, k, j);
        A(Ai, k, j) = A(Ai, r, j);
        A(Ai, r, j) = t;
      }
    }

    if (absfunc(A(Ai, k, k)) == 0.0) {
      fail = k + 1;
      return fail;
    }

    for (int i = k + 1; i < N; i++) {
      A(Ai, i, k) = A(Ai, i, k) / A(Ai, k, k);
    }

    for (int i = k + 1; i < N; i++) {
      for (int j = k + 1; j < N; j++) {
        A(Ai, i, j) -= A(Ai, i, k) * A(Ai, k, j);
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
        Ainv(i, k) -= A(Ai, i, j) * Ainv(j, k);
      }
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int j = i + 1; j < N; j++) {
        Ainv(i, k) -= A(Ai, i, j) * Ainv(j, k);
      }
      Ainv(i, k) = Ainv(i, k) / A(Ai, i, i);
    }
  }

  return fail;
}

}  // namespace A2D

#endif  //  A2D_SLICE_NUMERIC_H