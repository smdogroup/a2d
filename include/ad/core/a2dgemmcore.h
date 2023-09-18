#ifndef A2D_GEMMCORE_H
#define A2D_GEMMCORE_H

#include "a2ddefs.h"

namespace A2D {

template <bool B, int i, int j>
struct int_conditional {};

template <int i, int j>
struct int_conditional<true, i, j> {
  static constexpr int value = i;
};

template <int i, int j>
struct int_conditional<false, i, j> {
  static constexpr int value = j;
};

template <typename T, MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
KOKKOS_FUNCTION void MatMatMultCore3x3(const T A[], const T B[], T C[]) {
  if constexpr (opA == MatOp::NORMAL and opB == MatOp::NORMAL) {
    C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
    C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
    C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::NORMAL) {
    C[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
    C[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
    C[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];
    C[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
    C[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
    C[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];
    C[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
    C[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
    C[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
  } else if constexpr (opA == MatOp::NORMAL and opB == MatOp::TRANSPOSE) {
    C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
    C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
    C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
    C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
    C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
    C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
    C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
    C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::TRANSPOSE) {
    C[0] = A[0] * B[0] + A[3] * B[1] + A[6] * B[2];
    C[1] = A[0] * B[3] + A[3] * B[4] + A[6] * B[5];
    C[2] = A[0] * B[6] + A[3] * B[7] + A[6] * B[8];
    C[3] = A[1] * B[0] + A[4] * B[1] + A[7] * B[2];
    C[4] = A[1] * B[3] + A[4] * B[4] + A[7] * B[5];
    C[5] = A[1] * B[6] + A[4] * B[7] + A[7] * B[8];
    C[6] = A[2] * B[0] + A[5] * B[1] + A[8] * B[2];
    C[7] = A[2] * B[3] + A[5] * B[4] + A[8] * B[5];
    C[8] = A[2] * B[6] + A[5] * B[7] + A[8] * B[8];
  }
}

template <typename T, MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
KOKKOS_FUNCTION void MatMatMultCore3x3Scale(T scalar, const T A[], const T B[],
                                            T C[]) {
  if constexpr (opA == MatOp::NORMAL and opB == MatOp::NORMAL) {
    C[0] = scalar * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
    C[1] = scalar * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
    C[2] = scalar * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
    C[3] = scalar * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
    C[4] = scalar * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
    C[5] = scalar * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
    C[6] = scalar * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
    C[7] = scalar * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
    C[8] = scalar * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::NORMAL) {
    C[0] = scalar * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
    C[1] = scalar * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
    C[2] = scalar * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
    C[3] = scalar * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
    C[4] = scalar * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
    C[5] = scalar * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
    C[6] = scalar * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
    C[7] = scalar * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
    C[8] = scalar * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::NORMAL and opB == MatOp::TRANSPOSE) {
    C[0] = scalar * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
    C[1] = scalar * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
    C[2] = scalar * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
    C[3] = scalar * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
    C[4] = scalar * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
    C[5] = scalar * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
    C[6] = scalar * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
    C[7] = scalar * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
    C[8] = scalar * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::TRANSPOSE) {
    C[0] = scalar * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
    C[1] = scalar * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
    C[2] = scalar * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
    C[3] = scalar * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
    C[4] = scalar * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
    C[5] = scalar * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
    C[6] = scalar * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
    C[7] = scalar * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
    C[8] = scalar * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
  }
}

template <typename T, MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
KOKKOS_FUNCTION void MatMatMultCore3x3Add(const T A[], const T B[], T C[]) {
  if constexpr (opA == MatOp::NORMAL and opB == MatOp::NORMAL) {
    C[0] += A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    C[1] += A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    C[2] += A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
    C[3] += A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    C[4] += A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    C[5] += A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
    C[6] += A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    C[7] += A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    C[8] += A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::NORMAL) {
    C[0] += A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
    C[1] += A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
    C[2] += A[0] * B[2] + A[3] * B[5] + A[6] * B[8];
    C[3] += A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
    C[4] += A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
    C[5] += A[1] * B[2] + A[4] * B[5] + A[7] * B[8];
    C[6] += A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
    C[7] += A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
    C[8] += A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
  } else if constexpr (opA == MatOp::NORMAL and opB == MatOp::TRANSPOSE) {
    C[0] += A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    C[1] += A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
    C[2] += A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
    C[3] += A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
    C[4] += A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
    C[5] += A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
    C[6] += A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
    C[7] += A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
    C[8] += A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::TRANSPOSE) {
    C[0] += A[0] * B[0] + A[3] * B[1] + A[6] * B[2];
    C[1] += A[0] * B[3] + A[3] * B[4] + A[6] * B[5];
    C[2] += A[0] * B[6] + A[3] * B[7] + A[6] * B[8];
    C[3] += A[1] * B[0] + A[4] * B[1] + A[7] * B[2];
    C[4] += A[1] * B[3] + A[4] * B[4] + A[7] * B[5];
    C[5] += A[1] * B[6] + A[4] * B[7] + A[7] * B[8];
    C[6] += A[2] * B[0] + A[5] * B[1] + A[8] * B[2];
    C[7] += A[2] * B[3] + A[5] * B[4] + A[8] * B[5];
    C[8] += A[2] * B[6] + A[5] * B[7] + A[8] * B[8];
  }
}

template <typename T, MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
KOKKOS_FUNCTION void MatMatMultCore3x3ScaleAdd(T scalar, const T A[],
                                               const T B[], T C[]) {
  if constexpr (opA == MatOp::NORMAL and opB == MatOp::NORMAL) {
    C[0] += scalar * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
    C[1] += scalar * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
    C[2] += scalar * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
    C[3] += scalar * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
    C[4] += scalar * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
    C[5] += scalar * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
    C[6] += scalar * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
    C[7] += scalar * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
    C[8] += scalar * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::NORMAL) {
    C[0] += scalar * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
    C[1] += scalar * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
    C[2] += scalar * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
    C[3] += scalar * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
    C[4] += scalar * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
    C[5] += scalar * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
    C[6] += scalar * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
    C[7] += scalar * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
    C[8] += scalar * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::NORMAL and opB == MatOp::TRANSPOSE) {
    C[0] += scalar * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
    C[1] += scalar * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
    C[2] += scalar * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
    C[3] += scalar * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
    C[4] += scalar * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
    C[5] += scalar * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
    C[6] += scalar * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
    C[7] += scalar * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
    C[8] += scalar * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
  } else if constexpr (opA == MatOp::TRANSPOSE and opB == MatOp::TRANSPOSE) {
    C[0] += scalar * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
    C[1] += scalar * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
    C[2] += scalar * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
    C[3] += scalar * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
    C[4] += scalar * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
    C[5] += scalar * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
    C[6] += scalar * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
    C[7] += scalar * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
    C[8] += scalar * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
  }
}

/**
 * @brief mat-mat multiplication C = alpha * Op(A) * Op(B), where op is normal
 * or transpose
 *
 * Note: consistency of the dimensions are not checked so check them when
 * calling!
 *
 * @tparam T: scalar type
 * @tparam Anrows: number of rows of A
 * @tparam Ancols: number of cols of A
 * @tparam Bnrows: number of rows of B
 * @tparam Bncols: number of cols of B
 * @tparam Cnrows: number of rows of C
 * @tparam Cncols: number of cols of C
 * @tparam opA: transpose A or not
 * @tparam opB: transpose B or not
 * @tparam additive: true for increment (C += ..), false for assignment (C = ..)
 *
 * @param[in] A: Anrows-by-Ancols matrix
 * @param[in] B: Bnrows-by-Bncols matrix
 * @param[in, out] C: Cnrows-by-Cncols matrix
 * @param[in] alpha: the scalar
 */
template <typename T, int Anrows, int Ancols, int Bnrows, int Bncols,
          int Cnrows, int Cncols, MatOp opA = MatOp::NORMAL,
          MatOp opB = MatOp::NORMAL, bool additive = false>
KOKKOS_FUNCTION void MatMatMultCoreGeneral(const T A[], const T B[], T C[]) {
  // Op(A) is M-by-P, Op(B) is P-by-N, C is M-by-N
  constexpr int M = Cnrows;
  constexpr int N = Cncols;
  constexpr int P =
      int_conditional<opA == MatOp::NORMAL, Ancols, Anrows>::value;

  if constexpr (opA == MatOp::NORMAL) {
    if (opB == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[Ancols * i];
          const T* aend = a + Ancols;
          const T* b = &B[j];

          T value = T(0.0);
          for (; a < aend; a++, b += Bncols) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += value;
          } else {
            C[0] = value;
          }
        }
      }
    } else {  // opB == MatOp::TRANSPOSE
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[Ancols * i];
          const T* aend = a + Ancols;
          const T* b = &B[Bncols * j];

          T value = T(0.0);
          for (; a < aend; a++, b++) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += value;
          } else {
            C[0] = value;
          }
        }
      }
    }
  } else {  // opA == MatOp::TRANSPOSE
    if (opB == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[i];
          const T* b = &B[j];
          const T* bend = b + Bnrows * Bncols;

          T value = T(0.0);
          for (; b < bend; a += Ancols, b += Bncols) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += value;
          } else {
            C[0] = value;
          }
        }
      }
    } else {  // opB == MatOp::TRANSPOSE
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[i];
          const T* b = &B[Bncols * j];
          const T* bend = b + Bncols;

          T value = T(0.0);
          for (; b < bend; a += Ancols, b++) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += value;
          } else {
            C[0] = value;
          }
        }
      }
    }
  }
}

/**
 * @brief mat-mat multiplication C = alpha * Op(A) * Op(B), where op is
 * normal or transpose
 *
 * Note: consistency of the dimensions are not checked so check them when
 * calling!
 *
 * @tparam T: scalar type
 * @tparam Anrows: number of rows of A
 * @tparam Ancols: number of cols of A
 * @tparam Bnrows: number of rows of B
 * @tparam Bncols: number of cols of B
 * @tparam Cnrows: number of rows of C
 * @tparam Cncols: number of cols of C
 * @tparam opA: transpose A or not
 * @tparam opB: transpose B or not
 * @tparam additive: true for increment (C += ..), false for assignment (C =
 * ..)
 *
 * @param[in] A: Anrows-by-Ancols matrix
 * @param[in] B: Bnrows-by-Bncols matrix
 * @param[in, out] C: Cnrows-by-Cncols matrix
 * @param[in] alpha: the scalar
 */
template <typename T, int Anrows, int Ancols, int Bnrows, int Bncols,
          int Cnrows, int Cncols, MatOp opA = MatOp::NORMAL,
          MatOp opB = MatOp::NORMAL, bool additive = false>
KOKKOS_FUNCTION void MatMatMultScaleCoreGeneral(const T alpha, const T A[],
                                                const T B[], T C[]) {
  // Op(A) is M-by-P, Op(B) is P-by-N, C is M-by-N
  constexpr int M = Cnrows;
  constexpr int N = Cncols;
  constexpr int P =
      int_conditional<opA == MatOp::NORMAL, Ancols, Anrows>::value;

  if constexpr (opA == MatOp::NORMAL) {
    if (opB == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[Ancols * i];
          const T* aend = a + Ancols;
          const T* b = &B[j];

          T value = T(0.0);
          for (; a < aend; a++, b += Bncols) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += alpha * value;
          } else {
            C[0] = alpha * value;
          }
        }
      }
    } else {  // opB == MatOp::TRANSPOSE
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[Ancols * i];
          const T* aend = a + Ancols;
          const T* b = &B[Bncols * j];

          T value = T(0.0);
          for (; a < aend; a++, b++) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += alpha * value;
          } else {
            C[0] = alpha * value;
          }
        }
      }
    }
  } else {  // opA == MatOp::TRANSPOSE
    if (opB == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[i];
          const T* b = &B[j];
          const T* bend = b + Bnrows * Bncols;

          T value = T(0.0);
          for (; b < bend; a += Ancols, b += Bncols) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += alpha * value;
          } else {
            C[0] = alpha * value;
          }
        }
      }
    } else {  // opB == MatOp::TRANSPOSE
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++, C++) {
          const T* a = &A[i];
          const T* b = &B[Bncols * j];
          const T* bend = b + Bncols;

          T value = T(0.0);
          for (; b < bend; a += Ancols, b++) {
            value += a[0] * b[0];
          }

          if constexpr (additive) {
            C[0] += alpha * value;
          } else {
            C[0] = alpha * value;
          }
        }
      }
    }
  }
}

/**
 * @brief matrix-matrix multiplication C = alpha * Op(A) * Op(B), where op
 * is normal (nominal) or transpose
 *
 * @note This function acts as a dispatcher for specialized versions of the
 * matrix multiplication function, falling back on the general
 * implementation. This also verifies the dimensions of the inputs and
 * outputs are appropriate for the multiplication to be executed.
 *
 * @tparam T: scalar type
 * @tparam Anrows: number of rows of A
 * @tparam Ancols: number of cols of A
 * @tparam Bnrows: number of rows of B
 * @tparam Bncols: number of cols of B
 * @tparam Cnrows: number of rows of C
 * @tparam Cncols: number of cols of C
 * @tparam opA: transpose A or not
 * @tparam opB: transpose B or not
 * @tparam additive: true for increment (C += ..), false for assignment (C =
 * ..)
 * @tparam scale: true if result should be scaled by alpha before output
 *
 * @param[in] A: Anrows-by-Ancols matrix
 * @param[in] B: Bnrows-by-Bncols matrix
 * @param[in, out] C: Cnrows-by-Cncols matrix
 * @param[in] alpha: the scalar
 */
template <typename T, int Anrows, int Ancols, int Bnrows, int Bncols,
          int Cnrows, int Cncols, MatOp opA = MatOp::NORMAL,
          MatOp opB = MatOp::NORMAL, bool additive = false>
KOKKOS_FUNCTION void MatMatMultCore(const T A[], const T B[], T C[]) {
  // Check if shapes are consistent
  if constexpr (opA == MatOp::TRANSPOSE && opB == MatOp::TRANSPOSE) {
    static_assert(Anrows == Bncols && Ancols == Cnrows && Bnrows == Cncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opA == MatOp::TRANSPOSE) {
    static_assert(Anrows == Bnrows && Ancols == Cnrows && Bncols == Cncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opB == MatOp::TRANSPOSE) {
    static_assert(Ancols == Bncols && Anrows == Cnrows && Bnrows == Cncols,
                  "Matrix dimensions must agree.");
  } else {
    static_assert(Ancols == Bnrows && Anrows == Cnrows && Bncols == Cncols,
                  "Matrix dimensions must agree.");
  }

  if constexpr (Anrows == 3 && Ancols == 3 && Bnrows == 3 && Bncols == 3) {
    if constexpr (additive) {
      MatMatMultCore3x3Add<T, opA, opB>(A, B, C);
    } else {
      MatMatMultCore3x3<T, opA, opB>(A, B, C);
    }
  } else {  // The general fallback implmentation
    MatMatMultCoreGeneral<T, Anrows, Ancols, Bnrows, Bncols, Cnrows, Cncols,
                          opA, opB, additive>(A, B, C);
  }
}

template <typename T, int Anrows, int Ancols, int Bnrows, int Bncols,
          int Cnrows, int Cncols, MatOp opA = MatOp::NORMAL,
          MatOp opB = MatOp::NORMAL, bool additive = false>
inline void MatMatMultScaleCore(T alpha, const T A[], const T B[], T C[]) {
  // Check if shapes are consistent
  if constexpr (opA == MatOp::TRANSPOSE && opB == MatOp::TRANSPOSE) {
    static_assert(Anrows == Bncols && Ancols == Cnrows && Bnrows == Cncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opA == MatOp::TRANSPOSE) {
    static_assert(Anrows == Bnrows && Ancols == Cnrows && Bncols == Cncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opB == MatOp::TRANSPOSE) {
    static_assert(Ancols == Bncols && Anrows == Cnrows && Bnrows == Cncols,
                  "Matrix dimensions must agree.");
  } else {
    static_assert(Ancols == Bnrows && Anrows == Cnrows && Bncols == Cncols,
                  "Matrix dimensions must agree.");
  }

  if constexpr (Anrows == 3 && Ancols == 3 && Bnrows == 3 && Bncols == 3) {
    if constexpr (additive) {
      MatMatMultCore3x3ScaleAdd<T, opA, opB>(alpha, A, B, C);
    } else {
      MatMatMultCore3x3Scale<T, opA, opB>(alpha, A, B, C);
    }
  } else {  // The general fallback implmentation
    MatMatMultScaleCoreGeneral<T, Anrows, Ancols, Bnrows, Bncols, Cnrows,
                               Cncols, opA, opB, additive>(alpha, A, B, C);
  }
}

}  // namespace A2D
#endif  // A2D_GEMMCORE_H