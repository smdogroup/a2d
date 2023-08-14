#ifndef A2D_MAT_VEC_CORE_H
#define A2D_MAT_VEC_CORE_H

#include "a2denum.h"

namespace A2D {
/*
  Compute the matrix-vector products

  y = op(A) * x

  Or

  y = alpha * op(A) * x
*/
template <typename T, int M, int N, MatOp opA = MatOp::NORMAL,
          bool additive = false>
inline void MatVecCore(const T A[], const T x[], T y[]) noexcept {
  if constexpr (additive) {
    if consexpr (opA == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        T value = 0.0;
        for (int j = 0; j < N; j++, A++) {
          value += A[0] * x[j];
        }
        y[i] += value;
      }
    } else {
      for (int i = 0; i < M; i++) {
        const T value = x[i];
        for (int j = 0; j < N; j++, A++) {
          y[j] += A[0] * value;
        }
      }
    }
  } else {
    if consexpr (opA == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        T value = 0.0;
        for (int j = 0; j < N; j++, A++) {
          value += A[0] * x[j];
        }
        y[i] = value;
      }
    } else {
      for (int j = 0; j < N; j++) {
        y[j] = T(0.0);
      }

      for (int i = 0; i < M; i++) {
        const T value = x[i];

        for (int j = 0; j < N; j++, A++) {
          y[j] += A[0] * value;
        }
      }
    }
  }
}

template <typename T, int M, int N, MatOp opA = MatOp::NORMAL,
          bool additive = false>
inline void MatVecCoreScale(const T alpha, const T A[], const T x[],
                            T y[]) noexcept {
  if constexpr (additive) {
    if consexpr (opA == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        T value = 0.0;
        for (int j = 0; j < N; j++, A++) {
          value += A[0] * x[j];
        }

        y[i] += alpha * value;
      }
    } else {
      for (int i = 0; i < M; i++) {
        const T value = alpha * x[i];
        for (int j = 0; j < N; j++, A++) {
          y[j] += A[0] * value;
        }
      }
    }
  } else {
    if consexpr (opA == MatOp::NORMAL) {
      for (int i = 0; i < M; i++) {
        T value = 0.0;
        for (int j = 0; j < N; j++, A++) {
          value += A[0] * x[j];
        }

        y[i] = alpha * value;
      }
    } else {
      for (int j = 0; j < N; j++) {
        y[j] = T(0.0);
      }

      for (int i = 0; i < M; i++) {
        const T value = alpha * x[i];

        for (int j = 0; j < N; j++, A++) {
          y[j] += A[0] * value;
        }
      }
    }
  }
}

}  // namespace A2D

#endif  //  A2D_MAT_VEC_CORE_H