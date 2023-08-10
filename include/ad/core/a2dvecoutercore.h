#ifndef A2D_VEC_OUTER_CORE_H
#define A2D_VEC_OUTER_CORE_H

/*
  Compute the outer product of two vectors

  A = alpha * x * y^{T}

  Or

  A = A + alpha * x * y^{T}
*/
template <typename T, int M, int N, bool additive = false>
inline void VecOuterCore(const T alpha, const T x[], const T y[], T A[]) {
  if constexpr (additive) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] += alpha * x[i] * y[j];
      }
    }

  } else {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] += alpha * x[i] * y[j];
      }
    }
  }
}

#endif  //  A2D_VEC_OUTER_CORE_H