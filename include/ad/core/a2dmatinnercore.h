#ifndef A2D_MAT_INNER_CORE_H
#define A2D_MAT_INNER_CORE_H

namespace A2D {

template <typename T, int M, int N>
inline T MatInnerCore(const T A[], const T x[], const T y[]) noexcept {
  T value = 0.0;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      value += x[i] * A[0] * y[j];
    }
  }

  return value;
}

}  // namespace A2D

#endif  // A2D_MAT_INNER_CORE_H