#ifndef A2D_MAT_DET_CORE_H
#define A2D_MAT_DET_CORE_H

#include "a2denum.h"

namespace A2D {

template <typename T, int N>
T MatDetCore(const T A[]) {
  static_assert((N >= 1 && N <= 3), "MatDet not implemented for N >= 4");

  if constexpr (N == 1) {
    return A[0];
  } else if constexpr (N == 2) {
    return A[0] * A[3] - A[1] * A[2];
  } else {  // N == 3
    T det = (A[8] * (A[0] * A[4] - A[3] * A[1]) -
             A[7] * (A[0] * A[5] - A[3] * A[2]) +
             A[6] * (A[1] * A[5] - A[2] * A[4]));
    return det;
  }
}

template <typename T, int N>
T MatDetForwardCore(const T A[], const T Ad[]) {
  static_assert((N >= 1 && N <= 3),
                "MatDetForwardCore not implemented for N >= 4");

  if constexpr (N == 1) {
    return Ad[0];
  } else if constexpr (N == 2) {
    T detd = Ad[0] * A[3] + A[0] * Ad[3] - Ad[1] * A[2] - A[1] * Ad[2];
    return detd;
  } else {  // N == 3
    T detd = (Ad[0] * (A[8] * A[4] - A[7] * A[5]) +
              Ad[1] * (A[6] * A[5] - A[8] * A[3]) +
              Ad[2] * (A[7] * A[3] - A[6] * A[4]) +
              Ad[3] * (A[7] * A[2] - A[8] * A[1]) +
              Ad[4] * (A[8] * A[0] - A[6] * A[2]) +
              Ad[5] * (A[6] * A[1] - A[7] * A[0]) +
              Ad[6] * (A[1] * A[5] - A[2] * A[4]) +
              Ad[7] * (A[3] * A[2] - A[0] * A[5]) +
              Ad[8] * (A[0] * A[4] - A[3] * A[1]));
    return detd;
  }
}

template <typename T, int N>
void MatDetReverseCore(const T bdet, const T A[], T Ab[]) {
  static_assert((N >= 1 && N <= 3),
                "MatDetReverseCore not implemented for N >= 4");

  if constexpr (N == 1) {
    Ab[0] += bdet;
  } else if constexpr (N == 2) {
    Ab[0] += A[3] * bdet;
    Ab[1] += -A[2] * bdet;
    Ab[2] += -A[1] * bdet;
    Ab[3] += A[0] * bdet;
  } else {  // N == 3
    Ab[0] += (A[8] * A[4] - A[7] * A[5]) * bdet;
    Ab[1] += (A[6] * A[5] - A[8] * A[3]) * bdet;
    Ab[2] += (A[7] * A[3] - A[6] * A[4]) * bdet;
    Ab[3] += (A[7] * A[2] - A[8] * A[1]) * bdet;
    Ab[4] += (A[8] * A[0] - A[6] * A[2]) * bdet;
    Ab[5] += (A[6] * A[1] - A[7] * A[0]) * bdet;
    Ab[6] += (A[1] * A[5] - A[2] * A[4]) * bdet;
    Ab[7] += (A[3] * A[2] - A[0] * A[5]) * bdet;
    Ab[8] += (A[0] * A[4] - A[3] * A[1]) * bdet;
  }
}

template <typename T, int N>
void MatDetHReverseCore(const T bdet, const T hdet, const T A[], const T Ap[],
                        T Ah[]) {
  if constexpr (N == 1) {
    Ah[0] += hdet;
  } else if constexpr (N == 2) {
    Ah[0] += Ap[3] * bdet;
    Ah[1] += -Ap[2] * bdet;
    Ah[2] += -Ap[1] * bdet;
    Ah[3] += Ap[0] * bdet;

    Ah[0] += A[3] * hdet;
    Ah[1] += -A[2] * hdet;
    Ah[2] += -A[1] * hdet;
    Ah[3] += A[0] * hdet;
  } else {  // N == 3
    Ah[0] += (A[8] * Ap[4] - A[7] * Ap[5] + Ap[8] * A[4] - Ap[7] * A[5]) * bdet;
    Ah[1] += (A[6] * Ap[5] - A[8] * Ap[3] + Ap[6] * A[5] - Ap[8] * A[3]) * bdet;
    Ah[2] += (A[7] * Ap[3] - A[6] * Ap[4] + Ap[7] * A[3] - Ap[6] * A[4]) * bdet;
    Ah[3] += (A[7] * Ap[2] - A[8] * Ap[1] + Ap[7] * A[2] - Ap[8] * A[1]) * bdet;
    Ah[4] += (A[8] * Ap[0] - A[6] * Ap[2] + Ap[8] * A[0] - Ap[6] * A[2]) * bdet;
    Ah[5] += (A[6] * Ap[1] - A[7] * Ap[0] + Ap[6] * A[1] - Ap[7] * A[0]) * bdet;
    Ah[6] += (A[1] * Ap[5] - A[2] * Ap[4] + Ap[1] * A[5] - Ap[2] * A[4]) * bdet;
    Ah[7] += (A[3] * Ap[2] - A[0] * Ap[5] + Ap[3] * A[2] - Ap[0] * A[5]) * bdet;
    Ah[8] += (A[0] * Ap[4] - A[3] * Ap[1] + Ap[0] * A[4] - Ap[3] * A[1]) * bdet;

    Ah[0] += (A[8] * A[4] - A[7] * A[5]) * hdet;
    Ah[1] += (A[6] * A[5] - A[8] * A[3]) * hdet;
    Ah[2] += (A[7] * A[3] - A[6] * A[4]) * hdet;
    Ah[3] += (A[7] * A[2] - A[8] * A[1]) * hdet;
    Ah[4] += (A[8] * A[0] - A[6] * A[2]) * hdet;
    Ah[5] += (A[6] * A[1] - A[7] * A[0]) * hdet;
    Ah[6] += (A[1] * A[5] - A[2] * A[4]) * hdet;
    Ah[7] += (A[3] * A[2] - A[0] * A[5]) * hdet;
    Ah[8] += (A[0] * A[4] - A[3] * A[1]) * hdet;
  }
}

}  // namespace A2D

#endif  // A2D_MAT_DET_CORE_H