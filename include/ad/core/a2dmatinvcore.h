#ifndef A2D_MAT_INV_CORE_H
#define A2D_MAT_INV_CORE_H

template <typename T, int N>
void MatInvCore(const T A[], T Ainv[]) {
  static_assert((N >= 1 && N <= 3), "MatInvCore not implemented for N >= 4");

  if constexpr (N == 1) {
    Ainv[0] = 1.0 / A[0];
  } else if constexpr (N == 2) {
    T det = A[0] * A[3] - A[1] * A[2];
    T detinv = 1.0 / det;

    Ainv[0] = A[3] * detinv;
    Ainv[1] = -A[1] * detinv;
    Ainv[2] = -A[2] * detinv;
    Ainv[3] = A[0] * detinv;
  } else {  // N == 3
    T det = (A[8] * (A[0] * A[4] - A[3] * A[1]) -
             A[7] * (A[0] * A[5] - A[3] * A[2]) +
             A[6] * (A[1] * A[5] - A[2] * A[4]));
    T detinv = 1.0 / det;

    Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * detinv;
    Ainv[1] = -(A[1] * A[8] - A[2] * A[7]) * detinv;
    Ainv[2] = (A[1] * A[5] - A[2] * A[4]) * detinv;

    Ainv[3] = -(A[3] * A[8] - A[5] * A[6]) * detinv;
    Ainv[4] = (A[0] * A[8] - A[2] * A[6]) * detinv;
    Ainv[5] = -(A[0] * A[5] - A[2] * A[3]) * detinv;

    Ainv[6] = (A[3] * A[7] - A[4] * A[6]) * detinv;
    Ainv[7] = -(A[0] * A[7] - A[1] * A[6]) * detinv;
    Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * detinv;
  }
}

#endif  // A2D_MAT_INV_CORE_H