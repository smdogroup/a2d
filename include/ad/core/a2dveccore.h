#ifndef A2D_VEC_CORE_H
#define A2D_VEC_CORE_H

#include "a2dobjs.h"

namespace A2D {

template <typename T>
KOKKOS_FUNCTION void VecCrossCore(const T x[], const T y[], T v[]) {
  v[0] = x[1] * y[2] - x[2] * y[1];
  v[1] = x[2] * y[0] - x[0] * y[2];
  v[2] = x[0] * y[1] - x[1] * y[0];
}

template <typename T>
KOKKOS_FUNCTION void VecCrossCoreAdd(const T x[], const T y[], T v[]) {
  v[0] += x[1] * y[2] - x[2] * y[1];
  v[1] += x[2] * y[0] - x[0] * y[2];
  v[2] += x[0] * y[1] - x[1] * y[0];
}

template <typename T, int size>
KOKKOS_FUNCTION void VecZeroCore(T A[]) {
  for (int i = 0; i < size; i++) {
    A[0] = T(0.0);
    A++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecCopyCore(const T A[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] = A[0];
    C++, A++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecScaleCore(const T alpha, const T A[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] = alpha * A[0];
    C++, A++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecAddCore(const T A[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] += A[0];
    C++, A++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecAddCore(const T alpha, const T A[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] += alpha * A[0];
    C++, A++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION T VecDotCore(const T A[], const T B[]) {
  T dot = 0.0;
  for (int i = 0; i < size; i++) {
    dot += A[0] * B[0];
    A++, B++;
  }
  return dot;
}

template <typename T, int size>
KOKKOS_FUNCTION void VecSumCore(const T A[], const T B[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] = A[0] + B[0];
    C++, A++, B++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecSumCore(const T alpha, const T A[], const T beta,
                                const T B[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] = alpha * A[0] + beta * B[0];
    C++, A++, B++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecSumCoreAdd(const T A[], const T B[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] += A[0] + B[0];
    C++, A++, B++;
  }
}

template <typename T, int size>
KOKKOS_FUNCTION void VecSumCoreAdd(const T alpha, const T A[], const T beta,
                                   const T B[], T C[]) {
  for (int i = 0; i < size; i++) {
    C[0] += alpha * A[0] + beta * B[0];
    C++, A++, B++;
  }
}

}  // namespace A2D

#endif  // A2D_VEC_CORE_H