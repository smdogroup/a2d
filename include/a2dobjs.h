#ifndef A2D_OBJS_H
#define A2D_OBJS_H

#include <cmath>
#include <cstdint>

#include "Kokkos_Core.hpp"
#include "Kokkos_UnorderedMap.hpp"

#define A2D_LAMBDA KOKKOS_LAMBDA
#define A2D_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION

#ifdef KOKKOS_ENABLE_CUDA
#include "cuda/std/complex"
#include "thrust/fill.h"
template <typename T>
using A2D_complex_t = cuda::std::complex<T>;
#else
#include <algorithm>
#include <complex>
template <typename T>
using A2D_complex_t = std::complex<T>;
#endif

namespace A2D {
using index_t = uint32_t;  // TODO: size_t may be a better choice here
static constexpr index_t MAX_INDEX = std::numeric_limits<index_t>::max();
static constexpr index_t NO_INDEX = MAX_INDEX;

/**
 * @brief Check if a type is complex.
 *
 * Usage:
 *   given an arbitrary type T, if T is a complex type, then
 *     is_complex<T>::value == true, otherwise
 *     is_complex<T>::value == false
 */
template <typename T>
struct is_complex : public std::false_type {};

template <typename T>
struct is_complex<A2D_complex_t<T>> : public std::true_type {};

/*
 Convert scalar value to printf-able format
*/
template <typename T>
A2D_INLINE_FUNCTION double fmt(A2D_complex_t<T> val) {
  return val.real();
}

A2D_INLINE_FUNCTION double fmt(double val) { return val; }

#ifdef KOKKOS_ENABLE_CUDA
template <typename T>
A2D_INLINE_FUNCTION T sqrt(T val) {
  return cuda::std::sqrt(val);
}

template <typename T>
A2D_INLINE_FUNCTION T exp(T val) {
  return cuda::std::exp(val);
}

template <typename T>
A2D_INLINE_FUNCTION T log(T val) {
  return cuda::std::log(val);
}

template <class ForwardIt, class T>
A2D_INLINE_FUNCTION void fill(ForwardIt first, ForwardIt last, const T& value) {
  thrust::fill(first, last, value);
}
#else
template <typename T>
A2D_INLINE_FUNCTION T sqrt(T val) {
  return std::sqrt(val);
}

template <typename T>
A2D_INLINE_FUNCTION T exp(T val) {
  return std::exp(val);
}

template <typename T>
A2D_INLINE_FUNCTION T log(T val) {
  return std::log(val);
}
template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T& value) {
  std::fill(first, last, value);
}
#endif

template <typename T, int N>
class Vec {
 public:
  typedef T type;

  A2D_INLINE_FUNCTION Vec() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  A2D_INLINE_FUNCTION Vec(const T* vals) {
    for (int i = 0; i < N; i++) {
      x[i] = vals[i];
    }
  }
  template <class VecType>
  A2D_INLINE_FUNCTION Vec(const VecType& vec) {
    for (int i = 0; i < N; i++) {
      x[i] = vec(i);
    }
  }
  A2D_INLINE_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION T& operator()(const IdxType i) {
    return x[i];
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType i) const {
    return x[i];
  }

  T* data() { return x; }

  T x[N];
};

template <typename T, int M, int N>
class Mat {
 public:
  typedef T type;

  A2D_INLINE_FUNCTION Mat() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  A2D_INLINE_FUNCTION Mat(const T* vals) {
    for (int i = 0; i < M * N; i++) {
      A[i] = vals[i];
    }
  }
  template <class MatType>
  A2D_INLINE_FUNCTION Mat(const MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] = mat(i, j);
      }
    }
  }
  A2D_INLINE_FUNCTION void zero() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  template <class MatType>
  A2D_INLINE_FUNCTION void set(const MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] = mat(i, j);
      }
    }
  }
  template <class MatType>
  A2D_INLINE_FUNCTION void get(MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        mat(i, j) = A[N * i + j];
      }
    }
  }
  template <class IdxType1, class IdxType2>
  A2D_INLINE_FUNCTION T& operator()(const IdxType1 i, const IdxType2 j) {
    return A[N * i + j];
  }
  template <class IdxType1, class IdxType2>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType1 i,
                                          const IdxType2 j) const {
    return A[N * i + j];
  }

  T* data() { return A; }

  T A[M * N];
};

template <typename T, int N>
class SymmMat {
 public:
  typedef T type;
  static const int MAT_SIZE = (N * (N + 1)) / 2;

  A2D_INLINE_FUNCTION SymmMat() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  A2D_INLINE_FUNCTION SymmMat(const T* vals) {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = vals[i];
    }
  }
  A2D_INLINE_FUNCTION void zero() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class SymmMat>
  A2D_INLINE_FUNCTION void set(const SymmMat& mat) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        A[i + j * (j + 1) / 2] = mat(i, j);
      }
    }
  }
  template <class SymmMat>
  A2D_INLINE_FUNCTION void get(SymmMat& mat) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        mat(i, j) = A[i + j * (j + 1) / 2];
      }
    }
  }

  template <class IdxType1, class IdxType2>
  A2D_INLINE_FUNCTION T& operator()(const IdxType1 i, const IdxType2 j) {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }
  template <class IdxType1, class IdxType2>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType1 i,
                                          const IdxType2 j) const {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }

  T* data() { return A; }

  T A[MAT_SIZE];
};

/*
  Full 4-th order tensor without any symmetry.

  This tensor class is required for mixed second order derivatives of
  two matrices X and Y such that

  A(i, j, k, l) = d^2 f/dX(i, j) dY(k, l)
*/
template <typename T, int M, int N, int P, int Q>
class Tensor {
 public:
  typedef T type;
  static const int TENSOR_SIZE = M * N * P * Q;
  A2D_INLINE_FUNCTION Tensor() {
    for (int i = 0; i < TENSOR_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION T& operator()(const IdxType i, const IdxType j,
                                    const IdxType k, const IdxType l) {
    return A[l + Q * (k + P * (j + N * i))];
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType i, const IdxType j,
                                          const IdxType k,
                                          const IdxType l) const {
    return A[l + Q * (k + P * (j + N * i))];
  }

  T* data() { return A; }

  T A[TENSOR_SIZE];
};

/*
  Basic 4-th order symmetric tensor.

  This class stores a symmetric tensor found by taking the second order
  derivatives of a scalar function f with respect to a non-symmetric M-by-N
  matrix:

  A(i, j, k, l) = d^2 f/dX(i, j) dX(k, l)

  As a result:

  A(i, j, k, l) = A(k, l, i, j).
*/
template <typename T, int M, int N>
class SymmTensor {
 public:
  typedef T type;
  static const int TENSOR_SIZE = (M * N * (M * N + 1)) / 2;
  A2D_INLINE_FUNCTION SymmTensor() {
    for (int i = 0; i < TENSOR_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION T& operator()(const IdxType i, const IdxType j,
                                    const IdxType k, const IdxType l) {
    const int ii = N * i + j;
    const int jj = N * k + l;

    if (ii >= jj) {
      return A[jj + ii * (ii + 1) / 2];
    } else {
      return A[ii + jj * (jj + 1) / 2];
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType i, const IdxType j,
                                          const IdxType k,
                                          const IdxType l) const {
    const int ii = N * i + j;
    const int jj = N * k + l;

    if (ii >= jj) {
      return A[jj + ii * (ii + 1) / 2];
    } else {
      return A[ii + jj * (jj + 1) / 2];
    }
  }

  T* data() { return A; }

  T A[TENSOR_SIZE];
};

}  // namespace A2D

#endif  // A2D_OBJS_H
