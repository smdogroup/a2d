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
static constexpr index_t INDEX_NBITS = std::numeric_limits<index_t>::digits;
static constexpr index_t NO_INDEX = MAX_INDEX;

// Free heap memory and set pointer to nullptr
#define DELETE_ARRAY(array_ptr) \
  delete[] array_ptr;           \
  array_ptr = nullptr;

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

double absfunc(A2D_complex_t<double> a) {
  if (a.real() >= 0.0) {
    return a.real();
  } else {
    return -a.real();
  }
}

double absfunc(double a) {
  if (a >= 0.0) {
    return a;
  } else {
    return -a;
  }
}

double RealPart(double a) { return a; }

double RealPart(A2D_complex_t<double> a) { return a.real(); }

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
T* get_data(Mat<T, M, N>& A) {
  return A.A;
}

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

template <typename T, int M, int N>
T* get_data(Mat<T, M, N>& A) {
  return A.A;
}

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

template <typename T, int M>
T* get_data(SymmMat<T, M>& A) {
  return A.A;
}

}  // namespace A2D

#endif  // A2D_OBJS_H
