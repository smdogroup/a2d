#ifndef A2D_MAT_H
#define A2D_MAT_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

enum class MatSymType { NORMAL, SYMMETRIC };

template <typename T, int M, int N>
class Mat {
 public:
  typedef T type;

  static const index_t ncomp = M * N;
  static const int nrows = M;
  static const int ncols = N;

  KOKKOS_FUNCTION Mat() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  KOKKOS_FUNCTION Mat(const T* vals) {
    for (int i = 0; i < M * N; i++) {
      A[i] = vals[i];
    }
  }
  template <class MatType>
  KOKKOS_FUNCTION Mat(const MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] = mat(i, j);
      }
    }
  }
  KOKKOS_FUNCTION void zero() {
    for (int i = 0; i < M * N; i++) {
      A[i] = 0.0;
    }
  }
  template <class MatType>
  KOKKOS_FUNCTION void set(const MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        A[N * i + j] = mat(i, j);
      }
    }
  }
  template <class MatType>
  KOKKOS_FUNCTION void get(MatType& mat) {
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
        mat(i, j) = A[N * i + j];
      }
    }
  }
  template <class IdxType1, class IdxType2>
  KOKKOS_FUNCTION T& operator()(const IdxType1 i, const IdxType2 j) {
    return A[N * i + j];
  }
  template <class IdxType1, class IdxType2>
  KOKKOS_FUNCTION const T& operator()(const IdxType1 i,
                                      const IdxType2 j) const {
    return A[N * i + j];
  }

  T* get_data() { return A; }
  const T* get_data() const { return A; }

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return A[i];
  }

  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return A[i];
  }

 private:
  T A[M * N];
};

template <typename T, int N>
class SymMat {
 public:
  typedef T type;
  static const int MAT_SIZE = (N * (N + 1)) / 2;
  static constexpr int nrows = N;
  static constexpr int ncols = N;

  static const index_t ncomp = MAT_SIZE;

  KOKKOS_FUNCTION SymMat() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  KOKKOS_FUNCTION SymMat(const T* vals) {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = vals[i];
    }
  }
  KOKKOS_FUNCTION void zero() {
    for (int i = 0; i < MAT_SIZE; i++) {
      A[i] = 0.0;
    }
  }
  template <class SymMat>
  KOKKOS_FUNCTION void set(const SymMat& mat) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        A[i + j * (j + 1) / 2] = mat(i, j);
      }
    }
  }
  template <class SymMat>
  KOKKOS_FUNCTION void get(SymMat& mat) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= i; j++) {
        mat(i, j) = A[i + j * (j + 1) / 2];
      }
    }
  }

  template <class IdxType1, class IdxType2>
  KOKKOS_FUNCTION T& operator()(const IdxType1 i, const IdxType2 j) {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }
  template <class IdxType1, class IdxType2>
  KOKKOS_FUNCTION const T& operator()(const IdxType1 i,
                                      const IdxType2 j) const {
    if (i >= j) {
      return A[j + i * (i + 1) / 2];
    } else {
      return A[i + j * (j + 1) / 2];
    }
  }

  T* get_data() { return A; }
  const T* get_data() const { return A; }

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return A[i];
  }
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return A[i];
  }

 private:
  T A[MAT_SIZE];
};

}  // namespace A2D

#endif  // A2D_MAT_H