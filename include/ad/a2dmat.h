#ifndef A2D_MAT_H
#define A2D_MAT_H

#include "a2denum.h"
#include "a2dobjs.h"
#include "a2dscalar.h"
#include "a2dvec.h"

namespace A2D {

enum class MatSymType { NORMAL, SYMMETRIC };

template <typename T, int M, int N>
class Mat {
 public:
  typedef T type;

  static const index_t num_components = M * N;
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

  T* data() { return A; }

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return A[i];
  }
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return A[i];
  }

  T A[M * N];
};

template <typename T, int N>
class SymMat {
 public:
  typedef T type;
  static const int MAT_SIZE = (N * (N + 1)) / 2;
  static constexpr int nrows = N;
  static constexpr int ncols = N;

  static const index_t num_components = MAT_SIZE;

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

  T* data() { return A; }

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return A[i];
  }
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return A[i];
  }

  T A[MAT_SIZE];
};

template <class MatType>
class ADMat {
 public:
  KOKKOS_FUNCTION ADMat(MatType& A, MatType& Ab) : A(A), Ab(Ab) {}

  KOKKOS_FUNCTION MatType& value() { return A; }
  KOKKOS_FUNCTION const MatType& value() const { return A; }

  KOKKOS_FUNCTION MatType& bvalue() { return Ab; }
  KOKKOS_FUNCTION const MatType& bvalue() const { return Ab; }

  MatType& A;   // Matrix
  MatType& Ab;  // Reverse mode derivative value
};

template <class MatType>
class A2DMat {
 public:
  KOKKOS_FUNCTION A2DMat() {}
  KOKKOS_FUNCTION A2DMat(const MatType& A) : A(A) {}
  KOKKOS_FUNCTION A2DMat(const MatType& A, const MatType& Ab) : A(A), Ab(Ab) {}
  KOKKOS_FUNCTION A2DMat(const MatType& A, const MatType& Ab, const MatType& Ap)
      : A(A), Ab(Ab), Ap(Ap) {}
  KOKKOS_FUNCTION A2DMat(const MatType& A, const MatType& Ab, const MatType& Ap,
                         const MatType& Ah)
      : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  KOKKOS_FUNCTION MatType& value() { return A; }
  KOKKOS_FUNCTION const MatType& value() const { return A; }

  KOKKOS_FUNCTION void set_bvalue(const MatType& val) { Ab.set(val); }
  KOKKOS_FUNCTION void get_bvalue(MatType& val) { Ab.get(val); }
  KOKKOS_FUNCTION MatType& bvalue() { return Ab; }
  KOKKOS_FUNCTION const MatType& bvalue() const { return Ab; }

  KOKKOS_FUNCTION void set_pvalue(const MatType& val) { Ap.set(val); }
  KOKKOS_FUNCTION void get_pvalue(MatType& val) { Ap.get(val); }
  KOKKOS_FUNCTION MatType& pvalue() { return Ap; }
  KOKKOS_FUNCTION const MatType& pvalue() const { return Ap; }

  KOKKOS_FUNCTION void set_hvalue(const MatType& val) { Ah.set(val); }
  KOKKOS_FUNCTION void get_hvalue(MatType& val) { Ah.get(val); }
  KOKKOS_FUNCTION MatType& hvalue() { return Ah; }
  KOKKOS_FUNCTION const MatType& hvalue() const { return Ah; }

  MatType A;   // Matrix
  MatType Ab;  // Reverse mode derivative value
  MatType Ap;  // Projected second derivative value
  MatType Ah;  // Reverse mode second derivative
};

/**
 * @brief Get data pointers from Mat/ADMat/A2DMat objects
 */
template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(Mat<T, m, n>& mat) {
  return mat.A;
}

template <typename T, int m, int n>
KOKKOS_FUNCTION const T* get_data(const Mat<T, m, n>& mat) {
  return mat.A;
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(ADMat<Mat<T, m, n>>& mat) {
  return mat.A.A;
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(A2DMat<Mat<T, m, n>>& mat) {
  return mat.A.A;
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(SymMat<T, m>& mat) {
  return mat.A;
}

template <typename T, int m>
KOKKOS_FUNCTION const T* get_data(const SymMat<T, m>& mat) {
  return mat.A;
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(ADMat<SymMat<T, m>>& mat) {
  return mat.A.A;
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(A2DMat<SymMat<T, m>>& mat) {
  return mat.A.A;
}

/**
 * @brief Get pointer to seed data (bvalue, pvalue, hvalue) from
 * Mat/ADMat/A2DMat objects
 */
template <ADseed seed>
class GetSeed {
 public:
  template <typename T>
  static KOKKOS_FUNCTION T& get_data(ADScalar<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADScalar");
    return value.bvalue;
  }

  template <typename T>
  static KOKKOS_FUNCTION T& get_data(A2DScalar<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DScalar");
    if constexpr (seed == ADseed::b) {
      return value.bvalue;
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue;
    } else {  // seed == ADseed::h
      return value.hvalue;
    }
  }

  template <typename T>
  static KOKKOS_FUNCTION T& get_obj(ADScalar<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADScalar");
    return value.bvalue;
  }

  template <typename T>
  static KOKKOS_FUNCTION T& get_obj(A2DScalar<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DScalar");
    if constexpr (seed == ADseed::b) {
      return value.bvalue;
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue;
    } else {  // seed == ADseed::h
      return value.hvalue;
    }
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(ADVec<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADScalar");
    return value.Vb.V;
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(A2DVec<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DScalar");
    if constexpr (seed == ADseed::b) {
      return value.Vb.V;
    } else if constexpr (seed == ADseed::p) {
      return value.Vp.V;
    } else {  // seed == ADseed::h
      return value.Vh.V;
    }
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION Vec<T, N>& get_obj(ADVec<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADScalar");
    return value.Vb;
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION Vec<T, N>& get_obj(A2DVec<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DScalar");
    if constexpr (seed == ADseed::b) {
      return value.Vb;
    } else if constexpr (seed == ADseed::p) {
      return value.Vp;
    } else {  // seed == ADseed::h
      return value.Vh;
    }
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(ADMat<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADMat");
    return mat.Ab.A;
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(A2DMat<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DMat");
    if constexpr (seed == ADseed::b) {
      return mat.Ab.A;
    } else if constexpr (seed == ADseed::p) {
      return mat.Ap.A;
    } else {  // seed == ADseed::h
      return mat.Ah.A;
    }
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION Mat<T, m, n>& get_obj(ADMat<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADMat");
    return mat.Ab;
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION Mat<T, m, n>& get_obj(A2DMat<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DMat");
    if constexpr (seed == ADseed::b) {
      return mat.Ab;
    } else if constexpr (seed == ADseed::p) {
      return mat.Ap;
    } else {  // seed == ADseed::h
      return mat.Ah;
    }
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(ADMat<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADMat");
    return mat.Ab.A;
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(A2DMat<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DMat");
    if constexpr (seed == ADseed::b) {
      return mat.Ab.A;
    } else if constexpr (seed == ADseed::p) {
      return mat.Ap.A;
    } else {  // seed == ADseed::h
      return mat.Ah.A;
    }
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION SymMat<T, m>& get_obj(ADMat<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADMat");
    return mat.Ab;
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION SymMat<T, m>& get_obj(A2DMat<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DMat");
    if constexpr (seed == ADseed::b) {
      return mat.Ab;
    } else if constexpr (seed == ADseed::p) {
      return mat.Ap;
    } else {  // seed == ADseed::h
      return mat.Ah;
    }
  }
};

/**
 * @brief Select type based on whether the matrix is passive or active (can be
 * differentiated)
 *
 * For example, the following types are equivalent:
 *
 * Mat<...>    == ADMatType<ADiffType::PASSIVE, *,               Mat<...>>;
 * ADMat<...>  == ADMatType<ADiffType::ACTIVE,  ADorder::FIRST,  Mat<...>>;
 * A2DMat<...> == ADMatType<ADiffType::ACTIVE,  ADorder::SECOND, Mat<...>>;
 *
 * @tparam adiff_type passive or active
 * @tparam order first (AD) or second (A2D)
 * @tparam MatType the numeric type of the matrix
 */
template <ADiffType adiff_type, ADorder order, class MatType>
using ADMatType = typename std::conditional<
    adiff_type == ADiffType::ACTIVE,
    typename std::conditional<order == ADorder::FIRST, ADMat<MatType>,
                              A2DMat<MatType>>::type,
    const MatType>::type;
}  // namespace A2D

#endif  // A2D_MAT_H