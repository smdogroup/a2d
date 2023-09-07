#ifndef A2D_OBJECTS_H
#define A2D_OBJECTS_H

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dvec.h"

namespace A2D {

/*
  Remove the const-ness and references for a type
*/
template <class T>
struct remove_const_and_refs
    : std::remove_const<typename std::remove_reference<T>::type> {};

/*
  Check if a type is numeric or not
*/
template <class T>
struct __is_numeric_type : std::is_floating_point<T> {};

template <class T>
struct __is_numeric_type<std::complex<T>> : std::is_floating_point<T> {};

template <class T>
struct is_numeric_type
    : __is_numeric_type<typename remove_const_and_refs<T>::type> {};

/*
  Get the numeric type of the object
*/
template <class T>
struct __get_object_numeric_type {
  using type = typename T::type;
};

template <>
struct __get_object_numeric_type<double> {
  using type = double;
};

template <>
struct __get_object_numeric_type<std::complex<double>> {
  using type = std::complex<double>;
};

/*
  Get the numeric type of the underlying object.

  All A2D numeric objects must either be scalar values (float, double,
  std::complex) or must use a typedef statement to define the static scalar
  type.
*/
template <class T>
struct get_object_numeric_type
    : __get_object_numeric_type<typename remove_const_and_refs<T>::type> {};

/*
  Get the type of object
*/
template <class T>
struct __get_a2d_object_type {
  static constexpr ADObjType value = T::obj_type;
};

template <>
struct __get_a2d_object_type<double> {
  static constexpr ADObjType value = ADObjType::SCALAR;
};

template <>
struct __get_a2d_object_type<std::complex<double>> {
  static constexpr ADObjType value = ADObjType::SCALAR;
};

template <class T>
struct get_a2d_object_type
    : __get_a2d_object_type<typename std::remove_reference<T>::type> {};

/*
  The base first-order AD class
*/
template <class ObjType>
class ADObj {
 public:
  typedef typename get_object_numeric_type<ObjType>::type type;
  static constexpr ADObjType obj_type = get_a2d_object_type<ObjType>::value;

  KOKKOS_FUNCTION ADObj() {}
  KOKKOS_FUNCTION ADObj(const ObjType& A) : A(A) {}
  KOKKOS_FUNCTION ADObj(const ObjType& A, const ObjType& Ab) : A(A), Ab(Ab) {}
  KOKKOS_FUNCTION ObjType& value() { return A; }
  KOKKOS_FUNCTION const ObjType& value() const { return A; }
  KOKKOS_FUNCTION ObjType& bvalue() { return Ab; }
  KOKKOS_FUNCTION const ObjType& bvalue() const { return Ab; }

 private:
  ObjType A;   // Object
  ObjType Ab;  // Reverse mode derivative value
};

/*
  The base second-order AD class
*/
template <class ObjType>
class A2DObj {
 public:
  typedef typename get_object_numeric_type<ObjType>::type type;
  static constexpr ADObjType obj_type = get_a2d_object_type<ObjType>::value;

  KOKKOS_FUNCTION A2DObj() {}
  KOKKOS_FUNCTION A2DObj(const ObjType& A) : A(A) {}
  KOKKOS_FUNCTION A2DObj(const ObjType& A, const ObjType& Ab) : A(A), Ab(Ab) {}
  KOKKOS_FUNCTION A2DObj(const ObjType& A, const ObjType& Ab, const ObjType& Ap)
      : A(A), Ab(Ab), Ap(Ap) {}
  KOKKOS_FUNCTION A2DObj(const ObjType& A, const ObjType& Ab, const ObjType& Ap,
                         const ObjType& Ah)
      : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  KOKKOS_FUNCTION ObjType& value() { return A; }
  KOKKOS_FUNCTION const ObjType& value() const { return A; }
  KOKKOS_FUNCTION ObjType& bvalue() { return Ab; }
  KOKKOS_FUNCTION const ObjType& bvalue() const { return Ab; }
  KOKKOS_FUNCTION ObjType& pvalue() { return Ap; }
  KOKKOS_FUNCTION const ObjType& pvalue() const { return Ap; }
  KOKKOS_FUNCTION ObjType& hvalue() { return Ah; }
  KOKKOS_FUNCTION const ObjType& hvalue() const { return Ah; }

 private:
  ObjType A;   // Object
  ObjType Ab;  // Reverse mode derivative value
  ObjType Ap;  // Projected second derivative value
  ObjType Ah;  // Reverse mode second derivative
};

/*
  Remove the A2D template wrapper remove_a2dobj<T>::type provides the base type
*/
template <class T>
struct __remove_a2dobj {
  typedef T type;
};

template <class T>
struct __remove_a2dobj<ADObj<T>> {
  typedef typename std::remove_reference<T>::type type;
};

template <class T>
struct __remove_a2dobj<A2DObj<T>> {
  typedef typename std::remove_reference<T>::type type;
};

template <class T>
struct remove_a2dobj
    : __remove_a2dobj<typename remove_const_and_refs<T>::type> {};

/*
  Get whether the type is passive or not...
*/
template <class T>
struct get_diff_type {
  static constexpr ADiffType diff_type = ADiffType::PASSIVE;
};

template <class T>
struct get_diff_type<ADObj<T>> {
  static constexpr ADiffType diff_type = ADiffType::ACTIVE;
};

template <class T>
struct get_diff_type<A2DObj<T>> {
  static constexpr ADiffType diff_type = ADiffType::ACTIVE;
};

template <class T>
struct get_diff_type<ADObj<T>&> {
  static constexpr ADiffType diff_type = ADiffType::ACTIVE;
};

template <class T>
struct get_diff_type<A2DObj<T>&> {
  static constexpr ADiffType diff_type = ADiffType::ACTIVE;
};

/*
  Get whether the type is passive or not...
*/
template <class T>
struct get_diff_order {
  static constexpr ADorder order = ADorder::ZERO;
};

template <class T>
struct get_diff_order<ADObj<T>> {
  static constexpr ADorder order = ADorder::FIRST;
};

template <class T>
struct get_diff_order<A2DObj<T>> {
  static constexpr ADorder order = ADorder::SECOND;
};

/*
  Get the vector size
*/
template <class T>
struct __get_vec_size {
  static constexpr int size = 0;
};

// This will pick up vectors too...
template <template <typename, int> class Vec, typename T, int N>
struct __get_vec_size<Vec<T, N>> {
  static constexpr int size = N;
};

template <class T>
struct get_vec_size : __get_vec_size<typename remove_a2dobj<T>::type> {
  static_assert(get_a2d_object_type<T>::value == ADObjType::VECTOR,
                "get_symmatrix_size called on incorrect type");
};

/*
  Get the symmetric matrix size
*/
template <class T>
struct __get_symmatrix_size {
  static constexpr int size = 0;
};

// This will pick up vectors too...
template <template <typename, int> class SymMat, typename T, int N>
struct __get_symmatrix_size<SymMat<T, N>> {
  static constexpr int size = N;
};

template <class T>
struct get_symmatrix_size
    : __get_symmatrix_size<typename remove_a2dobj<T>::type> {
  static_assert(get_a2d_object_type<T>::value == ADObjType::SYMMAT,
                "get_symmatrix_size called on incorrect type");
};

/*
  Get the number of matrix rows
*/
template <class T>
struct __get_matrix_rows {
  static constexpr int size = 0;
};

template <template <typename, int, int> class Mat, typename T, int N, int M>
struct __get_matrix_rows<Mat<T, N, M>> {
  static constexpr int size = N;
};

template <class T>
struct get_matrix_rows : __get_matrix_rows<typename remove_a2dobj<T>::type> {
  static_assert(get_a2d_object_type<T>::value == ADObjType::MATRIX,
                "get_matrix_rows called on incorrect type");
};

/*
  Get the number of matrix columns
*/
template <class T>
struct __get_matrix_columns {
  static constexpr int size = 0;
};

template <template <typename, int, int> class Mat, typename T, int N, int M>
struct __get_matrix_columns<Mat<T, N, M>> {
  static constexpr int size = M;
};

template <class T>
struct get_matrix_columns
    : __get_matrix_columns<typename remove_a2dobj<T>::type> {
  static_assert(get_a2d_object_type<T>::value == ADObjType::MATRIX,
                "get_matrix_rows called on incorrect type");
};

/*
  Get the number of entries in the matrix -- symmetric or non-symmetric
*/
template <class T>
struct __get_num_matrix_entries {
  static constexpr int size = 0;
};

template <template <typename, int> class SymMat, typename T, int N>
struct __get_num_matrix_entries<SymMat<T, N>> {
  static constexpr int size = N * (N + 1) / 2;
};

template <template <typename, int, int> class Mat, typename T, int N, int M>
struct __get_num_matrix_entries<Mat<T, N, M>> {
  static constexpr int size = N * M;
};

template <class T>
struct get_num_matrix_entries
    : __get_num_matrix_entries<typename remove_a2dobj<T>::type> {
  static_assert((get_a2d_object_type<T>::value == ADObjType::MATRIX ||
                 get_a2d_object_type<T>::value == ADObjType::SYMMAT),
                "get_num_matrix_entries called on incorrect type");
};

/**
 * @brief Select type based on supplied ADiffType value
 *
 * For example, the following types are equivalent:
 *
 * Obj<...>         ==
 *    ADObjSelect<ADiffType::PASSIVE, *,               Obj<...>>;
 * ADObj<Obj<...>>  ==
 *    ADObjSelect<ADiffType::ACTIVE,  ADorder::FIRST,  Obj<...>>;
 * A2DObj<Obj<...>> ==
 *    ADObjSelect<ADiffType::ACTIVE,  ADorder::SECOND, Obj<...>>;
 *
 * @tparam adiff_type passive or active
 * @tparam order first (AD) or second (A2D)
 * @tparam Obj the numeric type of the matrix
 */
template <ADiffType adiff_type, ADorder order, class Obj>
using ADObjSelect = typename std::conditional<
    adiff_type == ADiffType::ACTIVE,
    typename std::conditional<order == ADorder::FIRST, ADObj<Obj>,
                              A2DObj<Obj>>::type,
    Obj>::type;

/**
 * @brief Get objects and pointers to seed data (bvalue(), pvalue(), hvalue)
 */
template <ADseed seed>
class GetSeed {
 public:
  template <typename T>
  static KOKKOS_FUNCTION T& get_obj(ADObj<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return value.bvalue();
  }

  template <typename T>
  static KOKKOS_FUNCTION T& get_obj(A2DObj<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return value.bvalue();
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue();
    } else {  // seed == ADseed::h
      return value.hvalue();
    }
  }

  template <typename T,
            std::enable_if_t<is_numeric_type<T>::value, bool> = true>
  static KOKKOS_FUNCTION T& get_data(ADObj<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return value.bvalue();
  }

  template <typename T,
            std::enable_if_t<is_numeric_type<T>::value, bool> = true>
  static KOKKOS_FUNCTION T& get_data(A2DObj<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return value.bvalue();
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue();
    } else {  // seed == ADseed::h
      return value.hvalue();
    }
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(ADObj<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return value.bvalue().get_data();
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(A2DObj<Vec<T, N>>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return value.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue().get_data();
    } else {  // seed == ADseed::h
      return value.hvalue().get_data();
    }
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(ADObj<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return mat.bvalue().get_data();
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(A2DObj<Mat<T, m, n>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return mat.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return mat.pvalue().get_data();
    } else {  // seed == ADseed::h
      return mat.hvalue().get_data();
    }
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(ADObj<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return mat.bvalue().get_data();
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(A2DObj<SymMat<T, m>>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return mat.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return mat.pvalue().get_data();
    } else {  // seed == ADseed::h
      return mat.hvalue().get_data();
    }
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(ADObj<Vec<T, N>&>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return value.bvalue().get_data();
  }

  template <typename T, int N>
  static KOKKOS_FUNCTION T* get_data(A2DObj<Vec<T, N>&>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return value.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue().get_data();
    } else {  // seed == ADseed::h
      return value.hvalue().get_data();
    }
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(ADObj<Mat<T, m, n>&>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return mat.bvalue().get_data();
  }

  template <typename T, int m, int n>
  static KOKKOS_FUNCTION T* get_data(A2DObj<Mat<T, m, n>&>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return mat.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return mat.pvalue().get_data();
    } else {  // seed == ADseed::h
      return mat.hvalue().get_data();
    }
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(ADObj<SymMat<T, m>&>& mat) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return mat.bvalue().get_data();
  }

  template <typename T, int m>
  static KOKKOS_FUNCTION T* get_data(A2DObj<SymMat<T, m>&>& mat) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return mat.bvalue().get_data();
    } else if constexpr (seed == ADseed::p) {
      return mat.pvalue().get_data();
    } else {  // seed == ADseed::h
      return mat.hvalue().get_data();
    }
  }
};

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
T& get_data(T& value) {
  return value;
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
const T& get_data(const T& value) {
  return value;
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
T& get_data(ADObj<T>& value) {
  return value.value();
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
const T& get_data(const ADObj<T>& value) {
  return value.value();
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
T& get_data(A2DObj<T>& value) {
  return value.value();
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
const T& get_data(const A2DObj<T>& value) {
  return value.value();
}

/**
 * @brief Get data pointers from objects
 */
template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(Mat<T, m, n>& mat) {
  return mat.get_data();
}

template <typename T, int m, int n>
KOKKOS_FUNCTION const T* get_data(const Mat<T, m, n>& mat) {
  return mat.get_data();
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(ADObj<Mat<T, m, n>>& mat) {
  return mat.value().get_data();
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(A2DObj<Mat<T, m, n>>& mat) {
  return mat.value().get_data();
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(ADObj<Mat<T, m, n>&>& mat) {
  return mat.value().get_data();
}

template <typename T, int m, int n>
KOKKOS_FUNCTION T* get_data(A2DObj<Mat<T, m, n>&>& mat) {
  return mat.value().get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(SymMat<T, m>& mat) {
  return mat.get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION const T* get_data(const SymMat<T, m>& mat) {
  return mat.get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(ADObj<SymMat<T, m>>& mat) {
  return mat.value().get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(A2DObj<SymMat<T, m>>& mat) {
  return mat.value().get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(ADObj<SymMat<T, m>&>& mat) {
  return mat.value().get_data();
}

template <typename T, int m>
KOKKOS_FUNCTION T* get_data(A2DObj<SymMat<T, m>&>& mat) {
  return mat.value().get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(Vec<T, n>& vec) {
  return vec.get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION const T* get_data(const Vec<T, n>& vec) {
  return vec.get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(ADObj<Vec<T, n>>& vec) {
  return vec.value().get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(A2DObj<Vec<T, n>>& vec) {
  return vec.value().get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(ADObj<Vec<T, n>&>& vec) {
  return vec.value().get_data();
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(A2DObj<Vec<T, n>&>& vec) {
  return vec.value().get_data();
}

}  // namespace A2D

#endif  // A2D_OBJECTS_H