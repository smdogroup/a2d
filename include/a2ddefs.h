#ifndef A2D_DEFS_H
#define A2D_DEFS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>

#ifndef __CUDACC__
template <typename T>
using A2D_complex_t = std::complex<T>;
#else
#include <thrust/fill.h>

#include <cuda/std/complex>
template <typename T>
using A2D_complex_t = cuda::std::complex<T>;
#endif

// CUDA headers
#ifndef A2D_FUNCTION
#ifdef __CUDACC__
#define A2D_FUNCTION \
  __host__ __device__  // A2D_FUNCTION does nothing in this scenario
#else                  // not __CUDACC__
#define A2D_FUNCTION
#endif
#endif

namespace A2D {

/**
 * @brief Whether a variable should be automatically differentiated or not
 */
enum class ADiffType { PASSIVE, ACTIVE };

/**
 * @brief Whether we compute first or second order derivatives
 */
enum class ADorder { ZERO, FIRST, SECOND };

/**
 * @brief The class of object
 */
enum class ADObjType { SCALAR, VECTOR, MATRIX, SYMMAT };

/**
 * @brief Is the matrix normal (not transposed) or transposed
 */
enum class MatOp { NORMAL, TRANSPOSE };

/**
 * @brief The symmetry type of the matrix
 */
enum class MatSymType { NORMAL, SYMMETRIC };

/**
 * @brief Specify bvalue, pvalue or hvalue
 */
enum class ADseed { b, p, h };

// Forward declarations of classes
template <typename T, int N>
class ADScalar;

template <class T>
class ADObj;

template <class T>
class A2DObj;

/**
 * @brief Check if a type is complex.
 *
 * Usage:
 *   given an arbitrary type T, if T is a complex type, then
 *     is_complex<T>::value == true, otherwise
 *     is_complex<T>::value == false
 */
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<A2D_complex_t<T>> : std::true_type {};

/*
 Convert scalar value to printf-able format
*/
template <typename T>
A2D_FUNCTION double fmt(A2D_complex_t<T> val) {
  return val.real();
}

A2D_FUNCTION inline double fmt(double val) { return val; }

A2D_FUNCTION inline double absfunc(A2D_complex_t<double> a) {
  if (a.real() >= 0.0) {
    return a.real();
  } else {
    return -a.real();
  }
}

A2D_FUNCTION inline double absfunc(double a) {
  if (a >= 0.0) {
    return a;
  } else {
    return -a;
  }
}

A2D_FUNCTION inline double RealPart(double a) { return a; }

A2D_FUNCTION inline double RealPart(A2D_complex_t<double> a) {
  return a.real();
}

A2D_FUNCTION inline double ImagPart(double a) { return 0.0; }

A2D_FUNCTION inline double ImagPart(A2D_complex_t<double> a) {
  return a.imag();
}

/**
 * @brief Remove the constant value and references
 */
template <typename T>
using remove_cvref_t = typename std::remove_cv_t<std::remove_reference_t<T>>;

/**
 * @brief Detect passive scalar types - e.g. float, double, complex
 */
template <typename T>
struct is_passive_scalar_impl : std::false_type {};

template <>
struct is_passive_scalar_impl<float> : std::true_type {};

template <>
struct is_passive_scalar_impl<double> : std::true_type {};

template <typename T>
struct is_passive_scalar_impl<A2D_complex_t<T>>
    : std::bool_constant<std::is_floating_point_v<T>> {};

template <typename T>
struct is_passive_scalar : is_passive_scalar_impl<remove_cvref_t<T>> {};

/**
 * @brief Detect scalar types, discarding const/ref qualifiers
 */
template <typename T>
inline constexpr bool is_passive_scalar_v = is_passive_scalar<T>::value;

/*
 * Detect passive scalars or integral types
 */
template <typename T>
struct is_passive_scalar_or_integral
    : std::bool_constant<is_passive_scalar<T>::value ||
                         std::is_integral<T>::value> {};

template <typename T>
inline constexpr bool is_passive_scalar_or_integral_v =
    is_passive_scalar_or_integral<T>::value;

/**
 * @brief Detect ADScalar types
 */
template <typename T>
struct is_adscalar_impl : std::false_type {};

template <typename T, int N>
struct is_adscalar_impl<ADScalar<T, N>> : std::true_type {};

template <typename T>
struct is_adscalar : is_adscalar_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr bool is_adscalar_v = is_adscalar<T>::value;

/**
 * @brief Detect scalar or ad types - numeric types
 */
template <typename T>
struct is_numeric_impl
    : std::bool_constant<is_passive_scalar<T>::value || is_adscalar<T>::value> {
};

template <typename T>
struct is_numeric : is_numeric_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr bool is_numeric_v = is_numeric<T>::value;

/**
 * Detect numeric or integral types
 */
template <typename T>
struct is_numeric_or_integral
    : std::bool_constant<is_numeric<T>::value || std::is_integral<T>::value> {};

template <typename T>
inline constexpr bool is_numeric_or_integral_v =
    is_numeric_or_integral<T>::value;

/**
 * @brief Get the scalar type used for the computation.
 *
 * This works for ADScalar, ADObj and complex types.
 */
template <typename T>
struct object_scalar_impl {
  using type = typename T::type;
};

template <>
struct object_scalar_impl<float> {
  using type = float;
};

template <>
struct object_scalar_impl<double> {
  using type = double;
};

template <>
struct object_scalar_impl<A2D_complex_t<float>> {
  using type = A2D_complex_t<float>;
};

template <>
struct object_scalar_impl<A2D_complex_t<double>> {
  using type = A2D_complex_t<double>;
};

template <typename T>
struct object_scalar : object_scalar_impl<remove_cvref_t<T>> {};

/**
 * @brief Get the underlying scalar computational type, e.g. double, float,
 * complex or ADScalar
 */
template <typename T>
using object_scalar_t = typename object_scalar<T>::type;

/*
  Get the type of object. Note that ADObj and ADScalar types have obj_type.
*/
template <class T>
struct a2d_object_type_impl {
  static constexpr ADObjType value = T::obj_type;
};

template <>
struct a2d_object_type_impl<double> {
  static constexpr ADObjType value = ADObjType::SCALAR;
};

template <>
struct a2d_object_type_impl<A2D_complex_t<double>> {
  static constexpr ADObjType value = ADObjType::SCALAR;
};

template <class T, int N>
struct a2d_object_type_impl<ADScalar<T, N>> {
  static constexpr ADObjType value = ADObjType::SCALAR;
};

template <class T>
struct a2d_obj_type : a2d_object_type_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr ADObjType a2d_obj_type_v = a2d_obj_type<T>::value;

/**
 * @brief compile-time conditional for non-type values, works like
 * std::conditional, but user gets non-type values instead.
 *
 * @tparam T type of the non-type value
 * @tparam B condition, true or false
 * @tparam TrueVal the value user gets if B is true
 * @tparam FalseVal the value user gets if B is false
 */
template <class T, bool B, T TrueVal, T FalseVal>
struct conditional_value {
  static constexpr T value = FalseVal;
};

template <class T, T TrueVal, T FalseVal>
struct conditional_value<T, true, TrueVal, FalseVal> {
  static constexpr T value = TrueVal;
};

/**
 * @brief Check if two types are the same when they're unqualified by
 * const/references
 */
template <class Ta, class Tb>
struct is_same_unqualified {
  const static bool value =
      std::is_same<remove_cvref_t<Ta>, remove_cvref_t<Tb>>::value;
};

template <typename Ta, typename Tb>
inline constexpr bool is_same_unqualified_v =
    is_same_unqualified<Ta, Tb>::value;

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

/**
 * Forward declaration of vector and matrix types
 */

template <typename T, int N>
class Vec;

template <typename T, int N>
class VecRef;

template <typename T, int M, int N>
class Mat;

template <typename T, int M, int N>
class MatRef;

template <typename T, int N>
class SymMat;

template <typename T, int N>
class SymMatRef;

/*
  Remove the A2D template wrapper remove_a2dobj<T>::type provides the base type
*/
template <typename T>
struct remove_a2dobj_impl {
  typedef T type;
};

template <typename T>
struct remove_a2dobj_impl<ADObj<T>> {
  typedef typename std::remove_reference<T>::type type;
};

template <typename T>
struct remove_a2dobj_impl<A2DObj<T>> {
  typedef typename std::remove_reference<T>::type type;
};

template <typename T>
struct remove_a2dobj : remove_a2dobj_impl<remove_cvref_t<T>> {};

template <typename T>
using remove_a2dobj_t = typename remove_a2dobj<T>::type;

/*
  Get the vector size
*/
template <typename T>
struct get_vec_size_impl {
  static constexpr int size = 0;
};

template <typename T, int N>
struct get_vec_size_impl<Vec<T, N>> {
  static constexpr int size = N;
};

template <typename T, int N>
struct get_vec_size_impl<VecRef<T, N>> {
  static constexpr int size = N;
};

template <class T>
struct get_vec_size : get_vec_size_impl<remove_a2dobj_t<T>> {
  static_assert(a2d_obj_type_v<T> == ADObjType::VECTOR,
                "get_symmatrix_size called on incorrect type");
};

/*
  Get the symmetric matrix size
*/
template <typename T>
struct get_symmatrix_size_impl {
  static constexpr int size = 0;
};

template <typename T, int N>
struct get_symmatrix_size_impl<SymMat<T, N>> {
  static constexpr int size = N;
};

template <typename T, int N>
struct get_symmatrix_size_impl<SymMatRef<T, N>> {
  static constexpr int size = N;
};

template <typename T>
struct get_symmatrix_size : get_symmatrix_size_impl<remove_a2dobj_t<T>> {
  static_assert(a2d_obj_type_v<T> == ADObjType::SYMMAT,
                "get_symmatrix_size called on incorrect type");
};

/*
  Get the number of matrix rows
*/
template <typename T>
struct get_matrix_rows_impl {
  static constexpr int size = 0;
};

template <typename T, int N, int M>
struct get_matrix_rows_impl<Mat<T, N, M>> {
  static constexpr int size = N;
};

template <typename T, int N, int M>
struct get_matrix_rows_impl<MatRef<T, N, M>> {
  static constexpr int size = N;
};

template <typename T, int N>
struct get_matrix_rows_impl<SymMat<T, N>> {
  static constexpr int size = N;
};

template <typename T, int N>
struct get_matrix_rows_impl<SymMatRef<T, N>> {
  static constexpr int size = N;
};

template <typename T>
struct get_matrix_rows : get_matrix_rows_impl<remove_a2dobj_t<T>> {
  static_assert((a2d_obj_type_v<T> == ADObjType::MATRIX ||
                 a2d_obj_type_v<T> == ADObjType::SYMMAT),
                "get_matrix_rows called on incorrect type");
};

/*
  Get the number of matrix columns
*/
template <typename T>
struct get_matrix_columns_impl {
  static constexpr int size = 0;
};

template <typename T, int N, int M>
struct get_matrix_columns_impl<Mat<T, N, M>> {
  static constexpr int size = M;
};

template <typename T, int N, int M>
struct get_matrix_columns_impl<MatRef<T, N, M>> {
  static constexpr int size = M;
};

template <typename T, int N>
struct get_matrix_columns_impl<SymMat<T, N>> {
  static constexpr int size = N;
};

template <typename T, int N>
struct get_matrix_columns_impl<SymMatRef<T, N>> {
  static constexpr int size = N;
};

template <class T>
struct get_matrix_columns : get_matrix_columns_impl<remove_a2dobj_t<T>> {
  static_assert((a2d_obj_type_v<T> == ADObjType::MATRIX ||
                 a2d_obj_type_v<T> == ADObjType::SYMMAT),
                "get_matrix_rows called on incorrect type");
};

/*
  Get the number of entries in the matrix -- symmetric or non-symmetric
*/
template <class T>
struct get_num_matrix_entries_impl {
  static constexpr int size = 0;
};

template <typename T, int N>
struct get_num_matrix_entries_impl<SymMat<T, N>> {
  static constexpr int size = N * (N + 1) / 2;
};

template <typename T, int N>
struct get_num_matrix_entries_impl<SymMatRef<T, N>> {
  static constexpr int size = N * (N + 1) / 2;
};

template <typename T, int N, int M>
struct get_num_matrix_entries_impl<Mat<T, N, M>> {
  static constexpr int size = N * M;
};

template <typename T, int N, int M>
struct get_num_matrix_entries_impl<MatRef<T, N, M>> {
  static constexpr int size = N * M;
};

template <class T>
struct get_num_matrix_entries
    : get_num_matrix_entries_impl<remove_a2dobj_t<T>> {
  static_assert((a2d_obj_type_v<T> == ADObjType::MATRIX ||
                 a2d_obj_type_v<T> == ADObjType::SYMMAT),
                "get_num_matrix_entries called on incorrect type");
};

/*
 * Vector class detection
 */
template <typename T>
struct is_a2d_vector_impl : std::false_type {};

template <typename U, int N>
struct is_a2d_vector_impl<Vec<U, N>> : std::true_type {};

template <typename U, int N>
struct is_a2d_vector_impl<VecRef<U, N>> : std::true_type {};

template <typename T>
struct is_a2d_vector : is_a2d_vector_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr bool is_a2d_vector_v = is_a2d_vector<T>::value;

/*
 * Matrix class detection
 */
template <typename T>
struct is_a2d_matrix_impl : std::false_type {};

template <typename U, int M, int N>
struct is_a2d_matrix_impl<Mat<U, M, N>> : std::true_type {};

template <typename U, int M, int N>
struct is_a2d_matrix_impl<MatRef<U, M, N>> : std::true_type {};

template <typename T>
struct is_a2d_matrix : is_a2d_matrix_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr bool is_a2d_matrix_v = is_a2d_matrix<T>::value;

/*
 * Symmetric matrix class detection
 */
template <typename T>
struct is_a2d_sym_matrix_impl : std::false_type {};

template <typename U, int N>
struct is_a2d_sym_matrix_impl<SymMat<U, N>> : std::true_type {};

template <typename U, int N>
struct is_a2d_sym_matrix_impl<SymMatRef<U, N>> : std::true_type {};

template <typename T>
struct is_a2d_sym_matrix : is_a2d_sym_matrix_impl<remove_cvref_t<T>> {};

template <typename T>
inline constexpr bool is_a2d_sym_matrix_v = is_a2d_sym_matrix<T>::value;

}  // namespace A2D

#endif  // A2D_DEFS_H
