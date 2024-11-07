#ifndef A2D_DEFS_H
#define A2D_DEFS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
template <typename T>
using A2D_complex_t = std::complex<T>;

#ifndef A2D_FUNCTION
#define A2D_FUNCTION  // A2D_FUNCTION does nothing in this scenario
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

/**
 * @brief The type of variable type
 *
 */
enum class FEVarType { DATA, GEOMETRY, STATE };

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
A2D_FUNCTION double fmt(A2D_complex_t<T> val) {
  return val.real();
}

A2D_FUNCTION double fmt(double val) { return val; }

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

/*
  Remove the const-ness and references for a type
*/
template <class T>
struct remove_const_and_refs
    : std::remove_const<typename std::remove_reference<T>::type> {};

/*
  Check if a type is numeric or not - this includes complex numbers
*/
template <class T>
struct __is_numeric_type : std::is_floating_point<T> {};

template <class T>
struct __is_numeric_type<std::complex<T>> : std::is_floating_point<T> {};

template <class T>
struct is_numeric_type
    : __is_numeric_type<typename remove_const_and_refs<T>::type> {};

/*
  Check if the type is a scalar - an arithmetic or complex type
*/
template <class T>
struct is_scalar_type {
  static const bool value =
      is_numeric_type<T>::value || std::is_arithmetic<T>::value;
};

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

/**
 * @brief compile-time conditional for type values
 */
template <bool B, class TrueType, class FalseType>
struct conditional_type {
  using type = TrueType;
};

template <class TrueType, class FalseType>
struct conditional_type<false, TrueType, FalseType> {
  using type = FalseType;
};

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

#ifdef KOKKOS_ENABLE_CUDA
template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T sqrt(T val) {
  return cuda::std::sqrt(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T exp(T val) {
  return cuda::std::exp(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T log(T val) {
  return cuda::std::log(val);
}

template <class ForwardIt, class T>
A2D_FUNCTION void fill(ForwardIt first, ForwardIt last, const T& value) {
  thrust::fill(first, last, value);
}
#else

template <typename T, typename R,
          std::enable_if_t<is_scalar_type<T>::value, bool> = true,
          std::enable_if_t<is_scalar_type<R>::value, bool> = true>
A2D_FUNCTION T pow(T val, R exponent) {
  return std::pow(val, exponent);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T fabs(T val) {
  return std::fabs(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T sqrt(T val) {
  return std::sqrt(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T exp(T val) {
  return std::exp(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T log(T val) {
  return std::log(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T sin(T val) {
  return std::sin(val);
}

template <typename T, std::enable_if_t<is_scalar_type<T>::value, bool> = true>
A2D_FUNCTION T cos(T val) {
  return std::cos(val);
}

template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T& value) {
  std::fill(first, last, value);
}
#endif

}  // namespace A2D

#endif  // A2D_DEFS_H
