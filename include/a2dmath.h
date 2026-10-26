#ifndef A2D_MATH_H
#define A2D_MATH_H

#include "a2ddefs.h"

namespace A2D {

template <typename T, typename R,
          std::enable_if_t<is_passive_scalar_v<T>, int> = 0,
          std::enable_if_t<is_passive_scalar_v<R>, int> = 0>
A2D_FUNCTION T pow(T val, R exponent) {
#ifndef __CUDACC__
  return std::pow(val, exponent);
#else
  return cuda::std::pow(val, exponent);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T fabs(T val) {
#ifndef __CUDACC__
  return std::fabs(val);
#else
  return cuda::std::fabs(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T fsgn(T val) {
  return std::copysign(1.0, val);
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T sqrt(T val) {
#ifndef __CUDACC__
  return std::sqrt(val);
#else
  return cuda::std::sqrt(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T exp(T val) {
#ifndef __CUDACC__
  return std::exp(val);
#else
  return cuda::std::exp(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T log(T val) {
#ifndef __CUDACC__
  return std::log(val);
#else
  return cuda::std::log(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T sin(T val) {
#ifndef __CUDACC__
  return std::sin(val);
#else
  return cuda::std::sin(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T cos(T val) {
#ifndef __CUDACC__
  return std::cos(val);
#else
  return cuda::std::cos(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T tan(T val) {
#ifndef __CUDACC__
  return std::tan(val);
#else
  return cuda::std::tan(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T asin(T val) {
#ifndef __CUDACC__
  return std::asin(val);
#else
  return cuda::std::asin(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T acos(T val) {
#ifndef __CUDACC__
  return std::acos(val);
#else
  return cuda::std::acos(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T atan(T val) {
#ifndef __CUDACC__
  return std::atan(val);
#else
  return cuda::std::atan(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T atan2(T y, T x) {
#ifndef __CUDACC__
  return std::atan2(y, x);
#else
  return cuda::std::atan2(y, x);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T sinh(T val) {
#ifndef __CUDACC__
  return std::sinh(val);
#else
  return cuda::std::sinh(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T cosh(T val) {
#ifndef __CUDACC__
  return std::cosh(val);
#else
  return cuda::std::cosh(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T tanh(T val) {
#ifndef __CUDACC__
  return std::tanh(val);
#else
  return cuda::std::tanh(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T asinh(T val) {
#ifndef __CUDACC__
  return std::asinh(val);
#else
  return cuda::std::asinh(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T acosh(T val) {
#ifndef __CUDACC__
  return std::acosh(val);
#else
  return cuda::std::acosh(val);
#endif
}

template <typename T, std::enable_if_t<is_passive_scalar_v<T>, int> = 0>
A2D_FUNCTION T atanh(T val) {
#ifndef __CUDACC__
  return std::atanh(val);
#else
  return cuda::std::atanh(val);
#endif
}

template <class ForwardIt, class T>
A2D_FUNCTION void fill(ForwardIt first, ForwardIt last, const T& value) {
#ifdef __CUDACC__
  thrust::fill(first, last, value);
#else
  std::fill(first, last, value);
#endif
}

// A2D_FUNCTION int A2D_rand() {
// #ifdef __CUDACC__
//   return 123456;
// #else
//   return rand();
// #endif
// }

}  // namespace A2D

#endif  // A2D_MATH_H