#ifndef A2D_OBJS_H
#define A2D_OBJS_H

#include <cmath>
#include <cstdint>

#include "Kokkos_Core.hpp"
#include "Kokkos_UnorderedMap.hpp"

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
KOKKOS_FUNCTION double fmt(A2D_complex_t<T> val) {
  return val.real();
}

KOKKOS_FUNCTION double fmt(double val) { return val; }

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
KOKKOS_FUNCTION T sqrt(T val) {
  return cuda::std::sqrt(val);
}

template <typename T>
KOKKOS_FUNCTION T exp(T val) {
  return cuda::std::exp(val);
}

template <typename T>
KOKKOS_FUNCTION T log(T val) {
  return cuda::std::log(val);
}

template <class ForwardIt, class T>
KOKKOS_FUNCTION void fill(ForwardIt first, ForwardIt last, const T& value) {
  thrust::fill(first, last, value);
}
#else
// template <typename T>
// KOKKOS_FUNCTION T sqrt(T val) {
//   return std::sqrt(val);
// }

// template <typename T>
// KOKKOS_FUNCTION T exp(T val) {
//   return std::exp(val);
// }

// template <typename T>
// KOKKOS_FUNCTION T log(T val) {
//   return std::log(val);
// }
template <class ForwardIt, class T>
void fill(ForwardIt first, ForwardIt last, const T& value) {
  std::fill(first, last, value);
}
#endif

}  // namespace A2D

#endif  // A2D_OBJS_H
