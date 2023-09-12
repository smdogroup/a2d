/**
 * The multi-dimensional array data structure and BLAS operations.
 */

#ifndef A2D_ARRAY_H
#define A2D_ARRAY_H

#include <numeric>

#include "Kokkos_Core.hpp"
#include "a2ddefs.h"

namespace A2D {

/**
 * Set the default memory layout.
 *
 * Options:
 *   Kokkos::LayoutRight right-most dimension increases fastest
 *   Kokkos::LayoutLeft left-most dimension increases fastest
 */
using DefaultLayout = Kokkos::LayoutRight;

/**
 * Set the default memory space.
 *
 * Options:
 *   Kokkos::HostSpace CPU memory space
 *   Kokkos::CudaSpace device space
 *   Kokkos::CudaUVMSpace (CPU and device) unified memory space
 */
#ifdef KOKKOS_ENABLE_CUDA
using DefaultMemSpace = Kokkos::CudaUVMSpace;
#else
using DefaultMemSpace = Kokkos::HostSpace;
#endif

/**
 * The alias of Kokkos' Multi-dimensional array.
 *
 * Example usage:
 *   To create a multiarray of shape (n, N1, N2), where n is runtime dimension,
 *   N1 and N2 are compile-time dimensions, and T is the scalar type:
 *   MultiArrayNew<T*[N1][N2]> array("label", n);
 */
template <class DataType>
using MultiArrayNew = Kokkos::View<DataType, DefaultLayout, DefaultMemSpace>;

/**
 * A collection of BLAS operations for the multiarray data type.
 *
 * TODO: use KokkosKernel instead
 */
namespace BLAS {

/*
  Zero all elements in the array
*/
template <class View>
KOKKOS_FUNCTION void zero(const View& x) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  T* data = x.data();
  A2D::fill(data, data + x.size(), T(0));
}

/*
  Fill all the values in the array with the specified value
*/
template <class View>
KOKKOS_FUNCTION void fill(const View& x,
                          typename View::const_value_type& value) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  T* data = x.data();
  A2D::fill(data, data + x.size(), value);
}

/*
  Copy elements from the source to this vector
*/
template <class DestView, class SrcView>
void copy(const DestView& dest, const SrcView& src) {
  Kokkos::deep_copy(dest, src);
}

/*
  Scale the array
*/
template <class View>
KOKKOS_FUNCTION void scale(const View& x,
                           typename View::const_value_type& alpha) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  T* data = x.data();
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    data[i] *= alpha;
  }
}

/*
  Set a random seed for array x
*/
template <class View>
void random(const View& x, typename View::const_value_type& lower = -1.0,
            typename View::const_value_type& upper = 1.0) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  T* data = x.data();
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    data[i] =
        lower + ((upper - lower) * (1.0 * std::rand())) / (1.0 * RAND_MAX);
  }
}

/*
  Take the dot product with the source vector data
*/
template <class View>
typename View::value_type dot(const View& x, const View& y) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  T* data = x.data();
  size_t len = x.size();
  return std::inner_product(data, data + len, y.data(), T(0));
}

/*
  Norm of the array
*/
template <class View>
typename View::value_type norm(const View& x) {
  using T = typename View::value_type;
  T* data = x.data();
  size_t len = x.size();
  return sqrt(std::inner_product(data, data + len, data, T(0)));
}

/*
  Axpy: y += alpha * x
*/
template <class View>
void axpy(View& y, typename View::const_value_type& alpha, const View& x) {
  assert(x.span_is_contiguous());
  assert(y.span_is_contiguous());
  using T = typename View::value_type;
  T* x_data = x.data();
  T* y_data = y.data();
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    y_data[i] += alpha * x_data[i];
  }
}

/*
  Axpby: y = alpha * x + beta * y
*/
template <class View>
void axpby(View& y, typename View::const_value_type& alpha,
           typename View::const_value_type& beta, const View& x) {
  using T = typename View::value_type;
  assert(x.span_is_contiguous());
  assert(y.span_is_contiguous());
  T* x_data = x.data();
  T* y_data = y.data();
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    y_data[i] = alpha * x_data[i] + beta * y_data[i];
  }
}

}  // namespace BLAS

// Array type shortcuts
using IdxArray1D_t = MultiArrayNew<index_t*>;
template <typename T>
using ValArray1D_t = MultiArrayNew<T*>;

}  // namespace A2D

#endif  // A2D_ARRAY_H