#ifndef A2D_SCALAR_H
#define A2D_SCALAR_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <class T>
class ADScalar {
 public:
  KOKKOS_FUNCTION ADScalar(T value = 0.0, T bvalue = 0.0)
      : value(value), bvalue(bvalue) {}

  T value;
  T bvalue;
};

template <class T>
class A2DScalar {
 public:
  KOKKOS_FUNCTION A2DScalar(T value = 0.0, T bvalue = 0.0, T pvalue = 0.0,
                                T hvalue = 0.0)
      : value(value), bvalue(bvalue), pvalue(pvalue), hvalue(hvalue) {}

  T value;
  T bvalue;
  T pvalue;
  T hvalue;
};

double& get_data(double& value) { return value; }

std::complex<double>& get_data(std::complex<double>& value) { return value; }

template <typename T>
T& get_data(ADScalar<T>& value) {
  return value.value;
}

template <typename T>
T& get_data(A2DScalar<T>& value) {
  return value.value;
}

/**
 * @brief Select type based on whether the scalar is passive or active (can be
 * differentiated)
 *
 * @tparam adiff_type passive or active
 * @tparam order first (AD) or second (A2D)
 * @tparam MatType the numeric type of the matrix
 */
template <ADiffType adiff_type, ADorder order, typename T>
using ADScalarType = typename std::conditional<
    adiff_type == ADiffType::ACTIVE,
    typename std::conditional<order == ADorder::FIRST, ADScalar<T>,
                              A2DScalar<T>>::type,
    T>::type;

}  // namespace A2D

#endif  // A2D_SCALAR_H