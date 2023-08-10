#ifndef A2D_SCALAR_H
#define A2D_SCALAR_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <class T>
class ADScalar {
 public:
  A2D_INLINE_FUNCTION ADScalar(T value = 0.0, T bvalue = 0.0)
      : value(value), bvalue(bvalue) {}

  T value;
  T bvalue;
};

template <class T>
class A2DScalar {
 public:
  A2D_INLINE_FUNCTION A2DScalar(T value = 0.0, T bvalue = 0.0, T pvalue = 0.0,
                                T hvalue = 0.0)
      : value(value), bvalue(bvalue), pvalue(pvalue), hvalue(hvalue) {}

  T value;
  T bvalue;
  T pvalue;
  T hvalue;
};

/**
 * @brief Select type based on whether the scalar is passive or active (can be
 * differentiated)
 *
 * @tparam adiff_type passive or active
 * @tparam MatType the numeric type of the matrix
 */
template <ADiffType adiff_type, typename T>
using ADScalarType = typename std::conditional<adiff_type == ADiffType::ACTIVE,
                                               ADScalar<T>, T>::type;

template <ADiffType adiff_type, typename T>
using A2DScalarType = typename std::conditional<adiff_type == ADiffType::ACTIVE,
                                                A2DScalar<T>, T>::type;

}  // namespace A2D
#endif  // A2D_SCALAR_H