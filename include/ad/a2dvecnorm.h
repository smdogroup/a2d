#ifndef A2D_VEC_NORM_H
#define A2D_VEC_NORM_H

#include "a2denum.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "a2dvec.h"
#include "ad/core/a2dveccore.h"

namespace A2D {

template <typename T, int N>
void VecNorm(const Vec<T, N> &x, T &alpha) {
  alpha = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
}

template <typename T, int N>
void VecNormalize(const Vec<T, N> &x, Vec<T, N> &y) {
  T alpha = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
  VecScaleCore<T, N>(1.0 / alpha, get_data(x), get_data(y));
}

template <typename T, int N>
void VecScale(const T alpha, const Vec<T, N> &x, Vec<T, N> &y) {
  VecScaleCore<T, N>(alpha, get_data(x), get_data(y));
}

template <typename T, int N>
void VecDot(const Vec<T, N> &x, const Vec<T, N> &y, T &alpha) {
  alpha = VecDotCore<T, N>(get_data(x), get_data(y));
}

namespace Test {}

}  // namespace A2D

#endif  //  A2D_VEC_NORM_H