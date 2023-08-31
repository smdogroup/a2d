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
void VecSum(const Vec<T, N> &x, const Vec<T, N> &y, Vec<T, N> &z) {
  VecSumCore<T, N>(get_data(x), get_data(y), get_data(z));
}

template <typename T, int N>
void VecSum(const Vec<T, N> &x, const Vec<T, N> &y, Vec<T, N> &z) {
  VecSumCore<T, N * M>(get_data(x), get_data(y), get_data(z));
}

namespace Test {}

}  // namespace A2D

#endif  // A2D_VEC_NORM_H