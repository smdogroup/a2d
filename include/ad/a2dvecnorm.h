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

template <typename T, int N, ADorder order>
class VecNormExpr {
 public:
  using vtype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION VecNormExpr(vtype &x, dtype &alpha) : x(x), alpha(alpha) {}

  KOKKOS_FUNCTION void eval() {
    get_data(alpha) = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {}
  KOKKOS_FUNCTION void reverse() {}
  KOKKOS_FUNCTION void hreverse() {}

  vtype &x;
  dtype &alpha;
};

template <typename T, int N>
void VecNormalize(const Vec<T, N> &x, Vec<T, N> &y) {
  T alpha = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
  VecScaleCore<T, N>(1.0 / alpha, get_data(x), get_data(y));
}

template <typename T, int N, ADorder order>
class VecNormalizeExpr {
 public:
  using vtype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION VecNormalizeExpr(vtype &x, vtype &y) : x(x), y(y) {}

  KOKKOS_FUNCTION void eval() {
    alpha = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
    VecScaleCore<T, N>(1.0 / alpha, get_data(x), get_data(y));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {}
  KOKKOS_FUNCTION void reverse() {}
  KOKKOS_FUNCTION void hreverse() {}

  T alpha;
  vtype &x;
  vtype &y;
};

template <typename T, int N>
void VecScale(const T alpha, const Vec<T, N> &x, Vec<T, N> &y) {
  VecScaleCore<T, N>(alpha, get_data(x), get_data(y));
}

template <typename T, int N, ADorder order>
class VecScaleExpr {
 public:
  using vtype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;
  using dtype = ADScalarInputType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION VecNormExpr(vtype &x, dtype &alpha) : x(x), alpha(alpha) {}

  KOKKOS_FUNCTION void eval() {
    get_data(alpha) = std::sqrt(VecDotCore<T, N>(get_data(x), get_data(x)));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {}
  KOKKOS_FUNCTION void reverse() {}
  KOKKOS_FUNCTION void hreverse() {}

  vtype &x;
  dtype &alpha;
};

template <typename T, int N>
void VecDot(const Vec<T, N> &x, const Vec<T, N> &y, T &alpha) {
  alpha = VecDotCore<T, N>(get_data(x), get_data(y));
}

template <typename T, int N, ADorder order>
class VecDotExpr {
 public:
  using xtype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;
  using ytype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION VecDotExpr(xtype &x, ytype &y, dtype &alpha)
      : x(x), y(y), alpha(alpha) {}

  KOKKOS_FUNCTION void eval() {
    get_data(alpha) = VecDotCore<T, N>(get_data(x), get_data(y));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {}
  KOKKOS_FUNCTION void reverse() {}
  KOKKOS_FUNCTION void hreverse() {}

  vtype &x;
  dtype &alpha;
};

namespace Test {}

}  // namespace A2D

#endif  //  A2D_VEC_NORM_H