#ifndef A2D_VEC_CROSS_H
#define A2D_VEC_CROSS_H

#include "a2ddefs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "a2dvec.h"
#include "ad/core/a2dveccore.h"

namespace A2D {

template <typename T>
void VecCross(const Vec<T, 3> &x, const Vec<T, 3> &y, Vec<T, 3> &z) {
  VecCrossCore<T>(get_data(x), get_data(y), get_data(z));
}

template <class xtype, class ytype, class ztype>
class VecCrossExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<ztype>::type T;

  // Extract the dimensions of the underlying vectors
  static constexpr int N = get_vec_size<xtype>::size;
  static constexpr int M = get_vec_size<ytype>::size;
  static constexpr int K = get_vec_size<ztype>::size;

  // Get the types of the matrices
  static constexpr ADiffType adx = get_diff_type<xtype>::diff_type;
  static constexpr ADiffType ady = get_diff_type<ytype>::diff_type;

  // Make sure the matrix dimensions are consistent
  static_assert((N == M && M == K && K == 3), "Vector dimensions must agree");

  KOKKOS_FUNCTION VecCrossExpr(xtype &x, ytype &y, ztype &z)
      : x(x), y(y), z(z) {}

  KOKKOS_FUNCTION void eval() {
    VecCrossCore<T>(get_data(x), get_data(y), get_data(z));
  }

  KOKKOS_FUNCTION void bzero() { z.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      VecCrossCore<T>(GetSeed<seed>::get_data(x), get_data(y),
                      GetSeed<seed>::get_data(z));
      VecCrossCoreAdd<T>(get_data(x), GetSeed<seed>::get_data(y),
                         GetSeed<seed>::get_data(z));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      VecCrossCore<T>(GetSeed<seed>::get_data(x), get_data(y),
                      GetSeed<seed>::get_data(z));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecCrossCore<T>(get_data(x), GetSeed<seed>::get_data(y),
                      GetSeed<seed>::get_data(z));
    }
  }
  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;

    if constexpr (adx == ADiffType::ACTIVE) {
      VecCrossCoreAdd<T>(get_data(y), GetSeed<seed>::get_data(z),
                         GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecCrossCoreAdd<T>(GetSeed<seed>::get_data(z), get_data(x),
                         GetSeed<seed>::get_data(y));
    }
  }

  KOKKOS_FUNCTION void hzero() { z.hzero(); }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (adx == ADiffType::ACTIVE) {
      VecCrossCoreAdd<T>(get_data(y), GetSeed<ADseed::h>::get_data(z),
                         GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecCrossCoreAdd<T>(GetSeed<ADseed::h>::get_data(z), get_data(x),
                         GetSeed<ADseed::h>::get_data(y));
    }
    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      VecCrossCoreAdd<T>(GetSeed<ADseed::p>::get_data(y),
                         GetSeed<ADseed::b>::get_data(z),
                         GetSeed<ADseed::h>::get_data(x));
      VecCrossCoreAdd<T>(GetSeed<ADseed::b>::get_data(z),
                         GetSeed<ADseed::p>::get_data(x),
                         GetSeed<ADseed::h>::get_data(y));
    }
  }

  xtype &x;
  ytype &y;
  ztype &z;
};

template <class xtype, class ytype, class ztype>
auto VecCross(ADObj<xtype> &x, ADObj<ytype> &y, ADObj<ztype> &z) {
  return VecCrossExpr<ADObj<xtype>, ADObj<ytype>, ADObj<ztype>>(x, y, z);
}

template <class xtype, class ytype, class ztype>
auto VecCross(const xtype &x, ADObj<ytype> &y, ADObj<ztype> &z) {
  return VecCrossExpr<const xtype, ADObj<ytype>, ADObj<ztype>>(x, y, z);
}

template <class xtype, class ytype, class ztype>
auto VecCross(ADObj<xtype> &x, const ytype &y, ADObj<ztype> &z) {
  return VecCrossExpr<ADObj<xtype>, const ytype, ADObj<ztype>>(x, y, z);
}

template <class xtype, class ytype, class ztype>
auto VecCross(A2DObj<xtype> &x, A2DObj<ytype> &y, A2DObj<ztype> &z) {
  return VecCrossExpr<A2DObj<xtype>, A2DObj<ytype>, A2DObj<ztype>>(x, y, z);
}

template <class xtype, class ytype, class ztype>
auto VecCross(const xtype &x, A2DObj<ytype> &y, A2DObj<ztype> &z) {
  return VecCrossExpr<const xtype, A2DObj<ytype>, A2DObj<ztype>>(x, y, z);
}

template <class xtype, class ytype, class ztype>
auto VecCross(A2DObj<xtype> &x, const ytype &y, A2DObj<ztype> &z) {
  return VecCrossExpr<A2DObj<xtype>, const ytype, A2DObj<ztype>>(x, y, z);
}

namespace Test {

template <typename T>
class VecCrossTest : public A2DTest<T, Vec<T, 3>, Vec<T, 3>, Vec<T, 3>> {
 public:
  using Input = VarTuple<T, Vec<T, 3>, Vec<T, 3>>;
  using Output = VarTuple<T, Vec<T, 3>>;

  // Assemble a string to describe the test
  std::string name() { return "VecCross"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input &X) {
    Vec<T, 3> x, y, z;
    X.get_values(x, y);
    VecCross(x, y, z);
    return MakeVarTuple<T>(z);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &X, Input &g) {
    ADObj<Vec<T, 3>> x, y, z;
    X.get_values(x.value(), y.value());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    stack.reverse();
    g.set_values(x.bvalue(), y.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &X,
             const Input &p, Input &h) {
    A2DObj<Vec<T, 3>> x, y, z;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    hval.get_values(z.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

bool VecCrossTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecCrossTest<Tc> test1;
  passed = passed && Run(test1, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_VEC_CROSS_H