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
    inv = 1.0 / get_data(alpha);
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(alpha) =
        inv * VecDotCore<T, N>(GetSeed<seed>::get_data(x), get_data(x));
  }
  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    VecAddCore<T, N>(inv * GetSeed<seed>::get_data(alpha), get_data(x),
                     GetSeed<seed>::get_data(x));
  }
  KOKKOS_FUNCTION void hreverse() {
    VecAddCore<T, N>(inv * GetSeed<ADseed::h>::get_data(alpha), get_data(x),
                     GetSeed<ADseed::h>::get_data(x));

    VecAddCore<T, N>(inv * GetSeed<ADseed::b>::get_data(alpha),
                     GetSeed<ADseed::p>::get_data(x),
                     GetSeed<ADseed::h>::get_data(x));

    T scale = -inv * inv * inv * GetSeed<ADseed::b>::get_data(alpha) *
              VecDotCore<T, N>(GetSeed<ADseed::p>::get_data(x), get_data(x));
    VecAddCore<T, N>(scale, get_data(x), GetSeed<ADseed::h>::get_data(x));
  }

  vtype &x;
  dtype &alpha;
  T inv;
};

template <typename T, int N>
auto VecNorm(ADVec<Vec<T, N>> &x, ADScalar<T> &alpha) {
  return VecNormExpr<T, N, ADorder::FIRST>(x, alpha);
}

template <typename T, int N>
auto VecNorm(A2DVec<Vec<T, N>> &x, A2DScalar<T> &alpha) {
  return VecNormExpr<T, N, ADorder::SECOND>(x, alpha);
}

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
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
  }
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

template <typename T, int N, ADorder order, ADiffType ada, ADiffType adx>
class VecScaleExpr {
 public:
  using dtype = ADScalarInputType<ada, order, T>;
  using xtype = ADVecType<adx, order, Vec<T, N>>;
  using ytype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;

  KOKKOS_FUNCTION VecScaleExpr(dtype alpha, xtype &x, ytype &y)
      : alpha(alpha), x(x), y(y) {}

  KOKKOS_FUNCTION void eval() {
    VecScaleCore<T, N>(get_data(alpha), get_data(x), get_data(y));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (ada == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      VecScaleCore<T, N>(GetSeed<seed>::get_data(alpha), get_data(x),
                         GetSeed<seed>::get_data(y));
      VecAddCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(x),
                       GetSeed<seed>::get_data(y));
    } else if constexpr (ada == ADiffType::ACTIVE) {
      VecScaleCore<T, N>(GetSeed<seed>::get_data(alpha), get_data(x),
                         GetSeed<seed>::get_data(y));
    } else if constexpr (adx == ADiffType::ACTIVE) {
      VecScaleCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(x),
                         GetSeed<seed>::get_data(y));
    }
  }
  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, N>(GetSeed<seed>::get_data(y), get_data(x));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(y),
                       GetSeed<seed>::get_data(x));
    }
  }
  KOKKOS_FUNCTION void hreverse() {
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(alpha) +=
          VecDotCore<T, N>(GetSeed<ADseed::h>::get_data(y), get_data(x));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(alpha), GetSeed<ADseed::h>::get_data(y),
                       GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (ada == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(alpha) += VecDotCore<T, N>(
          GetSeed<ADseed::b>::get_data(y), GetSeed<ADseed::p>::get_data(x));
      VecAddCore<T, N>(GetSeed<ADseed::p>::get_data(alpha),
                       GetSeed<ADseed::b>::get_data(y),
                       GetSeed<ADseed::h>::get_data(x));
    }
  }

  dtype alpha;
  xtype &x;
  ytype &y;
};

template <typename T, int N>
auto VecScale(ADScalar<T> &alpha, ADVec<Vec<T, N>> &x, ADVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(alpha, x, y);
}

template <typename T, int N>
auto VecScale(const T alpha, ADVec<Vec<T, N>> &x, ADVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::FIRST, ADiffType::PASSIVE,
                      ADiffType::ACTIVE>(alpha, x, y);
}

template <typename T, int N>
auto VecScale(ADScalar<T> alpha, const Vec<T, N> &x, ADVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                      ADiffType::PASSIVE>(alpha, x, y);
}

template <typename T, int N>
auto VecScale(A2DScalar<T> &alpha, A2DVec<Vec<T, N>> &x, A2DVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(alpha, x, y);
}

template <typename T, int N>
auto VecScale(const T alpha, A2DVec<Vec<T, N>> &x, A2DVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::SECOND, ADiffType::PASSIVE,
                      ADiffType::ACTIVE>(alpha, x, y);
}

template <typename T, int N>
auto VecScale(A2DScalar<T> &alpha, const Vec<T, N> &x, A2DVec<Vec<T, N>> &y) {
  return VecScaleExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                      ADiffType::PASSIVE>(alpha, x, y);
}

template <typename T, int N>
void VecDot(const Vec<T, N> &x, const Vec<T, N> &y, T &alpha) {
  alpha = VecDotCore<T, N>(get_data(x), get_data(y));
}

template <typename T, int N, ADorder order, ADiffType adx, ADiffType ady>
class VecDotExpr {
 public:
  using xtype = ADVecType<adx, order, Vec<T, N>>;
  using ytype = ADVecType<ady, order, Vec<T, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION VecDotExpr(xtype &x, ytype &y, dtype &alpha)
      : x(x), y(y), alpha(alpha) {}

  KOKKOS_FUNCTION void eval() {
    get_data(alpha) = VecDotCore<T, N>(get_data(x), get_data(y));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) =
          VecDotCore<T, N>(GetSeed<seed>::get_data(x), get_data(y)) +
          VecDotCore<T, N>(get_data(x), GetSeed<seed>::get_data(y));

    } else if constexpr (adx == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) =
          VecDotCore<T, N>(GetSeed<seed>::get_data(x), get_data(y));
    } else if constexpr (ady == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) =
          VecDotCore<T, N>(get_data(x), GetSeed<seed>::get_data(y));
    }
  }
  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(alpha), get_data(y),
                       GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(alpha), get_data(x),
                       GetSeed<seed>::get_data(y));
    }
  }
  KOKKOS_FUNCTION void hreverse() {
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<ADseed::h>::get_data(alpha), get_data(y),
                       GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<ADseed::h>::get_data(alpha), get_data(x),
                       GetSeed<ADseed::h>::get_data(y));
    }
    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<ADseed::b>::get_data(alpha),
                       GetSeed<ADseed::p>::get_data(y),
                       GetSeed<ADseed::h>::get_data(x));
      VecAddCore<T, N>(GetSeed<ADseed::b>::get_data(alpha),
                       GetSeed<ADseed::p>::get_data(x),
                       GetSeed<ADseed::h>::get_data(y));
    }
  }

  xtype &x;
  ytype &y;
  dtype &alpha;
};

template <typename T, int N>
auto VecDot(ADVec<Vec<T, N>> &x, ADVec<Vec<T, N>> &y, ADScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      x, y, alpha);
}

template <typename T, int N>
auto VecDot(const Vec<T, N> &x, ADVec<Vec<T, N>> &y, ADScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::FIRST, ADiffType::PASSIVE,
                    ADiffType::ACTIVE>(x, y, alpha);
}

template <typename T, int N>
auto VecDot(ADVec<Vec<T, N>> &x, const Vec<T, N> &y, ADScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::PASSIVE>(x, y, alpha);
}

template <typename T, int N>
auto VecDot(A2DVec<Vec<T, N>> &x, A2DVec<Vec<T, N>> &y, A2DScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::ACTIVE>(x, y, alpha);
}

template <typename T, int N>
auto VecDot(const Vec<T, N> &x, A2DVec<Vec<T, N>> &y, A2DScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::SECOND, ADiffType::PASSIVE,
                    ADiffType::ACTIVE>(x, y, alpha);
}

template <typename T, int N>
auto VecDot(A2DVec<Vec<T, N>> &x, const Vec<T, N> &y, A2DScalar<T> &alpha) {
  return VecDotExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::PASSIVE>(x, y, alpha);
}

namespace Test {

template <typename T, int N>
class VecNormTest : public A2DTest<T, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "VecNorm"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input &X) {
    T alpha;
    Vec<T, N> x;
    X.get_values(x);
    VecNorm(x, alpha);
    return MakeVarTuple<T>(alpha);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &X, Input &g) {
    ADScalar<T> alpha;
    Vec<T, N> x0, xb;
    ADVec<Vec<T, N>> x(x0, xb);
    X.get_values(x0);
    auto op = VecNorm(x, alpha);
    auto stack = MakeStack(op);
    seed.get_values(alpha.bvalue);
    stack.reverse();
    g.set_values(xb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &X,
             const Input &p, Input &h) {
    A2DScalar<T> alpha;
    A2DVec<Vec<T, N>> x;
    X.get_values(x.value());
    p.get_values(x.pvalue());
    auto op = VecNorm(x, alpha);
    auto stack = MakeStack(op);
    seed.get_values(alpha.bvalue);
    hval.get_values(alpha.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(x.hvalue());
  }
};

bool VecNormTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecNormTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecNormTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

template <typename T, int N>
class VecScaleTest : public A2DTest<T, Vec<T, N>, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, T, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() { return "VecScale"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input &X) {
    T alpha;
    Vec<T, N> x, y;
    X.get_values(alpha, x);
    VecScale(alpha, x, y);
    return MakeVarTuple<T>(y);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &X, Input &g) {
    ADScalar<T> alpha;
    Vec<T, N> x0, xb, y0, yb;
    ADVec<Vec<T, N>> x(x0, xb), y(y0, yb);
    X.get_values(alpha.value, x0);
    auto op = VecScale(alpha, x, y);
    auto stack = MakeStack(op);
    seed.get_values(yb);
    stack.reverse();
    g.set_values(alpha.bvalue, xb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &X,
             const Input &p, Input &h) {
    A2DScalar<T> alpha;
    A2DVec<Vec<T, N>> x, y;
    X.get_values(alpha.value, x.value());
    p.get_values(alpha.pvalue, x.pvalue());
    auto op = VecScale(alpha, x, y);
    auto stack = MakeStack(op);
    seed.get_values(y.bvalue());
    hval.get_values(y.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(alpha.hvalue, x.hvalue());
  }
};

bool VecScaleTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecScaleTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecScaleTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

template <typename T, int N>
class VecDotTest : public A2DTest<T, T, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "VecDot"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input &X) {
    T alpha;
    Vec<T, N> x, y;
    X.get_values(x, y);
    VecDot(x, y, alpha);
    return MakeVarTuple<T>(alpha);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &X, Input &g) {
    ADScalar<T> alpha;
    Vec<T, N> x0, xb, y0, yb;
    ADVec<Vec<T, N>> x(x0, xb), y(y0, yb);
    X.get_values(x0, y0);
    auto op = VecDot(x, y, alpha);
    auto stack = MakeStack(op);
    seed.get_values(alpha.bvalue);
    stack.reverse();
    g.set_values(xb, yb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &X,
             const Input &p, Input &h) {
    A2DScalar<T> alpha;
    A2DVec<Vec<T, N>> x, y;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto op = VecDot(x, y, alpha);
    auto stack = MakeStack(op);
    seed.get_values(alpha.bvalue);
    hval.get_values(alpha.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

bool VecDotTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecDotTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecDotTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  //  A2D_VEC_NORM_H