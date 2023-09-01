#ifndef A2D_VEC_SUM_H
#define A2D_VEC_SUM_H

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

template <typename T, int N, ADorder order, ADiffType adx, ADiffType ady>
class VecSumExpr {
 public:
  using xtype = ADVecType<adx, order, Vec<T, N>>;
  using ytype = ADVecType<ady, order, Vec<T, N>>;
  using ztype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;

  KOKKOS_FUNCTION
  VecSumExpr(xtype &x, ytype &y, ztype &z) : x(x), y(y), z(z) {}

  KOKKOS_FUNCTION void eval() {
    VecSumCore<T, N>(get_data(x), get_data(y), get_data(z));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      VecSumCore<T, N>(GetSeed<seed>::get_data(x), GetSeed<seed>::get_data(y),
                       GetSeed<seed>::get_data(z));
    } else if constexpr (adx == ADiffType::ACTIVE) {
      VecCopyCore<T, N>(GetSeed<seed>::get_data(x), GetSeed<seed>::get_data(z));
    } else if constexpr (ady == ADiffType::ACTIVE) {
      VecCopyCore<T, N>(GetSeed<seed>::get_data(y), GetSeed<seed>::get_data(z));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(z), GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(z), GetSeed<seed>::get_data(y));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(z), GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(GetSeed<seed>::get_data(z), GetSeed<seed>::get_data(y));
    }
  }

  xtype &x;
  ytype &y;
  ztype &z;
};

template <typename T, int N>
auto VecSum(ADVec<Vec<T, N>> &x, ADVec<Vec<T, N>> &y, ADVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      x, y, z);
}

template <typename T, int N>
auto VecSum(const Vec<T, N> &x, ADVec<Vec<T, N>> &y, ADVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::FIRST, ADiffType::PASSIVE,
                    ADiffType::ACTIVE>(x, y, z);
}

template <typename T, int N>
auto VecSum(ADVec<Vec<T, N>> &x, const Vec<T, N> &y, ADVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                    ADiffType::PASSIVE>(x, y, z);
}

template <typename T, int N>
auto VecSum(A2DVec<Vec<T, N>> &x, A2DVec<Vec<T, N>> &y, A2DVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::ACTIVE>(x, y, z);
}

template <typename T, int N>
auto VecSum(const Vec<T, N> &x, A2DVec<Vec<T, N>> &y, A2DVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::SECOND, ADiffType::PASSIVE,
                    ADiffType::ACTIVE>(x, y, z);
}

template <typename T, int N>
auto VecSum(A2DVec<Vec<T, N>> &x, const Vec<T, N> &y, A2DVec<Vec<T, N>> &z) {
  return VecSumExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                    ADiffType::PASSIVE>(x, y, z);
}

template <typename T, int N>
void VecSum(const T alpha, const Vec<T, N> &x, const T beta, const Vec<T, N> &y,
            Vec<T, N> &z) {
  VecSumCore<T, N>(alpha, get_data(x), beta, get_data(y), get_data(z));
}

template <typename T, int N, ADorder order, ADiffType ada, ADiffType adx,
          ADiffType adb, ADiffType ady>
class VecSumScaleExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using xtype = ADVecType<adx, order, Vec<T, N>>;
  using btype = ADScalarInputType<adb, order, T>;
  using ytype = ADVecType<ady, order, Vec<T, N>>;
  using ztype = ADVecType<ADiffType::ACTIVE, order, Vec<T, N>>;

  KOKKOS_FUNCTION
  VecSumScaleExpr(atype alpha, xtype &x, btype beta, ytype &y, ztype &z)
      : alpha(alpha), x(x), beta(beta), y(y), z(z) {}

  KOKKOS_FUNCTION void eval() {
    VecSumCore<T, N>(get_data(alpha), get_data(x), get_data(beta), get_data(y),
                     get_data(z));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (ada == ADiffType::ACTIVE && adx == ADiffType::ACTIVE &&
                  adb == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      VecSumCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(x),
                       get_data(beta), GetSeed<seed>::get_data(y),
                       GetSeed<seed>::get_data(z));
      VecSumCoreAdd<T, N>(GetSeed<seed>::get_data(alpha), get_data(x),
                          GetSeed<seed>::get_data(beta), get_data(y),
                          GetSeed<seed>::get_data(z));
    } else {
      VecZeroCore<T, N>(GetSeed<seed>::get_data(z));
      if constexpr (adx == ADiffType::ACTIVE) {
        VecAddCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(x),
                         GetSeed<seed>::get_data(z));
      }
      if constexpr (ady == ADiffType::ACTIVE) {
        VecAddCore<T, N>(get_data(beta), GetSeed<seed>::get_data(y),
                         GetSeed<seed>::get_data(z));
      }
      if constexpr (ada == ADiffType::ACTIVE) {
        VecAddCore<T, N>(GetSeed<seed>::get_data(alpha), get_data(x),
                         GetSeed<seed>::get_data(z));
      }
      if constexpr (adb == ADiffType::ACTIVE) {
        VecAddCore<T, N>(GetSeed<seed>::get_data(beta), get_data(y),
                         GetSeed<seed>::get_data(z));
      }
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(z),
                       GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(beta), GetSeed<seed>::get_data(z),
                       GetSeed<seed>::get_data(y));
    }
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, N>(GetSeed<seed>::get_data(z), get_data(x));
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) +=
          VecDotCore<T, N>(GetSeed<seed>::get_data(z), get_data(y));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    if constexpr (adx == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(alpha), GetSeed<seed>::get_data(z),
                       GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      VecAddCore<T, N>(get_data(beta), GetSeed<seed>::get_data(z),
                       GetSeed<seed>::get_data(y));
    }
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, N>(GetSeed<seed>::get_data(z), get_data(x));
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) +=
          VecDotCore<T, N>(GetSeed<seed>::get_data(z), get_data(y));
    }
    if constexpr (adx == ADiffType::ACTIVE && ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) += VecDotCore<T, N>(
          GetSeed<ADseed::b>::get_data(z), GetSeed<ADseed::p>::get_data(x));
      VecAddCore<T, N>(GetSeed<ADseed::p>::get_data(alpha),
                       GetSeed<ADseed::b>::get_data(z),
                       GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(beta) += VecDotCore<T, N>(
          GetSeed<ADseed::b>::get_data(z), GetSeed<ADseed::p>::get_data(y));
      VecAddCore<T, N>(GetSeed<ADseed::p>::get_data(beta),
                       GetSeed<ADseed::b>::get_data(z),
                       GetSeed<seed>::get_data(y));
    }
  }

  atype alpha;
  xtype &x;
  btype beta;
  ytype &y;
  ztype &z;
};

// First-order AD
template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(ADScalar<T> &alpha, ADVec<Vec<T, N>> &x,
                            ADScalar<T> &beta, ADVec<Vec<T, N>> &y,
                            ADVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE>(alpha, x, beta, y, z);
}

template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(const T alpha, ADVec<Vec<T, N>> &x, const T beta,
                            ADVec<Vec<T, N>> &y, ADVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::FIRST, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE>(alpha, x, beta, y, z);
}

template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(ADScalar<T> &alpha, const Vec<T, N> &x,
                            ADScalar<T> &beta, const Vec<T, N> &y,
                            ADVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::FIRST, ADiffType::ACTIVE,
                         ADiffType::PASSIVE, ADiffType::ACTIVE,
                         ADiffType::PASSIVE>(alpha, x, beta, y, z);
}

// Second-order AD
template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(A2DScalar<T> &alpha, A2DVec<Vec<T, N>> &x,
                            A2DScalar<T> &beta, A2DVec<Vec<T, N>> &y,
                            A2DVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                         ADiffType::ACTIVE, ADiffType::ACTIVE,
                         ADiffType::ACTIVE>(alpha, x, beta, y, z);
}

template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(const T alpha, A2DVec<Vec<T, N>> &x, const T beta,
                            A2DVec<Vec<T, N>> &y, A2DVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::SECOND, ADiffType::PASSIVE,
                         ADiffType::ACTIVE, ADiffType::PASSIVE,
                         ADiffType::ACTIVE>(alpha, x, beta, y, z);
}

template <typename T, int N>
KOKKOS_FUNCTION auto VecSum(A2DScalar<T> &alpha, const Vec<T, N> &x,
                            A2DScalar<T> &beta, const Vec<T, N> &y,
                            A2DVec<Vec<T, N>> &z) {
  return VecSumScaleExpr<T, N, ADorder::SECOND, ADiffType::ACTIVE,
                         ADiffType::PASSIVE, ADiffType::ACTIVE,
                         ADiffType::PASSIVE>(alpha, x, beta, y, z);
}

namespace Test {

template <typename T, int N>
class VecSumTest : public A2DTest<T, Vec<T, N>, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "VecSum<" << N << ">";
    return s.str();
  }

  // Evaluate the operation
  Output eval(const Input &x) {
    Vec<T, N> A, B, C;
    x.get_values(A, B);
    VecSum(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    Vec<T, N> A0, Ab, B0, Bb, C0, Cb;
    ADVec<Vec<T, N>> A(A0, Ab), B(B0, Bb), C(C0, Cb);
    x.get_values(A0, B0);
    auto op = VecSum(A, B, C);
    auto stack = MakeStack(op);
    seed.get_values(Cb);
    stack.reverse();
    g.set_values(Ab, Bb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    A2DVec<Vec<T, N>> A, B, C;
    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto op = VecSum(A, B, C);
    auto stack = MakeStack(op);
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

template <typename T, int N>
class VecSumScaleTest
    : public A2DTest<T, Vec<T, N>, T, Vec<T, N>, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, T, Vec<T, N>, T, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "VecSum<" << N << ">";
    return s.str();
  }

  // Evaluate the Vecrix-Vecrix product
  Output eval(const Input &x) {
    T alpha, beta;
    Vec<T, N> A, B, C;
    x.get_values(alpha, A, beta, B);
    VecSum(alpha, A, beta, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    ADScalar<T> alpha, beta;
    Vec<T, N> A0, Ab, B0, Bb, C0, Cb;
    ADVec<Vec<T, N>> A(A0, Ab), B(B0, Bb), C(C0, Cb);
    x.get_values(alpha.value, A0, beta.value, B0);
    auto op = VecSum(alpha, A, beta, B, C);
    auto stack = MakeStack(op);
    seed.get_values(Cb);
    stack.reverse();
    g.set_values(alpha.bvalue, Ab, beta.bvalue, Bb);
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    A2DScalar<T> alpha, beta;
    A2DVec<Vec<T, N>> A, B, C;
    x.get_values(alpha.value, A.value(), beta.value, B.value());
    p.get_values(alpha.pvalue, A.pvalue(), beta.pvalue, B.pvalue());
    auto op = VecSum(alpha, A, beta, B, C);
    auto stack = MakeStack(op);
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(alpha.hvalue, A.hvalue(), beta.hvalue, B.hvalue());
  }
};

bool VecSumTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecSumTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecSumTest<Tc, 5> test2;
  passed = passed && Run(test2, component, write_output);

  VecSumScaleTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_VEC_SUM_H