#ifndef A2D_VEC_OUTER_H
#define A2D_VEC_OUTER_H

#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dscalar.h"
#include "a2dvec.h"
#include "ad/core/a2dmatveccore.h"
#include "ad/core/a2dveccore.h"

namespace A2D {

template <typename T, int M, int N>
KOKKOS_FUNCTION void VecOuter(const Vec<T, M>& x, const Vec<T, N>& y,
                              Mat<T, M, N>& A) {
  VecOuterCore<T, M, N>(get_data(x), get_data(y), get_data(A));
}

template <typename T, int M, int N>
KOKKOS_FUNCTION void VecOuter(const T alpha, const Vec<T, M>& x,
                              const Vec<T, N>& y, Mat<T, M, N>& A) {
  VecOuterCore<T, M, N>(alpha, get_data(x), get_data(y), get_data(A));
}

/*
  Compute the vector outer product

  A = alpha * x * y^{T}

  where alpha is a numeric constant scalar parameter.

  The forward mode derivative is

  dot{A} = alpha * (dot{x} * y^{T} + x * dot{y}^T)

  Using the trace identity the reverse mode derivatives are

  bar{x} = alpha * bar{A} * y
  bar{y} = alpha * bar{A}^{T} * x
*/
template <typename T, int M, int N, ADorder order, ADiffType adx, ADiffType ady>
class VecOuterExpr {
 private:
  using xtype = ADVecType<adx, order, Vec<T, M>>;
  using ytype = ADVecType<ady, order, Vec<T, N>>;
  using Atype = ADMatType<ADiffType::ACTIVE, order, Mat<T, M, N>>;

 public:
  KOKKOS_FUNCTION VecOuterExpr(const T alpha, xtype& x, ytype& y, Atype& A)
      : alpha(alpha), x(x), y(y), A(A) {}

  KOKKOS_FUNCTION void eval() {
    VecOuterCore<T, M, N>(alpha, get_data(x), get_data(y), get_data(A));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      constexpr bool additive = true;
      VecOuterCore<T, M, N>(alpha, GetSeed<seed>::get_data(x), get_data(y),
                            GetSeed<seed>::get_data(A));
      VecOuterCore<T, M, N, additive>(alpha, get_data(x),
                                      GetSeed<seed>::get_data(y),
                                      GetSeed<seed>::get_data(A));
    } else if constexpr (adx == ADiffType::ACTIVE) {
      VecOuterCore<T, M, N>(alpha, GetSeed<seed>::get_data(x), get_data(y),
                            GetSeed<seed>::get_data(A));
    } else if constexpr (ady == ADiffType::ACTIVE) {
      VecOuterCore<T, M, N>(alpha, get_data(x), GetSeed<seed>::get_data(y),
                            GetSeed<seed>::get_data(A));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adx == ADiffType::ACTIVE) {
      constexpr bool additive = true;
      MatVecCoreScale<T, M, N, MatOp::NORMAL, additive>(
          alpha, GetSeed<seed>::get_data(A), get_data(y),
          GetSeed<seed>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      constexpr bool additive = true;
      MatVecCoreScale<T, M, N, MatOp::TRANSPOSE, additive>(
          alpha, GetSeed<seed>::get_data(A), get_data(x),
          GetSeed<seed>::get_data(y));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    constexpr bool additive = true;
    if constexpr (adx == ADiffType::ACTIVE) {
      MatVecCoreScale<T, M, N, MatOp::NORMAL, additive>(
          alpha, GetSeed<ADseed::h>::get_data(A), get_data(y),
          GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      MatVecCoreScale<T, M, N, MatOp::TRANSPOSE, additive>(
          alpha, GetSeed<ADseed::h>::get_data(A), get_data(x),
          GetSeed<ADseed::h>::get_data(y));
    }
    if constexpr (adx == ADiffType::ACTIVE && ady == ADiffType::ACTIVE) {
      MatVecCoreScale<T, M, N, MatOp::NORMAL, additive>(
          alpha, GetSeed<ADseed::b>::get_data(A),
          GetSeed<ADseed::p>::get_data(y), GetSeed<ADseed::h>::get_data(x));
      MatVecCoreScale<T, M, N, MatOp::TRANSPOSE, additive>(
          alpha, GetSeed<ADseed::b>::get_data(A),
          GetSeed<ADseed::p>::get_data(x), GetSeed<ADseed::h>::get_data(y));
    }
  }

 private:
  const T alpha;
  xtype& x;
  ytype& y;
  Atype& A;
};

template <typename T, int M, int N>
KOKKOS_FUNCTION auto VecOuter(ADVec<Vec<T, M>>& x, ADVec<Vec<T, N>>& y,
                              ADMat<Mat<T, M, N>>& A) {
  return VecOuterExpr<T, M, N, ADorder::FIRST, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(T(1.0), x, y, A);
}

template <typename T, int M, int N>
KOKKOS_FUNCTION auto VecOuter(const T alpha, ADVec<Vec<T, M>>& x,
                              ADVec<Vec<T, N>>& y, ADMat<Mat<T, M, N>>& A) {
  return VecOuterExpr<T, M, N, ADorder::FIRST, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(alpha, x, y, A);
}

template <typename T, int M, int N>
KOKKOS_FUNCTION auto VecOuter(A2DVec<Vec<T, M>>& x, A2DVec<Vec<T, N>>& y,
                              A2DMat<Mat<T, M, N>>& A) {
  return VecOuterExpr<T, M, N, ADorder::SECOND, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(T(1.0), x, y, A);
}

template <typename T, int M, int N>
KOKKOS_FUNCTION auto VecOuter(const T alpha, A2DVec<Vec<T, M>>& x,
                              A2DVec<Vec<T, N>>& y, A2DMat<Mat<T, M, N>>& A) {
  return VecOuterExpr<T, M, N, ADorder::SECOND, ADiffType::ACTIVE,
                      ADiffType::ACTIVE>(alpha, x, y, A);
}

namespace Test {

template <typename T, int N, int M>
class VecOuterTest : public A2DTest<T, Mat<T, N, M>, Vec<T, N>, Vec<T, M>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() { return "VecOuter"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    T alpha(0.3157);
    Vec<T, N> x;
    Vec<T, M> y;
    Mat<T, N, M> A;
    X.get_values(x, y);
    VecOuter(alpha, x, y, A);
    return MakeVarTuple<T>(A);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    T alpha(0.3157);
    Vec<T, N> x0, xb;
    Vec<T, M> y0, yb;
    Mat<T, N, M> A0, Ab;
    ADVec<Vec<T, N>> x(x0, xb);
    ADVec<Vec<T, M>> y(y0, yb);
    ADMat<Mat<T, N, M>> A(A0, Ab);
    X.get_values(x0, y0);
    auto op = VecOuter(alpha, x, y, A);
    auto stack = MakeStack(op);
    seed.get_values(Ab);
    stack.reverse();
    g.set_values(xb, yb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    T alpha(0.3157);
    A2DVec<Vec<T, N>> x;
    A2DVec<Vec<T, M>> y;
    A2DMat<Mat<T, N, M>> A;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto op = VecOuter(alpha, x, y, A);
    auto stack = MakeStack(op);
    seed.get_values(A.bvalue());
    hval.get_values(A.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

bool VecOuterTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  VecOuterTest<Tc, 3, 5> test1;
  passed = passed && Run(test1, component, write_output);
  VecOuterTest<Tc, 6, 4> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_VEC_OUTER_H