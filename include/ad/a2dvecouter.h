#ifndef A2D_VEC_OUTER_H
#define A2D_VEC_OUTER_H

#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dscalar.h"
#include "a2dvec.h"
#include "ad/core/a2dmatinnercore.h"
#include "ad/core/a2dmatveccore.h"
#include "ad/core/a2dvecoutercore.h"

namespace A2D {

template <typename T, int M, int N>
void VecOuter(T alpha, Vec<T, M>& x, Vec<T, N>& y, Mat<T, M, N>& A) {
  VecOuterCore<T, M, N>(alpha, get_data(x), get_data(y), get_data(A));
}

/*
  Compute the vector outer product

  A = alpha * x * y^{T}

  The forward mode derivative is

  dot{A} = (dot{alpha} * x * y^{T} +
        alpha * dot{x} * y^{T} +
        alpha * x * dot{y}^T

  Using the trace identity the reverse mode derivatives are

  bar{alpha} = x^{T} * bar{A} * y
  bar{x} = alpha * bar{A} * y
  bar{y} = alpha * bar{A}^{T} * x
*/
template <typename T, int M, int N, ADiffType adAlpha = ADiffType::ACTIVE,
          ADiffType adx = ADiffType::ACTIVE, ADiffType ady = ADiffType::ACTIVE,
          bool additive = false>
class ADVecOuterExpr {
 private:
  using alphaType = ADScalarType<adAlpha, T>;
  using xType = ADVecType<adx, Vec<T, M>>;
  using yType = ADVecType<ady, Vec<T, N>>;
  using AType = ADMat<Mat<T, M, N>>;

 public:
  KOKKOS_FUNCTION ADVecOuterExpr(alphaType& alpha, xType& x, yType& y, AType& A)
      : alpha(alpha), x(x), y(y), A(A) {
    VecOuterCore<T, M, N>(get_data(alpha), get_data(x), get_data(y),
                          get_data(A));
  }

  void forward() {
    if constexpr (adAlpha == ADiffType::ACTIVE) {
      constexpr bool additive = false;
      VecOuterCore<T, M, N, additive>(get_bdata(alpha), get_data(x),
                                      get_data(y), get_bdata(A));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      constexpr bool additive = (adAlpha == ADiffType::ACTIVE);
      VecOuterCore<T, M, N, additive>(get_data(alpha), get_bdata(x),
                                      get_data(y), get_bdata(A));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      constexpr bool additive =
          (adAlpha == ADiffType::ACTIVE or adx == ADiffType::ACTIVE);
      VecOuterCore<T, M, N, additive>(get_data(alpha), get_data(x),
                                      get_bdata(y), get_bdata(A));
    }
  }
  void reverse() {
    if constexpr (adAlpha == ADiffType::ACTIVE) {
      alpha.bvalue +=
          MatInnerCore<T, M, N>(get_bdata(A), get_data(x), get_data(y));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      MatVecCoreScale<T, M, N, MatOp::NORMAL, true>(
          get_data(alpha), get_bdata(A), get_data(y), get_bdata(x));
    }
    if constexpr (ady == ADiffType::ACTIVE) {
      MatVecCoreScale<T, M, N, MatOp::TRANSPOSE, true>(
          get_data(alpha), get_bdata(A), get_data(x), get_bdata(y));
    }
  }

 private:
  alphaType& alpha;
  xType& x;
  yType& y;
  ADMat<Mat<T, M, N>>& A;
};

}  // namespace A2D

#endif  // A2D_VEC_OUTER_H