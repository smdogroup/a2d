#ifndef A2D_VEC_OUTER_H
#define A2D_VEC_OUTER_H

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "ad/core/a2dmatvecinnercore.h"
#include "ad/core/a2dvecoutercore.h"
#include "ad/corea/2dmatveccore.h"

namespace A2D {

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

template <typename T, int M, int N, VarAd adAlpha = VarAD::AD,
          VarAD adx = VarAd::AD, VarAd ady = VarAD::AD, bool additive = false,
          bool scale = false>
class ADVecOuterExpr {
 private:
  using alphaType =
      typename std::conditional<adAlpha == VarAd::AD, ADScalar<T>, T>::type;
  using xType = typename std::conditional<adx == VarAd::AD, ADVec<Vec<T, M>>,
                                          Vec<T, M>>::type;
  using yType = typename std::conditional<ady == VarAd::AD, ADVec<Vec<T, N>>,
                                          Vec<T, N>>::type;

 public:
  A2D_INLINE_FUNCTION
  ADVecOuterExpr(alphaType& alpha, xType& x, yType& y, ADMat<Mat<T, M, N>>& A)
      : alpha(alpha), x(x), y(y), A(A) {
    VecOuterCore<T, M, N>(get_data(alpha), get_data(x), get_data(y),
                          get_data(A));
  }

  void forward() {
    if constexpr (adAlpha == VarAD::AD) {
      VecOuterCore<T, M, N>(get_bdata(alpha), get_data(x), get_data(y),
                            get_bdata(A));
    }
    if constexpr (xad == VarAD::AD) {
      VecOuterCore<T, M, N>(get_data(alpha), get_bdata(x), get_data(y),
                            get_bdata(A));
    }
    if constexpr (yad == VarAD::AD) {
      VecOuterCore<T, M, N>(get_data(alpha), get_data(x), get_bdata(y),
                            get_bdata(A));
    }
  }
  void reverse() {
    if constexpr (adAlpha == VarAD::AD) {
      alpha.bvalue +=
          MatInnerCore<T, M, N>(get_bdata(A), get_data(x), get_data(y));
    }
    if constexpr (xad == VarAD::AD) {
      MatVecCoreScale<T, M, N, MatOp::NORMAL, true>(
          get_data(alpha), get_bdata(A), get_data(y), get_bdata(x));
    }
    if constexpr (yad == VarAD::AD) {
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

template <typename T, int M, int N>
void VecOuter(T alpha, Vec<T, M>& x, Vec<T, N>& y, Mat<T, M, N>& A) {
  VecOuterCore<T, M, N>(alpha, get_data(x), get_data(y), get_data(A));
}

}  // namespace A2D

#endif  // A2D_VEC_OUTER_H