#ifndef A2D_MAT_TRACE_H
#define A2D_MAT_TRACE_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dscalar.h"

namespace A2D {

template <typename T, int M>
A2D_INLINE_FUNCTION T MatTraceCore(const T* A) {
  T trace = T(0.0);
  for (int i = 0; i < M; i++) {
    trace += A[0];
    A += M + 1;
  }
  return trace;
}

template <typename T, int M>
A2D_INLINE_FUNCTION void MatAddDiagCore(const T diag, T* A) {
  for (int i = 0; i < M; i++) {
    A[0] += diag;
    A += M + 1;
  }
}

template <typename T, int M>
A2D_INLINE_FUNCTION void MatTrace(Mat<T, M, M>& A, T& trace) {
  trace = MatTraceCore<T, M>(get_data(A));
}

template <typename T, int M, ADorder order, ADiffType adA>
class MatTraceExpr {
 private:
  using Atype = ADMatType<adA, order, Mat<T, M, M>>;
  using ScalarType = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  A2D_INLINE_FUNCTION MatTraceExpr(Atype& A, ScalarType& tr) : A(A), tr(tr) {
    get_data(tr) = MatTraceCore<T, M>(get_data(A));
  }

  template <ADorder forder>
  A2D_INLINE_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(tr) =
          MatTraceCore<T, M>(GetSeed<seed>::get_data(A));
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      MatAddDiagCore<T, M>(GetSeed<ADseed::b>::get_data(tr),
                           GetSeed<ADseed::b>::get_data(A));
    }
  }

  // A2D_INLINE_FUNCTION void hreverse() {
  // }

 private:
  Atype& A;
  ScalarType& tr;
};

template <typename T, int M>
A2D_INLINE_FUNCTION auto MatTrace(ADMat<Mat<T, M, M>>& A, ADScalar<T>& tr) {
  return MatTraceExpr<T, M, ADorder::FIRST, ADiffType::ACTIVE>(A, tr);
}
template <typename T, int M>
A2D_INLINE_FUNCTION auto MatTrace(A2DMat<Mat<T, M, M>>& A, A2DScalar<T>& tr) {
  return MatTraceExpr<T, M, ADorder::SECOND, ADiffType::ACTIVE>(A, tr);
}

}  // namespace A2D

#endif  // A2D_MAT_TRACE_H