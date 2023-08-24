#ifndef A2D_GEMM_H
#define A2D_GEMM_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "ad/core/a2dgemmcore.h"

namespace A2D {

/*
  Compute C = opA(A) * opB(B)

  Forward mode-derivative:

  dot{C} =  opA(dot{A}) * opB(B) + opA(A) * opB(dot{B})

  Derivation of the reverse-mode derivative:

  tr(bar{C}^{T} * dot{C})
  = tr(bar{C}^{T} * opA(dot{A}) * opB(B)) +
    tr(bar{C}^{T} * opA(A) * opB(dot{B}))
  = tr(opB(B) * bar{C}^{T} * opA(dot{A})) +
    tr(bar{C}^{T} * opA(A) * opB(dot{B}))

  In general:

  if opA == NORMAL then
    tr(bar{B}^{T} * dot{A})
    ==> bar{A} = bar{B}
  if opA == TRANSPOSE
    tr(bar{B}^{T} * dot{A^{T}}) = tr(dot{A} * bar{B}) = tr(bar{B} * dot{A})
    ==> bar{A} = bar{B}^{T}

  Therefore:

  if opA == NORMAL then:
    bar{A} = (opB(B) * bar{C}^{T})^{T} = bar{C} * not_opB(B)
  if opA == TRANSPOSE then:
    bar{A} = opB(B) * bar{C}^{T}

  if opB == NORMAL
    bar{B} = (bar{C}^{T} * opA(A))^{T} = not_opA(A) * bar{C}
  if opB == TRANSPOSE then:
    bar{B} = bar{C}^{T} * opA(A)
*/

// compute C = op(A) * op(B) and returns nothing, where A and B are all
// passive variables
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION void MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    Mat<T, P, Q>& C) {
  MatMatMultCore<T, N, M, K, L, P, Q, MatOp::NORMAL, MatOp::NORMAL>(
      get_data(A), get_data(B), get_data(C));
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION void MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    Mat<T, P, Q>& C) {
  MatMatMultCore<T, N, M, K, L, P, Q, opA, opB>(get_data(A), get_data(B),
                                                get_data(C));
}

template <typename T, int N, int M, int K, int L, int P, int Q, ADorder order,
          MatOp opA, MatOp opB, ADiffType adA, ADiffType adB>
class MatMatMultExpr {
 private:
  static constexpr MatOp not_opA =
      conditional_value<MatOp, opA == MatOp::NORMAL, MatOp::TRANSPOSE,
                        MatOp::NORMAL>::value;
  static constexpr MatOp not_opB =
      conditional_value<MatOp, opB == MatOp::NORMAL, MatOp::TRANSPOSE,
                        MatOp::NORMAL>::value;

  using Atype = ADMatType<adA, order, Mat<T, N, M>>;
  using Btype = ADMatType<adB, order, Mat<T, K, L>>;
  using Ctype =
      ADMatType<ADiffType::ACTIVE, order, Mat<T, P, Q>>;  // C is always active

 public:
  A2D_INLINE_FUNCTION MatMatMultExpr(Atype& A, Btype& B, Ctype& C)
      : A(A), B(B), C(C) {
    MatMatMultCore<T, N, M, K, L, P, Q, opA, opB>(get_data(A), get_data(B),
                                                  get_data(C));
  }

  template <ADorder forder>
  A2D_INLINE_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE) {
      MatMatMultCore<T, N, M, K, L, P, Q, opA, opB>(
          GetSeed<seed>::get_data(A), get_data(B), GetSeed<seed>::get_data(C));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      MatMatMultCore<T, N, M, K, L, P, Q, opA, opB, adA == ADiffType::ACTIVE>(
          get_data(A), GetSeed<seed>::get_data(B), GetSeed<seed>::get_data(C));
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      if constexpr (opA == MatOp::NORMAL) {
        // bar{A} += bar{C} * not_opB(B)
        MatMatMultCore<T, P, Q, K, L, N, M, MatOp::NORMAL, not_opB, true>(
            GetSeed<ADseed::b>::get_data(C), get_data(B),
            GetSeed<ADseed::b>::get_data(A));
      } else {
        // bar{A} += opB(B) * bar{C}^{T}
        MatMatMultCore<T, K, L, P, Q, N, M, opB, MatOp::TRANSPOSE, true>(
            get_data(B), GetSeed<ADseed::b>::get_data(C),
            GetSeed<ADseed::b>::get_data(A));
      }
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      if constexpr (opB == MatOp::NORMAL) {
        // bar{B} += not_opA(A) * bar{C}
        MatMatMultCore<T, N, M, P, Q, K, L, not_opA, MatOp::NORMAL, true>(
            get_data(A), GetSeed<ADseed::b>::get_data(C),
            GetSeed<ADseed::b>::get_data(B));
      } else {
        // bar{B} += bar{C}^{T} * opA(A)
        MatMatMultCore<T, P, Q, N, M, K, L, MatOp::TRANSPOSE, opA, true>(
            GetSeed<ADseed::b>::get_data(C), get_data(A),
            GetSeed<ADseed::b>::get_data(B));
      }
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    if constexpr (adA == ADiffType::ACTIVE) {
      if constexpr (opA == MatOp::NORMAL) {
        MatMatMultCore<T, P, Q, K, L, N, M, MatOp::NORMAL, not_opB, true>(
            GetSeed<ADseed::h>::get_data(C), get_data(B),
            GetSeed<ADseed::h>::get_data(A));
      } else {
        MatMatMultCore<T, K, L, P, Q, N, M, opB, MatOp::TRANSPOSE, true>(
            get_data(B), GetSeed<ADseed::h>::get_data(C),
            GetSeed<ADseed::h>::get_data(A));
      }
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      if constexpr (opB == MatOp::NORMAL) {
        MatMatMultCore<T, N, M, P, Q, K, L, not_opA, MatOp::NORMAL, true>(
            get_data(A), GetSeed<ADseed::h>::get_data(C),
            GetSeed<ADseed::h>::get_data(B));
      } else {
        MatMatMultCore<T, P, Q, N, M, K, L, MatOp::TRANSPOSE, opA, true>(
            GetSeed<ADseed::h>::get_data(C), get_data(A),
            GetSeed<ADseed::h>::get_data(B));
      }
    }
    if constexpr (adA == ADiffType::ACTIVE and adB == ADiffType::ACTIVE) {
      if constexpr (opA == MatOp::NORMAL) {
        MatMatMultCore<T, P, Q, K, L, N, M, MatOp::NORMAL, not_opB, true>(
            GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::p>::get_data(B),
            GetSeed<ADseed::h>::get_data(A));

      } else {
        MatMatMultCore<T, K, L, P, Q, N, M, opB, MatOp::TRANSPOSE, true>(
            GetSeed<ADseed::p>::get_data(B), GetSeed<ADseed::b>::get_data(C),
            GetSeed<ADseed::h>::get_data(A));
      }

      if constexpr (opB == MatOp::NORMAL) {
        MatMatMultCore<T, N, M, P, Q, K, L, not_opA, MatOp::NORMAL, true>(
            GetSeed<ADseed::p>::get_data(A), GetSeed<ADseed::b>::get_data(C),
            GetSeed<ADseed::h>::get_data(B));

      } else {
        MatMatMultCore<T, P, Q, N, M, K, L, MatOp::TRANSPOSE, opA, true>(
            GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::p>::get_data(A),
            GetSeed<ADseed::h>::get_data(B));
      }
    }
  }

 private:
  Atype& A;
  Btype& B;
  Ctype& C;
};

// compute C = op(A) * op(B) and return an expression, where A and B are all
// active variables
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(ADMat<Mat<T, N, M>>& A,
                                    ADMat<Mat<T, K, L>>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      A, B, C);
}
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(A2DMat<Mat<T, N, M>>& A,
                                    A2DMat<Mat<T, K, L>>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(ADMat<Mat<T, N, M>>& A,
                                    ADMat<Mat<T, K, L>>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, opA, opB,
                        ADiffType::ACTIVE, ADiffType::ACTIVE>(A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(A2DMat<Mat<T, N, M>>& A,
                                    A2DMat<Mat<T, K, L>>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, opA, opB,
                        ADiffType::ACTIVE, ADiffType::ACTIVE>(A, B, C);
}

// compute C = op(A) * op(B) and return an expression, where A is passive, B is
// active variables
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, ADMat<Mat<T, K, L>>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      A, B, C);
}
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, A2DMat<Mat<T, K, L>>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, ADMat<Mat<T, K, L>>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, opA, opB,
                        ADiffType::PASSIVE, ADiffType::ACTIVE>(A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, A2DMat<Mat<T, K, L>>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, opA, opB,
                        ADiffType::PASSIVE, ADiffType::ACTIVE>(A, B, C);
}

// compute C = op(A) * op(B) and return an expression, where A is active, B is
// passive variables
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(ADMat<Mat<T, N, M>>& A, Mat<T, K, L>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      A, B, C);
}
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(A2DMat<Mat<T, N, M>>& A, Mat<T, K, L>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(ADMat<Mat<T, N, M>>& A, Mat<T, K, L>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, opA, opB,
                        ADiffType::ACTIVE, ADiffType::PASSIVE>(A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(A2DMat<Mat<T, N, M>>& A, Mat<T, K, L>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, opA, opB,
                        ADiffType::ACTIVE, ADiffType::PASSIVE>(A, B, C);
}

// compute C = op(A) * op(B) and return an expression, where A and B are both
// passive variables, but C is active variable
// Note: this still returns an expression, but it's empty, i.e. the expression
// doesn't have meaningful implementation of forward(), reverse(), etc.
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::PASSIVE, ADiffType::PASSIVE>(
      A, B, C);
}
template <typename T, int N, int M, int K, int L, int P, int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, MatOp::NORMAL,
                        MatOp::NORMAL, ADiffType::PASSIVE, ADiffType::PASSIVE>(
      A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    ADMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::FIRST, opA, opB,
                        ADiffType::PASSIVE, ADiffType::PASSIVE>(A, B, C);
}
template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
A2D_INLINE_FUNCTION auto MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B,
                                    A2DMat<Mat<T, P, Q>>& C) {
  return MatMatMultExpr<T, N, M, K, L, P, Q, ADorder::SECOND, opA, opB,
                        ADiffType::PASSIVE, ADiffType::PASSIVE>(A, B, C);
}

}  // namespace A2D

#endif  // A2D_GEMM_H