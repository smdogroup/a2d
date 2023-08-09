#ifndef A2D_GEMM_H
#define A2D_GEMM_H

#include <type_traits>

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "ad/core/a2dgemmcore.h"

namespace A2D {

enum class VarAd { CONST, AD };

template <typename T, int m, int n>
A2D_INLINE_FUNCTION T* get_mat_data(Mat<T, m, n> mat) {
  return mat.data();
}

template <typename T, int m, int n>
A2D_INLINE_FUNCTION T* get_mat_data(ADMat<Mat<T, m, n>> mat) {
  return mat.value().data();
}

// template <typename T, int N, int M, int K, int L, int P, int Q,
//           MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
// A2D_INLINE_FUNCTION void MatMatMult(const Mat<T, N, M>& A,
//                                     const Mat<T, K, L>& B,
//                                     Mat<T, P, Q>& C) noexcept {
//   MatMatMultCore<T, N, M, K, L, P, Q, opA, opB>(A.data(), B.data(),
//   C.data());
// }

template <typename T, int N, int M, int K, int L, int P, int Q,
          MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL,
          VarAd adA = VarAd::AD, VarAd adB = VarAd::AD>
class ADMatMatMultExpr {
 private:
  static constexpr MatOp not_opA = negate_op<opA>::value;
  static constexpr MatOp not_opB = negate_op<opB>::value;
  // static constexpr VarAd adC = any_VarAd<adA, adB>::value;

  using Atype = typename std::conditional<adA == VarAd::AD, ADMat<Mat<T, N, M>>,
                                          Mat<T, N, M>>::type;
  using Btype = typename std::conditional<adB == VarAd::AD, ADMat<Mat<T, K, L>>,
                                          Mat<T, K, L>>::type;
  using Ctype =
      typename std::conditional<adA == VarAd::AD or adB == VarAd::AD,
                                ADMat<Mat<T, P, Q>>, Mat<T, P, Q>>::type;

 public:
  A2D_INLINE_FUNCTION ADMatMatMultExpr(Atype& A, Btype& B, Ctype& C)
      : A(A), B(B), C(C) {
    MatMatMultCore<T, N, M, K, L, P, Q, opA, opB>(
        get_mat_data(A), get_mat_data(B), get_mat_data(C));
  }

  A2D_INLINE_FUNCTION void forward() {
    if constexpr (adA == VarAd::AD) {
      MatMatMultCore<T, N, M, K, L, P, Q, opA, opB, false>(
          A.bvalue().data(), get_mat_data(B), C.bvalue().data());
    }
    if constexpr (adB == VarAd::AD) {
      MatMatMultCore<T, N, M, K, L, P, Q, opA, opB, adA == VarAd::AD>(
          get_mat_data(A), B.bvalue().data(), C.bvalue().data());
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    if constexpr (adA == VarAd::AD) {
      MatMatMultCore<T, N, M, K, L, P, Q, MatOp::NORMAL, not_opB, opA, true>(
          C.bvalue().data(), get_mat_data(B), A.bvalue().data());
    }
    if constexpr (adB == VarAd::AD) {
      MatMatMultCore<T, N, M, K, L, P, Q, not_opA, MatOp::NORMAL, opB, true>(
          get_mat_data(A), C.bvalue().data(), B.bvalue().data());
    }
  }

 private:
  Atype& A;
  Btype& B;
  Ctype& C;
};

template <typename T, int N, int M, int K, int L, int P, int Q,
          MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
auto MatMatMult(ADMat<Mat<T, N, M>>& A, ADMat<Mat<T, K, L>>& B,
                ADMat<Mat<T, P, Q>>& C) {
  return ADMatMatMultExpr<T, N, M, K, L, P, Q, opA, opB, VarAd::AD, VarAd::AD>(
      A, B, C);
}

template <typename T, int N, int M, int K, int L, int P, int Q,
          MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
auto MatMatMult(Mat<T, N, M>& A, ADMat<Mat<T, K, L>>& B,
                ADMat<Mat<T, P, Q>>& C) {
  return ADMatMatMultExpr<T, N, M, K, L, P, Q, opA, opB, VarAd::CONST,
                          VarAd::AD>(A, B, C);
}

template <typename T, int N, int M, int K, int L, int P, int Q,
          MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
auto MatMatMult(ADMat<Mat<T, N, M>>& A, Mat<T, K, L>& B,
                ADMat<Mat<T, P, Q>>& C) {
  return ADMatMatMultExpr<T, N, M, K, L, P, Q, opA, opB, VarAd::AD,
                          VarAd::CONST>(A, B, C);
}

template <typename T, int N, int M, int K, int L, int P, int Q,
          MatOp opA = MatOp::NORMAL, MatOp opB = MatOp::NORMAL>
auto MatMatMult(Mat<T, N, M>& A, Mat<T, K, L>& B, Mat<T, P, Q>& C) {
  return ADMatMatMultExpr<T, N, M, K, L, P, Q, opA, opB, VarAd::CONST,
                          VarAd::CONST>(A, B, C);
}

}  // namespace A2D
#endif  // A2D_GEMM_H