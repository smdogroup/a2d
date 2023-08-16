#ifndef A2D_MAT_DET_H
#define A2D_MAT_DET_H

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "ad/core/a2dmatdetcore.h"

namespace A2D {

template <typename T, int N>
A2D_INLINE_FUNCTION void MatDet(Mat<T, N, N>& A, T& det) {
  det = MatDetCore<T, N>(get_data(A));
}

// template <typename T, int N>
// A2D_INLINE_FUNCTION void MatDet(SymMat<T, N>& A, T& det) {
//   det = SymMatDetCore<T, N>(get_data(A));
// }

template <typename T, int N, ADorder order>
class MatDetExpr {
 private:
  using Atype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, N>>;
  using ScalarType = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  A2D_INLINE_FUNCTION MatDetExpr(Atype& A, ScalarType& det) : A(A), det(det) {
    get_data(det) = MatDetCore<T, N>(get_data(A));
  }

  template <ADorder forder>
  A2D_INLINE_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(det) =
        MatDetForwardCore<T, N>(get_data(A), GetSeed<seed>::get_data(A));
  }

  A2D_INLINE_FUNCTION void reverse() {
    MatDetReverseCore<T, N>(GetSeed<ADseed::b>::get_data(det), get_data(A),
                            GetSeed<ADseed::b>::get_data(A));
  }

  A2D_INLINE_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    MatDetHReverseCore<T, N>(GetSeed<ADseed::b>::get_data(det),
                             GetSeed<ADseed::h>::get_data(det), get_data(A),
                             GetSeed<ADseed::p>::get_data(A),
                             GetSeed<ADseed::h>::get_data(A));
  }

  Atype& A;
  ScalarType& det;
};

template <typename T, int N>
A2D_INLINE_FUNCTION auto MatDet(ADMat<Mat<T, N, N>>& A, ADScalar<T>& det) {
  return MatDetExpr<T, N, ADorder::FIRST>(A, det);
}

template <typename T, int N>
A2D_INLINE_FUNCTION auto MatDet(A2DMat<Mat<T, N, N>>& A, A2DScalar<T>& det) {
  return MatDetExpr<T, N, ADorder::SECOND>(A, det);
}

}  // namespace A2D

#endif  //  A2D_MAT_DET_H