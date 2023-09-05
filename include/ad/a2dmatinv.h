#ifndef A2D_MAT_INV_H
#define A2D_MAT_INV_H

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "ad/core/a2dgemmcore.h"
#include "ad/core/a2dmatinvcore.h"

namespace A2D {

/*
  Compute Ainv = A^{-1} for small matrices

  dot{Ainv} = - A^{-1} * dot{A} * A^{-1}

  hA = A^{-T} * Ap^{T} * A^{-T} * Ainvb * A^{-T} +
      A^{-T} * Ainvb * A^{-T} * Ap^{T} * A^{-T} =
    = - (A^{-T} * Ap^{T} * Ab + Ab * Ap^{T} * A^{-T})
*/

template <typename T, int N>
KOKKOS_FUNCTION void MatInv(const Mat<T, N, N>& A, Mat<T, N, N>& Ainv) {
  MatInvCore<T, N>(get_data(A), get_data(Ainv));
}

template <typename T, int N, ADorder order>
class MatInvExpr {
 private:
  using Atype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, N>>;

  static constexpr MatOp NORMAL = MatOp::NORMAL;
  static constexpr MatOp TRANSPOSE = MatOp::TRANSPOSE;

 public:
  KOKKOS_FUNCTION MatInvExpr(Atype& A, Atype& Ainv) : A(A), Ainv(Ainv) {}

  KOKKOS_FUNCTION void eval() { MatInvCore<T, N>(get_data(A), get_data(Ainv)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    T temp[N * N];
    MatMatMultCore<T, N, N, N, N, N, N, NORMAL, NORMAL>(
        get_data(Ainv), GetSeed<seed>::get_data(A), temp);
    MatMatMultScaleCore<T, N, N, N, N, N, N, NORMAL, NORMAL, true>(
        T(-1.0), temp, get_data(Ainv), GetSeed<seed>::get_data(Ainv));
  }

  KOKKOS_FUNCTION void reverse() {
    T temp[N * N];
    MatMatMultCore<T, N, N, N, N, N, N, TRANSPOSE, NORMAL>(
        get_data(Ainv), GetSeed<ADseed::b>::get_data(Ainv), temp);
    MatMatMultScaleCore<T, N, N, N, N, N, N, NORMAL, TRANSPOSE, true>(
        T(-1.0), temp, get_data(Ainv), GetSeed<ADseed::b>::get_data(A));
  }

  KOKKOS_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    T temp[N * N];

    // - A^{-T} * Ap^{T} * Ab
    MatMatMultCore<T, N, N, N, N, N, N, TRANSPOSE, TRANSPOSE>(
        get_data(Ainv), GetSeed<ADseed::p>::get_data(A), temp);
    MatMatMultScaleCore<T, N, N, N, N, N, N, NORMAL, NORMAL, true>(
        T(-1.0), temp, GetSeed<ADseed::b>::get_data(A),
        GetSeed<ADseed::h>::get_data(A));

    // - Ab * Ap^{T} * A^{-T}
    MatMatMultCore<T, N, N, N, N, N, N, NORMAL, TRANSPOSE>(
        GetSeed<ADseed::b>::get_data(A), GetSeed<ADseed::p>::get_data(A), temp);
    MatMatMultScaleCore<T, N, N, N, N, N, N, NORMAL, TRANSPOSE, true>(
        T(-1.0), temp, get_data(Ainv), GetSeed<ADseed::h>::get_data(A));

    // - A^{-T} * Ainvh * A^{-T}
    MatMatMultCore<T, N, N, N, N, N, N, TRANSPOSE, NORMAL>(
        get_data(Ainv), GetSeed<ADseed::h>::get_data(Ainv), temp);
    MatMatMultScaleCore<T, N, N, N, N, N, N, NORMAL, TRANSPOSE, true>(
        T(-1.0), temp, get_data(Ainv), GetSeed<ADseed::h>::get_data(A));
  }

  Atype& A;
  Atype& Ainv;
};

template <typename T, int N>
KOKKOS_FUNCTION auto MatInv(ADMat<Mat<T, N, N>>& A, ADMat<Mat<T, N, N>>& Ainv) {
  return MatInvExpr<T, N, ADorder::FIRST>(A, Ainv);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatInv(A2DMat<Mat<T, N, N>>& A,
                            A2DMat<Mat<T, N, N>>& Ainv) {
  return MatInvExpr<T, N, ADorder::SECOND>(A, Ainv);
}

namespace Test {

template <typename T, int N>
class MatInvTest : public A2DTest<T, Mat<T, N, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, Mat<T, N, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatInv<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> A;
    Mat<T, N, N> B;
    x.get_values(A);
    MatInv(A, B);
    return MakeVarTuple<T>(B);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    Mat<T, N, N> A0, Ab;
    Mat<T, N, N> B0, Bb;
    ADMat<Mat<T, N, N>> A(A0, Ab);
    ADMat<Mat<T, N, N>> B(B0, Bb);

    x.get_values(A0);
    auto op = MatInv(A, B);
    auto stack = MakeStack(op);
    seed.get_values(Bb);
    stack.reverse();
    g.set_values(Ab);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DMat<Mat<T, N, N>> A;
    A2DMat<Mat<T, N, N>> B;
    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto op = MatInv(A, B);
    auto stack = MakeStack(op);

    seed.get_values(B.bvalue());
    hval.get_values(B.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

bool MatInvTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatInvTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatInvTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_MAT_INV_H