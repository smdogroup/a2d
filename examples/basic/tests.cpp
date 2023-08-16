#include "a2denum.h"
#include "ad/a2dgemm.h"
#include "ad/a2dmatdet.h"
#include "ad/a2dmatinv.h"
#include "ad/a2dmattrace.h"
#include "ad/a2dstack.h"
#include "ad/a2dtest.h"

template <typename T, int N>
class TestMatDet : public A2D::Test::A2DTest<T, A2D::Mat<T, N, N>> {
 public:
  using TupleType = A2D::VarTuple<T, A2D::Mat<T, N, N>>;

  TestMatDet() {}

  T eval(const TupleType& x) {
    A2D::Mat<T, N, N> A;
    x.get_values(A);

    T det;
    A2D::MatDet(A, det);
    return det;
  }

  void deriv(const TupleType& x, TupleType& g) {
    A2D::Mat<T, N, N> A0, Ab;
    A2D::ADMat<A2D::Mat<T, N, N>> A(A0, Ab);
    x.get_values(A0);

    A2D::ADScalar<T> det;
    auto op1 = A2D::MatDet(A, det);
    det.bvalue = 1.0;
    op1.reverse();

    g.set_values(Ab);
  }

  void hprod(const TupleType& x, const TupleType& p, TupleType& h) {
    A2D::A2DScalar<T> det;
    A2D::A2DMat<A2D::Mat<T, N, N>> A;

    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto op1 = A2D::MatDet(A, det);
    auto stack = A2D::MakeStack(op1);

    det.bvalue = 1.0;
    stack.reverse();
    stack.hforward();
    stack.hreverse();

    h.set_values(A.hvalue());
  }
};

template <typename T, int N>
class TestMatInv : public A2D::Test::A2DTest<T, A2D::Mat<T, N, N>> {
 public:
  using TupleType = A2D::VarTuple<T, A2D::Mat<T, N, N>>;

  TestMatInv() {}

  T eval(const TupleType& x) {
    A2D::Mat<T, N, N> A, B;
    x.get_values(A);

    T det;
    A2D::MatInv(A, B);
    A2D::MatDet(B, det);

    return det;
  }

  void deriv(const TupleType& x, TupleType& g) {
    A2D::Mat<T, N, N> A0, Ab, B0, Bb;
    A2D::ADMat<A2D::Mat<T, N, N>> A(A0, Ab), B(B0, Bb);
    A2D::ADScalar<T> det;
    x.get_values(A0);

    auto op1 = A2D::MatInv(A, B);
    auto op2 = A2D::MatDet(B, det);
    auto stack = A2D::MakeStack(op1, op2);
    det.bvalue = 1.0;

    stack.reverse();

    g.set_values(Ab);
  }

  void hprod(const TupleType& x, const TupleType& p, TupleType& h) {
    A2D::A2DScalar<T> det;
    A2D::A2DMat<A2D::Mat<T, N, N>> A, B;

    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto op1 = A2D::MatInv(A, B);
    auto op2 = A2D::MatDet(B, det);
    auto stack = A2D::MakeStack(op1, op2);
    det.bvalue = 1.0;

    stack.reverse();
    stack.hforward();
    stack.hreverse();

    h.set_values(A.hvalue());
  }
};

template <typename T, int M, int N, int K>
class TestMatMult
    : public A2D::Test::A2DTest<T, A2D::Mat<T, M, N>, A2D::Mat<T, N, K>> {
 public:
  using TupleType = A2D::VarTuple<T, A2D::Mat<T, M, N>, A2D::Mat<T, N, K>>;

  TestMatMult() {}

  static constexpr A2D::MatOp NORMAL = A2D::MatOp::NORMAL;
  static constexpr A2D::MatOp TRANSPOSE = A2D::MatOp::TRANSPOSE;

  T eval(const TupleType& x) {
    A2D::Mat<T, M, N> A;
    A2D::Mat<T, N, K> B;
    A2D::Mat<T, M, K> C;
    A2D::Mat<T, K, K> D;
    A2D::Mat<T, K, M> E;
    A2D::Mat<T, K, K> F;

    x.get_values(A, B);

    // C = A * B
    A2D::MatMatMult(A, B, C);
    // D = C^{T} * C
    A2D::MatMatMult<TRANSPOSE, NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    A2D::MatMatMult<TRANSPOSE, TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    A2D::MatMatMult<NORMAL, TRANSPOSE>(E, E, F);

    T trace;
    A2D::MatTrace(D, trace);

    return trace;
  }

  void deriv(const TupleType& x, TupleType& g) {
    A2D::Mat<T, M, N> A0, Ab;
    A2D::Mat<T, N, K> B0, Bb;
    A2D::Mat<T, M, K> C0, Cb;
    A2D::Mat<T, K, K> D0, Db;
    A2D::Mat<T, K, M> E0, Eb;
    A2D::Mat<T, K, K> F0, Fb;

    A2D::ADMat<A2D::Mat<T, M, N>> A(A0, Ab);
    A2D::ADMat<A2D::Mat<T, N, K>> B(B0, Bb);
    A2D::ADMat<A2D::Mat<T, M, K>> C(C0, Cb);
    A2D::ADMat<A2D::Mat<T, K, K>> D(D0, Db);
    A2D::ADMat<A2D::Mat<T, K, M>> E(E0, Eb);
    A2D::ADMat<A2D::Mat<T, K, K>> F(F0, Fb);

    x.get_values(A0, B0);

    // C = A * B
    auto mult1 = A2D::MatMatMult(A, B, C);
    // D = C^{T} * C
    auto mult2 = A2D::MatMatMult<TRANSPOSE, NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    auto mult3 = A2D::MatMatMult<TRANSPOSE, TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    auto mult4 = A2D::MatMatMult<NORMAL, TRANSPOSE>(E, E, F);

    A2D::ADScalar<T> trace;
    auto tr1 = A2D::MatTrace(D, trace);

    auto stack = A2D::MakeStack(mult1, mult2, mult3, mult4, tr1);

    trace.bvalue = 1.0;
    stack.reverse();

    g.set_values(Ab, Bb);
  }

  void hprod(const TupleType& x, const TupleType& p, TupleType& h) {
    A2D::A2DMat<A2D::Mat<T, M, N>> A;
    A2D::A2DMat<A2D::Mat<T, N, K>> B;
    A2D::A2DMat<A2D::Mat<T, M, K>> C;
    A2D::A2DMat<A2D::Mat<T, K, K>> D;
    A2D::A2DMat<A2D::Mat<T, K, M>> E;
    A2D::A2DMat<A2D::Mat<T, K, K>> F;

    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());

    // C = A * B
    auto mult1 = A2D::MatMatMult(A, B, C);
    // D = C^{T} * C
    auto mult2 = A2D::MatMatMult<TRANSPOSE, NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    auto mult3 = A2D::MatMatMult<TRANSPOSE, TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    auto mult4 = A2D::MatMatMult<NORMAL, TRANSPOSE>(E, E, F);

    A2D::A2DScalar<T> trace;
    auto tr1 = A2D::MatTrace(D, trace);

    auto stack = A2D::MakeStack(mult1, mult2, mult3, mult4, tr1);

    trace.bvalue = 1.0;
    stack.reverse();
    stack.hforward();
    stack.hreverse();

    h.set_values(A.hvalue(), B.hvalue());
  }
};

int main() {
  TestMatDet<std::complex<double>, 3> test_det;
  A2D::Test::RunADTest(test_det);

  TestMatInv<std::complex<double>, 3> test_inv;
  A2D::Test::RunADTest(test_inv);

  TestMatMult<std::complex<double>, 3, 3, 3> test1;
  A2D::Test::RunADTest(test1);

  TestMatMult<std::complex<double>, 5, 4, 1> test2;
  A2D::Test::RunADTest(test2);

  return (0);
}