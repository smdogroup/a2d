#include "a2denum.h"
#include "ad/a2dgemm.h"
#include "ad/a2dmattrace.h"
#include "ad/a2dstack.h"
#include "ad/a2dtest.h"

template <typename T, int M, int N, int K>
class TestMatMult
    : public A2D::Test::A2DTest<T, A2D::Mat<T, M, N>, A2D::Mat<T, N, K>> {
 public:
  using TupleType = A2D::VarTuple<T, A2D::Mat<T, M, N>, A2D::Mat<T, N, K>>;

  TestMatMult() {}

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
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    A2D::MatMatMult<A2D::MatOp::NORMAL, A2D::MatOp::TRANSPOSE>(E, E, F);

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
    auto mult2 =
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    auto mult3 =
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    auto mult4 =
        A2D::MatMatMult<A2D::MatOp::NORMAL, A2D::MatOp::TRANSPOSE>(E, E, F);

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
    auto mult2 =
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(C, C, D);
    // E = D^{T} * C^{T}
    auto mult3 =
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::TRANSPOSE>(D, C, E);
    // F = E * E^{T}
    auto mult4 =
        A2D::MatMatMult<A2D::MatOp::NORMAL, A2D::MatOp::TRANSPOSE>(E, E, F);

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
  TestMatMult<std::complex<double>, 3, 3, 3> test1;
  A2D::Test::RunADTest(test1);

  TestMatMult<std::complex<double>, 3, 4, 5> test2;
  A2D::Test::RunADTest(test2);

  return (0);
}