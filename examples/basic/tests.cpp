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
    // T scalar;
    A2D::Mat<T, M, N> A;
    A2D::Mat<T, N, K> B;
    A2D::Mat<T, M, K> C;
    A2D::Mat<T, K, K> D;

    x.get_values(A, B);

    // Compute C = A * B
    A2D::MatMatMult(A, B, C);
    // Compute D = C^{T} * C
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(C, C, D);

    T trace;
    A2D::MatTrace(D, trace);

    return trace;
  }

  void deriv(const TupleType& x, TupleType& g) {
    // T scalar;
    A2D::Mat<T, M, N> A0, Ab;
    A2D::Mat<T, N, K> B0, Bb;
    A2D::Mat<T, M, K> C0, Cb;
    A2D::Mat<T, K, K> D0, Db;

    A2D::ADMat<A2D::Mat<T, M, N>> A(A0, Ab);
    A2D::ADMat<A2D::Mat<T, N, K>> B(B0, Bb);
    A2D::ADMat<A2D::Mat<T, M, K>> C(C0, Cb);
    A2D::ADMat<A2D::Mat<T, K, K>> D(D0, Db);

    x.get_values(A0, B0);

    // Compute C = A * B
    auto mult1 = A2D::MatMatMult(A, B, C);
    // Compute D = C^{T} * C
    auto mult2 =
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(C, C, D);

    A2D::ADScalar<T> trace;
    auto tr1 = A2D::MatTrace(D, trace);

    auto stack = A2D::MakeStack(mult1, mult2, tr1);

    trace.bvalue = 1.0;
    stack.reverse();

    g.set_values(Ab, Bb);
  }
};

int main() {
  TestMatMult<std::complex<double>, 3, 3, 3> test1;

  A2D::Test::RunADTest(test1);

  return (0);
}