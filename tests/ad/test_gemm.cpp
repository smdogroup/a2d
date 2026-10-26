#include <gtest/gtest.h>

#include "a2dcore.h"

using namespace A2D;

template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
          int Q>
class MatMatMultTest
    : public Test::A2DTest<T, Mat<T, P, Q>, Mat<T, N, M>, Mat<T, K, L>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>, Mat<T, K, L>>;
  using Output = VarTuple<T, Mat<T, P, Q>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatMatMult<";
    if (opA == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    if (opB == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    s << N << "," << M << "," << K << "," << L << "," << P << "," << Q << ">";

    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, M> A;
    Mat<T, K, L> B;
    Mat<T, P, Q> C;

    x.get_values(A, B);
    MatMatMult<opA, opB>(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, M>> A;
    ADObj<Mat<T, K, L>> B;
    ADObj<Mat<T, P, Q>> C;

    x.get_values(A.value(), B.value());
    auto stack = MakeStack(MatMatMult<opA, opB>(A, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, M>> A;
    A2DObj<Mat<T, K, L>> B;
    A2DObj<Mat<T, P, Q>> C;

    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto stack = MakeStack(MatMatMult<opA, opB>(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

template <typename T, int N, int M, int K>
bool MatMatMultTestHelper(bool component = false, bool write_output = true) {
  const MatOp NORMAL = MatOp::NORMAL;
  const MatOp TRANSPOSE = MatOp::TRANSPOSE;
  using Tc = A2D_complex_t<T>;

  bool passed = true;
  MatMatMultTest<NORMAL, NORMAL, Tc, N, M, M, K, N, K> test1;
  passed = passed && Run(test1, component, write_output);
  MatMatMultTest<NORMAL, TRANSPOSE, Tc, N, M, K, M, N, K> test2;
  passed = passed && Run(test2, component, write_output);
  MatMatMultTest<TRANSPOSE, NORMAL, Tc, N, M, N, K, M, K> test3;
  passed = passed && Run(test3, component, write_output);
  MatMatMultTest<TRANSPOSE, TRANSPOSE, Tc, N, M, K, N, M, K> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

inline bool MatMatMultTestAll(bool component = false,
                              bool write_output = true) {
  bool passed = true;
  passed =
      passed && MatMatMultTestHelper<double, 3, 3, 3>(component, write_output);
  passed =
      passed && MatMatMultTestHelper<double, 2, 3, 4>(component, write_output);
  passed =
      passed && MatMatMultTestHelper<double, 5, 4, 2>(component, write_output);

  return passed;
}

TEST(test_gemm, MatMatMultTestAll) { EXPECT_TRUE(MatMatMultTestAll()); }