#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <MatOp op, typename T, int N, int M, int P>
class SymMatRKTest : public Test::A2DTest<T, SymMat<T, P>, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>>;
  using Output = VarTuple<T, SymMat<T, P>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatRKTest<";
    if (op == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    s << N << "," << M << "," << P << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, M> A;
    SymMat<T, P> S;
    x.get_values(A);
    SymMatRK<op>(A, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, M>> A;
    ADObj<SymMat<T, P>> S;

    x.get_values(A.value());
    auto stack = MakeStack(SymMatRK<op>(A, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, M>> A;
    A2DObj<SymMat<T, P>> S;

    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(SymMatRK<op>(A, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

template <MatOp op, typename T, int N, int M, int P>
class SymMatRKScaleTest : public Test::A2DTest<T, SymMat<T, P>, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>>;
  using Output = VarTuple<T, SymMat<T, P>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatRKScaleTest<";
    if (op == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    s << N << "," << M << "," << P << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    const T alpha = 0.25;
    Mat<T, N, M> A;
    SymMat<T, P> S;
    x.get_values(A);
    SymMatRK<op>(alpha, A, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    const T alpha = 0.25;
    ADObj<Mat<T, N, M>> A;
    ADObj<SymMat<T, P>> S;

    x.get_values(A.value());
    auto stack = MakeStack(SymMatRK<op>(alpha, A, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    const T alpha = 0.25;
    A2DObj<Mat<T, N, M>> A;
    A2DObj<SymMat<T, P>> S;

    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(SymMatRK<op>(alpha, A, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N, int K>
bool SymMatRKTestHelper(bool component = false, bool write_output = true) {
  const MatOp NORMAL = MatOp::NORMAL;
  const MatOp TRANSPOSE = MatOp::TRANSPOSE;
  using Tc = A2D_complex_t<T>;

  bool passed = true;
  SymMatRKTest<NORMAL, Tc, N, K, N> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatRKTest<TRANSPOSE, Tc, N, K, K> test2;
  passed = passed && Run(test2, component, write_output);

  SymMatRKScaleTest<NORMAL, Tc, N, K, N> test3;
  passed = passed && Run(test3, component, write_output);
  SymMatRKScaleTest<TRANSPOSE, Tc, N, K, K> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

inline bool SymMatRKTestAll(bool component = false, bool write_output = true) {
  bool passed = true;
  passed = passed && SymMatRKTestHelper<double, 3, 3>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 2, 4>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 5, 2>(component, write_output);
  passed = passed && SymMatRKTestHelper<double, 7, 5>(component, write_output);

  return passed;
}

TEST(test_symrk, SymMatRKTestAll) { EXPECT_TRUE(SymMatRKTestAll()); }