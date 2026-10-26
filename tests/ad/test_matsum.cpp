#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N, int M>
class MatSumTest
    : public Test::A2DTest<T, Mat<T, N, M>, Mat<T, N, M>, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>, Mat<T, N, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatSum<" << N << "," << M << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, M> A, B, C;
    x.get_values(A, B);
    MatSum(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, M>> A, B, C;
    x.get_values(A.value(), B.value());
    auto stack = MakeStack(MatSum(A, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, M>> A, B, C;
    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto stack = MakeStack(MatSum(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

template <typename T, int N, int M>
class MatSumScaleTest
    : public Test::A2DTest<T, Mat<T, N, M>, T, Mat<T, N, M>, T, Mat<T, N, M>> {
 public:
  using Input = VarTuple<T, T, Mat<T, N, M>, T, Mat<T, N, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatSum<" << N << "," << M << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T alpha, beta;
    Mat<T, N, M> A, B, C;
    x.get_values(alpha, A, beta, B);
    MatSum(alpha, A, beta, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> alpha, beta;
    ADObj<Mat<T, N, M>> A, B, C;
    x.get_values(alpha.value(), A.value(), beta.value(), B.value());
    auto stack = MakeStack(MatSum(alpha, A, beta, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(alpha.bvalue(), A.bvalue(), beta.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> alpha, beta;
    A2DObj<Mat<T, N, M>> A, B, C;
    x.get_values(alpha.value(), A.value(), beta.value(), B.value());
    p.get_values(alpha.pvalue(), A.pvalue(), beta.pvalue(), B.pvalue());
    auto stack = MakeStack(MatSum(alpha, A, beta, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(alpha.hvalue(), A.hvalue(), beta.hvalue(), B.hvalue());
  }
};

inline bool MatSumTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  MatSumTest<Tc, 3, 4> test1;
  passed = passed && Run(test1, component, write_output);
  MatSumTest<Tc, 5, 3> test2;
  passed = passed && Run(test2, component, write_output);

  MatSumScaleTest<Tc, 3, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

TEST(test_matsum, MatSumTestAll) { EXPECT_TRUE(MatSumTestAll()); }