#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class SymMatSumTest : public Test::A2DTest<T, SymMat<T, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatSum<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> A;
    SymMat<T, N> S;
    x.get_values(A);
    SymMatSum(A, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, N>> A;
    ADObj<SymMat<T, N>> S;
    x.get_values(A.value());
    auto stack = MakeStack(SymMatSum(A, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> S;
    A2DObj<Mat<T, N, N>> A;
    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(SymMatSum(A, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N>
class SymMatSumScaleTest
    : public Test::A2DTest<T, SymMat<T, N>, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, T, Mat<T, N, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatSum<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T alpha;
    Mat<T, N, N> A;
    SymMat<T, N> S;
    x.get_values(alpha, A);
    SymMatSum(alpha, A, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> alpha;
    ADObj<Mat<T, N, N>> A;
    ADObj<SymMat<T, N>> S;
    x.get_values(alpha.value(), A.value());
    auto stack = MakeStack(SymMatSum(alpha, A, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(alpha.bvalue(), A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> alpha;
    A2DObj<SymMat<T, N>> S;
    A2DObj<Mat<T, N, N>> A;
    x.get_values(alpha.value(), A.value());
    p.get_values(alpha.pvalue(), A.pvalue());
    auto stack = MakeStack(SymMatSum(alpha, A, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(alpha.hvalue(), A.hvalue());
  }
};

inline bool SymMatSumTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  SymMatSumTest<Tc, 4> test1;
  passed = passed && Run(test1, component, write_output);

  SymMatSumScaleTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_symmatsum, SymMatSumTestAll) { EXPECT_TRUE(SymMatSumTestAll()); }