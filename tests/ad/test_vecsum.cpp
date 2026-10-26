#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class VecSumTest : public Test::A2DTest<T, Vec<T, N>, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "VecSum<" << N << ">";
    return s.str();
  }

  // Evaluate the operation
  Output eval(const Input& x) {
    Vec<T, N> A, B, C;
    x.get_values(A, B);
    VecSum(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Vec<T, N>> A, B, C;
    x.get_values(A.value(), B.value());
    auto stack = MakeStack(VecSum(A, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Vec<T, N>> A, B, C;
    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto stack = MakeStack(VecSum(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

template <typename T, int N>
class VecSumScaleTest
    : public Test::A2DTest<T, Vec<T, N>, T, Vec<T, N>, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, T, Vec<T, N>, T, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "VecSum<" << N << ">";
    return s.str();
  }

  // Evaluate the Vecrix-Vecrix product
  Output eval(const Input& x) {
    T alpha, beta;
    Vec<T, N> A, B, C;
    x.get_values(alpha, A, beta, B);
    VecSum(alpha, A, beta, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> alpha, beta;
    ADObj<Vec<T, N>> A, B, C;
    x.get_values(alpha.value(), A.value(), beta.value(), B.value());
    auto stack = MakeStack(VecSum(alpha, A, beta, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(alpha.bvalue(), A.bvalue(), beta.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> alpha, beta;
    A2DObj<Vec<T, N>> A, B, C;
    x.get_values(alpha.value(), A.value(), beta.value(), B.value());
    p.get_values(alpha.pvalue(), A.pvalue(), beta.pvalue(), B.pvalue());
    auto stack = MakeStack(VecSum(alpha, A, beta, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(alpha.hvalue(), A.hvalue(), beta.hvalue(), B.hvalue());
  }
};

inline bool VecSumTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecSumTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecSumTest<Tc, 5> test2;
  passed = passed && Run(test2, component, write_output);

  VecSumScaleTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

TEST(test_vecsum, VecSumTestAll) { EXPECT_TRUE(VecSumTestAll()); }