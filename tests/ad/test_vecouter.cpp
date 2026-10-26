#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N, int M>
class VecOuterTest
    : public Test::A2DTest<T, Mat<T, N, M>, Vec<T, N>, Vec<T, M>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, M>>;
  using Output = VarTuple<T, Mat<T, N, M>>;

  // Assemble a string to describe the test
  std::string name() { return "VecOuter"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    T alpha(0.3157);
    Vec<T, N> x;
    Vec<T, M> y;
    Mat<T, N, M> A;
    X.get_values(x, y);
    VecOuter(alpha, x, y, A);
    return MakeVarTuple<T>(A);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    T alpha(0.3157);
    ADObj<Vec<T, N>> x;
    ADObj<Vec<T, M>> y;
    ADObj<Mat<T, N, M>> A;
    X.get_values(x.value(), y.value());
    auto stack = MakeStack(VecOuter(alpha, x, y, A));
    seed.get_values(A.bvalue());
    stack.reverse();
    g.set_values(x.bvalue(), y.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    T alpha(0.3157);
    A2DObj<Vec<T, N>> x;
    A2DObj<Vec<T, M>> y;
    A2DObj<Mat<T, N, M>> A;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto stack = MakeStack(VecOuter(alpha, x, y, A));
    seed.get_values(A.bvalue());
    hval.get_values(A.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

inline bool VecOuterTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecOuterTest<Tc, 3, 5> test1;
  passed = passed && Run(test1, component, write_output);
  VecOuterTest<Tc, 6, 4> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_vecouter, VecOuterTestAll) { EXPECT_TRUE(VecOuterTestAll()); }