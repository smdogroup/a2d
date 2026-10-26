#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class MatInvTest : public Test::A2DTest<T, Mat<T, N, N>, Mat<T, N, N>> {
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
    ADObj<Mat<T, N, N>> A;
    ADObj<Mat<T, N, N>> B;

    x.get_values(A.value());
    auto stack = MakeStack(MatInv(A, B));
    seed.get_values(B.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, N>> A;
    A2DObj<Mat<T, N, N>> B;

    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatInv(A, B));
    seed.get_values(B.bvalue());
    hval.get_values(B.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

inline bool MatInvTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  MatInvTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatInvTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_matinv, MatInvTestAll) { EXPECT_TRUE(MatInvTestAll()); }