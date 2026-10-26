#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class VecHadamardTest
    : public Test::A2DTest<T, Vec<T, N>, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "VecHadamard<" << N << ">";
    return s.str();
  }

  // Evaluate the operation
  Output eval(const Input& x) {
    Vec<T, N> A, B, C;
    x.get_values(A, B);
    VecHadamard(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Vec<T, N>> A, B, C;
    x.get_values(A.value(), B.value());
    auto stack = MakeStack(VecHadamard(A, B, C));
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
    auto stack = MakeStack(VecHadamard(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

inline bool VecHadamardTestAll(bool component = false,
                               bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecHadamardTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);

  return passed;
}

TEST(test_hadamard, VecHadamardTestAll) { EXPECT_TRUE(VecHadamardTestAll()); }