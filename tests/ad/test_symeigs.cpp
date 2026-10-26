#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class SymEigsTest : public Test::A2DTest<T, Vec<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymEigsTest<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    SymMat<T, N> S;
    Vec<T, N> eigs;
    x.get_values(S);
    SymEigs(S, eigs);
    return MakeVarTuple<T>(eigs);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<SymMat<T, N>> S;
    ADObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    auto stack = MakeStack(SymEigs(S, eigs));
    seed.get_values(eigs.bvalue());
    stack.reverse();
    g.set_values(S.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> S;
    A2DObj<Vec<T, N>> eigs;
    x.get_values(S.value());
    p.get_values(S.pvalue());
    auto stack = MakeStack(SymEigs(S, eigs));
    seed.get_values(eigs.bvalue());
    hval.get_values(eigs.hvalue());
    stack.hproduct();
    h.set_values(S.hvalue());
  }
};

inline bool SymEigsTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 2> test1;
    passed = passed && Run(test1, component, write_output);
  }

  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 3> test1;
    passed = passed && Run(test1, component, write_output);
  }

  for (int i = 0; i < 10; i++) {
    SymEigsTest<Tc, 10> test1;
    test1.set_tolerances(1e-8, 1e-12);
    passed = passed && Run(test1, component, write_output);
  }

  return passed;
}

TEST(test_symeigs, SymEigsTestAll) { EXPECT_TRUE(SymEigsTestAll()); }