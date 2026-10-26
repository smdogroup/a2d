#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class SymMatMultTraceTest
    : public Test::A2DTest<T, T, SymMat<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatMultTrace<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T output;
    SymMat<T, N> S, E;
    x.get_values(S, E);
    SymMatMultTrace(S, E, output);
    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> output;
    ADObj<SymMat<T, N>> S, E;

    x.get_values(S.value(), E.value());
    auto stack = MakeStack(SymMatMultTrace(S, E, output));
    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(S.bvalue(), E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> output;
    A2DObj<SymMat<T, N>> S, E;

    x.get_values(S.value(), E.value());
    p.get_values(S.pvalue(), E.pvalue());
    auto stack = MakeStack(SymMatMultTrace(S, E, output));
    seed.get_values(output.bvalue());
    hval.get_values(output.hvalue());
    stack.hproduct();
    h.set_values(S.hvalue(), E.hvalue());
  }
};

inline bool SymMatMultTraceTestAll(bool component = false,
                                   bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  SymMatMultTraceTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatMultTraceTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  SymMatMultTraceTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

TEST(test_symmatmulttrace, SymMatMultTraceTestAll) {
  EXPECT_TRUE(SymMatMultTraceTestAll());
}