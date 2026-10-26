#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class SymIsotropicConstTest
    : public Test::A2DTest<T, SymMat<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  std::string name() {
    std::stringstream s;
    s << "SymIsotropic<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    SymMat<T, N> E, S;
    x.get_values(E);
    SymIsotropic(T(0.314), T(0.731), E, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<SymMat<T, N>> E, S;

    x.get_values(E.value());
    auto stack = MakeStack(SymIsotropic(T(0.314), T(0.731), E, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> E, S;

    x.get_values(E.value());
    p.get_values(E.pvalue());
    auto stack = MakeStack(SymIsotropic(T(0.314), T(0.731), E, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(E.hvalue());
  }
};

template <typename T, int N>
class SymIsotropicTest
    : public Test::A2DTest<T, SymMat<T, N>, T, T, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, T, T, SymMat<T, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  std::string name() {
    std::stringstream s;
    s << "SymIsotropic<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T mu, lambda;
    SymMat<T, N> E, S;
    x.get_values(mu, lambda, E);
    SymIsotropic(mu, lambda, E, S);
    return MakeVarTuple<T>(S);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> mu, lambda;
    ADObj<SymMat<T, N>> E, S;

    x.get_values(mu.value(), lambda.value(), E.value());
    auto stack = MakeStack(SymIsotropic(mu, lambda, E, S));
    seed.get_values(S.bvalue());
    stack.reverse();
    g.set_values(mu.bvalue(), lambda.bvalue(), E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> mu, lambda;
    A2DObj<SymMat<T, N>> E, S;

    x.get_values(mu.value(), lambda.value(), E.value());
    p.get_values(mu.pvalue(), lambda.pvalue(), E.pvalue());
    auto stack = MakeStack(SymIsotropic(mu, lambda, E, S));
    seed.get_values(S.bvalue());
    hval.get_values(S.hvalue());
    stack.hproduct();
    h.set_values(mu.hvalue(), lambda.hvalue(), E.hvalue());
  }
};

inline bool SymIsotropicTestAll(bool component = false,
                                bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  SymIsotropicConstTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymIsotropicConstTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  SymIsotropicTest<Tc, 2> test3;
  passed = passed && Run(test3, component, write_output);
  SymIsotropicTest<Tc, 3> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

TEST(test_isotropic, SymIsotropicTestAll) {
  EXPECT_TRUE(SymIsotropicTestAll());
}