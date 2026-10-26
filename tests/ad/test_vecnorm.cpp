#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class VecNormTest : public Test::A2DTest<T, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "VecNorm"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    T alpha;
    Vec<T, N> x;
    X.get_values(x);
    VecNorm(x, alpha);
    return MakeVarTuple<T>(alpha);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<T> alpha;
    ADObj<Vec<T, N>> x;
    X.get_values(x.value());
    auto stack = MakeStack(VecNorm(x, alpha));
    seed.get_values(alpha.bvalue());
    stack.reverse();
    g.set_values(x.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<T> alpha;
    A2DObj<Vec<T, N>> x;
    X.get_values(x.value());
    p.get_values(x.pvalue());
    auto stack = MakeStack(VecNorm(x, alpha));
    seed.get_values(alpha.bvalue());
    hval.get_values(alpha.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue());
  }
};

inline bool VecNormTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecNormTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecNormTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_vecnorm, VecNormTestAll) { EXPECT_TRUE(VecNormTestAll()); }

template <typename T, int N>
class VecScaleTest : public Test::A2DTest<T, Vec<T, N>, T, Vec<T, N>> {
 public:
  using Input = VarTuple<T, T, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() { return "VecScale"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    T alpha;
    Vec<T, N> x, y;
    X.get_values(alpha, x);
    VecScale(alpha, x, y);
    return MakeVarTuple<T>(y);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<T> alpha;
    ADObj<Vec<T, N>> x, y;
    X.get_values(alpha.value(), x.value());
    auto stack = MakeStack(VecScale(alpha, x, y));
    seed.get_values(y.bvalue());
    stack.reverse();
    g.set_values(alpha.bvalue(), x.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<T> alpha;
    A2DObj<Vec<T, N>> x, y;
    X.get_values(alpha.value(), x.value());
    p.get_values(alpha.pvalue(), x.pvalue());
    auto stack = MakeStack(VecScale(alpha, x, y));
    seed.get_values(y.bvalue());
    hval.get_values(y.hvalue());
    stack.hproduct();
    h.set_values(alpha.hvalue(), x.hvalue());
  }
};

inline bool VecScaleTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecScaleTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecScaleTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_vecscale, VecScaleTestAll) { EXPECT_TRUE(VecScaleTestAll()); }

template <typename T, int N>
class VecNormalizeTest : public Test::A2DTest<T, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() { return "VecNormalize"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Vec<T, N> x, y;
    X.get_values(x);
    VecNormalize(x, y);
    return MakeVarTuple<T>(y);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Vec<T, N>> x, y;
    X.get_values(x.value());
    auto stack = MakeStack(VecNormalize(x, y));
    seed.get_values(y.bvalue());
    stack.reverse();
    g.set_values(x.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Vec<T, N>> x, y;
    X.get_values(x.value());
    p.get_values(x.pvalue());
    auto stack = MakeStack(VecNormalize(x, y));
    seed.get_values(y.bvalue());
    hval.get_values(y.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue());
  }
};

inline bool VecNormalizeTestAll(bool component = false,
                                bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecNormalizeTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecNormalizeTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_vecnormalize, VecNormalizeTestAll) {
  EXPECT_TRUE(VecNormalizeTestAll());
}

template <typename T, int N>
class VecDotTest : public Test::A2DTest<T, T, Vec<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, Vec<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "VecDot"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    T alpha;
    Vec<T, N> x, y;
    X.get_values(x, y);
    VecDot(x, y, alpha);
    return MakeVarTuple<T>(alpha);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<T> alpha;
    ADObj<Vec<T, N>> x, y;
    X.get_values(x.value(), y.value());
    auto stack = MakeStack(VecDot(x, y, alpha));
    seed.get_values(alpha.bvalue());
    stack.reverse();
    g.set_values(x.bvalue(), y.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<T> alpha;
    A2DObj<Vec<T, N>> x, y;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto stack = MakeStack(VecDot(x, y, alpha));
    seed.get_values(alpha.bvalue());
    hval.get_values(alpha.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

inline bool VecDotTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecDotTest<Tc, 3> test1;
  passed = passed && Run(test1, component, write_output);
  VecDotTest<Tc, 6> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_vecdot, VecDotTestAll) { EXPECT_TRUE(VecDotTestAll()); }
