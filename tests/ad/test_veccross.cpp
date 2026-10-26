#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T>
class VecCross3DTest
    : public Test::A2DTest<T, Vec<T, 3>, Vec<T, 3>, Vec<T, 3>> {
 public:
  using Input = VarTuple<T, Vec<T, 3>, Vec<T, 3>>;
  using Output = VarTuple<T, Vec<T, 3>>;

  // Assemble a string to describe the test
  std::string name() { return "VecCross"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Vec<T, 3> x, y, z;
    X.get_values(x, y);
    VecCross(x, y, z);
    return MakeVarTuple<T>(z);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Vec<T, 3>> x, y, z;
    X.get_values(x.value(), y.value());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    stack.reverse();
    g.set_values(x.bvalue(), y.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Vec<T, 3>> x, y, z;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    hval.get_values(z.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

template <typename T>
class VecCross2DTest : public Test::A2DTest<T, Vec<T, 2>, Vec<T, 2>, T> {
 public:
  using Input = VarTuple<T, Vec<T, 2>, T>;
  using Output = VarTuple<T, Vec<T, 2>>;

  // Assemble a string to describe the test
  std::string name() { return "VecCross"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Vec<T, 2> x, z;
    T y;
    X.get_values(x, y);
    VecCross(x, y, z);
    return MakeVarTuple<T>(z);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Vec<T, 2>> x, z;
    ADObj<T> y;
    X.get_values(x.value(), y.value());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    stack.reverse();
    g.set_values(x.bvalue(), y.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Vec<T, 2>> x, z;
    A2DObj<T> y;
    X.get_values(x.value(), y.value());
    p.get_values(x.pvalue(), y.pvalue());
    auto stack = MakeStack(VecCross(x, y, z));
    seed.get_values(z.bvalue());
    hval.get_values(z.hvalue());
    stack.hproduct();
    h.set_values(x.hvalue(), y.hvalue());
  }
};

inline bool VecCrossTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  VecCross3DTest<Tc> test1;
  passed = passed && Run(test1, component, write_output);

  VecCross2DTest<Tc> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}
