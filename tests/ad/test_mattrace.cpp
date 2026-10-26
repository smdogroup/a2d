#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class MatTraceTest : public Test::A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatTrace<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T trace;
    Mat<T, N, N> A;
    x.get_values(A);
    MatTrace(A, trace);
    return MakeVarTuple<T>(trace);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> trace;
    ADObj<Mat<T, N, N>> A;

    x.get_values(A.value());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> trace;
    A2DObj<Mat<T, N, N>> A;
    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    hval.get_values(trace.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N>
class SymTraceTest : public Test::A2DTest<T, T, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatTrace<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T trace;
    SymMat<T, N> A;
    x.get_values(A);
    MatTrace(A, trace);
    return MakeVarTuple<T>(trace);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> trace;
    ADObj<SymMat<T, N>> A;

    x.get_values(A.value());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> trace;
    A2DObj<SymMat<T, N>> A;

    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    hval.get_values(trace.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

inline bool MatTraceTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  MatTraceTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatTraceTest<Tc, 4> test2;
  passed = passed && Run(test2, component, write_output);

  SymTraceTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

TEST(test_mattrace, MatTraceTestAll) { EXPECT_TRUE(MatTraceTestAll()); }