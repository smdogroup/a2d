#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <MatOp op, typename T, int N, int M, int K, int P>
class MatVecMultTest
    : public Test::A2DTest<T, Vec<T, P>, Mat<T, N, M>, Vec<T, K>> {
 public:
  using Input = VarTuple<T, Mat<T, N, M>, Vec<T, K>>;
  using Output = VarTuple<T, Vec<T, P>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatVecMult<";
    if (op == MatOp::NORMAL) {
      s << "N,";
    } else {
      s << "T,";
    }
    s << N << "," << M << "," << K << "," << P << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Mat<T, N, M> A;
    Vec<T, K> x;
    Vec<T, P> y;

    X.get_values(A, x);
    MatVecMult<op>(A, x, y);
    return MakeVarTuple<T>(y);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Mat<T, N, M>> A;
    ADObj<Vec<T, K>> x;
    ADObj<Vec<T, P>> y;

    X.get_values(A.value(), x.value());
    auto stack = MakeStack(MatVecMult<op>(A, x, y));
    seed.get_values(y.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), x.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, M>> A;
    A2DObj<Vec<T, K>> x;
    A2DObj<Vec<T, P>> y;

    X.get_values(A.value(), x.value());
    p.get_values(A.pvalue(), x.pvalue());
    auto stack = MakeStack(MatVecMult<op>(A, x, y));
    seed.get_values(y.bvalue());
    hval.get_values(y.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), x.hvalue());
  }
};

template <typename T, int N>
class SymMatVecMultTest
    : public Test::A2DTest<T, Vec<T, N>, SymMat<T, N>, Vec<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>, Vec<T, N>>;
  using Output = VarTuple<T, Vec<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatVecMult<";
    s << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    SymMat<T, N> A;
    Vec<T, N> x;
    Vec<T, N> y;

    X.get_values(A, x);
    MatVecMult(A, x, y);
    return MakeVarTuple<T>(y);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<SymMat<T, N>> A;
    ADObj<Vec<T, N>> x;
    ADObj<Vec<T, N>> y;

    X.get_values(A.value(), x.value());
    auto stack = MakeStack(MatVecMult(A, x, y));
    seed.get_values(y.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), x.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<SymMat<T, N>> A;
    A2DObj<Vec<T, N>> x;
    A2DObj<Vec<T, N>> y;

    X.get_values(A.value(), x.value());
    p.get_values(A.pvalue(), x.pvalue());
    auto stack = MakeStack(MatVecMult(A, x, y));
    seed.get_values(y.bvalue());
    hval.get_values(y.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), x.hvalue());
  }
};

template <typename T, int N>
bool SymMatVecMultTestHelper(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<T>;

  bool passed = true;
  SymMatVecMultTest<Tc, N> test1;
  passed = passed && Run(test1, component, write_output);

  return passed;
}

template <typename T, int N, int M>
bool MatVecMultTestHelper(bool component = false, bool write_output = true) {
  const MatOp NORMAL = MatOp::NORMAL;
  const MatOp TRANSPOSE = MatOp::TRANSPOSE;
  using Tc = A2D_complex_t<T>;

  bool passed = true;
  MatVecMultTest<NORMAL, Tc, N, M, M, N> test1;
  passed = passed && Run(test1, component, write_output);

  MatVecMultTest<TRANSPOSE, Tc, M, N, M, N> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

inline bool MatVecMultTestAll(bool component = false,
                              bool write_output = true) {
  bool passed = true;
  passed =
      passed && MatVecMultTestHelper<double, 3, 3>(component, write_output);
  passed =
      passed && MatVecMultTestHelper<double, 2, 4>(component, write_output);
  passed =
      passed && MatVecMultTestHelper<double, 5, 3>(component, write_output);
  passed =
      passed && SymMatVecMultTestHelper<double, 1>(component, write_output);
  passed =
      passed && SymMatVecMultTestHelper<double, 2>(component, write_output);
  passed =
      passed && SymMatVecMultTestHelper<double, 3>(component, write_output);
  passed =
      passed && SymMatVecMultTestHelper<double, 4>(component, write_output);

  return passed;
}

TEST(test_matvecmult, MatVecMultTestAll) { EXPECT_TRUE(MatVecMultTestAll()); }