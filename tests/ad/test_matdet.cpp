#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T, int N>
class MatDetTest : public Test::A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatDet<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T det;
    Mat<T, N, N> A, Ainv;
    x.get_values(A);
    MatDet(A, det);
    MatInv(A, Ainv);
    return MakeVarTuple<T>(det);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> det;
    ADObj<Mat<T, N, N>> A, Ainv;

    x.get_values(A.value());
    auto stack = MakeStack(MatDet(A, det), MatInv(A, Ainv));
    seed.get_values(det.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> det;
    A2DObj<Mat<T, N, N>> A, Ainv;
    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatDet(A, det), MatInv(A, Ainv));
    seed.get_values(det.bvalue());
    hval.get_values(det.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N>
class SymMatDetTest : public Test::A2DTest<T, T, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatDet<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T det;
    SymMat<T, N> S, Sinv;
    x.get_values(S);
    MatDet(S, det);
    return MakeVarTuple<T>(det);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> det;
    ADObj<SymMat<T, N>> S, Sinv;

    x.get_values(S.value());
    auto stack = MakeStack(MatDet(S, det));
    seed.get_values(det.bvalue());
    stack.reverse();
    g.set_values(S.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> det;
    A2DObj<SymMat<T, N>> S, Sinv;
    x.get_values(S.value());
    p.get_values(S.pvalue());
    auto stack = MakeStack(MatDet(S, det));
    seed.get_values(det.bvalue());
    hval.get_values(det.hvalue());
    stack.hproduct();
    h.set_values(S.hvalue());
  }
};

inline bool MatDetTestAll(bool component = false, bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  MatDetTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatDetTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  SymMatDetTest<Tc, 2> test3;
  passed = passed && Run(test3, component, write_output);
  SymMatDetTest<Tc, 3> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

TEST(test_matdet, MatDetTestAll) { EXPECT_TRUE(MatDetTestAll()); }