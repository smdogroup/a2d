#include <gtest/gtest.h>

#include "a2dcore.h"
#include "a2dtest.h"

using namespace A2D;

template <typename T>
class QuaternionMatrixTest : public Test::A2DTest<T, Mat<T, 3, 3>, Vec<T, 4>> {
 public:
  using Input = VarTuple<T, Vec<T, 4>>;
  using Output = VarTuple<T, Mat<T, 3, 3>>;

  // Assemble a string to describe the test
  std::string name() { return "QuaternionMatrix"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Vec<T, 4> q;
    Mat<T, 3, 3> C;
    X.get_values(q);
    QuaternionMatrix(q, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Vec<T, 4>> q;
    ADObj<Mat<T, 3, 3>> C;
    X.get_values(q.value());
    auto stack = MakeStack(QuaternionMatrix(q, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(q.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Vec<T, 4>> q;
    A2DObj<Mat<T, 3, 3>> C;
    X.get_values(q.value());
    p.get_values(q.pvalue());
    auto stack = MakeStack(QuaternionMatrix(q, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(q.hvalue());
  }
};

template <typename T>
class QuaternionAngularVelocityTest
    : public Test::A2DTest<T, Vec<T, 3>, Vec<T, 4>, Vec<T, 4>> {
 public:
  using Input = VarTuple<T, Vec<T, 4>, Vec<T, 4>>;
  using Output = VarTuple<T, Vec<T, 3>>;

  // Assemble a string to describe the test
  std::string name() { return "QuaternionMatrix"; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& X) {
    Vec<T, 4> q, qdot;
    Vec<T, 3> omega;
    X.get_values(q, qdot);
    QuaternionAngularVelocity(q, qdot, omega);
    return MakeVarTuple<T>(omega);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& X, Input& g) {
    ADObj<Vec<T, 4>> q, qdot;
    ADObj<Vec<T, 3>> omega;
    X.get_values(q.value(), qdot.value());
    auto stack = MakeStack(QuaternionAngularVelocity(q, qdot, omega));
    seed.get_values(omega.bvalue());
    stack.reverse();
    g.set_values(q.bvalue(), qdot.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& X,
             const Input& p, Input& h) {
    A2DObj<Vec<T, 4>> q, qdot;
    A2DObj<Vec<T, 3>> omega;
    X.get_values(q.value(), qdot.value());
    p.get_values(q.pvalue(), qdot.pvalue());
    auto stack = MakeStack(QuaternionAngularVelocity(q, qdot, omega));
    seed.get_values(omega.bvalue());
    hval.get_values(omega.hvalue());
    stack.hproduct();
    h.set_values(q.hvalue(), qdot.hvalue());
  }
};

template <typename T>
bool TestQuaternions(QuaternionAngularVelocityTest<T> test,
                     bool write_output = false, const double dh = 1e-30,
                     const double rtol = 1e-14) {
  using Tc = A2D_complex_t<double>;

  Vec<Tc, 4> q, qdot;
  Mat<Tc, 3, 3> C, Cdot, A;
  Vec<Tc, 3> omega;
  Mat<Tc, 3, 3> OmegaX, negOmegaX;

  Tc qnrm = 0.0;
  for (int i = 0; i < 4; i++) {
    q(i) = -1.0 + 2.0 * (static_cast<double>(rand()) / RAND_MAX);
    qnrm += q(i) * q(i);
  }
  qnrm = sqrt(qnrm);

  for (int i = 0; i < 4; i++) {
    q(i) = q(i) / qnrm;
  }

  // Compute a qdot such that q^{T} * qdot = 0
  for (int i = 1; i < 4; i++) {
    qdot(i) = -1.0 + 2.0 * (static_cast<double>(rand()) / RAND_MAX);
    qdot(0) -= q(i) * qdot(i) / q(0);
  }

  for (int i = 0; i < 4; i++) {
    q(i) = Tc(RealPart(q(i)), dh * RealPart(qdot(i)));
  }

  // Compute the quaternion matrix and the angular rate
  QuaternionMatrix(q, C);
  QuaternionAngularVelocity(q, qdot, omega);

  // Set Cdot
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Cdot(i, j) = ImagPart(C(i, j)) / dh;
    }
  }

  // Set the output matrix
  OmegaX(0, 1) = -omega(2);
  OmegaX(0, 2) = omega(1);

  OmegaX(1, 0) = omega(2);
  OmegaX(1, 2) = -omega(0);

  OmegaX(2, 0) = -omega(1);
  OmegaX(2, 1) = omega(0);

  double err = 0.0;

  // Make sure C^{T} * C = I
  MatMatMult<MatOp::TRANSPOSE, MatOp::NORMAL>(C, C, A);

  // Compute the error for C^{T} * C - I
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        err += std::fabs(RealPart(A(i, j)) - 1.0);
      } else {
        err += std::fabs(RealPart(A(i, j)));
      }

      if (write_output) {
        if (i == j) {
          test.write_result("C^{T} * C", std::cout, A(i, j), 1.0);
        } else {
          test.write_result("C^{T} * C", std::cout, A(i, j), 0.0);
        }
      }
    }
  }

  // Omega^{x} = - dot{C} * C^{T}
  MatMatMult<MatOp::NORMAL, MatOp::TRANSPOSE>(Cdot, C, negOmegaX);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Tc ans = OmegaX(i, j) + negOmegaX(i, j);
      err += std::fabs(RealPart(ans));

      if (write_output) {
        test.write_result("omega^{x} + dot{C} * C^{T}", std::cout, ans, 0.0);
      }
    }
  }

  bool passed = (std::fabs(err) < rtol);

  return passed;
}

inline bool QuaternionTestAll(bool component = false,
                              bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  QuaternionMatrixTest<Tc> test1;
  bool passed = Run(test1, component, write_output);

  QuaternionAngularVelocityTest<Tc> test2;
  test1.set_tolerances(1e-10, 1e-14);
  passed = passed && TestQuaternions(test2, write_output, 1e-30, 1e-12);
  passed = passed && Run(test2, component, write_output);

  return passed;
}

TEST(test_quaternion, QuaternionTestAll) { EXPECT_TRUE(QuaternionTestAll()); }