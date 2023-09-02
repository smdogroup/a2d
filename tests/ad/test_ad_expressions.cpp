#include <vector>

#include "ad/a2dgemm.h"
#include "ad/a2dgreenstrain.h"
#include "ad/a2disotropic.h"
#include "ad/a2dmatdet.h"
#include "ad/a2dmatinv.h"
#include "ad/a2dmatsum.h"
#include "ad/a2dmattrace.h"
#include "ad/a2dscalarops.h"
#include "ad/a2dsymrk.h"
#include "ad/a2dsymsum.h"
#include "ad/a2dsymtrace.h"
#include "ad/a2dveccross.h"
#include "ad/a2dvecnorm.h"
#include "ad/a2dvecouter.h"
#include "ad/a2dvecsum.h"

using namespace A2D;

template <typename T, int N>
class StrainTest : public A2D::Test::A2DTest<T, T, Mat<T, N, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "StrainTest<" << N << ">";
    return s.str();
  }

  // Evaluate the function
  Output eval(const Input& x) {
    Mat<T, N, N> Uxi, J;
    x.get_values(Uxi, J);
    T output;

    // Set the intermediary values
    Mat<T, N, N> Jinv, Ux;
    SymMat<T, N> E1, E2, E, S;

    MatInv(J, Jinv);                       // Jinv = J^{-1}
    MatMatMult(Uxi, Jinv, Ux);             // Ux = Uxi * Jinv
    SymMatSum(T(0.5), Ux, E1);             // E1 = 0.5 * (Ux + Ux^{T})
    SymMatRK<MatOp::TRANSPOSE>(Ux, E2);    // E2 = Ux^{T} * Ux
    MatSum(T(1.0), E1, T(0.5), E2, E);     // E = E1 + 0.5 * E2
    SymIsotropic(T(0.35), T(0.51), E, S);  // S = 2 * mu * E + lam * tr(E) * I
    SymMatTrace(E, S, output);             // output = tr(E * S)

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    // Input
    Mat<T, N, N> Uxi0, Uxib, J0, Jb;
    x.get_values(Uxi0, J0);
    ADScalar<T> output;

    // Set the intermediary values
    Mat<T, N, N> Jinv0, Jinvb, Ux0, Uxb;
    SymMat<T, N> E10, E1b, E20, E2b, E0, Eb, S0, Sb;

    // The AD objects
    ADMat<Mat<T, N, N>> Uxi(Uxi0, Uxib), J(J0, Jb);
    ADMat<Mat<T, N, N>> Jinv(Jinv0, Jinvb), Ux(Ux0, Uxb);
    ADMat<SymMat<T, N>> E1(E10, E1b), E2(E20, E2b), E(E0, Eb), S(S0, Sb);

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  SymMatSum(T(0.5), Ux, E1),  // E1 = 0.5 * (Ux + Ux^{T})
                  SymMatRK<MatOp::TRANSPOSE>(Ux, E2),    // E2 = Ux^{T} * Ux
                  MatSum(T(1.0), E1, T(0.5), E2, E),     // E = E1 + 0.5 * E2
                  SymIsotropic(T(0.35), T(0.51), E, S),  // S = S(E)
                  SymMatTrace(E, S, output));

    seed.get_values(output.bvalue);
    stack.reverse();
    g.set_values(Uxib, Jb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    // The AD objects
    A2DMat<Mat<T, N, N>> Uxi, J;
    A2DScalar<T> output;
    A2DMat<Mat<T, N, N>> Jinv, Ux;
    A2DMat<SymMat<T, N>> E1, E2, E, S;

    x.get_values(Uxi.value(), J.value());
    p.get_values(Uxi.pvalue(), J.pvalue());

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  SymMatSum(T(0.5), Ux, E1),  // E1 = 0.5 * (Ux + Ux^{T})
                  SymMatRK<MatOp::TRANSPOSE>(Ux, E2),    // E2 = Ux^{T} * Ux
                  MatSum(T(1.0), E1, T(0.5), E2, E),     // E = E1 + 0.5 * E2
                  SymIsotropic(T(0.35), T(0.51), E, S),  // S = S(E)
                  SymMatTrace(E, S, output));

    seed.get_values(output.bvalue);
    hval.get_values(output.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(Uxi.hvalue(), J.hvalue());
  }
};

template <typename T>
class MooneyRivlin
    : public A2D::Test::A2DTest<T, T, Mat<T, 3, 3>, Mat<T, 3, 3>> {
 public:
  using Input = VarTuple<T, Mat<T, 3, 3>, Mat<T, 3, 3>>;
  using Output = VarTuple<T, T>;
  static const int N = 3;

  // Assemble a string to describe the test
  std::string name() { return std::string("MooneyRivlin"); }

  void get_point(Input &x) {
    x.set_rand();
    for (int i = 0; i < 9; i++) {
      x[i] *= 0.05;
    }
  }

  Output eval(const Input& x) {
    Mat<T, N, N> Uxi, J;
    x.get_values(Uxi, J);
    T W;

    // The intermediary values
    Mat<T, N, N> Jinv, Ux, F;
    SymMat<T, N> B;
    T detF, trB, trB2, I2, I1bar, I2bar;
    T t0, inv, t1;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    const T C1(0.1), C2(0.23);

    MatInv(J, Jinv);                     // Jinv = J^{-1}
    MatMatMult(Uxi, Jinv, Ux);           // Ux = Uxi * Jinv
    MatSum(Id, Ux, F);                   // F = I + Ux
    MatDet(F, detF);                     // detF = det(F)
    SymMatRK(F, B);                      // B = F * F^{T}
    MatTrace(B, trB);                    // trB = tr(B)
    SymMatTrace(B, B, trB2);             // trB2 = tr(B * B)
    Mult(trB, trB, t0);                  // t0 = trB * trB
    Sum(T(0.5), t0, T(-0.5), trB2, I2);  // I2 = 0.5 * (trB * trB - tr(B * B))
    Pow(detF, T(-2.0 / 3.0), inv);       // inv = (detF)^{-2/3}
    Mult(inv, trB, I1bar);               // I1bar = inv * tr(B)
    Mult(inv, I2, t1);                   // t1 = inv * I2
    Mult(inv, t1, I2bar);                // I2 = inv * t1 = inv * inv * I2
    Sum(C1, I1bar, C2, I2bar, W);        // W = C1 * I1bar + C2 * I2bar

    return MakeVarTuple<T>(W);
  }

  void deriv(const Output& seed, const Input& x, Input& g) {
    // Input
    Mat<T, N, N> Uxi0, Uxib, J0, Jb;
    x.get_values(Uxi0, J0);
    ADScalar<T> W;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    // Set the intermediary values
    Mat<T, N, N> Jinv0, Jinvb, Ux0, Uxb, F0, Fb;
    SymMat<T, N> B0, Bb;

    // The AD objects
    ADMat<Mat<T, N, N>> Uxi(Uxi0, Uxib), J(J0, Jb);
    ADMat<Mat<T, N, N>> Jinv(Jinv0, Jinvb), Ux(Ux0, Uxb), F(F0, Fb);
    ADMat<SymMat<T, N>> B(B0, Bb);

    ADScalar<T> detF, trB, trB2, I2, I1bar, I2bar;
    ADScalar<T> t0, inv, t1;

    const T C1(0.1), C2(0.23);

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  MatSum(Id, Ux, F),          // F = I + Ux
                  MatDet(F, detF),            // detF = det(F)
                  SymMatRK(F, B),             // B = F * F^{T}
                  MatTrace(B, trB),           // trB = tr(B)
                  SymMatTrace(B, B, trB2),    // trB2 = tr(B * B)
                  Mult(trB, trB, t0),         // t0 = trB * trB
                  Sum(T(0.5), t0, T(-0.5), trB2,
                      I2),  // I2 = 0.5 * (trB * trB - tr(B * B))
                  Pow(detF, T(-2.0 / 3.0), inv),  // inv = (detF)^{-2/3}
                  Mult(inv, trB, I1bar),          // I1bar = inv * tr(B)
                  Mult(inv, I2, t1),              // t1 = inv * I2
                  Mult(inv, t1, I2bar),  // I2 = inv * t1 = inv * inv * I2
                  Sum(C1, I1bar, C2, I2bar, W)  // W = C1 * I1bar + C2 * I2bar
        );

    seed.get_values(W.bvalue);
    stack.reverse();
    g.set_values(Uxib, Jb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    // The AD objects
    A2DMat<Mat<T, N, N>> Uxi, J;
    A2DScalar<T> W;

    // Intermediate values
    A2DMat<Mat<T, N, N>> Jinv, Ux, F;
    A2DMat<SymMat<T, N>> B;
    A2DScalar<T> detF, trB, trB2, I2, I1bar, I2bar;
    A2DScalar<T> t0, inv, t1;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    const T C1(0.1), C2(0.23);

    x.get_values(Uxi.value(), J.value());
    p.get_values(Uxi.pvalue(), J.pvalue());

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  MatSum(Id, Ux, F),          // F = I + Ux
                  MatDet(F, detF),            // detF = det(F)
                  SymMatRK(F, B),             // B = F * F^{T}
                  MatTrace(B, trB),           // trB = tr(B)
                  SymMatTrace(B, B, trB2),    // trB2 = tr(B * B)
                  Mult(trB, trB, t0),         // t0 = trB * trB
                  Sum(T(0.5), t0, T(-0.5), trB2,
                      I2),  // I2 = 0.5 * (trB * trB - tr(B * B))
                  Pow(detF, T(-2.0 / 3.0), inv),  // inv = (detF)^{-2/3}
                  Mult(inv, trB, I1bar),          // I1bar = inv * tr(B)
                  Mult(inv, I2, t1),              // t1 = inv * I2
                  Mult(inv, t1, I2bar),  // I2 = inv * t1 = inv * inv * I2
                  Sum(C1, I1bar, C2, I2bar, W)  // W = C1 * I1bar + C2 * I2bar
        );

    seed.get_values(W.bvalue);
    hval.get_values(W.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(Uxi.hvalue(), J.hvalue());
  }
};

template <typename T, int N>
class DefGradTest
    : public A2D::Test::A2DTest<T, T, Mat<T, N, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "DefGradTest<" << N << ">";
    return s.str();
  }

  // Evaluate the function
  Output eval(const Input& x) {
    Mat<T, N, N> Uxi, J;
    x.get_values(Uxi, J);
    T output;

    // Set the intermediary values
    Mat<T, N, N> Jinv, Ux, F;
    SymMat<T, N> E, S;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    MatInv(J, Jinv);                           // Jinv = J^{-1}
    MatMatMult(Uxi, Jinv, Ux);                 // Ux = Uxi * Jinv
    MatSum(Ux, Id, F);                         // F = I + Ux
    SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E);  // E = 0.5 * F^{T} * F
    SymIsotropic(T(0.35), T(0.51), E, S);  // S = 2 * mu * E + lam * tr(E) * I
    SymMatTrace(E, S, output);             // output = tr(E * S)

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    // Input
    Mat<T, N, N> Uxi0, Uxib, J0, Jb;
    x.get_values(Uxi0, J0);
    ADScalar<T> output;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    // Set the intermediary values
    Mat<T, N, N> Jinv0, Jinvb, Ux0, Uxb, F0, Fb;
    SymMat<T, N> E0, Eb, S0, Sb;

    // The AD objects
    ADMat<Mat<T, N, N>> Uxi(Uxi0, Uxib), J(J0, Jb);
    ADMat<Mat<T, N, N>> Jinv(Jinv0, Jinvb), Ux(Ux0, Uxb), F(F0, Fb);
    ADMat<SymMat<T, N>> E(E0, Eb), S(S0, Sb);

    auto stack = MakeStack(
        MatInv(J, Jinv),                           // Jinv = J^{-1}
        MatMatMult(Uxi, Jinv, Ux),                 // Ux = Uxi * Jinv
        MatSum(Ux, Id, F),                         // F = I + Ux
        SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E),  // E = 0.5 * F^{T} * F
        SymIsotropic(T(0.35), T(0.51), E, S),      // S = S(E)
        SymMatTrace(E, S, output));                // output = tr(E * S)

    seed.get_values(output.bvalue);
    stack.reverse();
    g.set_values(Uxib, Jb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    // The AD objects
    A2DMat<Mat<T, N, N>> Uxi, J;
    A2DScalar<T> output;
    A2DMat<Mat<T, N, N>> Jinv, Ux, F;
    A2DMat<SymMat<T, N>> E, S;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    x.get_values(Uxi.value(), J.value());
    p.get_values(Uxi.pvalue(), J.pvalue());

    auto stack = MakeStack(
        MatInv(J, Jinv),                           // Jinv = J^{-1}
        MatMatMult(Uxi, Jinv, Ux),                 // Ux = Uxi * Jinv
        MatSum(Ux, Id, F),                         // F = I + Ux
        SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E),  // E = 0.5 * F^{T} * F
        SymIsotropic(T(0.35), T(0.51), E, S),      // S = S(E)
        SymMatTrace(E, S, output));                // output = tr(E * S)

    seed.get_values(output.bvalue);
    hval.get_values(output.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(Uxi.hvalue(), J.hvalue());
  }
};

bool MatIntegrationTests(bool component, bool write_output) {
  bool passed = true;

  StrainTest<std::complex<double>, 3> test1;
  passed = passed && A2D::Test::Run(test1, component, write_output);

  DefGradTest<std::complex<double>, 3> test2;
  passed = passed && A2D::Test::Run(test2, component, write_output);

  MooneyRivlin<std::complex<double>> test3;
  passed = passed && A2D::Test::Run(test3, component, write_output);

  return passed;
}

int main() {
  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

  typedef std::function<bool(bool, bool)> TestFunc;
  std::vector<TestFunc> tests;

  tests.push_back(MatIntegrationTests);
  tests.push_back(A2D::Test::MatMatMultTestAll);
  tests.push_back(A2D::Test::MatDetTestAll);
  tests.push_back(A2D::Test::MatInvTestAll);
  tests.push_back(A2D::Test::MatTraceTestAll);
  tests.push_back(A2D::Test::MatGreenStrainTestAll);
  tests.push_back(A2D::Test::SymMatTraceTestAll);
  tests.push_back(A2D::Test::SymIsotropicTestAll);
  tests.push_back(A2D::Test::MatSumTestAll);
  tests.push_back(A2D::Test::SymMatRKTestAll);
  tests.push_back(A2D::Test::SymMatSumTestAll);
  tests.push_back(A2D::Test::VecCrossTestAll);
  tests.push_back(A2D::Test::VecNormTestAll);
  tests.push_back(A2D::Test::VecDotTestAll);
  tests.push_back(A2D::Test::VecScaleTestAll);
  tests.push_back(A2D::Test::VecNormalizeTestAll);
  tests.push_back(A2D::Test::VecSumTestAll);
  tests.push_back(A2D::Test::VecOuterTestAll);
  tests.push_back(A2D::Test::ScalarTestAll);

  bool passed = true;
  for (int i = 0; i < tests.size(); i++) {
    bool test_passed = tests[i](component, write_output);

    // If the test fails, perform a component test and print the results
    // to screen
    if (!test_passed) {
      bool comp = true;
      bool write = true;
      tests[i](comp, write);
    }

    passed = test_passed && passed;
  }

  // Convert to a fail flag
  int fail = 0;
  if (!passed) {
    fail = 1;
  }
  return fail;
}
