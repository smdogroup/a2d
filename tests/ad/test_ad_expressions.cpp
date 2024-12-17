#include <functional>
#include <vector>

#include "a2dcore.h"

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
  Output eval(const Input &x) {
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
    SymMatMultTrace(E, S, output);         // output = tr(E * S)

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    // The AD objects
    ADObj<T> output;
    ADObj<Mat<T, N, N>> Uxi, J, Jinv, Ux;
    ADObj<SymMat<T, N>> E1, E2, E, S;

    x.get_values(Uxi.value(), J.value());

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  SymMatSum(T(0.5), Ux, E1),  // E1 = 0.5 * (Ux + Ux^{T})
                  SymMatRK<MatOp::TRANSPOSE>(Ux, E2),    // E2 = Ux^{T} * Ux
                  MatSum(T(1.0), E1, T(0.5), E2, E),     // E = E1 + 0.5 * E2
                  SymIsotropic(T(0.35), T(0.51), E, S),  // S = S(E)
                  SymMatMultTrace(E, S, output));

    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(Uxi.bvalue(), J.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // The AD objects
    A2DObj<Mat<T, N, N>> Uxi, J;
    A2DObj<T> output;
    A2DObj<Mat<T, N, N>> Jinv, Ux;
    A2DObj<SymMat<T, N>> E1, E2, E, S;

    x.get_values(Uxi.value(), J.value());
    p.get_values(Uxi.pvalue(), J.pvalue());

    auto stack =
        MakeStack(MatInv(J, Jinv),            // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),  // Ux = Uxi * Jinv
                  SymMatSum(T(0.5), Ux, E1),  // E1 = 0.5 * (Ux + Ux^{T})
                  SymMatRK<MatOp::TRANSPOSE>(Ux, E2),    // E2 = Ux^{T} * Ux
                  MatSum(T(1.0), E1, T(0.5), E2, E),     // E = E1 + 0.5 * E2
                  SymIsotropic(T(0.35), T(0.51), E, S),  // S = S(E)
                  SymMatMultTrace(E, S, output));

    seed.get_values(output.bvalue());
    hval.get_values(output.hvalue());
    stack.hproduct();
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

  Output eval(const Input &x) {
    Mat<T, N, N> Uxi, J;
    x.get_values(Uxi, J);
    T W;

    // The intermediary values
    Mat<T, N, N> Jinv, Ux, F;
    SymMat<T, N> B;
    T detF, trB, trB2;
    T inv, I2;
    T I1bar, I2bar;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    const T C1(0.1), C2(0.23);

    MatInv(J, Jinv);                // Jinv = J^{-1}
    MatMatMult(Uxi, Jinv, Ux);      // Ux = Uxi * Jinv
    MatSum(Id, Ux, F);              // F = I + Ux
    MatDet(F, detF);                // detF = det(F)
    SymMatRK(F, B);                 // B = F * F^{T}
    MatTrace(B, trB);               // trB = tr(B)
    SymMatMultTrace(B, B, trB2);    // trB2 = tr(B * B)
    inv = pow(detF, -2.0 / 3.0);    // inv = detF^{-2/3}
    I2 = 0.5 * (trB * trB - trB2);  // I2 = 0.5 * (trB * trB - tr(B * B))
    I1bar = inv * trB;              // I1bar = inv * I1 = inv * tr(B)
    I2bar = inv * inv * I2;         // I2bar = inv^2 * I2
    W = C1 * (I1bar - 3.0) + C2 * (I2bar - 3.0);

    return MakeVarTuple<T>(W);
  }

  void deriv(const Output &seed, const Input &x, Input &g) {
    // The AD objects
    ADObj<Mat<T, N, N>> Uxi, J;
    ADObj<T> W;

    // Intermediate values
    ADObj<Mat<T, N, N>> Jinv, Ux, F;
    ADObj<SymMat<T, N>> B;
    ADObj<T> detF, trB, trB2;
    ADObj<T> inv, I2;
    ADObj<T> I1bar, I2bar;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    const T C1(0.1), C2(0.23);

    x.get_values(Uxi.value(), J.value());

    auto stack =
        MakeStack(MatInv(J, Jinv),                   // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),         // Ux = Uxi * Jinv
                  MatSum(Id, Ux, F),                 // F = I + Ux
                  MatDet(F, detF),                   // detF = det(F)
                  SymMatRK(F, B),                    // B = F * F^{T}
                  MatTrace(B, trB),                  // trB = tr(B)
                  SymMatMultTrace(B, B, trB2),       // trB2 = tr(B * B)
                  Eval(pow(detF, -2.0 / 3.0), inv),  // inv = detF^{-2/3}
                  Eval(0.5 * (trB * trB - trB2),
                       I2),                // I2 = 0.5 * (trB * trB - tr(B * B))
                  Eval(inv * trB, I1bar),  // I1bar = inv * I1 = inv * tr(B)
                  Eval(inv * inv * I2, I2bar),  // I2bar = inv^2 * I2
                  Eval(C1 * (I1bar - 3.0) + C2 * (I2bar - 3.0), W));

    seed.get_values(W.bvalue());
    stack.reverse();
    g.set_values(Uxi.bvalue(), J.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // The AD objects
    A2DObj<Mat<T, N, N>> Uxi, J;
    A2DObj<T> W;

    // Intermediate values
    A2DObj<Mat<T, N, N>> Jinv, Ux, F;
    A2DObj<SymMat<T, N>> B;
    A2DObj<T> detF, trB, trB2;
    A2DObj<T> inv, I2;
    A2DObj<T> I1bar, I2bar;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    const T C1(0.1), C2(0.23);

    x.get_values(Uxi.value(), J.value());
    p.get_values(Uxi.pvalue(), J.pvalue());

    auto stack =
        MakeStack(MatInv(J, Jinv),                   // Jinv = J^{-1}
                  MatMatMult(Uxi, Jinv, Ux),         // Ux = Uxi * Jinv
                  MatSum(Id, Ux, F),                 // F = I + Ux
                  MatDet(F, detF),                   // detF = det(F)
                  SymMatRK(F, B),                    // B = F * F^{T}
                  MatTrace(B, trB),                  // trB = tr(B)
                  SymMatMultTrace(B, B, trB2),       // trB2 = tr(B * B)
                  Eval(pow(detF, -2.0 / 3.0), inv),  // inv = detF^{-2/3}
                  Eval(0.5 * (trB * trB - trB2),
                       I2),                // I2 = 0.5 * (trB * trB - tr(B * B))
                  Eval(inv * trB, I1bar),  // I1bar = inv * I1 = inv * tr(B)
                  Eval(inv * inv * I2, I2bar),  // I2bar = inv^2 * I2
                  Eval(C1 * (I1bar - 3.0) + C2 * (I2bar - 3.0), W));

    seed.get_values(W.bvalue());
    hval.get_values(W.hvalue());
    stack.hproduct();
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
  Output eval(const Input &x) {
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
    SymMatMultTrace(E, S, output);         // output = tr(E * S)

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    // The AD objects
    ADObj<Mat<T, N, N>> Uxi, J;
    ADObj<T> output;
    ADObj<Mat<T, N, N>> Jinv, Ux, F;
    ADObj<SymMat<T, N>> E, S;

    // Set the entries of the identity matrix
    Mat<T, N, N> Id;
    for (int k = 0; k < N; k++) {
      Id(k, k) = 1.0;
    }

    x.get_values(Uxi.value(), J.value());

    auto stack = MakeStack(
        MatInv(J, Jinv),                           // Jinv = J^{-1}
        MatMatMult(Uxi, Jinv, Ux),                 // Ux = Uxi * Jinv
        MatSum(Ux, Id, F),                         // F = I + Ux
        SymMatRK<MatOp::TRANSPOSE>(T(0.5), F, E),  // E = 0.5 * F^{T} * F
        SymIsotropic(T(0.35), T(0.51), E, S),      // S = S(E)
        SymMatMultTrace(E, S, output));            // output = tr(E * S)

    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(Uxi.bvalue(), J.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // The AD objects
    A2DObj<Mat<T, N, N>> Uxi, J;
    A2DObj<T> output;
    A2DObj<Mat<T, N, N>> Jinv, Ux, F;
    A2DObj<SymMat<T, N>> E, S;

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
        SymMatMultTrace(E, S, output));            // output = tr(E * S)

    seed.get_values(output.bvalue());
    hval.get_values(output.hvalue());
    stack.hproduct();
    h.set_values(Uxi.hvalue(), J.hvalue());
  }
};

template <typename T, int N>
class HExtractTest : public A2D::Test::A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "HExtract<" << N << ">";
    return s.str();
  }

  // Evaluate the function
  Output eval(const Input &x) {
    const T mu(0.197), lambda(0.839);
    Mat<T, N, N> Ux;
    x.get_values(Ux);
    T output;

    // Symmetric matrices
    SymMat<T, N> E, S;

    MatGreenStrain<GreenStrainType::NONLINEAR>(Ux, E);  // E = E(Ux)
    SymIsotropic(mu, lambda, E, S);                     // S = S(E)
    SymMatMultTrace(E, S, output);                      // output = tr(E * S)

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    // Input
    const T mu(0.197), lambda(0.839);
    ADObj<Mat<T, N, N>> Ux;
    x.get_values(Ux.value());
    ADObj<T> output;

    // Symmetric matrices
    ADObj<SymMat<T, N>> E, S;

    auto stack = MakeStack(
        MatGreenStrain<GreenStrainType::NONLINEAR>(Ux, E),  // E = E(Ux)
        SymIsotropic(mu, lambda, E, S),                     // S = S(E)
        SymMatMultTrace(E, S, output));  // output = tr(E * S)

    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(Ux.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // Input
    const T mu(0.197), lambda(0.839);
    A2DObj<Mat<T, N, N>> Ux;
    x.get_values(Ux.value());
    A2DObj<T> output;

    // Symmetric matrices
    A2DObj<SymMat<T, N>> E, S;

    auto stack = MakeStack(
        MatGreenStrain<GreenStrainType::NONLINEAR>(Ux, E),  // E = E(Ux)
        SymIsotropic(mu, lambda, E, S),                     // S = S(E)
        SymMatMultTrace(E, S, output));  // output = tr(E * S)

    // Number of components in the derivative
    constexpr index_t ncomp = N * N;

    output.bvalue() = hval[0];
    stack.reverse();

    auto g = MakeTieTuple<T, ADseed::b>(Ux);
    for (int i = 0; i < ncomp; i++) {
      h[i] = g[i];
    }

    auto temps = MakeTieTuple<T, ADseed::b>(Ux, S, E);
    temps.zero();

    // Set the seeds for the second-order part
    output.hvalue() = 0.0;
    output.bvalue() = seed[0];

    // Create data for extracting the Hessian-vector product
    auto in = MakeTieTuple<T, ADseed::p>(Ux);
    auto out = MakeTieTuple<T, ADseed::h>(Ux);

    // Extract the Hessian matrix
    Mat<T, ncomp, ncomp> jac;  // Symmetric only if hval = 0.0
    stack.hextract(in, out, jac);

    // Mupltiply the outputs
    for (int i = 0; i < ncomp; i++) {
      for (int j = 0; j < ncomp; j++) {
        h[i] += jac(i, j) * p[j];
      }
    }
  }
};

template <typename T>
class VonMisesPenaltyTest : public A2D::Test::A2DTest<T, T, T, Mat<T, 3, 3>> {
 public:
  using Input = VarTuple<T, T, Mat<T, 3, 3>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return "VonMisesPenalty"; }

  // Evaluate the function
  Output eval(const Input &x) {
    // Set constants
    double q = 5.0;
    double design_stress = 135.0;
    T mu = 3.374, lambda = 9.173;

    // Get the variables
    T rho;
    Mat<T, 3, 3> Ux;
    x.get_values(rho, Ux);

    // Intermediaries
    SymMat<T, 3> E, S;
    T trS, trSS, trS2;

    // Compute the strain and stress
    MatGreenStrain<GreenStrainType::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);

    // Compute the von Mises stress = sqrt(1.5 * tr(S * S) - 0.5 * tr(S)**2)
    MatTrace(S, trS);
    SymMatMultTrace(S, S, trSS);
    T vm = sqrt(1.5 * trSS - 0.5 * trS * trS);
    T relaxed_stress = (vm * ((q + 1.0) / (q * rho + 1.0)));
    T failure_index = relaxed_stress / design_stress;

    return MakeVarTuple<T>(failure_index);
  }

  void deriv(const Output &seed, const Input &x, Input &g) {
    // Set constants
    double q = 5.0;
    double design_stress = 135.0;
    T mu = 3.374, lambda = 9.173;

    // Get the variables
    ADObj<T> rho;
    ADObj<Mat<T, 3, 3>> Ux;
    x.get_values(rho.value(), Ux.value());

    // Intermediaries
    ADObj<SymMat<T, 3>> E, S;
    ADObj<T> trS, trSS;
    ADObj<T> vm, relaxed_stress, failure_index;

    auto stack = MakeStack(
        // Compute the strain and stress
        MatGreenStrain<GreenStrainType::LINEAR>(Ux, E),
        SymIsotropic(mu, lambda, E, S), MatTrace(S, trS),
        SymMatMultTrace(S, S, trSS),

        // Evaluate the von Mises stress output
        Eval(sqrt(1.5 * trSS - 0.5 * trS * trS), vm),
        Eval((vm * ((q + 1.0) / (q * rho + 1.0))), relaxed_stress),
        Eval(relaxed_stress / design_stress, failure_index));

    seed.get_values(failure_index.bvalue());
    stack.reverse();
    g.set_values(rho.bvalue(), Ux.bvalue());
  }

  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // Set constants
    double q = 5.0;
    double design_stress = 135.0;
    T mu = 3.374, lambda = 9.173;

    // Get the variables
    A2DObj<T> rho;
    A2DObj<Mat<T, 3, 3>> Ux;
    x.get_values(rho.value(), Ux.value());
    p.get_values(rho.pvalue(), Ux.pvalue());

    // Intermediaries
    A2DObj<SymMat<T, 3>> E, S;
    A2DObj<T> trS, trSS;
    A2DObj<T> vm, relaxed_stress, failure_index;

    auto stack = MakeStack(
        // Compute the strain and stress
        MatGreenStrain<GreenStrainType::LINEAR>(Ux, E),
        SymIsotropic(mu, lambda, E, S), MatTrace(S, trS),
        SymMatMultTrace(S, S, trSS),

        // Evaluate the von Mises stress output
        Eval(sqrt(1.5 * trSS - 0.5 * trS * trS), vm),
        Eval((vm * ((q + 1.0) / (q * rho + 1.0))), relaxed_stress),
        Eval(relaxed_stress / design_stress, failure_index));

    seed.get_values(failure_index.bvalue());
    hval.get_values(failure_index.hvalue());
    stack.hproduct();
    h.set_values(rho.hvalue(), Ux.hvalue());
  }
};

/* Test AD and A2D for a diamond-like computational graph */
template <typename T, int N>
class DiamondGraphTest : public A2D::Test::A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "DiamondGraphTest<" << N << ">";
    return s.str();
  }

  // Evaluate the function
  Output eval(const Input &x) {
    Mat<T, N, N> J;
    x.get_values(J);
    T Jdet, tr, output;

    // Set the intermediary values
    Mat<T, N, N> Jinv, Ux;
    SymMat<T, N> E1, E2, E, S;

    MatDet(J, Jdet);     // Jdet = det(J)
    MatInv(J, Jinv);     // Jinv = J^{-1}
    MatTrace(Jinv, tr);  // tr = tr(Jinv)
    output = tr + Jdet;  // output = tr + Jdet

    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    // The AD objects
    ADObj<Mat<T, N, N>> J, Jinv;
    ADObj<T> Jdet, tr, output;

    x.get_values(J.value());

    auto stack = MakeStack(MatDet(J, Jdet),           // Jdet = det(J)
                           MatInv(J, Jinv),           // Jinv = J^{-1}
                           MatTrace(Jinv, tr),        // tr = tr(Jinv)
                           Eval(tr + Jdet, output));  // output = tr + Jdet

    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(J.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    // The AD objects
    A2DObj<Mat<T, N, N>> J, Jinv;
    A2DObj<T> Jdet, tr, output;

    x.get_values(J.value());
    p.get_values(J.pvalue());

    auto stack = MakeStack(MatDet(J, Jdet),           // Jdet = det(J)
                           MatInv(J, Jinv),           // Jinv = J^{-1}
                           MatTrace(Jinv, tr),        // tr = tr(Jinv)
                           Eval(tr + Jdet, output));  // output = tr + Jdet

    seed.get_values(output.bvalue());
    hval.get_values(output.hvalue());
    stack.hproduct();
    h.set_values(J.hvalue());
  }
};

bool MatIntegrationTests(bool component, bool write_output) {
  bool passed = true;

  StrainTest<A2D_complex_t<double>, 3> test1;
  passed = passed && A2D::Test::Run(test1, component, write_output);

  DefGradTest<A2D_complex_t<double>, 3> test2;
  passed = passed && A2D::Test::Run(test2, component, write_output);

  MooneyRivlin<A2D_complex_t<double>> test3;
  passed = passed && A2D::Test::Run(test3, component, write_output);

  HExtractTest<A2D_complex_t<double>, 3> test4;
  passed = passed && A2D::Test::Run(test4, component, write_output);

  VonMisesPenaltyTest<A2D_complex_t<double>> test5;
  passed = passed && A2D::Test::Run(test5, component, write_output);

  DiamondGraphTest<A2D_complex_t<double>, 3> test6;
  passed = passed && A2D::Test::Run(test6, component, write_output);

  return passed;
}

int main(int argc, char *argv[]) {
  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

  // Check for the write_output flag
  for (int i = 0; i < argc; i++) {
    std::string str(argv[i]);
    if (str.compare("--write_output") == 0) {
      write_output = true;
    }
    if (str.compare("--component") == 0) {
      component = true;
    }
  }

  typedef std::function<bool(bool, bool)> TestFunc;
  std::vector<TestFunc> tests;

  tests.push_back(MatIntegrationTests);
  tests.push_back(A2D::Test::MatMatMultTestAll);
  tests.push_back(A2D::Test::MatVecMultTestAll);
  tests.push_back(A2D::Test::MatDetTestAll);
  tests.push_back(A2D::Test::MatInvTestAll);
  tests.push_back(A2D::Test::MatTraceTestAll);
  tests.push_back(A2D::Test::MatGreenStrainTestAll);
  tests.push_back(A2D::Test::SymMatMultTraceTestAll);
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
  tests.push_back(A2D::Test::SymEigsTestAll);
  tests.push_back(A2D::Test::QuaternionMatrixTestAll);
  tests.push_back(A2D::Test::QuaternionAngularVelocityTestAll);
  tests.push_back(A2D::Test::VecHadamardTestAll);

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
