#include "a2dtmp.h"

namespace A2D {

template <class MatType, class... MixedTypes>
class A2HMat {
 public:
  A2HMat(MatType& A, MatType& Ab, MixedTypes... mixed) {}

  MatType& A;
  MatType& Ab;
};

template <class A>
class A2HExpression {
 public:
  A& cast() { return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void reverse() { cast().reverse(); }
  void hreverse() { cast().hreverse(); }
};

template <class MatType, class EMatType>
class A2HMat3x3GreenStrainExpr
    : public A2HExpression<A2HMat3x3GreenStrainExpr<MatType, EMatType> > {
 public:
  A2HMat3x3GreenStrainExpr(A2HMat<MatType>& Ux, A2HMat<EMatType>& E)
      : Ux(Ux), E(E) {
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0) +
                                Ux(2, 0) * Ux(2, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1) +
                                Ux(2, 1) * Ux(2, 1));
    E(2, 2) = Ux(2, 2) + 0.5 * (Ux(0, 2) * Ux(0, 2) + Ux(1, 2) * Ux(1, 2) +
                                Ux(2, 2) * Ux(2, 2));

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) +
                     Ux(1, 0) * Ux(1, 1) + Ux(2, 0) * Ux(2, 1));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0) + Ux(0, 0) * Ux(0, 2) +
                     Ux(1, 0) * Ux(1, 2) + Ux(2, 0) * Ux(2, 2));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1) + Ux(0, 1) * Ux(0, 2) +
                     Ux(1, 1) * Ux(1, 2) + Ux(2, 1) * Ux(2, 2));
  }

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) +
                 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1) +
                 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +
                 Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) +
                 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1) +
                 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +
                 Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) += Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) + Ux(2, 1) * Eb(1, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +
                 (Ux(2, 2) + 1.0) * Eb(2, 2);
  }

  // df/dy * (d^2y/dx^2) + d^2f/dy^2 * (dy/dx) * (dy/dx)

  void hreverse() {
    // Compute the required second derivatives
    SymmMat& Uh;

    Uh(0, 0, 0, 0) += Eb(0, 0);
    Uh(0, 0, 0, 1) += 0.5 * Eb(0, 1);
    Uh(0, 0, 0, 2) += 0.5 * Eb(0, 2);

    Uh(0, 1, 0, 0) += 0.5 * Eb(0, 1);
    Uh(0, 1, 0, 1) += Eb(1, 1);
    Uh(0, 1, 0, 2) += 0.5 * Eb(1, 2);

    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) +
                 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1) +
                 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +
                 Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) +
                 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1) +
                 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +
                 Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) += Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) + Ux(2, 1) * Eb(1, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +
                 (Ux(2, 2) + 1.0) * Eb(2, 2);

    Uh(i, j);

    // Compute the other partial derivative terms - if any are provided
  }
};

/*
  This code is designed to directly compute the Hessian of the
  following expression:

  E = 0.5*(U + U^{T} + U^{T} * U)
  S = 2 * mu * E + I tr(E) * I
  tr(S * E)
*/
int main(int argc, char* argv[]) {
  typedef std::complex<double> ScalarType;
  typedef A2D::Mat<ScalarType, 3, 3> Mat3x3;
  typedef A2D::SymmMat<ScalarType, 3> SymmMat3x3;
  typedef A2D::SymmTensor<ScalarType, 3, 3> SymmTensor3x3;

  Mat3x3 U0, Ub;
  SymmTensor3x3 Uh;  // The Hessian matrix

  Symm3x3 E0, Eb;
  SymmSymmTensor3x3 Eh;  // Hessian of the strain

  Symm3x3 S0, Sb;
  SymmSymmTensor3x3 Sh;  // Hessian of the stress

  // Mixed partial of the stress/strain Hessian
  SymmTensor3x3 ESh;

  ScalarType mu(0.2533), lambda(0.71236);

  A2D::A2HMat<Mat3x3> U(U0, Ub, Uh);
  A2D::A2HMat<SymmMat3x3> S(S0, Sb, ESh);
  A2D::A2HMat<SymmMat3x3> E(E0, Eb, ESh);
  A2D::A2DScalar<ScalarType> output;

  auto strain = A2D::Mat3x3GreenStrain(U, E);
  auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
  auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

  output.bvalue = 1.0;

  // Compute the derivative
  trace.reverse();
  constitutive.reverse();
  strain.reverse();

  // Compute the full Hessian matrices
  trace.hreverse();
  constitutive.hreverse();
  strain.hreverse();

  return (0);
}