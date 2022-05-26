#include <complex>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"

/*
  The first derivative of a function f(y(x)) is

  df/dx = df/dy * dy/dx

  The second derivative of the same function is

  d^2f/dx^2 * px = df/dy * (d^2y/dx^2 * p) + d^2f/dy^2 * (dy/dx * p) * (dy/dx)

  # Definition of partial derivatives bar{x} = bx
  bx = df/dx
  by = df/dy

  # Definition of directional derivatives p_x, p_y
  px = p
  py = dy/dx * p

  # Definition of projected second derivative
  hx = d^2f/dx^2 * px
  hy = d^2f/dy^2 * py

  # The projected second derivative requires the computation
  hx = by * (d^2y/dx^2 * px) + hy * (dy/dx)
*/
int main(int argc, char* argv[]) {
  // typedef int32_t IndexType;
  typedef std::complex<double> ScalarType;
  typedef A2D::Mat<ScalarType, 3, 3> Mat3x3;
  typedef A2D::SymmMat<ScalarType, 3> SymmMat3x3;
  typedef A2D::SymmTensor<ScalarType, 3, 3> SymmTensor3x3;

  double dh = 1e-30;

  // Compute all components of the Jacobian matrix at once
  const int N = 9;
  Mat3x3 J0, Jb;
  Mat3x3 Jinv0, Jinvb;
  Mat3x3 Uxi0, Uxib;
  Mat3x3 Ux0, Uxb;
  SymmMat3x3 E0, Eb;
  SymmMat3x3 S0, Sb;

  A2D::A2DMat<N, Mat3x3> J(J0, Jb);
  A2D::A2DMat<N, Mat3x3> Jinv(Jinv0, Jinvb);
  A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
  A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
  A2D::A2DMat<N, SymmMat3x3> S(S0, Sb);
  A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
  A2D::A2DScalar<N, ScalarType> output;

  ScalarType mu(0.2533), lambda(0.71236);

  // Set random values
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      Uxi0(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      J0(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
    }
  }

  Mat3x3 P, result;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      P(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      Uxi0(i, j) = Uxi0(i, j) + ScalarType(0.0, dh) * P(i, j);
    }
  }

  // Set up the perturb value
  for (int k = 0; k < N; k++) {
    Mat3x3& Up = Uxi.pvalue(k);
    Up(k / 3, k % 3) = 1.0;
  }

  auto jinv = A2D::Mat3x3Inverse(J, Jinv);
  auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
  auto strain = A2D::Mat3x3GreenStrain(Ux, E);
  auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
  auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

  output.bvalue = 1.0;

  trace.reverse();
  constitutive.reverse();
  strain.reverse();
  mult.reverse();
  jinv.reverse();

  jinv.hforward();
  mult.hforward();
  strain.hforward();
  constitutive.hforward();
  trace.hforward();
  trace.hreverse();
  constitutive.hreverse();
  strain.hreverse();
  mult.hreverse();
  jinv.hreverse();

  // Compute the product of the Hessian with P
  for (int k = 0; k < N; k++) {
    const Mat3x3& Uh = Uxi.hvalue(k);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        result(i, j) += Uh(i, j) * P(k / 3, k % 3);
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double res = result(i, j).real();
      double fd = Uxib(i, j).imag() / dh;
      double error = (res - fd) / fd;

      std::cout << i << ", " << j << " result: " << std::setw(20) << res
                << " fd: " << std::setw(20) << fd << " error: " << std::setw(20)
                << error << std::endl;
    }
  }

  return (0);
}
