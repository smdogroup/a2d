#include "a2dtmp.h"
#include <iostream>
#include <iomanip>
#include <complex>

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
int main( int argc, char *argv[] ){
  // typedef int32_t IndexType;
  typedef std::complex<double> ScalarType;
  typedef A2D::Mat<ScalarType, 3, 3> Mat3x3;
  typedef A2D::SymmMat<ScalarType, 3> SymmMat3x3;
  typedef A2D::SymmTensor<ScalarType, 3, 3> SymmTensor3x3;

  double dh = 1e-30;

  Mat3x3 J0, Jb, Jp, Jh;
  Mat3x3 Jinv0, Jinvb, Jinvp, Jinvh;
  Mat3x3 Uxi0, Uxib, Uxip, Uxih;
  Mat3x3 Ux0, Uxb, Uxp, Uxh;
  SymmMat3x3 E0, Eb, Ep, Eh;
  SymmMat3x3 S0, Sb, Sp, Sh;

  A2D::A2DMat<Mat3x3> J(J0, Jb, Jp, Jh);
  A2D::A2DMat<Mat3x3> Jinv(Jinv0, Jinvb, Jinvp, Jinvh);
  A2D::A2DMat<Mat3x3> Uxi(Uxi0, Uxib, Uxip, Uxih);
  A2D::A2DMat<Mat3x3> Ux(Ux0, Uxb, Uxp, Uxh);
  A2D::A2DMat<SymmMat3x3> S(S0, Sb, Sp, Sh);
  A2D::A2DMat<SymmMat3x3> E(E0, Eb, Ep, Eh);
  A2D::A2DScalar<ScalarType> output;

  ScalarType mu(0.2533), lambda(0.71236);

  // Set random values
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Uxi0(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
      J0(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;

      Uxip(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
    }
  }

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Uxi0(i, j) = Uxi0(i, j) + ScalarType(0.0, dh) * Uxip(i, j);
    }
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

  ScalarType result = 0.0;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      result += Uxib(i, j) * Uxip(i, j) + Jb(i, j) * Jp(i, j);
    }
  }

  // double forward = output.valueb.real();
  double res = result.real();
  double fd = output.value.imag()/dh;
  double error = (res - fd)/fd;
  std::cout << "result: " << std::setw(20) << res << " fd: " << std::setw(20) << fd
    << " error: " << std::setw(20) << error << std::endl;

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

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      double res = Uxih(i, j).real();
      double fd = Uxib(i, j).imag()/dh;
      double error = (res - fd)/fd;

      std::cout << i << ", " << j << " result: " << std::setw(20) << res <<
        " fd: " << std::setw(20) << fd << " error: " << std::setw(20) <<  error << std::endl;
    }
  }

  return (0);
}
