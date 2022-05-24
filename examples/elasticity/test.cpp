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
  typedef A2D::Mat2ndDeriv<ScalarType, 3, 3> Mat2ndDeriv;

  double dh = 1e-30;

  /*
  Mat3x3 Jd, Uxid;

  Mat3x3 J, Jb;
  Mat3x3 Jinv, Jinvb;
  Mat3x3 Uxi, Uxib;
  Mat3x3 Ux, Uxb;
  SymmMat3x3 E, Eb;
  SymmMat3x3 S, Sb;

  A2D::ADMat<Mat3x3> JObj(J, Jb);
  A2D::ADMat<Mat3x3> JinvObj(Jinv, Jinvb);
  A2D::ADMat<Mat3x3> UxiObj(Uxi, Uxib);
  A2D::ADMat<Mat3x3> UxObj(Ux, Uxb);
  A2D::ADMat<SymmMat3x3> EObj(E, Eb);
  A2D::ADMat<SymmMat3x3> SObj(S, Sb);

  A2D::ADScalar<ScalarType> output;

  ScalarType mu(0.2533), lambda(0.71236);

  // Set random values
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Uxi(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
      J(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;

      Uxid(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
      Jd(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
    }
  }

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Uxi(i, j) = Uxi(i, j) + ScalarType(0.0, dh) * Uxid(i, j);
      J(i, j) = J(i, j) + ScalarType(0.0, dh) * Jd(i, j);
    }
  }

  auto jinv = A2D::Mat3x3Inverse(JObj, JinvObj);
  auto mult = A2D::Mat3x3MatMult(UxiObj, JinvObj, UxObj);
  auto strain = A2D::Mat3x3GreenStrain(UxObj, EObj);
  auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, EObj, SObj);
  auto trace = A2D::Symm3x3SymmMultTrace(SObj, EObj, output);

  output.bvalue = 1.0;

  trace.reverse();
  constitutive.reverse();
  strain.reverse();
  mult.reverse();
  jinv.reverse();

  ScalarType result = 0.0;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      result += Uxib(i, j) * Uxid(i, j) + Jb(i, j) * Jd(i, j);
    }
  }

  // double forward = output.valueb.real();
  double res = result.real();
  double fd = output.value.imag()/dh;
  double error = (res - fd)/fd;
  std::cout << "result: " << std::setw(20) << res << " fd: " << std::setw(20) << fd
    << " error: " << std::setw(20) << error << std::endl;

  // */
  Mat3x3 Ux, Uxb, Uxd, Uxh;
  Mat3x3 U, Ub, Ud, Uh;
  SymmMat3x3 E, Eb, Ed, Eh;
  SymmMat3x3 S, Sb, Sd, Sh;

  A2D::A2DScalar<ScalarType> output;
  A2D::A2DMat<Mat3x3> UxObj(Ux, Uxb, Uxd, Uxh);
  A2D::A2DMat<Mat3x3> UObj(U, Ub, Ud, Uh);
  A2D::A2DMat<SymmMat3x3> EObj(E, Eb, Ed, Eh);
  A2D::A2DMat<SymmMat3x3> SObj(S, Sb, Sd, Sh);

  ScalarType mu(0.2533), lambda(0.71236);

  // Set random values for Ux
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Ux(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
      Uxd(i, j) = -1.0 + 2.0 * rand()/RAND_MAX;
    }
  }

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      Ux(i, j) = Ux(i, j) + ScalarType(0.0, dh) * Uxd(i, j);
    }
  }

  auto inv = A2D::Mat3x3Inverse(UxObj, UObj);
  auto strain = A2D::Mat3x3GreenStrain(UObj, EObj);
  auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, EObj, SObj);
  auto trace = A2D::Symm3x3SymmMultTrace(SObj, EObj, output);

  output.bvalue = 1.0;

  trace.reverse();
  constitutive.reverse();
  strain.reverse();
  inv.reverse();

  ScalarType result = 0.0;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      result += Uxb(i, j) * Uxd(i, j);
    }
  }

  // double forward = output.valueb.real();
  double res = result.real();
  double fd = output.value.imag()/dh;
  double error = (res - fd)/fd;
  std::cout << "result: " << std::setw(20) << res << " fd: " << std::setw(20) << fd
    << " error: " << std::setw(20) << error << std::endl;

  inv.hforward();
  strain.hforward();
  constitutive.hforward();
  trace.hforward();
  trace.hreverse();
  constitutive.hreverse();
  strain.hreverse();
  inv.hreverse();

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      double res = Uxh(i, j).real();
      double fd = Uxb(i, j).imag()/dh;
      double error = (res - fd)/fd;

      std::cout << i << ", " << j << " result: " << std::setw(20) << res <<
        " fd: " << std::setw(20) << fd << " error: " << std::setw(20) <<  error << std::endl;
    }
  }

  // */

  return (0);
}
