#include "a2dtmp.h"
#include <iostream>
#include <iomanip>
#include <complex>

int main( int argc, char *argv[] ){
  // typedef int32_t IndexType;
  typedef std::complex<double> ScalarType;
  typedef A2D::Mat<ScalarType, 3, 3> Mat3x3;
  typedef A2D::SymmMat<ScalarType, 3> SymmMat3x3;
  typedef A2D::Mat2ndDeriv<ScalarType, 3, 3> Mat2ndDeriv;

  double dh = 1e-30;

  Mat3x3 Ux, Uxb, Uxd, Uxh;
  SymmMat3x3 E, Eb, Ed, Eh;
  SymmMat3x3 S, Sb, Sd, Sh;

  A2D::A2DScalarType<ScalarType> output;
  A2D::A2DMat<Mat3x3> UxObj(&Ux, &Uxb, &Uxd, &Uxh);
  A2D::A2DMat<SymmMat3x3> EObj(&E, &Eb, &Ed, &Eh);
  A2D::A2DMat<SymmMat3x3> SObj(&S, &Sb, &Sd, &Sh);

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

  A2D::Mat3x3GreenStrain<Mat3x3, SymmMat3x3> strain(UxObj, EObj);
  A2D::Symm3x3IsotropicConstitutive<ScalarType, SymmMat3x3, SymmMat3x3> constitutive(mu, lambda, EObj, SObj);
  A2D::Symm3x3SymmMultTrace<SymmMat3x3, SymmMat3x3, ScalarType> trace(SObj, EObj, output);

  output.bvalue = 1.0;

  trace.reverse();
  constitutive.reverse();
  strain.reverse();

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

  strain.hforward();
  constitutive.hforward();
  trace.hproduct();
  constitutive.hreverse();
  strain.hreverse();

  strain.hproduct();

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      double res = Uxh(i, j).real();
      double fd = Uxb(i, j).imag()/dh;
      double error = (res - fd)/fd;

      std::cout << i << ", " << j << " result: " << std::setw(20) << res <<
        " fd: " << std::setw(20) << fd << " error: " << std::setw(20) <<  error << std::endl;
    }
  }

  return (0);
}