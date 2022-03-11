#ifndef A2D_MAT_OPS_H
#define A2D_MAT_OPS_H

#include "a2dobjs.h"
#include "a2dmatcore.h"

/*
  Matrix trace operation
*/
class Symm3x3Trace {
public:
  Symm3x3Trace( Symm3x3& A, Scalar& alpha ){
    alpha.value = A.A[0] + A.A[3] + A.A[5];
  }
};

class ADSymm3x3Trace {
public:
  ADSymm3x3Trace( ADSymm3x3& A, ADScalar& alpha ) : A(A), alpha(alpha) {
    alpha.value = A.A[0] + A.A[3] + A.A[5];
  }
  void forward(){
    alpha.valued = A.Ad[0] + A.Ad[3] + A.Ad[5];
  }
  void reverse(){
    A.Ad[0] += alpha.valued;
    A.Ad[3] += alpha.valued;
    A.Ad[5] += alpha.valued;
  }

  ADSymm3x3& A;
  ADScalar& alpha;
};

class Mat3x3Trace {
public:
  Mat3x3Trace( Mat3x3& A, Scalar& alpha ){
    alpha.value = A.A[0] + A.A[4] + A.A[8];
  }
};

class ADMat3x3Trace {
public:
  ADMat3x3Trace( ADMat3x3& A, ADScalar& alpha ) : A(A), alpha(alpha) {
    alpha.value = A.A[0] + A.A[4] + A.A[8];
  }
  void forward(){
    alpha.valued = A.Ad[0] + A.Ad[4] + A.Ad[8];
  }
  void reverse(){
    A.Ad[0] += alpha.valued;
    A.Ad[4] += alpha.valued;
    A.Ad[8] += alpha.valued;
  }

  ADMat3x3& A;
  ADScalar& alpha;
};

/*
  Matrix determinant operations
*/
class Symm3x3Det {
public:
  Symm3x3Det( Symm3x3& A, Scalar& alpha ){
    alpha.value = Symm3x3DetCore(A.A);
  }
};

class ADSymm3x3Det {
public:
  ADSymm3x3Det( ADSymm3x3& A, ADScalar& alpha ) : A(A), alpha(alpha) {
    alpha.value = Symm3x3DetCore(A.A);
  }
  void forward(){
    alpha.valued = Symm3x3DetDerivForwardCore(A.A, A.Ad);
  }
  void reverse(){
    Symm3x3DetDerivReverseCore(alpha.valued, A.A, A.Ad);
  }
  ADSymm3x3& A;
  ADScalar& alpha;
};

class Mat3x3Det {
public:
  Mat3x3Det( Mat3x3& A, Scalar& alpha ){
    alpha.value = Mat3x3DetCore(A.A);
  }
};

class ADMat3x3Det {
public:
  ADMat3x3Det( ADMat3x3& A, ADScalar& alpha ) : A(A), alpha(alpha) {
    alpha.value = Mat3x3DetCore(A.A);
  }
  void forward(){
    alpha.valued = Mat3x3DetDerivForwardCore(A.A, A.Ad);
  }
  void reverse(){
    Mat3x3DetDerivReverseCore(alpha.valued, A.A, A.Ad);
  }

  ADMat3x3& A;
  ADScalar& alpha;
};

class Mat3x3Inverse {
public:
  Mat3x3Inverse( Mat3x3& A, Mat3x3& B ){
    Mat3x3InverseCore(A.A, B.A);
  }
};

class ADMat3x3Inverse {
public:
  ADMat3x3Inverse( ADMat3x3& A, ADMat3x3 &B ) : A(A), B(B) {
    Mat3x3InverseCore(A.A, B.A);
  }
  void forward(){
    Mat3x3InverseDerivForwardCore(B.A, A.Ad, B.Ad);
  }
  void reverse(){
    Mat3x3InverseDerivReverseCore(B.A, B.Ad, A.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
};

#endif // A2D_MAT_OPS_H
