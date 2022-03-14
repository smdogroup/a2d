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

/*
  Matrix inverse
*/
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

/*
  Matrix-matrix products C = A * B
*/
class Mat3x3MatMult {
public:
  Mat3x3MatMult( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    Mat3x3MatMultCore(A.A, B.A, C.A);
  }
};

class ADMat3x3MatMult {
public:
  ADMat3x3MatMult( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(C.Ad, B.A, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class Mat3x3ADMatMult {
public:
  Mat3x3ADMatMult( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMat3x3ADMatMult {
public:
  ADMat3x3ADMatMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultCore(A.Ad, B.A, C.Ad);
    Mat3x3MatMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(C.Ad, B.A, A.Ad);
    MatTrans3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

// /*
//   Matrix-matrix products C = A^{T} * B
// */
// class MatTrans3x3MatMult {
// public:
//   MatTrans3x3MatMult( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
//     MatTrans3x3MatMult(A.A, B.A, C.A);
//   }
// };

// class ADMatTrans3x3MatMult {
// public:
//   ADMat3x3MatMult( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ){
//     Mat3x3MatMult(A.A, B.A, C.A);
//   }
// };

// class Mat3x3ADMatMult {
// public:
//   Mat3x3ADMatMult( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ){
//     Mat3x3MatMult(A.A, B.A, C.A);
//   }
// };

// class ADMat3x3ADMatMult {
// public:
//   ADMat3x3ADMatMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ){
//     Mat3x3MatMult(A.A, B.A, C.A);
//   }
// };

/*
  Specific implementations that are handy for elasticity
*/
class Mat3x3LinearGreenStrain {
public:
  Mat3x3LinearGreenStrain( Mat3x3& Ux, Symm3x3& E ){
    Mat3x3LinearGreenStrainCore(Ux.A, E.A);
  }
};

class ADMat3x3LinearGreenStrain {
public:
  ADMat3x3LinearGreenStrain( ADMat3x3& Ux, ADSymm3x3& E ) : Ux(Ux), E(E) {
    Mat3x3LinearGreenStrainCore(Ux.A, E.A);
  }
  void forward(){
    Mat3x3LinearGreenStrainCore(Ux.Ad, E.Ad);
  }
  void reverse(){
    Mat3x3LinearGreenStrainReverseCore(E.Ad, Ux.Ad);
  }

  ADMat3x3& Ux;
  ADSymm3x3& E;
};

class Mat3x3GreenStrain {
public:
  Mat3x3GreenStrain( Mat3x3& Ux, Symm3x3& E ){
    Mat3x3GreenStrainCore(Ux.A, E.A);
  }
};

class ADMat3x3GreenStrain {
public:
  ADMat3x3GreenStrain( ADMat3x3& Ux, ADSymm3x3& E ) : Ux(Ux), E(E) {
    Mat3x3GreenStrainCore(Ux.A, E.A);
  }
  void forward(){
    Mat3x3GreenStrainForwardCore(Ux.A, Ux.Ad, E.Ad);
  }
  void reverse(){
    Mat3x3GreenStrainReverseCore(Ux.A, E.Ad, Ux.Ad);
  }

  ADMat3x3& Ux;
  ADSymm3x3& E;
};

class Symm3x3SymmMultTrace {
public:
  Symm3x3SymmMultTrace( Symm3x3& S, Symm3x3& T, Scalar& alpha ){
    alpha.value = Symm3x3MatMultTraceCore(S.A, T.A);
  }
};

class ADSymm3x3ADSymmMultTrace {
public:
  ADSymm3x3ADSymmMultTrace( ADSymm3x3& S, ADSymm3x3& T, ADScalar& alpha ) : S(S), T(T), alpha(alpha) {
    alpha.value = Symm3x3MatMultTraceCore(S.A, T.A);
  }
  void forward(){
    alpha.valued = Symm3x3MatMultTraceCore(S.Ad, T.A);
    alpha.valued += Symm3x3MatMultTraceCore(S.A, T.Ad);
  }
  void reverse(){
    Symm3x3MatMultTraceReverseCore(alpha.valued, S.A, T.Ad);
    Symm3x3MatMultTraceReverseCore(alpha.valued, T.A, S.Ad);
  }

  ADSymm3x3& S;
  ADSymm3x3& T;
  ADScalar& alpha;
};

class Symm3x3IsotropicConstitutive {
public:
  Symm3x3IsotropicConstitutive( Scalar& mu, Scalar& lambda, Symm3x3& E, Symm3x3& S ){
    Symm3x3IsotropicConstitutiveCore(mu.value, lambda.value, E.A, S.A);
  }
};

class ADSymm3x3IsotropicConstitutive {
public:
  ADSymm3x3IsotropicConstitutive( Scalar& mu, Scalar& lambda, ADSymm3x3& E, ADSymm3x3& S ) :
    mu(mu), lambda(lambda), E(E), S(S) {
    Symm3x3IsotropicConstitutiveCore(mu.value, lambda.value, E.A, S.A);
  }
  void forward(){
    Symm3x3IsotropicConstitutiveCore(mu.value, lambda.value, E.Ad, S.Ad);
  }
  void reverse(){
    Symm3x3IsotropicConstitutiveReverseCore(mu.value, lambda.value, S.Ad, E.Ad);
  }

  Scalar& mu;
  Scalar& lambda;
  ADSymm3x3& E;
  ADSymm3x3& S;
};

#endif // A2D_MAT_OPS_H
