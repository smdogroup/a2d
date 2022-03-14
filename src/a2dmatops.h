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
  Assemble 3x3 matrices from 3x2 matrices and vectors
*/
class Mat3x3FromMat3x2 {
public:
  Mat3x3FromMat3x2( Mat3x2& A, Mat3x3& B ){
    B.A[0] = A.A[0];
    B.A[1] = A.A[1];
    B.A[2] = 0.0;

    B.A[3] = A.A[2];
    B.A[4] = A.A[3];
    B.A[5] = 0.0;

    B.A[6] = A.A[4];
    B.A[7] = A.A[5];
    B.A[8] = 0.0;
  }
};

class Mat3x3FromMat3x2AndVec3 {
public:
  Mat3x3FromMat3x2AndVec3( Mat3x2& A, Vec3& B, Mat3x3& C ){
    C.A[0] = A.A[0];
    C.A[1] = A.A[1];
    C.A[2] = B.x[0];

    C.A[3] = A.A[2];
    C.A[4] = A.A[3];
    C.A[5] = B.x[1];

    C.A[6] = A.A[4];
    C.A[7] = A.A[5];
    C.A[8] = B.x[2];
  }
};

class Mat3x3FromVec3 {
public:
  Mat3x3FromVec3( Vec3&x, Vec3& y, Vec3&z, Mat3x3& C ){
    C.A[0] = x.x[0];
    C.A[3] = x.x[1];
    C.A[6] = x.x[2];

    C.A[1] = y.x[0];
    C.A[4] = y.x[1];
    C.A[7] = y.x[2];

    C.A[2] = z.x[0];
    C.A[5] = z.x[1];
    C.A[8] = z.x[2];
  }
};

class ADMat3x3FromADMat3x2 {
public:
  ADMat3x3FromADMat3x2( ADMat3x2& A, ADMat3x3& B ) : A(A), B(B) {
    B.A[0] = A.A[0];
    B.A[1] = A.A[1];
    B.A[2] = 0.0;

    B.A[3] = A.A[2];
    B.A[4] = A.A[3];
    B.A[5] = 0.0;

    B.A[6] = A.A[4];
    B.A[7] = A.A[5];
    B.A[8] = 0.0;
  }
  void forward(){
    B.Ad[0] = A.Ad[0];
    B.Ad[1] = A.Ad[1];
    B.Ad[2] = 0.0;

    B.Ad[3] = A.Ad[2];
    B.Ad[4] = A.Ad[3];
    B.Ad[5] = 0.0;

    B.Ad[6] = A.Ad[4];
    B.Ad[7] = A.Ad[5];
    B.Ad[8] = 0.0;
  }
  void reverse(){
    A.Ad[0] += B.Ad[0];
    A.Ad[1] += B.Ad[1];

    A.Ad[2] += B.Ad[3];
    A.Ad[3] += B.Ad[4];

    A.Ad[4] += B.Ad[6];
    A.Ad[5] += B.Ad[7];
  }

  ADMat3x2& A;
  ADMat3x3& B;
};

class ADMat3x3FromADMat3x2AndADVec3 {
public:
  ADMat3x3FromADMat3x2AndADVec3( ADMat3x2& A, ADVec3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    C.A[0] = A.A[0];
    C.A[1] = A.A[1];
    C.A[2] = B.x[0];

    C.A[3] = A.A[2];
    C.A[4] = A.A[3];
    C.A[5] = B.x[1];

    C.A[6] = A.A[4];
    C.A[7] = A.A[5];
    C.A[8] = B.x[2];
  }
  void forward(){
    C.Ad[0] = A.Ad[0];
    C.Ad[1] = A.Ad[1];
    C.Ad[2] = B.xd[0];

    C.Ad[3] = A.Ad[2];
    C.Ad[4] = A.Ad[3];
    C.Ad[5] = B.xd[1];

    C.Ad[6] = A.Ad[4];
    C.Ad[7] = A.Ad[5];
    C.Ad[8] = B.xd[2];
  }
  void reverse(){
    A.Ad[0] += C.Ad[0];
    A.Ad[1] += C.Ad[1];

    A.Ad[2] += C.Ad[3];
    A.Ad[3] += C.Ad[4];

    A.Ad[4] += C.Ad[6];
    A.Ad[5] += C.Ad[7];

    B.xd[0] += C.Ad[2];
    B.xd[1] += C.Ad[5];
    B.xd[2] += C.Ad[8];
  }

  ADMat3x2& A;
  ADVec3& B;
  ADMat3x3& C;
};

class ADMat3x3FromADVec3 {
public:
  ADMat3x3FromADVec3( ADVec3&x, ADVec3& y, ADVec3&z, ADMat3x3& C ) :
    x(x), y(y), z(z), C(C) {
    C.A[0] = x.x[0];
    C.A[3] = x.x[1];
    C.A[6] = x.x[2];

    C.A[1] = y.x[0];
    C.A[4] = y.x[1];
    C.A[7] = y.x[2];

    C.A[2] = z.x[0];
    C.A[5] = z.x[1];
    C.A[8] = z.x[2];
  }
  void forward(){
    C.Ad[0] = x.xd[0];
    C.Ad[3] = x.xd[1];
    C.Ad[6] = x.xd[2];

    C.Ad[1] = y.xd[0];
    C.Ad[4] = y.xd[1];
    C.Ad[7] = y.xd[2];

    C.Ad[2] = z.xd[0];
    C.Ad[5] = z.xd[1];
    C.Ad[8] = z.xd[2];
  }
  void reverse(){
    x.xd[0] += C.Ad[0];
    x.xd[1] += C.Ad[3];
    x.xd[2] += C.Ad[6];

    y.xd[0] += C.Ad[1];
    y.xd[1] += C.Ad[4];
    y.xd[2] += C.Ad[7];

    z.xd[0] += C.Ad[2];
    z.xd[1] += C.Ad[5];
    z.xd[2] += C.Ad[8];
  }

  ADVec3& x;
  ADVec3& y;
  ADVec3& z;
  ADMat3x3& C;
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

class Mat3x3MatMultAdd {
public:
  Mat3x3MatMultAdd( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    Mat3x3MatMultAddCore(A.A, B.A, C.A);
  }
};

class ADMat3x3MatMultAdd {
public:
  ADMat3x3MatMultAdd( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultAddCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(C.Ad, B.A, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class Mat3x3ADMatMultAdd {
public:
  Mat3x3ADMatMultAdd( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMat3x3ADMatMultAdd {
public:
  ADMat3x3ADMatMultAdd( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatMultAddCore(A.Ad, B.A, C.Ad);
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

/*
  Matrix-matrix products C = A^{T} * B
*/
class MatTrans3x3MatMult {
public:
  MatTrans3x3MatMult( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    MatTrans3x3MatMultCore(A.A, B.A, C.A);
  }
};

class ADMatTrans3x3MatMult {
public:
  ADMatTrans3x3MatMult( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3ADMatMult {
public:
  MatTrans3x3ADMatMult( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMatTrans3x3ADMatMult {
public:
  ADMatTrans3x3ADMatMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultCore(A.Ad, B.A, C.Ad);
    MatTrans3x3MatMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
    Mat3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3MatMultAdd {
public:
  MatTrans3x3MatMultAdd( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    MatTrans3x3MatMultAddCore(A.A, B.A, C.A);
  }
};

class ADMatTrans3x3MatMultAdd {
public:
  ADMatTrans3x3MatMultAdd( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultAddCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3ADMatMultAdd {
public:
  MatTrans3x3ADMatMultAdd( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMatTrans3x3ADMatMultAdd {
public:
  ADMatTrans3x3ADMatMultAdd( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatMultAddCore(A.Ad, B.A, C.Ad);
    MatTrans3x3MatMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
    Mat3x3MatMultAddCore(A.A, C.Ad, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

/*
  Matrix-matrix products C = A * B^{T}
*/
class Mat3x3MatTransMult {
public:
  Mat3x3MatTransMult( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    Mat3x3MatTransMultCore(A.A, B.A, C.A);
  }
};

class ADMat3x3MatTransMult {
public:
  ADMat3x3MatTransMult( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatTransMultCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(C.Ad, B.A, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class Mat3x3ADMatTransMult {
public:
  Mat3x3ADMatTransMult( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultCore(A.A, B.A, C.A);
  }
    void forward(){
    Mat3x3MatTransMultCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatMultAddCore(C.Ad, A.A, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMat3x3ADMatTransMult {
public:
  ADMat3x3ADMatTransMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatTransMultCore(A.Ad, B.A, C.Ad);
    Mat3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(C.Ad, B.A, A.Ad);
    MatTrans3x3MatMultAddCore(C.Ad, A.A, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class Mat3x3MatTransMultAdd {
public:
  Mat3x3MatTransMultAdd( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    Mat3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
};

class ADMat3x3MatTransMultAdd {
public:
  ADMat3x3MatTransMultAdd( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatTransMultAddCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(C.Ad, B.A, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class Mat3x3ADMatTransMultAdd {
public:
  Mat3x3ADMatTransMultAdd( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
    void forward(){
    Mat3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatMultAddCore(C.Ad, A.A, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMat3x3ADMatTransMultAdd {
public:
  ADMat3x3ADMatTransMultAdd( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    Mat3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    Mat3x3MatTransMultAddCore(A.Ad, B.A, C.Ad);
    Mat3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    Mat3x3MatMultAddCore(C.Ad, B.A, A.Ad);
    MatTrans3x3MatMultAddCore(C.Ad, A.A, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

/*
  Matrix-matrix products C = A^{T} * B^{T}
*/
class MatTrans3x3MatTransMult {
public:
  MatTrans3x3MatTransMult( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    MatTrans3x3MatTransMultCore(A.A, B.A, C.A);
  }
};

class ADMatTrans3x3MatTransMult {
public:
  ADMatTrans3x3MatTransMult( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3ADMatTransMult {
public:
  MatTrans3x3ADMatTransMult( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(C.Ad, A.A, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMatTrans3x3ADMatTransMult {
public:
  ADMatTrans3x3ADMatTransMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultCore(A.Ad, B.A, C.Ad);
    MatTrans3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
    MatTrans3x3MatTransMultAddCore(C.Ad, A.A, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3MatTransMultAdd {
public:
  MatTrans3x3MatTransMultAdd( Mat3x3& A, Mat3x3& B, Mat3x3& C ){
    MatTrans3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
};

class ADMatTrans3x3MatTransMultAdd {
public:
  ADMatTrans3x3MatTransMultAdd( ADMat3x3& A, Mat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultAddCore(A.Ad, B.A, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
  }

  ADMat3x3& A;
  Mat3x3& B;
  ADMat3x3& C;
};

class MatTrans3x3ADMatTransMultAdd {
public:
  MatTrans3x3ADMatTransMultAdd( Mat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(C.Ad, A.A, B.Ad);
  }

  Mat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

class ADMatTrans3x3ADMatTransMultAdd {
public:
  ADMatTrans3x3ADMatTransMultAdd( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    MatTrans3x3MatTransMultAddCore(A.A, B.A, C.A);
  }
  void forward(){
    MatTrans3x3MatTransMultAddCore(A.Ad, B.A, C.Ad);
    MatTrans3x3MatTransMultAddCore(A.A, B.Ad, C.Ad);
  }
  void reverse(){
    MatTrans3x3MatTransMultAddCore(B.A, C.Ad, A.Ad);
    MatTrans3x3MatTransMultAddCore(C.Ad, A.A, B.Ad);
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

/*
  Linear and nonlinear strain from the displacement gradient
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

/*
  Isotropic constitutive relationships
*/
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
