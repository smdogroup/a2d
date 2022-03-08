#include "a2d.h"
#include "matops.h"

/*
  Compute the inverse of a 3x3 matrix

  input:
  A:          a 3x3 matrix in row major order

  output:
  Ainv:       the inverse of the 3x3 matrix

  returns:    the determinant of A
*/
static inline TacsScalar inv3x3( const TacsScalar A[],
                                 TacsScalar Ainv[] ){
  TacsScalar det = (A[8]*(A[0]*A[4] - A[3]*A[1]) -
                    A[7]*(A[0]*A[5] - A[3]*A[2]) +
                    A[6]*(A[1]*A[5] - A[2]*A[4]));
  TacsScalar detinv = 1.0/det;

  Ainv[0] = (A[4]*A[8] - A[5]*A[7])*detinv;
  Ainv[1] =-(A[1]*A[8] - A[2]*A[7])*detinv;
  Ainv[2] = (A[1]*A[5] - A[2]*A[4])*detinv;

  Ainv[3] =-(A[3]*A[8] - A[5]*A[6])*detinv;
  Ainv[4] = (A[0]*A[8] - A[2]*A[6])*detinv;
  Ainv[5] =-(A[0]*A[5] - A[2]*A[3])*detinv;

  Ainv[6] = (A[3]*A[7] - A[4]*A[6])*detinv;
  Ainv[7] =-(A[0]*A[7] - A[1]*A[6])*detinv;
  Ainv[8] = (A[0]*A[4] - A[1]*A[3])*detinv;

  return det;
}

/*
  Compute the sensitivity of the 3x3 inverse matrix

  input:
  Ainv:   The 3x3 inverse of the matrix
  Ainvd:  The derivative of the 3x3 inverse

  output:
  Ainvd:      derivative of the inverse of the 3x3 matrix
*/
static inline void inv3x3Sens( const TacsScalar Ainv[],
                               const TacsScalar Ainvd[],
                               TacsScalar Ad[] ){
  // d(Ainv_{kl})/d(A_{ij})
  //  = -Ainv_{kn}*delta_{ni}*delta{mj}*Ainv_{ml}
  //  = -Ainv_{ki}*Ainv_{jl}

  // Ad_{ij}
  //  = d(Ainv_{kl})/d(A_{ij})*Ainvd_{kl}
  //  = -Ainv_{ki}*Ainv_{jl}*Ainvd_{kl}

  // Ad = -Ainv^{T}*Ainvd*Ainv^{T}
  TacsScalar t[9];
  MatTrans3x3Mat(Ainv, Ainvd, t);
  Mat3x3MatTransAdd(t, Ainv, Ad);

  Ad[0] = -Ad[0];
  Ad[1] = -Ad[1];
  Ad[2] = -Ad[2];
  Ad[3] = -Ad[3];
  Ad[4] = -Ad[4];
  Ad[5] = -Ad[5];
  Ad[6] = -Ad[6];
  Ad[7] = -Ad[7];
  Ad[8] = -Ad[8];
}

/*
  Compute the determinant of a 3x3 matrix

  input:
  A:        a 3x3 matrix in row-major order

  returns:  the determinant of A
*/
static inline TacsScalar Mat3x3Det( const TacsScalar A[] ){
  return (A[8]*(A[0]*A[4] - A[3]*A[1]) -
          A[7]*(A[0]*A[5] - A[3]*A[2]) +
          A[6]*(A[1]*A[5] - A[2]*A[4]));
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void Mat3x3DetDeriv( const TacsScalar A[],
                                   TacsScalar Ad[] ){
  Ad[0] = A[8]*A[4] - A[7]*A[5];
  Ad[1] = A[6]*A[5] - A[8]*A[3];
  Ad[2] = A[7]*A[3] - A[6]*A[4];
  Ad[3] = A[7]*A[2] - A[8]*A[1];
  Ad[4] = A[8]*A[0] - A[6]*A[2];
  Ad[5] = A[6]*A[1] - A[7]*A[0];
  Ad[6] = A[1]*A[5] - A[2]*A[4];
  Ad[7] = A[3]*A[2] - A[0]*A[5];
  Ad[8] = A[0]*A[4] - A[3]*A[1];
}


/*
  Inverse matrix operation
*/
template<class Atype, class Btype>
class Mat3x3InverseOper : public ADOp {
public:
  Mat3x3InverseOper( ADMat3x3 &A, ADMat3x3 &B ) : A(A), B(B) {
    B.op = this;
  }

  void computeDeriv(){
    // given B.Ad compute A.Ad
    inv3x3Sens(B.A, B.Ad, A.Ad);

    if (A.op){
      A::computeDeriv();
    }
  }

  Atype& A;
  Btype& B;
};

template<class Atype, class Btype>
Mat3x3InverseOper<Atype, Btype> Mat3x3Inverse( ADMat3x3& A, ADMat3x3& B ){
  inv3x3(A.A, B.A);
  return Mat3x3InverseOper<Atype, Btype>(A, B);
}

/*
  Matrix-matrix product operation
*/
class Mat3x3MatMultOper : public ADOp {
public:
  Mat3x3MatMultOper( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ) : A(A), B(B), C(C) {
    C.op = this;
  }

  // C = A * B
  void computeDeriv(){
    // dfdA = dfdC * B^{T}
    mat3x3MatTransMultAdd(C.Ad, B.A, A.Ad);

    // dfdB = A^{T} * dfdC
    mat3x3TransMatMultAdd(A.A, C.Ad, B.Ad);

    if (A.op && B.op && A.op == B.op){
      B.op->computeDeriv();
    }
    else if (B.op){
      B.op->computeDeriv();
    }
    else if (A.op){
      A.op->computeDeriv();
    }
  }

  ADMat3x3& A;
  ADMat3x3& B;
  ADMat3x3& C;
};

Mat3x3MatMultOper Mat3x3MatMult( ADMat3x3& A, ADMat3x3& B, ADMat3x3& C ){
  mat3x3MatMult(A.A, B.A, C.A);

  return Mat3x3MatMultOper(A, B, C);
}

/*
  Matrix-trace operation
*/
class Mat3x3TraceOper : public ADOp {
public:
  Mat3x3TraceOper( ADMat3x3& A, ADScalar& B ): A(A), B(B) {
    B.op = this;
  }

  void computeDeriv(){
    A.Ad[0] += B.adjoint;
    A.Ad[4] += B.adjoint;
    A.Ad[8] += B.adjoint;

    if (B.op){
      A.op->computeDeriv();
    }
  }

  ADMat3x3& A;
  ADScalar& B;
};

Mat3x3TraceOper Mat3x3Trace( ADMat3x3& A, ADScalar& B ){
  B.value = A.A[0] + A.A[4] + A.A[8];

  return Mat3x3TraceOper(A, B);
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void det3x3Sens( const TacsScalar s,
                               const TacsScalar A[],
                               TacsScalar Ad[] ){
  Ad[0] += s*(A[8]*A[4] - A[7]*A[5]);
  Ad[1] += s*(A[6]*A[5] - A[8]*A[3]);
  Ad[2] += s*(A[7]*A[3] - A[6]*A[4]);

  Ad[3] += s*(A[7]*A[2] - A[8]*A[1]);
  Ad[4] += s*(A[8]*A[0] - A[6]*A[2]);
  Ad[5] += s*(A[6]*A[1] - A[7]*A[0]);

  Ad[6] += s*(A[1]*A[5] - A[2]*A[4]);
  Ad[7] += s*(A[3]*A[2] - A[0]*A[5]);
  Ad[8] += s*(A[0]*A[4] - A[3]*A[1]);
}


class Mat3x3DetOper : public ADOp {
public:
  Mat3x3DetOper( ADMat3x3& A, ADScalar& B ): A(A), B(B) {
    B.op = this;
  }

  void computeDeriv(){
    det3x3Sens(B.adjoint, A.A, A.Ad);

    if (B.op){
      B.op->computeDeriv();
    }
  }
};

Mat3x3DetOper Mat3x3Det( ADMat3x3& A, ADScalar& B ){
  B.value = det3x3(A.A);
  return Mat3x3DetOper(A, B);
}

/*
  Input:
  Uxi

  Passive:
  J, lambda, mu

  Derivatives:
  Compute df/dUxi
  Compute d2f/d2Uxi

  Ux = Uxi * J^{-1}
  Ux^{T} = J^{-T} * Uxi^{T}

  E = 0.5*(Ux + Ux^{T} + Ux^{T}*Ux)

  S = lambda * tr(E) * I + 2 * mu * E

  energy += detJ * tr(E * S)
*/
int main( int argc, const char *argv[] ){

  TacsScalar Uvals[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  TacsScalar Jvals[] = {1, 0.5, 0.2, 1, 6, 2, 3, 7, 9};

  ADMat3x3 Uxi(Uvals), J(Jvals);

  TacsScalar lambda = 5.3, mu = 0.35;

  ADMat3x3 Jinv; // Passive values
  ADMat3x3 Ux, E, S, T;
  ADScalar energy;

  ADOp op1 = Mat3x3Inverse(J, Jinv);
  ADOp op2 = Mat3x3MatMult(Uxi, Jinv, Ux);
  // ADOp op3 = Mat3x3ToGreenStrain(Ux, E);
  // ADOp op4 = Mat3x3IsoStress(lambda, mu, E, S);
  // ADOp op5 = Mat3x3MatMult(E, S, T);
  ADOp op3 = Mat3x3MatMult(Ux, Ux, T);
  ADOp op6 = Mat3x3Trace(T, energy);

  energy.adjoint = 1.0;
  energy.op->computeDeriv();

  ADMat3x3 dfdUxi(Uxi.Ad);

  TacsScalar fd[9];
  double dh = 1e-6;

  for ( int k = 0; k < 9; k++ ){
    ADMat3x3 Jinv1; // Passive values
    ADMat3x3 Ux1, E1, S1, T1;
    ADScalar energy1;

    ADMat3x3 Uxi1(Uvals), J1(Jvals);
    Uxi1.A[k] += dh;

    Mat3x3Inverse(J1, Jinv1);
    Mat3x3MatMult(Uxi1, Jinv1, Ux1);
    // ADOp op3 = Mat3x3ToGreenStrain(Ux, E);
    // ADOp op4 = Mat3x3IsoStress(lambda, mu, E, S);
    // ADOp op5 = Mat3x3MatMult(E, S, T);
    Mat3x3MatMult(Ux1, Ux1, T1);
    Mat3x3Trace(T1, energy1);

    fd[k] = (energy1.value - energy.value)/dh;
  }

  for ( int i = 0; i < 9; i++ ){
    printf("dfdUxi[%d] = %25.15e  %25.15e  %10.5e\n", i, dfdUxi.A[i], fd[i], (fd[i] - dfdUxi.A[i])/fd[i]);
  }

  return 0;
}
