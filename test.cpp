#include "a2d.h"

TacsScalar test_vector_scalar( const TacsScalar beta,
                               const TacsScalar avals[],
                               const TacsScalar bvals[] ){
  Scalar beta_scalar(beta);
  Vec3 a(avals), b(bvals);
  Scalar output;

  Scalar abdot;
  Vec3 c, d, e;
  Vec3Dot dot(a, b, abdot);
  Vec3Scale scale(abdot, b, c);
  Vec3CrossProduct cross(c, a, d);
  Vec3Axpy axpy(beta_scalar, d, a, e);
  Vec3Norm norm(e, output);

  return output.value;
}

TacsScalar test_vector_forward( const TacsScalar beta,
                                const TacsScalar dbeta,
                                const TacsScalar avals[],
                                const TacsScalar davals[],
                                const TacsScalar bvals[],
                                const TacsScalar dbvals[] ){
  ADScalar beta_scalar(beta, dbeta);
  ADVec3 a(avals, davals), b(bvals, dbvals);
  ADScalar output;

  ADScalar abdot;
  ADVec3 c, d, e;
  ADVec3Dot dot(a, b, abdot);
  dot.forward();
  ADVec3Scale scale(abdot, b, c);
  scale.forward();
  ADVec3CrossProduct cross(c, a, d);
  cross.forward();
  ADVec3Axpy axpy(beta_scalar, d, a, e);
  axpy.forward();
  ADVec3Norm norm(e, output);
  norm.forward();

  return output.valued;
}

TacsScalar test_vector_reverse( const TacsScalar beta,
                                const TacsScalar dbeta,
                                const TacsScalar avals[],
                                const TacsScalar davals[],
                                const TacsScalar bvals[],
                                const TacsScalar dbvals[] ){
  ADScalar beta_scalar(beta);
  ADVec3 a(avals), b(bvals);
  ADScalar output;

  ADScalar abdot;
  ADVec3 c, d, e;
  ADVec3Dot dot(a, b, abdot);
  ADVec3Scale scale(abdot, b, c);
  ADVec3CrossProduct cross(c, a, d);
  ADVec3Axpy axpy(beta_scalar, d, a, e);
  ADVec3Norm norm(e, output);

  output.valued = 1.0;
  norm.reverse();
  axpy.reverse();
  cross.reverse();
  scale.reverse();
  dot.reverse();

  return (beta_scalar.valued * dbeta +
          a.xd[0] * davals[0] + a.xd[1] * davals[1] + a.xd[2] * davals[2] +
          b.xd[0] * dbvals[0] + b.xd[1] * dbvals[1] + b.xd[2] * dbvals[2]);
}

void test_vector(){
  TacsScalar beta = 0.794;
  TacsScalar avals[] = {1.0, -0.25, 0.333};
  TacsScalar bvals[] = {-0.2, 1.0, -0.4};

  TacsScalar dbeta = -0.145;
  TacsScalar davals[] = {0.74, 0.23, 0.89};
  TacsScalar dbvals[] = {0.19, 0.56, 0.32};

  TacsScalar f0 = test_vector_scalar(beta, avals, bvals);
  TacsScalar df = test_vector_forward(beta, dbeta,
                                      avals, davals, bvals, dbvals);
  TacsScalar dfr = test_vector_reverse(beta, dbeta,
                                       avals, davals, bvals, dbvals);

  TacsScalar dh = 1e-6;
  beta += dh * dbeta;
  for ( int i = 0; i < 3; i++ ){
    avals[i] += dh * davals[i];
    bvals[i] += dh * dbvals[i];
  }

  TacsScalar f1 = test_vector_scalar(beta, avals, bvals);
  TacsScalar fd = (f1 - f0)/dh;
  printf("Vector test\n");
  printf("Finite-difference: %15.8e\n", fd);
  printf("Forward mode:      %15.8e\n", df);
  printf("Reverse mode:      %15.8e\n", dfr);
}

TacsScalar test_mat_scalar( const TacsScalar avals[] ){
  Mat3x3 A(avals);
  Scalar output;

  Mat3x3 B;
  Mat3x3Inverse inv(A, B);
  Mat3x3Det det(B, output);

  return output.value;
}

TacsScalar test_mat_forward( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADMat3x3 A(avals, davals);
  ADScalar output;

  ADMat3x3 B;
  ADMat3x3Inverse inv(A, B);
  inv.forward();
  ADMat3x3Det det(B, output);
  det.forward();

  return output.valued;
}

TacsScalar test_mat_reverse( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADMat3x3 A(avals);
  ADScalar output;

  ADMat3x3 B;
  ADMat3x3Inverse inv(A, B);
  ADMat3x3Det det(B, output);

  output.valued = 1.0;
  det.reverse();
  inv.reverse();

  TacsScalar value = 0.0;
  for ( int i = 0; i < 9; i++ ){
    value += A.Ad[i] * davals[i];
  }
  return value;
}

void test_mat(){
  TacsScalar avals[] = {1.0, -0.25, 0.333,
                        -0.25, 0.23, 0.89,
                        0.333, 0.89, 2.34};
  TacsScalar davals[] = {1.0, 2.0, 3.0,
                         4.0, 5.0, 6.0,
                         7.0, 8.0, 9.0};

  TacsScalar f0 = test_mat_scalar(avals);
  TacsScalar df = test_mat_forward(avals, davals);
  TacsScalar dfr = test_mat_reverse(avals, davals);

  TacsScalar dh = 1e-6;
  for ( int i = 0; i < 9; i++ ){
    avals[i] += dh * davals[i];
  }

  TacsScalar f1 = test_mat_scalar(avals);
  TacsScalar fd = (f1 - f0)/dh;
  printf("3x3 matrix test\n");
  printf("Function value:    %15.8e\n", f0);
  printf("Finite-difference: %15.8e\n", fd);
  printf("Forward mode:      %15.8e\n", df);
  printf("Reverse mode:      %15.8e\n", dfr);
}

TacsScalar test_symm_scalar( const TacsScalar avals[] ){
  Symm3x3 A(avals);
  Scalar output;
  Symm3x3Det det(A, output);

  return output.value;
}

TacsScalar test_symm_forward( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADSymm3x3 A(avals, davals);
  ADScalar output;
  ADSymm3x3Det det(A, output);
  det.forward();

  return output.valued;
}

TacsScalar test_symm_reverse( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADSymm3x3 A(avals);
  ADScalar output;
  ADSymm3x3Det det(A, output);
  output.valued = 1.0;
  det.reverse();

  TacsScalar value = 0.0;
  for ( int i = 0; i < 6; i++ ){
    value += A.Ad[i] * davals[i];
  }
  return value;
}

void test_symm(){
  TacsScalar avals[] = {1.0, -0.25, 0.333,
                        0.23, 0.89,
                        2.34};
  TacsScalar davals[] = {1.0, 2.0, 3.0,
                         4.0, 5.0, 6.0};

  TacsScalar f0 = test_symm_scalar(avals);
  TacsScalar df = test_symm_forward(avals, davals);
  TacsScalar dfr = test_symm_reverse(avals, davals);

  TacsScalar dh = 1e-6;
  for ( int i = 0; i < 6; i++ ){
    avals[i] += dh * davals[i];
  }

  TacsScalar f1 = test_symm_scalar(avals);
  TacsScalar fd = (f1 - f0)/dh;
  printf("Symmetric 3x3 matrix test\n");
  printf("Function value:    %15.8e\n", f0);
  printf("Finite-difference: %15.8e\n", fd);
  printf("Forward mode:      %15.8e\n", df);
  printf("Reverse mode:      %15.8e\n", dfr);
}

int main( int argc, const char *argv[] ){
  test_vector();
  test_mat();
  test_symm();

  return (0);
}

// /*
//   Input:
//   Uxi

//   Passive:
//   J, lambda, mu

//   Derivatives:
//   Compute df/dUxi
//   Compute d2f/d2Uxi

//   Ux = Uxi * J^{-1}
//   Ux^{T} = J^{-T} * Uxi^{T}

//   E = 0.5*(Ux + Ux^{T} + Ux^{T}*Ux)

//   S = lambda * tr(E) * I + 2 * mu * E

//   energy += detJ * tr(E * S)
// */
// int main( int argc, const char *argv[] ){

//   TacsScalar Uvals[] = {1,   2,   3,
//                         4,   5,   6,
//                         7,   8,   9};
//   TacsScalar Jvals[] = {1, 0.5, 0.2,
//                         1, 6, 2,
//                         3, 7, 9};

//   ADMat3x3 Uxi(Uvals), J(Jvals);

//   TacsScalar lambda = 5.3, mu = 0.35;

//   ADMat3x3 Jinv; // Passive values
//   ADMat3x3 Ux, E, S, T;
//   ADScalar energy;

//   Mat3x3InverseOp op1(J, Jinv);
//   Mat3x3MatMultOp op2(Uxi, Jinv, Ux);
//   // ADOp op3 = Mat3x3ToGreenStrain(Ux, E);
//   // ADOp op4 = Mat3x3IsoStress(lambda, mu, E, S);
//   // ADOp op5 = Mat3x3MatMult(E, S, T);
//   Mat3x3MatMultOp op3(Ux, Ux, T);
//   Mat3x3TraceOp op4(T, energy);

//   energy.adjoint = 1.0;
//   op4.computeDeriv();
//   op3.computeDeriv();
//   op2.computeDeriv();
//   op1.computeDeriv();

//   ADMat3x3 dfdUxi(Uxi.Ad);

//   TacsScalar fd[9];
//   double dh = 1e-6;

//   for ( int k = 0; k < 9; k++ ){
//     ADMat3x3 Jinv1; // Passive values
//     ADMat3x3 Ux1, E1, S1, T1;
//     ADScalar energy1;

//     ADMat3x3 Uxi1(Uvals), J1(Jvals);
//     Uxi1.A[k] += dh;

//     Mat3x3InverseOp op1(J1, Jinv1);
//     Mat3x3MatMultOp op2(Uxi1, Jinv1, Ux1);
//     // ADOp op3 = Mat3x3ToGreenStrain(Ux, E);
//     // ADOp op4 = Mat3x3IsoStress(lambda, mu, E, S);
//     // ADOp op5 = Mat3x3MatMult(E, S, T);
//     Mat3x3MatMultOp op3(Ux1, Ux1, T1);
//     Mat3x3TraceOp op4(T1, energy1);

//     fd[k] = (energy1.value - energy.value)/dh;
//   }

//   for ( int i = 0; i < 9; i++ ){
//     printf("dfdUxi[%d] = %25.15e  %25.15e  %10.5e\n", i, dfdUxi.A[i], fd[i], (fd[i] - dfdUxi.A[i])/fd[i]);
//   }

//   return 0;
// }
