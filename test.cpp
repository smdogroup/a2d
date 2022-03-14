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
TacsScalar test_elasticity_scalar( const TacsScalar avals[] ){
  Mat3x3 Ux(avals);
  Scalar output;

  Scalar mu(0.5), lambda(0.75);
  Symm3x3 E, S;
  Mat3x3GreenStrain strain(Ux, E);
  Symm3x3IsotropicConstitutive stress(mu, lambda, E, S);
  Symm3x3SymmMultTrace trace(E, S, output);

  return output.value;
}

TacsScalar test_elasticity_forward( const TacsScalar avals[],
                                    const TacsScalar davals[] ){
  ADMat3x3 Ux(avals, davals);
  ADScalar output;

  Scalar mu(0.5), lambda(0.75);
  ADSymm3x3 E, S;
  ADMat3x3GreenStrain strain(Ux, E);
  strain.forward();
  ADSymm3x3IsotropicConstitutive stress(mu, lambda, E, S);
  stress.forward();
  ADSymm3x3ADSymmMultTrace trace(E, S, output);
  trace.forward();

  return output.valued;
}

TacsScalar test_elasticity_reverse( const TacsScalar avals[],
                                    const TacsScalar davals[] ){
  ADMat3x3 Ux(avals);
  ADScalar output;

  Scalar mu(0.5), lambda(0.75);
  ADSymm3x3 E, S;
  ADMat3x3GreenStrain strain(Ux, E);
  ADSymm3x3IsotropicConstitutive stress(mu, lambda, E, S);
  ADSymm3x3ADSymmMultTrace trace(E, S, output);

  output.valued = 1.0;
  trace.reverse();
  stress.reverse();
  strain.reverse();

  TacsScalar value = 0.0;
  for ( int i = 0; i < 9; i++ ){
    value += Ux.Ad[i] * davals[i];
  }

  return value;
}

void test_elasticity(){
  TacsScalar avals[] = {1.0, -0.25, 0.333,
                        -0.25, 0.23, 0.89,
                        0.333, 0.89, 2.34};
  TacsScalar davals[] = {1.0, 2.0, 3.0,
                         4.0, 5.0, 6.0,
                         7.0, 8.0, 9.0};

  TacsScalar f0 = test_elasticity_scalar(avals);
  TacsScalar df = test_elasticity_forward(avals, davals);
  TacsScalar dfr = test_elasticity_reverse(avals, davals);

  TacsScalar dh = 1e-6;
  for ( int i = 0; i < 9; i++ ){
    avals[i] += dh * davals[i];
  }

  TacsScalar f1 = test_elasticity_scalar(avals);
  TacsScalar fd = (f1 - f0)/dh;
  printf("Elasticity test\n");
  printf("Function value:    %15.8e\n", f0);
  printf("Finite-difference: %15.8e\n", fd);
  printf("Forward mode:      %15.8e\n", df);
  printf("Reverse mode:      %15.8e\n", dfr);
}

int main( int argc, const char *argv[] ){
  test_vector();
  test_mat();
  test_symm();
  test_elasticity();

  return 0;
}
