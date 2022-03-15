#include "a2d.h"
#include <stdlib.h>

void generate_random_array( int size, TacsScalar array[],
                            double lower=-1.0, double upper=1.0 ){
  for ( int i = 0; i < size; i++ ){
    array[i] = (upper - lower)*(rand()/((double)RAND_MAX + 1.0)) + lower;
  }
}

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

TacsScalar test_shell_scalar( const TacsScalar n0_data[],
                              const TacsScalar n0xi_data[],
                              const TacsScalar X0xi_data[],
                              const TacsScalar u0xi_data[],
                              const TacsScalar d0_data[],
                              const TacsScalar d0xi_data[] ){
  Vec3 n0(n0_data); // Interpolated normal at this point
  Mat3x2 n0xi(n0xi_data); // Derivative of the normal along the xi and eta directions
  Mat3x2 X0xi(X0xi_data); // Derivative of the mid-surface coordinate

  // Displacements and directors
  Mat3x2 u0xi(u0xi_data);
  Vec3 d0(d0_data);
  Mat3x2 d0xi(d0xi_data);

  // Get the in-plane tangent directions
  Vec3 t1_dir, t2_dir;
  Mat3x2ToVec3 X0tovecs(X0xi, t1_dir, t2_dir);

  // Compute the normal direction
  Vec3 n1_dir, n1;
  Vec3CrossProduct cross_normal(t1_dir, t2_dir, n1_dir);
  Vec3Normalize normalize_normal(n1_dir, n1);

  // Normalize the t1 direction
  Vec3 t1;
  Vec3Normalize normalize_t1(t1_dir, t1);

  // Find the t2 direction
  Vec3 t2;
  Vec3CrossProduct cross_t2(t1, n1, t2);

  // Form the transformation matrix
  Mat3x3 T;
  Mat3x3FromVec3 assembleT(t1, t2, n1, T);

  // Assemble the matrix Xd = [X0,xi | n0]
  Mat3x3 Xd, Xdz;
  Mat3x3FromMat3x2AndVec3 assembleXd(X0xi, n0, Xd);
  Mat3x3FromMat3x2 assembleXdz(n0xi, Xdz);

  // Assemble the zero-th and first order terms
  Mat3x3 u0d, u1d;
  Mat3x3FromMat3x2AndVec3 assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2 assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  Mat3x3 Xdinv;
  Scalar detXd;
  Mat3x3Det computedetXd(Xd, detXd);
  Mat3x3Inverse invXd(Xd, Xdinv);

  // Compute XdinvT = Xdinv * T
  Mat3x3 XdinvT;
  Mat3x3MatMult multXinvT(Xdinv, T, XdinvT);

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  Mat3x3 XdinvzT, XdzXdinvT;
  Mat3x3MatMult multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  Mat3x3MatMult multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  Mat3x3 u1dXdinvT, u1x;
  Mat3x3MatMult multu1d(u1d, XdinvT, u1dXdinvT);
  Mat3x3MatMultAdd multu1dadd(u0d, XdinvzT, u1dXdinvT);
  MatTrans3x3MatMult multu1x(T, u1dXdinvT, u1x);

  // Compute u0x = T^{T} * u0d * XdinvT
  Mat3x3 u0dXdinvT, u0x;
  Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
  MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

  // Compute the Green strain
  Symm3x3 e0x, e1x;
  Mat3x3GreenStrain strain0(u0x, e0x);
  Mat3x3GreenStrain strain1(u1x, e1x);

  // This isn't really the energy, but is okay for testing..
  Scalar energy;
  Symm3x3SymmMultTrace trace(e0x, e1x, energy);

  return energy.value;
}

TacsScalar test_shell_forward( const TacsScalar n0_data[],
                               const TacsScalar n0_ddata[],
                               const TacsScalar n0xi_data[],
                               const TacsScalar n0xi_ddata[],
                               const TacsScalar X0xi_data[],
                               const TacsScalar X0xi_ddata[],
                               const TacsScalar u0xi_data[],
                               const TacsScalar d0_data[],
                               const TacsScalar d0xi_data[] ){
  ADVec3 n0(n0_data, n0_ddata);
  ADMat3x2 n0xi(n0xi_data, n0xi_ddata);
  ADMat3x2 X0xi(X0xi_data, X0xi_ddata);

  // Displacements and directors
  Mat3x2 u0xi(u0xi_data);
  Vec3 d0(d0_data);
  Mat3x2 d0xi(d0xi_data);

  // Get the in-plane tangent directions
  ADVec3 t1_dir, t2_dir;
  ADMat3x2ToADVec3 X0tovecs(X0xi, t1_dir, t2_dir);
  X0tovecs.forward();

  // Compute the normal direction
  ADVec3 n1_dir, n1;
  ADVec3CrossProduct cross_normal(t1_dir, t2_dir, n1_dir);
  cross_normal.forward();
  ADVec3Normalize normalize_normal(n1_dir, n1);
  normalize_normal.forward();

  // Normalize the t1 direction
  ADVec3 t1;
  ADVec3Normalize normalize_t1(t1_dir, t1);
  normalize_t1.forward();

  // Find the t2 direction
  ADVec3 t2;
  ADVec3CrossProduct cross_t2(t1, n1, t2);
  cross_t2.forward();

  // Form the transformation matrix
  ADMat3x3 T;
  ADMat3x3FromADVec3 assembleT(t1, t2, n1, T);
  assembleT.forward();

  // Assemble the matrix Xd = [X0,xi | n0]
  ADMat3x3 Xd, Xdz;
  ADMat3x3FromADMat3x2AndADVec3 assembleXd(X0xi, n0, Xd);
  assembleXd.forward();
  ADMat3x3FromADMat3x2 assembleXdz(n0xi, Xdz);
  assembleXdz.forward();

  // Assemble the zero-th and first order terms
  Mat3x3 u0d, u1d;
  Mat3x3FromMat3x2AndVec3 assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2 assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  ADMat3x3 Xdinv;
  ADScalar detXd;
  ADMat3x3Det computedetXd(Xd, detXd);
  computedetXd.forward();
  ADMat3x3Inverse invXd(Xd, Xdinv);
  invXd.forward();

  // Compute XdinvT = Xdinv * T
  ADMat3x3 XdinvT;
  ADMat3x3ADMatMult multXinvT(Xdinv, T, XdinvT);
  multXinvT.forward();

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  ADMat3x3 XdinvzT, XdzXdinvT;
  ADMat3x3ADMatMult multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  multXdzXdinvT.forward();
  ADMat3x3ADMatMult multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);
  multXdinvzT.forward();

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  ADMat3x3 u1dXdinvT, u1x;
  Mat3x3ADMatMult multu1d(u1d, XdinvT, u1dXdinvT);
  multu1d.forward();
  Mat3x3ADMatMultAdd multu1dadd(u0d, XdinvzT, u1dXdinvT);
  multu1dadd.forward();
  ADMatTrans3x3ADMatMult multu1x(T, u1dXdinvT, u1x);
  multu1x.forward();

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3 u0dXdinvT, u0x;
  Mat3x3ADMatMult multu0d(u0d, XdinvT, u0dXdinvT);
  multu0d.forward();
  ADMatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);
  multu0x.forward();

  // Compute the Green strain
  ADSymm3x3 e0x, e1x;
  ADMat3x3GreenStrain strain0(u0x, e0x);
  strain0.forward();
  ADMat3x3GreenStrain strain1(u1x, e1x);
  strain1.forward();

  // This isn't really the energy, but is okay for testing..
  ADScalar energy;
  ADSymm3x3ADSymmMultTrace trace(e0x, e1x, energy);
  trace.forward();

  return energy.valued;
}

TacsScalar test_shell_reverse( const TacsScalar n0_data[],
                               const TacsScalar n0_ddata[],
                               const TacsScalar n0xi_data[],
                               const TacsScalar n0xi_ddata[],
                               const TacsScalar X0xi_data[],
                               const TacsScalar X0xi_ddata[],
                               const TacsScalar u0xi_data[],
                               const TacsScalar d0_data[],
                               const TacsScalar d0xi_data[] ){
  ADVec3 n0(n0_data);
  ADMat3x2 n0xi(n0xi_data);
  ADMat3x2 X0xi(X0xi_data);

  // Displacements and directors
  Mat3x2 u0xi(u0xi_data);
  Vec3 d0(d0_data);
  Mat3x2 d0xi(d0xi_data);

  // Get the in-plane tangent directions
  ADVec3 t1_dir, t2_dir;
  ADMat3x2ToADVec3 X0tovecs(X0xi, t1_dir, t2_dir);

  // Compute the normal direction
  ADVec3 n1_dir, n1;
  ADVec3CrossProduct cross_normal(t1_dir, t2_dir, n1_dir);
  ADVec3Normalize normalize_normal(n1_dir, n1);

  // Normalize the t1 direction
  ADVec3 t1;
  ADVec3Normalize normalize_t1(t1_dir, t1);

  // Find the t2 direction
  ADVec3 t2;
  ADVec3CrossProduct cross_t2(t1, n1, t2);

  // Form the transformation matrix
  ADMat3x3 T;
  ADMat3x3FromADVec3 assembleT(t1, t2, n1, T);

  // Assemble the matrix Xd = [X0,xi | n0]
  ADMat3x3 Xd, Xdz;
  ADMat3x3FromADMat3x2AndADVec3 assembleXd(X0xi, n0, Xd);
  ADMat3x3FromADMat3x2 assembleXdz(n0xi, Xdz);

  // Assemble the zero-th and first order terms
  Mat3x3 u0d, u1d;
  Mat3x3FromMat3x2AndVec3 assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2 assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  ADMat3x3 Xdinv;
  ADScalar detXd;
  ADMat3x3Det computedetXd(Xd, detXd);
  ADMat3x3Inverse invXd(Xd, Xdinv);

  // Compute XdinvT = Xdinv * T
  ADMat3x3 XdinvT;
  ADMat3x3ADMatMult multXinvT(Xdinv, T, XdinvT);

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  ADMat3x3 XdinvzT, XdzXdinvT;
  ADMat3x3ADMatMult multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  ADMat3x3ADMatMult multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  ADMat3x3 u1dXdinvT, u1x;
  Mat3x3ADMatMult multu1d(u1d, XdinvT, u1dXdinvT);
  Mat3x3ADMatMultAdd multu1dadd(u0d, XdinvzT, u1dXdinvT);
  ADMatTrans3x3ADMatMult multu1x(T, u1dXdinvT, u1x);

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3 u0dXdinvT, u0x;
  Mat3x3ADMatMult multu0d(u0d, XdinvT, u0dXdinvT);
  ADMatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

  // Compute the Green strain
  ADSymm3x3 e0x, e1x;
  ADMat3x3GreenStrain strain0(u0x, e0x);
  ADMat3x3GreenStrain strain1(u1x, e1x);

  // This isn't really the energy, but is okay for testing..
  ADScalar energy;
  ADSymm3x3ADSymmMultTrace trace(e0x, e1x, energy);

  // Reverse mode
  energy.valued = 1.0;

  trace.reverse();
  strain1.reverse();
  strain0.reverse();
  multu0x.reverse();
  multu0d.reverse();
  multu1x.reverse();
  multu1dadd.reverse();
  multu1d.reverse();
  multXdinvzT.reverse();
  multXdzXdinvT.reverse();
  multXinvT.reverse();
  invXd.reverse();
  computedetXd.reverse();
  assembleXdz.reverse();
  assembleXd.reverse();
  assembleT.reverse();
  cross_t2.reverse();
  normalize_t1.reverse();
  normalize_normal.reverse();
  cross_normal.reverse();
  X0tovecs.reverse();

  TacsScalar deriv = 0.0;
  for ( int i = 0; i < 3; i++ ){
    deriv += n0.xd[i] * n0_ddata[i];
  }

  for ( int i = 0; i < 6; i++ ){
    deriv += n0xi.Ad[i] * n0xi_ddata[i];
    deriv += X0xi.Ad[i] * X0xi_ddata[i];
  }

  return deriv;
}

void test_shell(){
  TacsScalar n0[3], n0d[3];
  TacsScalar n0xi[6], n0xid[6];
  TacsScalar X0xi[6], X0xid[6];
  TacsScalar u0xi[6];
  TacsScalar d0[3];
  TacsScalar d0xi[6];

  // Fill in the data
  generate_random_array(3, n0);
  generate_random_array(3, n0d);
  generate_random_array(6, n0xi);
  generate_random_array(6, n0xid);
  generate_random_array(6, X0xi);
  generate_random_array(6, X0xid);
  generate_random_array(6, u0xi);
  generate_random_array(3, d0);
  generate_random_array(6, d0xi);

  TacsScalar f0 = test_shell_scalar(n0, n0xi, X0xi, u0xi, d0, d0xi);
  TacsScalar df = test_shell_forward(n0, n0d, n0xi, n0xid,
                                     X0xi, X0xid,
                                     u0xi, d0, d0xi);
  TacsScalar dfr = test_shell_reverse(n0, n0d, n0xi, n0xid,
                                      X0xi, X0xid,
                                      u0xi, d0, d0xi);

  TacsScalar dh = 1e-7;
  for ( int i = 0; i < 3; i++ ){
    n0[i] += dh * n0d[i];
  }
  for ( int i = 0; i < 6; i++ ){
    n0xi[i] += dh * n0xid[i];
    X0xi[i] += dh * X0xid[i];
  }

  TacsScalar f1 = test_shell_scalar(n0, n0xi, X0xi, u0xi, d0, d0xi);
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
  test_shell();

  return 0;
}
