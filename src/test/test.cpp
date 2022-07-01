#include "a2d.h"
#include <stdlib.h>
#include <complex>

using namespace A2D;

/*
  Use the cplx type for TacsComplex
*/
typedef std::complex<double> TacsComplex;
typedef double TacsReal;

/*
  Define the basic scalar type TacsScalar
*/
typedef TacsComplex TacsScalar;

// Define the real part function for the complex data type
inline double TacsRealPart( const std::complex<double>& c ){
  return real(c);
}

// Define the imaginary part function for the complex data type
inline double TacsImagPart( const std::complex<double>& c ){
  return imag(c);
}

// Dummy function for real part
inline double TacsRealPart( const double& r ){
  return r;
}

void generate_random_array( int size, TacsScalar array[],
                            double lower=-1.0, double upper=1.0 ){
  for ( int i = 0; i < size; i++ ){
    array[i] = (upper - lower)*(rand()/((double)RAND_MAX + 1.0)) + lower;
  }
}

TacsScalar test_vector_scalar( const TacsScalar beta,
                               const TacsScalar avals[],
                               const TacsScalar bvals[] ){
  Scalar<TacsScalar> beta_scalar(beta);
  Vec3<TacsScalar> a(avals), b(bvals);
  Scalar<TacsScalar> output;

  Scalar<TacsScalar> abdot;
  Vec3<TacsScalar> c, d, e;
  Vec3Dot<TacsScalar> dot(a, b, abdot);
  Vec3Scale<TacsScalar> scale(abdot, b, c);
  Vec3CrossProduct<TacsScalar> cross(c, a, d);
  Vec3Axpy<TacsScalar> axpy(beta_scalar, d, a, e);
  Vec3Norm<TacsScalar> norm(e, output);

  return output.value;
}

TacsScalar test_vector_forward( const TacsScalar beta,
                                const TacsScalar dbeta,
                                const TacsScalar avals[],
                                const TacsScalar davals[],
                                const TacsScalar bvals[],
                                const TacsScalar dbvals[] ){
  ADScalar<TacsScalar> beta_scalar(beta, dbeta);
  ADVec3<TacsScalar> a(avals, davals), b(bvals, dbvals);
  ADScalar<TacsScalar> output;

  ADScalar<TacsScalar> abdot;
  ADVec3<TacsScalar> c, d, e;
  ADVec3Dot<TacsScalar> dot(a, b, abdot);
  dot.forward();
  ADVec3Scale<TacsScalar> scale(abdot, b, c);
  scale.forward();
  ADVec3CrossProduct<TacsScalar> cross(c, a, d);
  cross.forward();
  ADVec3Axpy<TacsScalar> axpy(beta_scalar, d, a, e);
  axpy.forward();
  ADVec3Norm<TacsScalar> norm(e, output);
  norm.forward();

  return output.valued;
}

TacsScalar test_vector_reverse( const TacsScalar beta,
                                const TacsScalar dbeta,
                                const TacsScalar avals[],
                                const TacsScalar davals[],
                                const TacsScalar bvals[],
                                const TacsScalar dbvals[] ){
  ADScalar<TacsScalar> beta_scalar(beta);
  ADVec3<TacsScalar> a(avals), b(bvals);
  ADScalar<TacsScalar> output;

  ADScalar<TacsScalar> abdot;
  ADVec3<TacsScalar> c, d, e;
  ADVec3Dot<TacsScalar> dot(a, b, abdot);
  ADVec3Scale<TacsScalar> scale(abdot, b, c);
  ADVec3CrossProduct<TacsScalar> cross(c, a, d);
  ADVec3Axpy<TacsScalar> axpy(beta_scalar, d, a, e);
  ADVec3Norm<TacsScalar> norm(e, output);

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

  TacsReal dh = 1e-30;
  beta += TacsScalar(0.0, dh) * dbeta;
  for ( int i = 0; i < 3; i++ ){
    avals[i] += TacsScalar(0.0, dh) * davals[i];
    bvals[i] += TacsScalar(0.0, dh) * dbvals[i];
  }

  TacsScalar f1 = test_vector_scalar(beta, avals, bvals);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("Vector test\n");
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
}

TacsScalar test_mat_scalar( const TacsScalar avals[] ){
  Mat3x3<TacsScalar> A(avals);
  Scalar<TacsScalar> output;

  Mat3x3<TacsScalar> B;
  Mat3x3Inverse<TacsScalar> inv(A, B);
  Mat3x3Det<TacsScalar> det(B, output);

  return output.value;
}

TacsScalar test_mat_forward( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADMat3x3<TacsScalar> A(avals, davals);
  ADScalar<TacsScalar> output;

  ADMat3x3<TacsScalar> B;
  ADMat3x3Inverse<TacsScalar> inv(A, B);
  inv.forward();
  ADMat3x3Det<TacsScalar> det(B, output);
  det.forward();

  return output.valued;
}

TacsScalar test_mat_reverse( const TacsScalar avals[],
                             const TacsScalar davals[] ){
  ADMat3x3<TacsScalar> A(avals);
  ADScalar<TacsScalar> output;

  ADMat3x3<TacsScalar> B;
  ADMat3x3Inverse<TacsScalar> inv(A, B);
  ADMat3x3Det<TacsScalar> det(B, output);

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

  TacsReal dh = 1e-30;
  for ( int i = 0; i < 9; i++ ){
    avals[i] += TacsScalar(0.0, dh) * davals[i];
  }

  TacsScalar f1 = test_mat_scalar(avals);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("3x3 matrix test\n");
  printf("Function value:    %25.16e\n", TacsRealPart(f0));
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
}

TacsScalar test_symm_scalar( const TacsScalar avals[] ){
  Symm3x3<TacsScalar> A(avals);
  Scalar<TacsScalar> output;
  Symm3x3Det<TacsScalar> det(A, output);

  return output.value;
}

TacsScalar test_symm_forward( const TacsScalar avals[],
                              const TacsScalar davals[] ){
  ADSymm3x3<TacsScalar> A(avals, davals);
  ADScalar<TacsScalar> output;
  ADSymm3x3Det<TacsScalar> det(A, output);
  det.forward();

  return output.valued;
}

TacsScalar test_symm_reverse( const TacsScalar avals[],
                              const TacsScalar davals[] ){
  ADSymm3x3<TacsScalar> A(avals);
  ADScalar<TacsScalar> output;
  ADSymm3x3Det<TacsScalar> det(A, output);
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

  TacsReal dh = 1e-30;
  for ( int i = 0; i < 6; i++ ){
    avals[i] += TacsScalar(0.0, dh) * davals[i];
  }

  TacsScalar f1 = test_symm_scalar(avals);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("Symmetric 3x3 matrix test\n");
  printf("Function value:    %25.16e\n", TacsRealPart(f0));
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
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
  Mat3x3<TacsScalar> Ux(avals);
  Scalar<TacsScalar> output;

  Scalar<TacsScalar> mu(0.5), lambda(0.75);
  Symm3x3<TacsScalar> E, S;
  Mat3x3GreenStrain<TacsScalar> strain(Ux, E);
  Symm3x3IsotropicConstitutive<TacsScalar> stress(mu, lambda, E, S);
  Symm3x3SymmMultTrace<TacsScalar> trace(E, S, output);

  return output.value;
}

TacsScalar test_elasticity_forward( const TacsScalar avals[],
                                    const TacsScalar davals[] ){
  ADMat3x3<TacsScalar> Ux(avals, davals);
  ADScalar<TacsScalar> output;

  Scalar<TacsScalar> mu(0.5), lambda(0.75);
  ADSymm3x3<TacsScalar> E, S;
  ADMat3x3GreenStrain<TacsScalar> strain(Ux, E);
  strain.forward();
  ADSymm3x3IsotropicConstitutive<TacsScalar> stress(mu, lambda, E, S);
  stress.forward();
  ADSymm3x3ADSymmMultTrace<TacsScalar> trace(E, S, output);
  trace.forward();

  return output.valued;
}

TacsScalar test_elasticity_reverse( const TacsScalar avals[],
                                    const TacsScalar davals[] ){
  ADMat3x3<TacsScalar> Ux(avals);
  ADScalar<TacsScalar> output;

  Scalar<TacsScalar> mu(0.5), lambda(0.75);
  ADSymm3x3<TacsScalar> E, S;
  ADMat3x3GreenStrain<TacsScalar> strain(Ux, E);
  ADSymm3x3IsotropicConstitutive<TacsScalar> stress(mu, lambda, E, S);
  ADSymm3x3ADSymmMultTrace<TacsScalar> trace(E, S, output);

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

  TacsReal dh = 1e-30;
  for ( int i = 0; i < 9; i++ ){
    avals[i] += TacsScalar(0.0, dh) * davals[i];
  }

  TacsScalar f1 = test_elasticity_scalar(avals);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("Elasticity test\n");
  printf("Function value:    %25.16e\n", TacsRealPart(f0));
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
}

TacsScalar test_shell_scalar( const TacsScalar n0_data[],
                              const TacsScalar n0xi_data[],
                              const TacsScalar X0xi_data[],
                              const TacsScalar u0xi_data[],
                              const TacsScalar d0_data[],
                              const TacsScalar d0xi_data[] ){
  Vec3<TacsScalar> n0(n0_data); // Interpolated normal at this point
  Mat3x2<TacsScalar> n0xi(n0xi_data); // Derivative of the normal along the xi and eta directions
  Mat3x2<TacsScalar> X0xi(X0xi_data); // Derivative of the mid-surface coordinate

  // Displacements and directors
  Mat3x2<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d0(d0_data);
  Mat3x2<TacsScalar> d0xi(d0xi_data);

  // Get the in-plane tangent directions
  Vec3<TacsScalar> t1_dir, t2_dir;
  Mat3x2ToVec3<TacsScalar> X0tovecs(X0xi, t1_dir, t2_dir);

  // Compute the normal direction
  Vec3<TacsScalar> n1_dir, n1;
  Vec3CrossProduct<TacsScalar> cross_normal(t1_dir, t2_dir, n1_dir);
  Vec3Normalize<TacsScalar> normalize_normal(n1_dir, n1);

  // Normalize the t1 direction
  Vec3<TacsScalar> t1;
  Vec3Normalize<TacsScalar> normalize_t1(t1_dir, t1);

  // Find the t2 direction
  Vec3<TacsScalar> t2;
  Vec3CrossProduct<TacsScalar> cross_t2(t1, n1, t2);

  // Form the transformation matrix
  Mat3x3<TacsScalar> T;
  Mat3x3FromThreeVec3<TacsScalar> assembleT(t1, t2, n1, T);

  // Assemble the matrix Xd = [X0,xi | n0]
  Mat3x3<TacsScalar> Xd, Xdz;
  Mat3x3FromMat3x2AndVec3<TacsScalar> assembleXd(X0xi, n0, Xd);
  Mat3x3FromMat3x2<TacsScalar> assembleXdz(n0xi, Xdz);

  // Assemble the zero-th and first order terms
  Mat3x3<TacsScalar> u0d, u1d;
  Mat3x3FromMat3x2AndVec3<TacsScalar> assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2<TacsScalar> assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  Mat3x3<TacsScalar> Xdinv;
  Scalar<TacsScalar> detXd;
  Mat3x3Det<TacsScalar> computedetXd(Xd, detXd);
  Mat3x3Inverse<TacsScalar> invXd(Xd, Xdinv);

  // Compute XdinvT = Xdinv * T
  Mat3x3<TacsScalar> XdinvT;
  Mat3x3MatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  Mat3x3<TacsScalar> XdinvzT, XdzXdinvT;
  Mat3x3MatMult<TacsScalar> multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  Mat3x3MatMult<TacsScalar> multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  Mat3x3<TacsScalar> u1dXdinvT, u1x;
  Mat3x3MatMult<TacsScalar> multu1d(u1d, XdinvT, u1dXdinvT);
  Mat3x3MatMultAdd<TacsScalar> multu1dadd(u0d, XdinvzT, u1dXdinvT);
  MatTrans3x3MatMult<TacsScalar> multu1x(T, u1dXdinvT, u1x);

  // Compute u0x = T^{T} * u0d * XdinvT
  Mat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3MatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  MatTrans3x3MatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);

  // Compute the Green strain
  Symm3x3<TacsScalar> e0x, e1x;
  Mat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);
  Mat3x3GreenStrain<TacsScalar> strain1(u1x, e1x);

  // This isn't really the energy, but is okay for testing..
  Scalar<TacsScalar> energy;
  Symm3x3SymmMultTraceScale<TacsScalar> trace(detXd, e0x, e1x, energy);

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
  ADVec3<TacsScalar> n0(n0_data, n0_ddata);
  ADMat3x2<TacsScalar> n0xi(n0xi_data, n0xi_ddata);
  ADMat3x2<TacsScalar> X0xi(X0xi_data, X0xi_ddata);

  // Displacements and directors
  Mat3x2<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d0(d0_data);
  Mat3x2<TacsScalar> d0xi(d0xi_data);

  // Get the in-plane tangent directions
  ADVec3<TacsScalar> t1_dir, t2_dir;
  ADMat3x2ToADVec3<TacsScalar> X0tovecs(X0xi, t1_dir, t2_dir);
  X0tovecs.forward();

  // Compute the normal direction
  ADVec3<TacsScalar> n1_dir, n1;
  ADVec3CrossProduct<TacsScalar> cross_normal(t1_dir, t2_dir, n1_dir);
  cross_normal.forward();
  ADVec3Normalize<TacsScalar> normalize_normal(n1_dir, n1);
  normalize_normal.forward();

  // Normalize the t1 direction
  ADVec3<TacsScalar> t1;
  ADVec3Normalize<TacsScalar> normalize_t1(t1_dir, t1);
  normalize_t1.forward();

  // Find the t2 direction
  ADVec3<TacsScalar> t2;
  ADVec3CrossProduct<TacsScalar> cross_t2(t1, n1, t2);
  cross_t2.forward();

  // Form the transformation matrix
  ADMat3x3<TacsScalar> T;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleT(t1, t2, n1, T);
  assembleT.forward();

  // Assemble the matrix Xd = [X0,xi | n0]
  ADMat3x3<TacsScalar> Xd, Xdz;
  ADMat3x3FromADMat3x2AndADVec3<TacsScalar> assembleXd(X0xi, n0, Xd);
  assembleXd.forward();
  ADMat3x3FromADMat3x2<TacsScalar> assembleXdz(n0xi, Xdz);
  assembleXdz.forward();

  // Assemble the zero-th and first order terms
  Mat3x3<TacsScalar> u0d, u1d;
  Mat3x3FromMat3x2AndVec3<TacsScalar> assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2<TacsScalar> assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  ADMat3x3<TacsScalar> Xdinv;
  ADScalar<TacsScalar> detXd;
  ADMat3x3Det<TacsScalar> computedetXd(Xd, detXd);
  computedetXd.forward();
  ADMat3x3Inverse <TacsScalar>invXd(Xd, Xdinv);
  invXd.forward();

  // Compute XdinvT = Xdinv * T
  ADMat3x3<TacsScalar> XdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);
  multXinvT.forward();

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  ADMat3x3<TacsScalar> XdinvzT, XdzXdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  multXdzXdinvT.forward();
  ADMat3x3ADMatMult<TacsScalar> multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);
  multXdinvzT.forward();

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  ADMat3x3<TacsScalar> u1dXdinvT, u1x;
  Mat3x3ADMatMult<TacsScalar> multu1d(u1d, XdinvT, u1dXdinvT);
  multu1d.forward();
  Mat3x3ADMatMultAdd<TacsScalar> multu1dadd(u0d, XdinvzT, u1dXdinvT);
  multu1dadd.forward();
  ADMatTrans3x3ADMatMult<TacsScalar> multu1x(T, u1dXdinvT, u1x);
  multu1x.forward();

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3ADMatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  multu0d.forward();
  ADMatTrans3x3ADMatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);
  multu0x.forward();

  // Compute the Green strain
  ADSymm3x3<TacsScalar> e0x, e1x;
  ADMat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);
  strain0.forward();
  ADMat3x3GreenStrain<TacsScalar> strain1(u1x, e1x);
  strain1.forward();

  // This isn't really the energy, but is okay for testing..
  ADScalar<TacsScalar> energy;
  ADSymm3x3ADSymmMultTraceADScale<TacsScalar> trace(detXd, e0x, e1x, energy);
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
  ADVec3<TacsScalar> n0(n0_data);
  ADMat3x2<TacsScalar> n0xi(n0xi_data);
  ADMat3x2<TacsScalar> X0xi(X0xi_data);

  // Displacements and directors
  Mat3x2<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d0(d0_data);
  Mat3x2<TacsScalar> d0xi(d0xi_data);

  // Get the in-plane tangent directions
  ADVec3<TacsScalar> t1_dir, t2_dir;
  ADMat3x2ToADVec3<TacsScalar> X0tovecs(X0xi, t1_dir, t2_dir);

  // Compute the normal direction
  ADVec3<TacsScalar> n1_dir, n1;
  ADVec3CrossProduct<TacsScalar> cross_normal(t1_dir, t2_dir, n1_dir);
  ADVec3Normalize<TacsScalar> normalize_normal(n1_dir, n1);

  // Normalize the t1 direction
  ADVec3<TacsScalar> t1;
  ADVec3Normalize<TacsScalar> normalize_t1(t1_dir, t1);

  // Find the t2 direction
  ADVec3<TacsScalar> t2;
  ADVec3CrossProduct<TacsScalar> cross_t2(t1, n1, t2);

  // Form the transformation matrix
  ADMat3x3<TacsScalar> T;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleT(t1, t2, n1, T);

  // Assemble the matrix Xd = [X0,xi | n0]
  ADMat3x3<TacsScalar> Xd, Xdz;
  ADMat3x3FromADMat3x2AndADVec3<TacsScalar> assembleXd(X0xi, n0, Xd);
  ADMat3x3FromADMat3x2<TacsScalar> assembleXdz(n0xi, Xdz);

  // Assemble the zero-th and first order terms
  Mat3x3<TacsScalar> u0d, u1d;
  Mat3x3FromMat3x2AndVec3<TacsScalar> assembleU0xi(u0xi, d0, u0d);
  Mat3x3FromMat3x2<TacsScalar> assembleU1xi(d0xi, u1d);

  // Compute the inverse matrix
  ADMat3x3<TacsScalar> Xdinv;
  ADScalar<TacsScalar> detXd;
  ADMat3x3Det<TacsScalar> computedetXd(Xd, detXd);
  ADMat3x3Inverse<TacsScalar> invXd(Xd, Xdinv);

  // Compute XdinvT = Xdinv * T
  ADMat3x3<TacsScalar> XdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);

  // Compute XdinvzT = - Xdinv * Xdz * Xdinv * T
  ADMat3x3<TacsScalar> XdinvzT, XdzXdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXdzXdinvT(Xdz, XdinvT, XdzXdinvT);
  ADMat3x3ADMatMult<TacsScalar> multXdinvzT(-1.0, Xdinv, XdzXdinvT, XdinvzT);

  // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  ADMat3x3<TacsScalar> u1dXdinvT, u1x;
  Mat3x3ADMatMult<TacsScalar> multu1d(u1d, XdinvT, u1dXdinvT);
  Mat3x3ADMatMultAdd<TacsScalar> multu1dadd(u0d, XdinvzT, u1dXdinvT);
  ADMatTrans3x3ADMatMult<TacsScalar> multu1x(T, u1dXdinvT, u1x);

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3ADMatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  ADMatTrans3x3ADMatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);

  // Compute the Green strain
  ADSymm3x3<TacsScalar> e0x, e1x;
  ADMat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);
  ADMat3x3GreenStrain<TacsScalar> strain1(u1x, e1x);

  // This isn't really the energy, but is okay for testing..
  ADScalar<TacsScalar> energy;
  ADSymm3x3ADSymmMultTraceADScale<TacsScalar> trace(detXd, e0x, e1x, energy);

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

  TacsReal dh = 1e-30;
  for ( int i = 0; i < 3; i++ ){
    n0[i] += TacsScalar(0.0, dh) * n0d[i];
  }
  for ( int i = 0; i < 6; i++ ){
    n0xi[i] += TacsScalar(0.0, dh) * n0xid[i];
    X0xi[i] += TacsScalar(0.0, dh) * X0xid[i];
  }

  TacsScalar f1 = test_shell_scalar(n0, n0xi, X0xi, u0xi, d0, d0xi);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("Shell test\n");
  printf("Function value:    %25.16e\n", TacsRealPart(f0));
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
}

/*
  Build up the terms required for the beam element
*/
TacsScalar test_beam_scalar( const TacsScalar X0xi_data[],
                             const TacsScalar axis_data[],
                             const TacsScalar n1_data[],
                             const TacsScalar n1xi_data[],
                             const TacsScalar n2_data[],
                             const TacsScalar n2xi_data[],
                             const TacsScalar u0xi_data[],
                             const TacsScalar d1_data[],
                             const TacsScalar d1xi_data[],
                             const TacsScalar d2_data[],
                             const TacsScalar d2xi_data[] ){
  Vec3<TacsScalar> X0xi(X0xi_data);
  Vec3<TacsScalar> axis(axis_data);
  Vec3<TacsScalar> n1(n1_data), n1xi(n1xi_data);
  Vec3<TacsScalar> n2(n2_data), n2xi(n2xi_data);

  Vec3<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d1(d1_data), d1xi(d1xi_data);
  Vec3<TacsScalar> d2(d2_data), d2xi(d2xi_data);

  // Compute the transformation to the local coordiantes.
  // Normalize the first direction.
  Vec3<TacsScalar> t1;
  Vec3Normalize<TacsScalar> normalizet1(X0xi, t1);

  // t2_dir = axis - dot(t1, axis) * t1
  Vec3<TacsScalar> t2_dir;
  Scalar<TacsScalar> dot;
  Vec3Dot<TacsScalar> dott1(t1, axis, dot);
  Vec3Axpy<TacsScalar> axpy(-1.0, dot, t1, axis, t2_dir);

  // Compute the n1 direction
  Vec3<TacsScalar> t2;
  Vec3Normalize<TacsScalar> normalizet2(t2_dir, t2);

  // Compute the n2 direction
  Vec3<TacsScalar> t3;
  Vec3CrossProduct<TacsScalar> cross(t1, t2, t3);

  // Assemble the referece frame
  Mat3x3<TacsScalar> T;
  Mat3x3FromThreeVec3<TacsScalar> assembleT(t1, t2, t3, T);

  // Compute the inverse
  Mat3x3<TacsScalar> Xd, Xdinv;
  Mat3x3FromThreeVec3<TacsScalar> assembleXd(X0xi, n1, n2, Xd);
  Mat3x3Inverse<TacsScalar> invXd(Xd, Xdinv);

  // Compute the determinant of the transform
  Scalar<TacsScalar> detXd;
  Mat3x3Det<TacsScalar> computedetXd(Xd, detXd);

  // Compute XdinvT = Xdinv * T
  Mat3x3<TacsScalar> XdinvT;
  Mat3x3MatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);

  // Assemble u0d
  Mat3x3<TacsScalar> u0d;
  Mat3x3FromThreeVec3<TacsScalar> assembleu0d(u0xi, d1, d2, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  Mat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3MatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  MatTrans3x3MatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);

  // Compute s0, sz1 and sz2
  Scalar<TacsScalar> s0, sz1, sz2;
  Vec3<TacsScalar> e1(1.0, 0.0, 0.0);
  Mat3x3VecVecInnerProduct<TacsScalar> inners0(XdinvT, e1, e1, s0);
  Mat3x3VecVecInnerProduct<TacsScalar> innersz1(Xdinv, e1, n1xi, sz1);
  Mat3x3VecVecInnerProduct<TacsScalar> innersz2(Xdinv, e1, n2xi, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  Vec3<TacsScalar> d1t, d1x;
  Vec3Axpy<TacsScalar> axpyd1t(-1.0, sz1, u0xi, d1xi, d1t);
  MatTrans3x3VecMultScale<TacsScalar> matmultd1x(s0, T, d1t, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  Vec3<TacsScalar> d2t, d2x;
  Vec3Axpy<TacsScalar> axpyd2t(-1.0, sz2, u0xi, d2xi, d2t);
  MatTrans3x3VecMultScale<TacsScalar> matmultd2x(s0, T, d2t, d2x);

  // Compute the Green strain
  Symm3x3<TacsScalar> e0x;
  Mat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);

  // This isn't really the energy, but is okay for testing..
  Scalar<TacsScalar> product;
  Vec3Dot<TacsScalar> dotd1d2(d1x, d2x, product);
  Scalar<TacsScalar> energy;
  Symm3x3SymmMultTraceScale<TacsScalar> trace(product, e0x, e0x, energy);

  return energy.value;
}

TacsScalar test_beam_forward( const TacsScalar X0xi_data[],
                              const TacsScalar X0xi_ddata[],
                              const TacsScalar axis_data[],
                              const TacsScalar axis_ddata[],
                              const TacsScalar n1_data[],
                              const TacsScalar n1_ddata[],
                              const TacsScalar n1xi_data[],
                              const TacsScalar n1xi_ddata[],
                              const TacsScalar n2_data[],
                              const TacsScalar n2_ddata[],
                              const TacsScalar n2xi_data[],
                              const TacsScalar n2xi_ddata[],
                              const TacsScalar u0xi_data[],
                              const TacsScalar d1_data[],
                              const TacsScalar d1xi_data[],
                              const TacsScalar d2_data[],
                              const TacsScalar d2xi_data[] ){
  ADVec3<TacsScalar> X0xi(X0xi_data, X0xi_ddata);
  ADVec3<TacsScalar> axis(axis_data, axis_ddata);
  ADVec3<TacsScalar> n1(n1_data, n1_ddata), n1xi(n1xi_data, n1xi_ddata);
  ADVec3<TacsScalar> n2(n2_data, n2_ddata), n2xi(n2xi_data, n2xi_ddata);

  Vec3<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d1(d1_data), d1xi(d1xi_data);
  Vec3<TacsScalar> d2(d2_data), d2xi(d2xi_data);

  // Compute the transformation to the local coordiantes.
  // Normalize the first direction.
  ADVec3<TacsScalar> t1;
  ADVec3Normalize<TacsScalar> normalizet1(X0xi, t1);
  normalizet1.forward();

  // t2_dir = axis - dot(t1, axis) * t1
  ADVec3<TacsScalar> t2_dir;
  ADScalar<TacsScalar> dot;
  ADVec3Dot<TacsScalar> dott1(t1, axis, dot);
  dott1.forward();
  ADVec3Axpy<TacsScalar> axpy(-1.0, dot, t1, axis, t2_dir);
  axpy.forward();

  // Compute the n1 direction
  ADVec3<TacsScalar> t2;
  ADVec3Normalize<TacsScalar> normalizet2(t2_dir, t2);
  normalizet2.forward();

  // Compute the n2 direction
  ADVec3<TacsScalar> t3;
  ADVec3CrossProduct<TacsScalar> cross(t1, t2, t3);
  cross.forward();

  // Assemble the referece frame
  ADMat3x3<TacsScalar> T;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleT(t1, t2, t3, T);
  assembleT.forward();

  // Compute the inverse
  ADMat3x3<TacsScalar> Xd, Xdinv;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleXd(X0xi, n1, n2, Xd);
  assembleXd.forward();
  ADMat3x3Inverse<TacsScalar> invXd(Xd, Xdinv);
  invXd.forward();

  // Compute the determinant of the transform
  ADScalar<TacsScalar> detXd;
  ADMat3x3Det<TacsScalar> computedetXd(Xd, detXd);
  computedetXd.forward();

  // Compute XdinvT = Xdinv * T
  ADMat3x3<TacsScalar> XdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);
  multXinvT.forward();

  // Assemble u0d
  Mat3x3<TacsScalar> u0d;
  Mat3x3FromThreeVec3<TacsScalar> assembleu0d(u0xi, d1, d2, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3ADMatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  multu0d.forward();
  ADMatTrans3x3ADMatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);
  multu0x.forward();

  // Compute s0, sz1 and sz2
  ADScalar<TacsScalar> s0, sz1, sz2;
  Vec3<TacsScalar> e1(1.0, 0.0, 0.0);
  ADMat3x3VecVecInnerProduct<TacsScalar> inners0(XdinvT, e1, e1, s0);
  inners0.forward();
  ADMat3x3VecADVecInnerProduct<TacsScalar> innersz1(Xdinv, e1, n1xi, sz1);
  innersz1.forward();
  ADMat3x3VecADVecInnerProduct<TacsScalar> innersz2(Xdinv, e1, n2xi, sz2);
  innersz2.forward();

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  ADVec3<TacsScalar> d1t, d1x;
  Vec3VecADScalarAxpy<TacsScalar> axpyd1t(-1.0, sz1, u0xi, d1xi, d1t);
  axpyd1t.forward();
  ADMatTrans3x3ADVecMultADScale<TacsScalar> matmultd1x(s0, T, d1t, d1x);
  matmultd1x.forward();

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  ADVec3<TacsScalar> d2t, d2x;
  Vec3VecADScalarAxpy<TacsScalar> axpyd2t(-1.0, sz2, u0xi, d2xi, d2t);
  axpyd2t.forward();
  ADMatTrans3x3ADVecMultADScale<TacsScalar> matmultd2x(s0, T, d2t, d2x);
  matmultd2x.forward();

  // Compute the Green strain
  ADSymm3x3<TacsScalar> e0x;
  ADMat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);
  strain0.forward();

  // This isn't really the energy, but is okay for testing..
  ADScalar<TacsScalar> product;
  ADVec3Dot<TacsScalar> dotd1d2(d1x, d2x, product);
  dotd1d2.forward();
  ADScalar<TacsScalar> energy;
  ADSymm3x3ADSymmMultTraceADScale<TacsScalar> trace(product, e0x, e0x, energy);
  trace.forward();

  return energy.valued;
}

TacsScalar test_beam_reverse( const TacsScalar X0xi_data[],
                              const TacsScalar X0xi_ddata[],
                              const TacsScalar axis_data[],
                              const TacsScalar axis_ddata[],
                              const TacsScalar n1_data[],
                              const TacsScalar n1_ddata[],
                              const TacsScalar n1xi_data[],
                              const TacsScalar n1xi_ddata[],
                              const TacsScalar n2_data[],
                              const TacsScalar n2_ddata[],
                              const TacsScalar n2xi_data[],
                              const TacsScalar n2xi_ddata[],
                              const TacsScalar u0xi_data[],
                              const TacsScalar d1_data[],
                              const TacsScalar d1xi_data[],
                              const TacsScalar d2_data[],
                              const TacsScalar d2xi_data[] ){
  ADVec3<TacsScalar> X0xi(X0xi_data);
  ADVec3<TacsScalar> axis(axis_data);
  ADVec3<TacsScalar> n1(n1_data), n1xi(n1xi_data);
  ADVec3<TacsScalar> n2(n2_data), n2xi(n2xi_data);

  Vec3<TacsScalar> u0xi(u0xi_data);
  Vec3<TacsScalar> d1(d1_data), d1xi(d1xi_data);
  Vec3<TacsScalar> d2(d2_data), d2xi(d2xi_data);


  // Compute the transformation to the local coordiantes.
  // Normalize the first direction.
  ADVec3<TacsScalar> t1;
  ADVec3Normalize<TacsScalar> normalizet1(X0xi, t1);

  // t2_dir = axis - dot(t1, axis) * t1
  ADVec3<TacsScalar> t2_dir;
  ADScalar<TacsScalar> dot;
  ADVec3Dot<TacsScalar> dott1(t1, axis, dot);
  ADVec3Axpy<TacsScalar> axpy(-1.0, dot, t1, axis, t2_dir);

  // Compute the n1 direction
  ADVec3<TacsScalar> t2;
  ADVec3Normalize<TacsScalar> normalizet2(t2_dir, t2);

  // Compute the n2 direction
  ADVec3<TacsScalar> t3;
  ADVec3CrossProduct<TacsScalar> cross(t1, t2, t3);

  // Assemble the referece frame
  ADMat3x3<TacsScalar> T;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleT(t1, t2, t3, T);

  // Compute the inverse
  ADMat3x3<TacsScalar> Xd, Xdinv;
  ADMat3x3FromThreeADVec3<TacsScalar> assembleXd(X0xi, n1, n2, Xd);
  ADMat3x3Inverse<TacsScalar> invXd(Xd, Xdinv);

  // Compute the determinant of the transform
  ADScalar<TacsScalar> detXd;
  ADMat3x3Det<TacsScalar> computedetXd(Xd, detXd);

  // Compute XdinvT = Xdinv * T
  ADMat3x3<TacsScalar> XdinvT;
  ADMat3x3ADMatMult<TacsScalar> multXinvT(Xdinv, T, XdinvT);

  // Assemble u0d
  Mat3x3<TacsScalar> u0d;
  Mat3x3FromThreeVec3<TacsScalar> assembleu0d(u0xi, d1, d2, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  ADMat3x3<TacsScalar> u0dXdinvT, u0x;
  Mat3x3ADMatMult<TacsScalar> multu0d(u0d, XdinvT, u0dXdinvT);
  ADMatTrans3x3ADMatMult<TacsScalar> multu0x(T, u0dXdinvT, u0x);

  // Compute s0, sz1 and sz2
  ADScalar<TacsScalar> s0, sz1, sz2;
  Vec3<TacsScalar> e1(1.0, 0.0, 0.0);
  ADMat3x3VecVecInnerProduct<TacsScalar> inners0(XdinvT, e1, e1, s0);
  ADMat3x3VecADVecInnerProduct<TacsScalar> innersz1(Xdinv, e1, n1xi, sz1);
  ADMat3x3VecADVecInnerProduct<TacsScalar> innersz2(Xdinv, e1, n2xi, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  ADVec3<TacsScalar> d1t, d1x;
  Vec3VecADScalarAxpy<TacsScalar> axpyd1t(-1.0, sz1, u0xi, d1xi, d1t);
  ADMatTrans3x3ADVecMultADScale<TacsScalar> matmultd1x(s0, T, d1t, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  ADVec3<TacsScalar> d2t, d2x;
  Vec3VecADScalarAxpy<TacsScalar> axpyd2t(-1.0, sz2, u0xi, d2xi, d2t);
  ADMatTrans3x3ADVecMultADScale<TacsScalar> matmultd2x(s0, T, d2t, d2x);

  // Compute the Green strain
  ADSymm3x3<TacsScalar> e0x;
  ADMat3x3GreenStrain<TacsScalar> strain0(u0x, e0x);

  // This isn't really the energy, but is okay for testing..
  ADScalar<TacsScalar> product;
  ADVec3Dot<TacsScalar> dotd1d2(d1x, d2x, product);
  ADScalar<TacsScalar> energy;
  ADSymm3x3ADSymmMultTraceADScale<TacsScalar> trace(product, e0x, e0x, energy);

  energy.valued = 1.0;

  trace.reverse();
  dotd1d2.reverse();
  strain0.reverse();
  matmultd2x.reverse();
  axpyd2t.reverse();
  matmultd1x.reverse();
  axpyd1t.reverse();
  innersz2.reverse();
  innersz1.reverse();
  inners0.reverse();
  multu0x.reverse();
  multu0d.reverse();
  multXinvT.reverse();
  computedetXd.reverse();
  invXd.reverse();
  assembleXd.reverse();
  assembleT.reverse();
  cross.reverse();
  normalizet2.reverse();
  axpy.reverse();
  dott1.reverse();
  normalizet1.reverse();

  TacsScalar deriv = 0.0;
  for ( int i = 0; i < 3; i++ ){
    deriv += X0xi.xd[i] * X0xi_ddata[i];
    deriv += axis.xd[i] * axis_ddata[i];
    deriv += n1.xd[i] * n1_ddata[i];
    deriv += n1xi.xd[i] * n1xi_ddata[i];
    deriv += n2.xd[i] * n2_ddata[i];
    deriv += n2xi.xd[i] * n2xi_ddata[i];
  }

  return deriv;
}

void test_beam(){
  TacsScalar X0xi[3], axis[3], n1[3], n1xi[3], n2[3], n2xi[3];
  TacsScalar X0xid[3], axisd[3], n1d[3], n1xid[3], n2d[3], n2xid[3];
  TacsScalar u0xi[3], d1[3], d1xi[3], d2[3], d2xi[3];

  // Fill in the data
  generate_random_array(3, X0xi);
  generate_random_array(3, axis);
  generate_random_array(3, n1);
  generate_random_array(3, n1xi);
  generate_random_array(3, n2);
  generate_random_array(3, n2xi);

  generate_random_array(3, X0xid);
  generate_random_array(3, axisd);
  generate_random_array(3, n1d);
  generate_random_array(3, n1xid);
  generate_random_array(3, n2d);
  generate_random_array(3, n2xid);

  generate_random_array(3, u0xi);
  generate_random_array(3, d1);
  generate_random_array(3, d1xi);
  generate_random_array(3, d2);
  generate_random_array(3, d2xi);

  TacsScalar f0 = test_beam_scalar(X0xi, axis, n1, n1xi, n2, n2xi,
                                   u0xi, d1, d1xi, d2, d2xi);
  TacsScalar df = test_beam_forward(X0xi, X0xid, axis, axisd,
                                    n1, n1d, n1xi, n1xid, n2, n2d, n2xi, n2xid,
                                    u0xi, d1, d1xi, d2, d2xi);
  TacsScalar dfr = test_beam_reverse(X0xi, X0xid, axis, axisd,
                                     n1, n1d, n1xi, n1xid, n2, n2d, n2xi, n2xid,
                                     u0xi, d1, d1xi, d2, d2xi);

  TacsReal dh = 1e-30;
  for ( int i = 0; i < 3; i++ ){
    X0xi[i] += TacsScalar(0.0, dh) * X0xid[i];
    axis[i] += TacsScalar(0.0, dh) * axisd[i];
    n1[i] += TacsScalar(0.0, dh) * n1d[i];
    n1xi[i] += TacsScalar(0.0, dh) * n1xid[i];
    n2[i] += TacsScalar(0.0, dh) * n2d[i];
    n2xi[i] += TacsScalar(0.0, dh) * n2xid[i];
  }

  TacsScalar f1 = test_beam_scalar(X0xi, axis, n1, n1xi, n2, n2xi,
                                   u0xi, d1, d1xi, d2, d2xi);
  TacsScalar fd = TacsImagPart(f1)/dh;
  printf("Beam test\n");
  printf("Function value:    %25.16e\n", TacsRealPart(f0));
  printf("Finite-difference: %25.16e\n", TacsRealPart(fd));
  printf("Forward mode:      %25.16e\n", TacsRealPart(df));
  printf("Reverse mode:      %25.16e\n", TacsRealPart(dfr));
}

int main( int argc, const char *argv[] ){
  test_vector();
  test_mat();
  test_symm();
  test_elasticity();
  test_shell();
  test_beam();

  return 0;
}
