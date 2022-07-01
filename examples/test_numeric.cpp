#include <iomanip>
#include <iostream>

#include "a2dtmp3d.h"
#include "a2dtypes.h"
#include "block_numeric.h"
#include "sparse_matrix.h"

template <typename T>
void reisdual(T mu, T lambda, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
              A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib) {
  A2D::Mat<T, 3, 3> Ux0, Uxb;
  A2D::SymmMat<T, 3> E0, Eb;

  A2D::ADMat<A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
  A2D::ADMat<A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
  A2D::ADMat<A2D::SymmMat<T, 3>> E(E0, Eb);
  A2D::ADScalar<T> output;

  auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
  auto strain = A2D::MatLinearGreenStrain(Ux, E);
  auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

  output.bvalue = wdetJ;

  energy.reverse();
  strain.reverse();
  mult.reverse();
}

template <typename T>
void adjoint_product(T mu0, T lambda0, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                     A2D::Mat<T, 3, 3> Uxi0, A2D::Mat<T, 3, 3> Pxi0, T& dmu,
                     T& dlambda) {
  A2D::Mat<T, 3, 3> Uxib, Ux0, Uxb;
  A2D::SymmMat<T, 3> E0, Eb;

  const int N = 1;
  A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
  A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
  A2D::A2DMat<N, A2D::SymmMat<T, 3>> E(E0, Eb);
  A2D::A2DScalar<N, T> output;
  A2D::A2DScalar<N, T> mu(mu0), lambda(lambda0);

  // Set the seed values
  A2D::Mat<T, 3, 3>& Psi = Uxi.pvalue(0);
  for (int k2 = 0; k2 < 3; k2++) {
    for (int k1 = 0; k1 < 3; k1++) {
      Psi(k1, k2) = Pxi0(k1, k2);
    }
  }

  auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
  auto strain = A2D::MatLinearGreenStrain(Ux, E);
  auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

  output.bvalue = wdetJ;

  energy.reverse();
  strain.reverse();
  mult.reverse();

  mult.hforward();
  strain.hforward();
  energy.hreverse();

  dmu = mu.hvalue[0];
  dlambda = lambda.hvalue[0];
}

int main(int argc, char* argv[]) {
  typedef std::complex<double> T;
  const int N = 10;

  A2D::Vec<A2D::index_t, N> ipiv;
  A2D::Mat<T, N, N> A;
  A2D::Mat<T, N, N> Acopy;
  A2D::Mat<T, N, N> Ainv;

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      Acopy(i, j) = A(i, j);
    }
  }

  A2D::blockInverse<T, N>(A, Ainv, ipiv);
  A2D::blockGemm<T, N, N, N>(Acopy, Ainv, A);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << "(" << i << ", " << j << "): " << A(i, j) << std::endl;
    }
  }

  T mu = 0.347;
  T lambda = 1.758;
  T wdetJ = 0.919;
  T dmu, dlambda;
  A2D::Mat<T, 3, 3> Jinv, Uxi, Pxi;
  A2D::Mat<T, 3, 3> Uxib;

  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      Jinv(ii, jj) = -1.0 + 2.0 * std::rand() / RAND_MAX;
      Uxi(ii, jj) = -1.0 + 2.0 * std::rand() / RAND_MAX;
      Pxi(ii, jj) = -1.0 + 2.0 * std::rand() / RAND_MAX;
    }
  }

  adjoint_product(mu, lambda, wdetJ, Jinv, Uxi, Pxi, dmu, dlambda);

  double dh = 1e-30;
  Uxib.zero();
  reisdual(mu + T(0.0, dh), lambda, wdetJ, Jinv, Uxi, Uxib);
  T fd = 0.0;
  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      fd += Uxib(ii, jj) * Pxi(ii, jj);
    }
  }
  fd = fd.imag() / dh;

  std::cout << "dmu: " << dmu.real() << " fd: " << fd.real() << std::endl;

  Uxib.zero();
  reisdual(mu, lambda + T(0.0, dh), wdetJ, Jinv, Uxi, Uxib);
  fd = 0.0;
  for (int ii = 0; ii < 3; ii++) {
    for (int jj = 0; jj < 3; jj++) {
      fd += Uxib(ii, jj) * Pxi(ii, jj);
    }
  }
  fd = fd.imag() / dh;

  std::cout << "dlambda: " << dlambda.real() << " fd: " << fd.real()
            << std::endl;

  return (0);
}