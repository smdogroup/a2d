#include <complex>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"
#include "a2dtypes.h"
#include "block_numeric.h"

int main(int argc, char* argv[]) {
  typedef double T;
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

  blockInverse<T, N>(A, Ainv, ipiv);
  blockGemm<T, N, N, N>(Acopy, Ainv, A);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      std::cout << "(" << i << ", " << j << "): " << A(i, j) << std::endl;
    }
  }

  return (0);
}