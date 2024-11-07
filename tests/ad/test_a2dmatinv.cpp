#include <gtest/gtest.h>

#include "a2ddefs.h"
#include "ad/a2dmat.h"
#include "ad/a2dmatinv.h"
#include "test_commons.h"

using namespace A2D;

template <typename T, int N>
void test_sym_mat_inv() {
  SymMat<T, N> S, Sinv;
  Mat<T, N, N> A, Ainv;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
    }
  }

  MatInv(A, Ainv);
  MatInv(S, Sinv);

  for (int i = 0; i < A.nrows; i++) {
    for (int j = 0; j < A.ncols; j++) {
      EXPECT_DOUBLE_EQ(Ainv(i, j), Sinv(i, j));
    }
  }
}

TEST(test_a2dmatinv, MatInv1x1) {
  using T = double;
  constexpr int N = 1;
  test_sym_mat_inv<T, N>();
}

TEST(test_a2dmatinv, MatInv2x2) {
  using T = double;
  constexpr int N = 2;
  test_sym_mat_inv<T, N>();
}

TEST(test_a2dmatinv, MatInv3x3) {
  using T = double;
  constexpr int N = 3;
  test_sym_mat_inv<T, N>();
}
