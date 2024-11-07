#include <gtest/gtest.h>

#include "a2ddefs.h"
#include "ad/a2dmat.h"
#include "ad/a2dmatdet.h"
#include "test_commons.h"

using namespace A2D;

template <int N>
void test_sym_mat_det() {
  using T = double;
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

  T detA, detS;
  MatDet(A, detA);
  MatDet(S, detS);
  EXPECT_DOUBLE_EQ(detA, detS);
}

TEST(test_a2dmatinv, MatDet) {
  test_sym_mat_det<1>();
  test_sym_mat_det<2>();
  test_sym_mat_det<3>();
}
