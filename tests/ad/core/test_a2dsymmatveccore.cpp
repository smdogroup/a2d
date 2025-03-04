#include <gtest/gtest.h>

#include "ad/a2dmat.h"
#include "ad/a2dobj.h"
#include "ad/core/a2dmatveccore.h"
#include "ad/core/a2dsymmatveccore.h"
#include "test_commons.h"

using namespace A2D;

template <int N>
void test_symmat_vec_core() {
  using T = double;

  using T = double;
  Mat<T, N, N> A;
  SymMat<T, N> S;
  Vec<T, N> x, y1, y2;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    x[i] = static_cast<T>(rand()) / RAND_MAX;
    y1[i] = static_cast<T>(rand()) / RAND_MAX;
    y2[i] = y1[i];
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
    }
  }

  // y = Ax
  SymMatVecCore<T, N, false>(get_data(S), get_data(x), get_data(y1));
  MatVecCore<T, N, N, MatOp::NORMAL, false>(get_data(A), get_data(x),
                                            get_data(y2));
  EXPECT_VEC_EQ(N, y1, y2);

  // y += Ax
  SymMatVecCore<T, N, true>(get_data(S), get_data(x), get_data(y1));
  MatVecCore<T, N, N, MatOp::NORMAL, true>(get_data(A), get_data(x),
                                           get_data(y2));
  EXPECT_VEC_EQ(N, y1, y2);
}

TEST(test_a2dsymmatveccore, SymMatVecCore) {
  test_symmat_vec_core<1>();
  test_symmat_vec_core<2>();
  test_symmat_vec_core<3>();
  test_symmat_vec_core<4>();
}
