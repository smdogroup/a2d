#include <gtest/gtest.h>

#include "ad/a2dmat.h"
#include "ad/a2dobj.h"
#include "ad/core/a2dmatdetcore.h"
#include "test_commons.h"

using namespace A2D;

template <int N>
void test_mat_det_core() {
  using T = double;
  Mat<T, N, N> A;
  SymMat<T, N> S;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
    }
  }

  T detA = MatDetCore<T, N>(get_data(A));
  T detS = SymMatDetCore<T, N>(get_data(S));

  EXPECT_DOUBLE_EQ(detA, detS);
}

template <int N>
void test_mat_det_forward_core() {
  using T = double;
  Mat<T, N, N> A, Ad;
  SymMat<T, N> S, Sd;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
    Sd[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
      Ad(i, j) = Sd(i, j);
    }
  }

  T detdA = MatDetForwardCore<T, N>(get_data(A), get_data(Ad));
  T detdS = SymMatDetForwardCore<T, N>(get_data(S), get_data(Sd));

  EXPECT_DOUBLE_EQ(detdA, detdS);
}

template <int N>
void test_mat_det_reverse_core() {
  using T = double;
  Mat<T, N, N> A, Ab;
  SymMat<T, N> S, Sb;
  T bdet = 1.2345;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
    Sb[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
      Ab(i, j) = Sb(i, j);
    }
  }

  MatDetReverseCore<T, N>(bdet, get_data(A), get_data(Ab));
  SymMatDetReverseCore<T, N>(bdet, get_data(S), get_data(Sb));

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      EXPECT_DOUBLE_EQ(Ab(i, j), Sb(i, j));
    }
  }
}

template <int N>
void test_mat_det_hreverse_core() {
  using T = double;
  Mat<T, N, N> A, Ap, Ah;
  SymMat<T, N> S, Sp, Sh;
  T bdet = 1.2345, hdet = -5.6678;

  for (int i = 0; i < S.ncomp; i++) {
    S[i] = static_cast<T>(rand()) / RAND_MAX;
    Sp[i] = static_cast<T>(rand()) / RAND_MAX;
    Sh[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A(i, j) = S(i, j);
      Ap(i, j) = Sp(i, j);
      Ah(i, j) = Sh(i, j);
    }
  }

  MatDetHReverseCore<T, N>(bdet, hdet, get_data(A), get_data(Ap), get_data(Ah));
  SymMatDetHReverseCore<T, N>(bdet, hdet, get_data(S), get_data(Sp),
                              get_data(Sh));

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      EXPECT_DOUBLE_EQ(Ah(i, j), Sh(i, j));
    }
  }
}

TEST(test_a2dmatdetcore, MatDetCore) {
  test_mat_det_core<1>();
  test_mat_det_core<2>();
  test_mat_det_core<3>();
}

TEST(test_a2dmatdetcore, MatDetForwardCore) {
  test_mat_det_forward_core<1>();
  test_mat_det_forward_core<2>();
  test_mat_det_forward_core<3>();
}

TEST(test_a2dmatdetcore, MatDetReverseCore) {
  test_mat_det_reverse_core<1>();
  test_mat_det_reverse_core<2>();
  test_mat_det_reverse_core<3>();
}

TEST(test_a2dmatdetcore, MatDetHReverseCore) {
  test_mat_det_hreverse_core<1>();
  test_mat_det_hreverse_core<2>();
  test_mat_det_hreverse_core<3>();
}
