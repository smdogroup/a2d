#include <gtest/gtest.h>

#include <complex>

#include "ad/a2dmat.h"
#include "test_commons.h"

using namespace A2D;

TEST(test_a2dmat, SymMat_CopyConstructor) {
  using T = double;
  using T2 = std::complex<T>;
  constexpr int N = 4;

  SymMat<T, N> S1;
  for (int i = 0; i < S1.ncomp; i++) {
    S1[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  SymMat<T2, N> S2(S1);
  for (int i = 0; i < S1.ncomp; i++) {
    EXPECT_DOUBLE_EQ(S1[i], S2[i].real());
    EXPECT_DOUBLE_EQ(0.0, S2[i].imag());
  }
}

TEST(test_a2dmat, SymMat_copy) {
  using T = double;
  using T2 = std::complex<T>;
  constexpr int N = 4;

  SymMat<T, N> S1;
  for (int i = 0; i < S1.ncomp; i++) {
    S1[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  SymMat<T2, N> S2;
  S2.copy(S1);
  for (int i = 0; i < S1.ncomp; i++) {
    EXPECT_DOUBLE_EQ(S1[i], S2[i].real());
    EXPECT_DOUBLE_EQ(0.0, S2[i].imag());
  }
}

TEST(test_a2dmat, SymMat_get) {
  using T = double;
  using T2 = std::complex<T>;
  constexpr int N = 4;

  SymMat<T, N> S1;
  for (int i = 0; i < S1.ncomp; i++) {
    S1[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  SymMat<T2, N> S2;
  S2.copy(S1);
  S1.get(S2);
  for (int i = 0; i < S1.ncomp; i++) {
    EXPECT_DOUBLE_EQ(S1[i], S2[i].real());
    EXPECT_DOUBLE_EQ(0.0, S2[i].imag());
  }
}

TEST(test_a2dmat, SymMat_operator) {
  using T = double;
  constexpr int N = 3;

  T data[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  T data_expect[] = {0.0, 1.0, 3.0, 1.0, 2.0, 4.0, 3.0, 4.0, 5.0};

  SymMat<T, N> S(data);
  for (int i = 0, index = 0; i < N; i++) {
    for (int j = 0; j < N; j++, index++) {
      EXPECT_DOUBLE_EQ(S(i, j), data_expect[index]);
    }
  }
}

TEST(test_a2dmat, SymMat_const_operator) {
  using T = double;
  constexpr int N = 3;

  T data[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  T data_expect[] = {0.0, 1.0, 3.0, 1.0, 2.0, 4.0, 3.0, 4.0, 5.0};

  const SymMat<T, N> S(data);
  for (int i = 0, index = 0; i < N; i++) {
    for (int j = 0; j < N; j++, index++) {
      EXPECT_DOUBLE_EQ(S(i, j), data_expect[index]);
    }
  }
}
