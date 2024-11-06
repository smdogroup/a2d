#include <gtest/gtest.h>

#include "a2ddefs.h"
#include "ad/a2dgemm.h"
#include "ad/a2dmat.h"
#include "test_commons.h"

using namespace A2D;

template <class Atype, class Btype, class Ctype>
void test_gemm() {
  Atype A;
  Btype B;
  Ctype C, C_expect;

  static_assert(
      std::is_same_v<typename get_object_numeric_type<Atype>::type, double>,
      "only support double type");
  static_assert(
      std::is_same_v<typename get_object_numeric_type<Btype>::type, double>,
      "only support double type");
  static_assert(
      std::is_same_v<typename get_object_numeric_type<Ctype>::type, double>,
      "only support double type");
  static_assert(A.ncols == B.nrows, "ncols of A not equal to nrows of B");
  static_assert(A.nrows == C.nrows, "nrows of A and C not equal");
  static_assert(B.ncols == C.ncols, "ncols of A and C not equal");

  for (int i = 0; i < A.ncomp; i++) {
    A[i] = static_cast<T>(rand()) / RAND_MAX;
  }
  for (int i = 0; i < B.ncomp; i++) {
    B[i] = static_cast<T>(rand()) / RAND_MAX;
  }

  for (int i = 0; i < A.nrows; i++) {
    for (int j = 0; j < A.ncols; j++) {
      for (int k = 0; k < B.ncols; k++) {
        C_expect(i, k) += A(i, j) * B(j, k);
      }
    }
  }

  MatMatMult(A, B, C);

  for (int i = 0; i < A.nrows; i++) {
    for (int j = 0; j < A.ncols; j++) {
      EXPECT_DOUBLE_EQ(C_expect(i, j), C(i, j));
    }
  }
}

TEST(test_a2dgemm, MatMatMult) {
  using T = double;
  using A = Mat<T, 3, 4>;
  using B = Mat<T, 4, 5>;
  using C = Mat<T, 3, 5>;
  test_gemm<A, B, C>();
}

TEST(test_a2dgemm, MatMatMult3x3) {
  using T = double;
  using A = Mat<T, 3, 3>;
  using B = Mat<T, 3, 3>;
  using C = Mat<T, 3, 3>;
  test_gemm<A, B, C>();
}

TEST(test_a2dgemm, SMatSMatMult) {
  using T = double;
  using A = SymMat<T, 4>;
  using B = SymMat<T, 4>;
  using C = Mat<T, 4, 4>;
  test_gemm<A, B, C>();
}

TEST(test_a2dgemm, SMatSMatMult2x2) {
  using T = double;
  using A = SymMat<T, 2>;
  using B = SymMat<T, 2>;
  using C = Mat<T, 2, 2>;
  test_gemm<A, B, C>();
}

TEST(test_a2dgemm, SMatSMatMult3x3) {
  using T = double;
  using A = SymMat<T, 3>;
  using B = SymMat<T, 3>;
  using C = Mat<T, 3, 3>;
  test_gemm<A, B, C>();
}
