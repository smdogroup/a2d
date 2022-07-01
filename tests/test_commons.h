#ifndef TEST_COMMONS_H
#define TEST_COMMONS_H

#include <gtest/gtest.h>

#include <iostream>

// Global typenames
typedef double T;
typedef int I;

/**
 * @brief Expect matrix == values
 *
 * @tparam m first dimension
 * @tparam n second dimension
 * @tparam MatType matrix type
 *
 * @param mat matrix
 * @param vals array of values (in row major order)
 * @param abs_err maximum allowed entry-wise absolute error
 */
template <I m, I n, class MatType>
void expect_mat_eq(const MatType& mat, const T vals[], T abs_err = 1e-15) {
  for (I i = 0; i < m; i++) {
    for (I j = 0; j < n; j++) {
      EXPECT_NEAR(mat(i, j), vals[n * i + j], abs_err);
    }
  }
}

template <typename ScalarType>
void expect_val_eq(const ScalarType& val1, const T val2, T abs_err = 1e-15) {
  EXPECT_NEAR(val1, val2, abs_err);
}

template <I m, I n, class MatType>
void print_mat(const MatType& mat) {
  for (I i = 0; i < m; i++) {
    for (I j = 0; j < n; j++) {
      std::cout << std::setw(15) << mat(i, j);
    }
    std::cout << std::endl;
  }
}

#endif  // TEST_COMMONS_H