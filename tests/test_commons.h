#ifndef TEST_COMMONS_H
#define TEST_COMMONS_H

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <string>

// Global typenames
typedef double T;
typedef int I;

// /**
//  * @brief Expect matrix == values
//  *
//  * @tparam m first dimension
//  * @tparam n second dimension
//  * @tparam MatType matrix type
//  *
//  * @param mat matrix
//  * @param vals array of values (in row major order)
//  * @param abs_err maximum allowed entry-wise absolute error
//  */
// template <I m, I n, class MatType>
// void expect_mat_eq(const MatType& mat, const T vals[], T abs_err = 1e-15) {
//   for (I i = 0; i < m; i++) {
//     for (I j = 0; j < n; j++) {
//       EXPECT_NEAR(mat(i, j), vals[n * i + j], abs_err);
//     }
//   }
// }

// template <I m, class VecType>
// void expect_vec_eq(const VecType& vec, const T vals[], T abs_err = 1e-15) {
//   for (I i = 0; i < m; i++) {
//     EXPECT_NEAR(vec(i), vals[i], abs_err);
//   }
// }

// template <typename ScalarType>
// void expect_val_eq(const ScalarType& val1, const T val2, T abs_err = 1e-15) {
//   EXPECT_NEAR(val1, val2, abs_err);
// }

template <I n, class VecType>
void print_vec(const VecType& vec) {
  for (I j = 0; j < n; j++) {
    std::cout << std::setw(5) << vec[j];
  }
  std::cout << std::endl;
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

template <I m, I n, class MatType>
void print_mat(std::string name, const MatType& mat) {
  std::cout << name << ":\n";
  print_mat<m, n>(mat);
}

// Helper macros
#define _EXPECT_VAL_NEAR(val1, val2) EXPECT_NEAR(val1, val2, 1e-15)

#define _EXPECT_VAL_NEAR_TOL(val1, val2, abs_err) \
  EXPECT_NEAR(val1, val2, abs_err)

#define _EXPECT_VEC_NEAR(m, vec, vals)   \
  for (I i = 0; i < m; i++) {            \
    EXPECT_NEAR(vec(i), vals[i], 1e-15); \
  }

#define _EXPECT_VEC_NEAR_TOL(m, vec, vals, abs_err) \
  for (I i = 0; i < m; i++) {                       \
    EXPECT_NEAR(vec(i), vals[i], abs_err);          \
  }

#define _EXPECT_MAT_NEAR(m, n, mat, vals)             \
  for (I i = 0; i < m; i++) {                         \
    for (I j = 0; j < n; j++) {                       \
      EXPECT_NEAR(mat(i, j), vals[n * i + j], 1e-15); \
    }                                                 \
  }

#define _EXPECT_MAT_NEAR_TOL(m, n, mat, vals, abs_err)  \
  for (I i = 0; i < m; i++) {                           \
    for (I j = 0; j < n; j++) {                         \
      EXPECT_NEAR(mat(i, j), vals[n * i + j], abs_err); \
    }                                                   \
  }

#define _GET_EXPECT_VAL_MACRO(_1, _2, _3, FUNC, ...) FUNC
#define _GET_EXPECT_VEC_MACRO(_1, _2, _3, _4, FUNC, ...) FUNC
#define _GET_EXPECT_MAT_MACRO(_1, _2, _3, _4, _5, FUNC, ...) FUNC

// Usage:
// - EXPECT_VAL_NEAR(val1, val2), or
// - EXPECT_VAL_NEAR(val1, val2, abs_err)
#define EXPECT_VAL_NEAR(...)                                                 \
  _GET_EXPECT_VAL_MACRO(__VA_ARGS__, _EXPECT_VAL_NEAR_TOL, _EXPECT_VAL_NEAR) \
  (__VA_ARGS__)

// Usage:
// - EXPECT_VEC_NEAR(m, vec, vals), or
// - EXPECT_VEC_NEAR(m, vec, vals, abs_err)
#define EXPECT_VEC_NEAR(...)                                                 \
  _GET_EXPECT_VEC_MACRO(__VA_ARGS__, _EXPECT_VEC_NEAR_TOL, _EXPECT_VEC_NEAR) \
  (__VA_ARGS__)

// Usage:
// - EXPECT_MAT_NEAR(m, n, mat, vals), or
// - EXPECT_MAT_NEAR(m, n, mat, vals, abs_err)
#define EXPECT_MAT_NEAR(...)                                                 \
  _GET_EXPECT_MAT_MACRO(__VA_ARGS__, _EXPECT_MAT_NEAR_TOL, _EXPECT_MAT_NEAR) \
  (__VA_ARGS__)

// Usage:
// - EXPECT_VEC_EQ(m, vec, vals)
#define EXPECT_VEC_EQ(m, vec, vals) \
  for (I i = 0; i < m; i++) {       \
    EXPECT_EQ(vec[i], vals[i]);     \
  }

#endif  // TEST_COMMONS_H
