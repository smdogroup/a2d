/*
  This is a set of unit tests for a2dmatcore2d.h using Google Test framework.
 */

#include <gtest/gtest.h>

#include "a2dmatcore2d.h"
#include "a2dobjs.h"
#include "test_commons.h"

// Global typenames
typedef double T;
typedef int I;
typedef A2D::Mat<T, 2, 2> Mat2x2;

class MatCoreTest : public ::testing::Test {
 protected:
  // Set up square matrices and results
  const T A2x2_data[4] = {0.96155402, 0.02855176, 0.95787560, 0.45439794};

  const T B2x2_data[4] = {0.80766462, 0.60212270, 0.86418474, 0.65304149};

  const T Z2x2_data[4] = {0.85512289, 0.39849899, 0.40246850, 0.35602724};

  const T AB_data[4] = {0.8012871574649149, 0.5976189866107764,
                        1.1663259981167076, 0.8734993503266507};

  const T ABT_data[4] = {0.7938048249937244, 0.8496057946621772,
                         1.0472455469885100, 1.1245221841288746};

  const T ATB_data[4] = {1.6043946385111165, 1.2045060117768980,
                         0.4157440120261668, 0.3139323706114826};

  const T ATBT_data[4] = {1.3533718047088925, 1.4564928198282989,
                          0.2966635608979692, 0.3214147030826730};
};

/**
 * Symm2x2SymmMatMult...Core
 */

TEST_F(MatCoreTest, Symm2x2SymmMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Symm2x2SymmMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, AB_data);
}

TEST_F(MatCoreTest, Symm2x2SymmMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Symm2x2SymmMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Symm2x2SymmMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2SymmMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2SymmMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2SymmMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2SymmMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Symm2x2SymmMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Symm2x2MatMult...Core
 */

TEST_F(MatCoreTest, Symm2x2MatMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Symm2x2MatMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, AB_data);
}

TEST_F(MatCoreTest, Symm2x2MatMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Symm2x2MatMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Symm2x2MatMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2MatMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2MatMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2MatMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2MatMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Symm2x2MatMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Symm2x2MatTransMult...Core
 */

TEST_F(MatCoreTest, Symm2x2MatTransMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Symm2x2MatTransMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, ABT_data);
}

TEST_F(MatCoreTest, Symm2x2MatTransMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Symm2x2MatTransMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * ABT_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Symm2x2MatTransMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2MatTransMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2MatTransMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Symm2x2MatTransMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Symm2x2MatTransMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Symm2x2MatTransMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B,
                                                                  C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Mat2x2SymmMult...Core
 */

TEST_F(MatCoreTest, Mat2x2SymmMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Mat2x2SymmMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, AB_data);
}

TEST_F(MatCoreTest, Mat2x2SymmMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Mat2x2SymmMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Mat2x2SymmMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2SymmMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2SymmMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2SymmMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2SymmMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Mat2x2SymmMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * MatTrans2x2SymmMult...Core
 */

TEST_F(MatCoreTest, MatTrans2x2SymmMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::MatTrans2x2SymmMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, ATB_data);
}

TEST_F(MatCoreTest, MatTrans2x2SymmMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::MatTrans2x2SymmMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * ATB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, MatTrans2x2SymmMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2SymmMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2SymmMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2SymmMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2SymmMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::MatTrans2x2SymmMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B,
                                                                  C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Mat2x2MatMult...Core
 */

TEST_F(MatCoreTest, Mat2x2MatMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Mat2x2MatMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, AB_data);
}

TEST_F(MatCoreTest, Mat2x2MatMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Mat2x2MatMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Mat2x2MatMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2MatMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2MatMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2MatMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2MatMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Mat2x2MatMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Mat2x2MatTransMult...Core
 */

TEST_F(MatCoreTest, Mat2x2MatTransMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::Mat2x2MatTransMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, ABT_data);
}

TEST_F(MatCoreTest, Mat2x2MatTransMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::Mat2x2MatTransMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * ABT_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, Mat2x2MatTransMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2MatTransMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2MatTransMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::Mat2x2MatTransMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, Mat2x2MatTransMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::Mat2x2MatTransMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B,
                                                                 C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * ABT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Mat3Transx3MatMult...Core
 */

TEST_F(MatCoreTest, MatTrans2x2MatMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::MatTrans2x2MatMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, ATB_data);
}

TEST_F(MatCoreTest, MatTrans2x2MatMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::MatTrans2x2MatMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B, AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * ATB_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2MatMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2MatMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::MatTrans2x2MatMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B,
                                                                 C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * ATB_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

/**
 * Mat3Transx3MatTransMult...Core
 */

TEST_F(MatCoreTest, MatTrans2x2MatTransMultCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  A2D::MatTrans2x2MatTransMultCore<Mat2x2, Mat2x2, Mat2x2>(A, B, AB);
  EXPECT_MAT_EQ(2, 2, AB, ATBT_data);
}

TEST_F(MatCoreTest, MatTrans2x2MatTransMultScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 AB;
  T scale = 2.34;
  A2D::MatTrans2x2MatTransMultScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A, B,
                                                                   AB);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = scale * ATBT_data[i];
  EXPECT_MAT_EQ(2, 2, AB, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatTransMultAddCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2MatTransMultAddCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + ATBT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatTransMultSubCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  A2D::MatTrans2x2MatTransMultSubCore<Mat2x2, Mat2x2, Mat2x2>(A, B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] - ATBT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}

TEST_F(MatCoreTest, MatTrans2x2MatTransMultAddScaleCore) {
  const Mat2x2 A(A2x2_data);
  const Mat2x2 B(B2x2_data);
  Mat2x2 C(Z2x2_data);
  T scale = 2.34;
  A2D::MatTrans2x2MatTransMultAddScaleCore<T, Mat2x2, Mat2x2, Mat2x2>(scale, A,
                                                                      B, C);
  T res[4];
  for (I i = 0; i < 4; i++) res[i] = Z2x2_data[i] + scale * ATBT_data[i];
  EXPECT_MAT_EQ(2, 2, C, res);
}