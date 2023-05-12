/*
  This is a set of unit tests for a2dmatcore.h using Google Test framework.
 */

#include "a2dmatcore3d.h"
#include "a2dobjs.h"
#include "test_commons.h"

typedef A2D::Mat<T, 3, 3> Mat3x3;
typedef A2D::Mat<T, 2, 3> Mat2x3;
typedef A2D::Mat<T, 3, 2> Mat3x2;
typedef A2D::Mat<T, 3, 5> Mat3x5;
typedef A2D::Mat<T, 5, 3> Mat5x3;
typedef A2D::Mat<T, 2, 5> Mat2x5;

class MatCoreTest : public ::testing::Test {
 protected:
  // Set up square matrices and results
  const T A3x3_data[9] = {0.54881350, 0.71518937, 0.60276338,
                          0.54488318, 0.42365480, 0.64589411,
                          0.43758721, 0.89177300, 0.96366276};

  const T B3x3_data[9] = {0.63263687, 0.19519233, 0.13896359,
                          0.83668126, 0.50298017, 0.20476354,
                          0.57524785, 0.83518588, 0.58878908};

  const T Z3x3_data[9] = {0.80726967, 0.41220247, 0.29267381,
                          0.63578301, 0.87809199, 0.29354763,
                          0.81275772, 0.20562505, 0.67229196};

  const T AB_data[9] = {
      1.2923235364876842, 0.9702697206623223, 0.5776102973399252,
      1.0707264194850581, 0.8588886214544922, 0.5427633782137270,
      1.5773084909754789, 1.3387953324733382, 0.8108053958867646};

  const T ABT_data[9] = {
      0.5705612975806112, 0.9423320050609680, 1.2679203452439609,
      0.5171631213239855, 0.8012390733492721, 1.0475687835547058,
      0.5848150892664311, 1.0119881514678650, 1.5639114293958993};

  const T ATB_data[9] = {
      1.0548143021889502, 0.7466562793675904, 0.4454838738193890,
      1.3199096973011699, 1.0974856606120882, 0.7112005432098703,
      1.4760827666702652, 1.2473642480645029, 0.7834118375247443};

  const T ATBT_data[9] = {
      0.5143653619754384, 0.8228493113698739, 1.0284290949071402,
      0.6590733095668260, 0.9940781029119421, 1.2903078586284185,
      0.6413179510169054, 1.0265157516138277, 1.4535740889415605};

  // Set up non-square matrices and results
  const T C2x3_data[6] = {0.23077168, 0.92204766, 0.37922783,
                          0.06597756, 0.09363517, 0.92310621};

  const T C3x2_data[6] = {0.23077168, 0.06597756, 0.92204766,
                          0.09363517, 0.37922783, 0.92310621};

  const T D3x5_data[15] = {0.68077325, 0.00794715, 0.79517780, 0.46285641,
                           0.09928354, 0.66088525, 0.73030111, 0.83340312,
                           0.83636412, 0.08635391, 0.86285983, 0.56180810,
                           0.46534083, 0.66092477, 0.93333843};

  const T D5x3_data[15] = {0.68077325, 0.66088525, 0.86285983, 0.00794715,
                           0.73030111, 0.56180810, 0.79517780, 0.83340312,
                           0.46534083, 0.46285641, 0.83636412, 0.66092477,
                           0.09928354, 0.08635391, 0.93333843};

  const T CD_data[10] = {1.0936913458176438, 0.8882596733670376,
                         1.1284121066087023, 1.1286227974087772,
                         0.4564821574340047, 0.9033091281150568,
                         0.5875147480902937, 0.5600587437694527,
                         0.7189549926500818, 0.8762067495250275};
};

TEST_F(MatCoreTest, IsRowMajor) {
  Mat2x3 C(C2x3_data);
  EXPECT_MAT_EQ(2, 3, C, C2x3_data);
}

TEST_F(MatCoreTest, MatMatMultCore) {
  const Mat2x3 C(C2x3_data);
  const Mat3x5 D(D3x5_data);
  Mat2x5 CD;
  A2D::MatMatMultCore<2, 3, 5, Mat2x3, Mat3x5, Mat2x5>(C, D, CD);
  EXPECT_MAT_EQ(2, 5, CD, CD_data);
}

TEST_F(MatCoreTest, MatTransMatMultCore) {
  const Mat3x2 CT(C3x2_data);
  const Mat3x5 D(D3x5_data);
  Mat2x5 CD;
  A2D::MatTransMatMultCore<2, 3, 5, Mat3x2, Mat3x5, Mat2x5>(CT, D, CD);
  EXPECT_MAT_EQ(2, 5, CD, CD_data);
}

TEST_F(MatCoreTest, MatMatTransMultCore) {
  const Mat2x3 C(C2x3_data);
  const Mat5x3 DT(D5x3_data);
  Mat2x5 CD;
  A2D::MatMatTransMultCore<2, 3, 5, Mat2x3, Mat5x3, Mat2x5>(C, DT, CD);
  EXPECT_MAT_EQ(2, 5, CD, CD_data);
}

TEST_F(MatCoreTest, MatTransMatTransMultCore) {
  const Mat3x2 CT(C3x2_data);
  const Mat5x3 DT(D5x3_data);
  Mat2x5 CD;
  A2D::MatTransMatTransMultCore<2, 3, 5, Mat3x2, Mat5x3, Mat2x5>(CT, DT, CD);
  EXPECT_MAT_EQ(2, 5, CD, CD_data);
}

/**
 * Symm3x3SymmMatMult...Core
 */

TEST_F(MatCoreTest, Symm3x3SymmMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Symm3x3SymmMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, AB_data);
}

TEST_F(MatCoreTest, Symm3x3SymmMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Symm3x3SymmMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Symm3x3SymmMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3SymmMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3SymmMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3SymmMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3SymmMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Symm3x3SymmMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Symm3x3MatMult...Core
 */

TEST_F(MatCoreTest, Symm3x3MatMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Symm3x3MatMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, AB_data);
}

TEST_F(MatCoreTest, Symm3x3MatMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Symm3x3MatMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Symm3x3MatMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3MatMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3MatMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3MatMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3MatMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Symm3x3MatMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Symm3x3MatTransMult...Core
 */

TEST_F(MatCoreTest, Symm3x3MatTransMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Symm3x3MatTransMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, ABT_data);
}

TEST_F(MatCoreTest, Symm3x3MatTransMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Symm3x3MatTransMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * ABT_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Symm3x3MatTransMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3MatTransMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3MatTransMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Symm3x3MatTransMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Symm3x3MatTransMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Symm3x3MatTransMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B,
                                                                  C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Mat3x3SymmMult...Core
 */

TEST_F(MatCoreTest, Mat3x3SymmMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Mat3x3SymmMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, AB_data);
}

TEST_F(MatCoreTest, Mat3x3SymmMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Mat3x3SymmMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Mat3x3SymmMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3SymmMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3SymmMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3SymmMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3SymmMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Mat3x3SymmMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * MatTrans3x3SymmMult...Core
 */

TEST_F(MatCoreTest, MatTrans3x3SymmMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::MatTrans3x3SymmMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, ATB_data);
}

TEST_F(MatCoreTest, MatTrans3x3SymmMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::MatTrans3x3SymmMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * ATB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, MatTrans3x3SymmMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3SymmMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3SymmMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3SymmMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3SymmMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::MatTrans3x3SymmMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B,
                                                                  C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Mat3x3MatMult...Core
 */

TEST_F(MatCoreTest, Mat3x3MatMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Mat3x3MatMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, AB_data);
}

TEST_F(MatCoreTest, Mat3x3MatMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Mat3x3MatMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Mat3x3MatMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3MatMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3MatMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3MatMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3MatMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Mat3x3MatMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * AB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Mat3x3MatTransMult...Core
 */

TEST_F(MatCoreTest, Mat3x3MatTransMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::Mat3x3MatTransMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, ABT_data);
}

TEST_F(MatCoreTest, Mat3x3MatTransMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::Mat3x3MatTransMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * ABT_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, Mat3x3MatTransMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3MatTransMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3MatTransMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::Mat3x3MatTransMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, Mat3x3MatTransMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::Mat3x3MatTransMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B,
                                                                 C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * ABT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Mat3Transx3MatMult...Core
 */

TEST_F(MatCoreTest, MatTrans3x3MatMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::MatTrans3x3MatMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, ATB_data);
}

TEST_F(MatCoreTest, MatTrans3x3MatMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::MatTrans3x3MatMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B, AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * ATB_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3MatMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3MatMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::MatTrans3x3MatMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B,
                                                                 C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * ATB_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

/**
 * Mat3Transx3MatTransMult...Core
 */

TEST_F(MatCoreTest, MatTrans3x3MatTransMultCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  A2D::MatTrans3x3MatTransMultCore<Mat3x3, Mat3x3, Mat3x3>(A, B, AB);
  EXPECT_MAT_EQ(3, 3, AB, ATBT_data);
}

TEST_F(MatCoreTest, MatTrans3x3MatTransMultScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 AB;
  T scale = 2.34;
  A2D::MatTrans3x3MatTransMultScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A, B,
                                                                   AB);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = scale * ATBT_data[i];
  EXPECT_MAT_EQ(3, 3, AB, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatTransMultAddCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3MatTransMultAddCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + ATBT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatTransMultSubCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  A2D::MatTrans3x3MatTransMultSubCore<Mat3x3, Mat3x3, Mat3x3>(A, B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] - ATBT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}

TEST_F(MatCoreTest, MatTrans3x3MatTransMultAddScaleCore) {
  const Mat3x3 A(A3x3_data);
  const Mat3x3 B(B3x3_data);
  Mat3x3 C(Z3x3_data);
  T scale = 2.34;
  A2D::MatTrans3x3MatTransMultAddScaleCore<T, Mat3x3, Mat3x3, Mat3x3>(scale, A,
                                                                      B, C);
  T res[9];
  for (I i = 0; i < 9; i++) res[i] = Z3x3_data[i] + scale * ATBT_data[i];
  EXPECT_MAT_EQ(3, 3, C, res);
}