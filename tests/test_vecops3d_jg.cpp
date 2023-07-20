/*
    This is a set of automatically generated unit tests for a2dvecops3d.h using
    Google Test framework.  These tests were written on 2022-10-18 using the
    A2DTestConstructor package.
*/

#include <gtest/gtest.h>

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2dvecops3d.h"
#include "test_commons.h"

using T = double;                    /*UNQ_TC_TD_01*/
using ADScalar_t = A2D::ADScalar<T>; /*UNQ_TC_TD_01*/
using I = int;                       /*UNQ_TC_TD_01*/
using Mat_t = A2D::Mat<T, 3, 3>;     /*UNQ_TC_TD_01*/
using ADMat_t = A2D::ADMat<Mat_t>;   /*UNQ_TC_TD_01*/
using Vec_t = A2D::Vec<T, 3>;        /*UNQ_TC_TD_01*/
using ADVec_t = A2D::ADVec<Vec_t>;   /*UNQ_TC_TD_01*/

class vecops3d_jgTest : public ::testing::Test {
 protected:
  // Input Options:
  const T x_data[3] = {0.20970639762168464, 0.10308079207850707,
                       0.9173999579053238};
  const T a_data = 0.6896026553439034;
  const T v_data[3] = {0.4242524109372384, 0.28692392009681067,
                       0.9492874123274755};
  const T y_data[3] = {0.6854063555079019, 0.7322747860081501,
                       0.3571801526014854};
  const T S_data[9] = {
      0.562738676354188,  0.7541181002046244,    0.5766762373185819,
      0.7541181002046244, 0.0011522173917457579, 1.4376110744805093,
      0.5766762373185819, 1.4376110744805093,    0.5936808838281609};
  const T xb_data[3] = {0.9028664628429675, 0.503710207783976,
                        0.715447313858641};
  const T ab_data = 0.2434746529642985;
  const T vb_data[3] = {0.7366391012272651, 0.6723153242502288,
                        0.22166250342375682};
  const T yb_data[3] = {0.9537280081498468, 0.7756965964838682,
                        0.7088390653667941};
  const T Sb_data[9] = {
      1.7739061006868846, 1.046007355988558,   0.9127544973325203,
      1.046007355988558,  0.7739829458714851,  0.23883618201612167,
      0.9127544973325203, 0.23883618201612167, 0.8576426401387423};
}; /*UNQ_TC_TIC_01*/

class Vec3Norm : public vecops3d_jgTest {
 protected:
}; /*UNQ_TC_TF_01*/

class Vec3Norm_Vec : public Vec3Norm {
 protected:
  // Results
  const T a_out = 0.9466916634594847;
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Norm_Vec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data); /*UNQ_T1F_TFP_01*/
  T a;             /*UNQ_T1F_TFP_02*/
  // Evaluations:
  A2D::Vec3Norm(x, a); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VAL_NEAR(a, a_out); /*UNQ_T1F_TFP_06*/
}

class Vec3Norm_ADVec : public Vec3Norm {
 protected:
  // Results
  const T a_out = 0.9466916634594847;
  const T ab_out = 0.9481556571332816;
  const T xb_out[3] = {0.053933286154386834, 0.02651080710575372,
                       0.23594127317465283};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Norm_ADVec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data);       /*UNQ_T1F_TFP_01*/
  Vec_t xb;              /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(0, 0); /*UNQ_T1F_TFP_03*/
  ADVec_t x_ad(x, xb);   /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Norm(x_ad, a_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Norm_ADVec, forward) {
  // Declarations and Initializations:
  Vec_t x(x_data), xb(xb_data); /*UNQ_T1F_TFF_01*/
  ADScalar_t a_ad(0, 0);        /*UNQ_T1F_TFF_03*/
  ADVec_t x_ad(x, xb);          /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Norm(x_ad, a_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                        /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Norm_ADVec, reverse) {
  // Declarations and Initializations:
  T ab(ab_data);          /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data);        /*UNQ_T1F_TFR_01*/
  Vec_t xb;               /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(0, ab); /*UNQ_T1F_TFR_03*/
  ADVec_t x_ad(x, xb);    /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Norm(x_ad, a_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                        /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);        /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
}

class Vec3Scale : public vecops3d_jgTest {
 protected:
}; /*UNQ_TC_TF_01*/

class Vec3Scale_VecScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.14461408864251815, 0.07108478793229128,
                      0.6326414469838965};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Scale_VecScalar, passive) {
  // Declarations and Initializations:
  T a(a_data);     /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data); /*UNQ_T1F_TFP_01*/
  Vec_t v;         /*UNQ_T1F_TFP_02*/
  // Evaluations:
  A2D::Vec3Scale(x, a, v); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v, v_out); /*UNQ_T1F_TFP_06*/
}

class Vec3Scale_ADVecScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.14461408864251815, 0.07108478793229128,
                      0.6326414469838965};
  const T vb_out[3] = {0.6226191101974682, 0.3473598968116592,
                       0.493374367395582};
  const T xb_out[3] = {0.5079882802364685, 0.4636304328313552,
                       0.15285905095119978};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Scale_ADVecScalar, passive) {
  // Declarations and Initializations:
  T a(a_data);                      /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);                  /*UNQ_T1F_TFP_01*/
  Vec_t v, vb, xb;                  /*UNQ_T1F_TFP_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Scale(x_ad, a, v_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Scale_ADVecScalar, forward) {
  // Declarations and Initializations:
  T a(a_data);                      /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data), xb(xb_data);     /*UNQ_T1F_TFF_01*/
  Vec_t v, vb;                      /*UNQ_T1F_TFF_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x_ad, a, v_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                            /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VEC_NEAR(3, v_ad.bvalue(), vb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Scale_ADVecScalar, reverse) {
  // Declarations and Initializations:
  T a(a_data);                      /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data), vb(vb_data);     /*UNQ_T1F_TFR_01*/
  Vec_t v, xb;                      /*UNQ_T1F_TFR_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x_ad, a, v_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                            /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
}

class Vec3Scale_VecADScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.14461408864251815, 0.07108478793229128,
                      0.6326414469838965};
  const T vb_out[3] = {0.05105819238533285, 0.025097560078599514,
                       0.22336363638046075};
  const T ab_out = 0.42713389972602034;
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Scale_VecADScalar, passive) {
  // Declarations and Initializations:
  T a(a_data);           /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);       /*UNQ_T1F_TFP_01*/
  Vec_t v, vb;           /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(a, 0); /*UNQ_T1F_TFP_03*/
  ADVec_t v_ad(v, vb);   /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Scale(x, a_ad, v_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Scale_VecADScalar, forward) {
  // Declarations and Initializations:
  T a(a_data), ab(ab_data); /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data);          /*UNQ_T1F_TFF_01*/
  Vec_t v, vb;              /*UNQ_T1F_TFF_02*/
  ADScalar_t a_ad(a, ab);   /*UNQ_T1F_TFF_03*/
  ADVec_t v_ad(v, vb);      /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x, a_ad, v_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                            /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VEC_NEAR(3, v_ad.bvalue(), vb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Scale_VecADScalar, reverse) {
  // Declarations and Initializations:
  T a(a_data);                  /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data), vb(vb_data); /*UNQ_T1F_TFR_01*/
  Vec_t v;                      /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(a, 0);        /*UNQ_T1F_TFR_03*/
  ADVec_t v_ad(v, vb);          /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x, a_ad, v_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                            /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out); /*UNQ_T1F_TFR_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out);    /*UNQ_T1F_TFR_07*/
}

class Vec3Scale_ADVecADScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.14461408864251815, 0.07108478793229128,
                      0.6326414469838965};
  const T vb_out[3] = {0.673677302582801, 0.3724574568902587,
                       0.7167380037760427};
  const T xb_out[3] = {0.5079882802364685, 0.4636304328313552,
                       0.15285905095119978};
  const T ab_out = 0.42713389972602034;
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Scale_ADVecADScalar, passive) {
  // Declarations and Initializations:
  T a(a_data);                      /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);                  /*UNQ_T1F_TFP_01*/
  Vec_t v, vb, xb;                  /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(a, 0);            /*UNQ_T1F_TFP_03*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Scale(x_ad, a_ad, v_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Scale_ADVecADScalar, forward) {
  // Declarations and Initializations:
  T a(a_data), ab(ab_data);         /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data), xb(xb_data);     /*UNQ_T1F_TFF_01*/
  Vec_t v, vb;                      /*UNQ_T1F_TFF_02*/
  ADScalar_t a_ad(a, ab);           /*UNQ_T1F_TFF_03*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x_ad, a_ad, v_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                               /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VEC_NEAR(3, v_ad.bvalue(), vb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Scale_ADVecADScalar, reverse) {
  // Declarations and Initializations:
  T a(a_data);                      /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data), vb(vb_data);     /*UNQ_T1F_TFR_01*/
  Vec_t v, xb;                      /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(a, 0);            /*UNQ_T1F_TFR_03*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Scale(x_ad, a_ad, v_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                               /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out);      /*UNQ_T1F_TFR_07*/
}

class Vec3Dot : public vecops3d_jgTest {
 protected:
}; /*UNQ_TC_TF_01*/

class Vec3Dot_VecVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.546894619642629;
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Dot_VecVec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data), y(y_data); /*UNQ_T1F_TFP_01*/
  T a;                        /*UNQ_T1F_TFP_02*/
  // Evaluations:
  A2D::Vec3Dot(x, y, a); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VAL_NEAR(a, a_out); /*UNQ_T1F_TFP_06*/
}

class Vec3Dot_ADVecVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.546894619642629;
  const T ab_out = 1.243228277164993;
  const T xb_out[3] = {0.166879074546811, 0.17829034939784028,
                       0.08696431370038184};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Dot_ADVecVec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data), y(y_data); /*UNQ_T1F_TFP_01*/
  Vec_t xb;                   /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(0, 0);      /*UNQ_T1F_TFP_03*/
  ADVec_t x_ad(x, xb);        /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Dot(x_ad, y, a_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Dot_ADVecVec, forward) {
  // Declarations and Initializations:
  Vec_t x(x_data), y(y_data), xb(xb_data); /*UNQ_T1F_TFF_01*/
  ADScalar_t a_ad(0, 0);                   /*UNQ_T1F_TFF_03*/
  ADVec_t x_ad(x, xb);                     /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Dot(x_ad, y, a_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                          /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Dot_ADVecVec, reverse) {
  // Declarations and Initializations:
  T ab(ab_data);              /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data), y(y_data); /*UNQ_T1F_TFR_01*/
  Vec_t xb;                   /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(0, ab);     /*UNQ_T1F_TFR_03*/
  ADVec_t x_ad(x, xb);        /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Dot(x_ad, y, a_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                          /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);        /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
}

class Vec3Dot_ADVecADVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.546894619642629;
  const T ab_out = 2.173479490372307;
  const T xb_out[3] = {0.166879074546811, 0.17829034939784028,
                       0.08696431370038184};
  const T yb_out[3] = {0.05105819238533286, 0.025097560078599517,
                       0.22336363638046075};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Dot_ADVecADVec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data), y(y_data);       /*UNQ_T1F_TFP_01*/
  Vec_t xb, yb;                     /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(0, 0);            /*UNQ_T1F_TFP_03*/
  ADVec_t x_ad(x, xb), y_ad(y, yb); /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Dot(x_ad, y_ad, a_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Dot_ADVecADVec, forward) {
  // Declarations and Initializations:
  Vec_t x(x_data), y(y_data), xb(xb_data), yb(yb_data); /*UNQ_T1F_TFF_01*/
  ADScalar_t a_ad(0, 0);                                /*UNQ_T1F_TFF_03*/
  ADVec_t x_ad(x, xb), y_ad(y, yb);                     /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Dot(x_ad, y_ad, a_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                             /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Dot_ADVecADVec, reverse) {
  // Declarations and Initializations:
  T ab(ab_data);                    /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data), y(y_data);       /*UNQ_T1F_TFR_01*/
  Vec_t xb, yb;                     /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(0, ab);           /*UNQ_T1F_TFR_03*/
  ADVec_t x_ad(x, xb), y_ad(y, yb); /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Dot(x_ad, y_ad, a_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                             /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VAL_NEAR(a_ad.value, a_out);        /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
  EXPECT_VEC_NEAR(3, y_ad.bvalue(), yb_out); /*UNQ_T1F_TFR_07*/
}

class Vec3Normalize : public vecops3d_jgTest {
 protected:
}; /*UNQ_TC_TF_01*/

class Vec3Normalize_Vec : public Vec3Normalize {
 protected:
  // Results
  const T v_out[3] = {0.22151499344079667, 0.10888528552350413,
                      0.9690588745155729};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Normalize_Vec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data); /*UNQ_T1F_TFP_01*/
  Vec_t v;         /*UNQ_T1F_TFP_02*/
  // Evaluations:
  A2D::Vec3Normalize(x, v); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v, v_out); /*UNQ_T1F_TFP_06*/
}

class Vec3Normalize_ADVec : public Vec3Normalize {
 protected:
  // Results
  const T v_out[3] = {0.22151499344079667, 0.10888528552350413,
                      0.9690588745155729};
  const T vb_out[3] = {0.7318494451936044, 0.42302052906313753,
                       -0.21482320797600826};
  const T xb_out[3] = {0.67254703009868, 0.6582796133772945,
                       -0.22770155706806752};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3Normalize_ADVec, passive) {
  // Declarations and Initializations:
  Vec_t x(x_data);                  /*UNQ_T1F_TFP_01*/
  Vec_t v, vb, xb;                  /*UNQ_T1F_TFP_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3Normalize(x_ad, v_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3Normalize_ADVec, forward) {
  // Declarations and Initializations:
  Vec_t x(x_data), xb(xb_data);     /*UNQ_T1F_TFF_01*/
  Vec_t v, vb;                      /*UNQ_T1F_TFF_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr = A2D::Vec3Normalize(x_ad, v_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                             /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_VEC_NEAR(3, v_ad.bvalue(), vb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3Normalize_ADVec, reverse) {
  // Declarations and Initializations:
  Vec_t x(x_data), vb(vb_data);     /*UNQ_T1F_TFR_01*/
  Vec_t v, xb;                      /*UNQ_T1F_TFR_02*/
  ADVec_t x_ad(x, xb), v_ad(v, vb); /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr = A2D::Vec3Normalize(x_ad, v_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                             /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_VEC_NEAR(3, v_ad.value(), v_out);   /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out); /*UNQ_T1F_TFR_07*/
}

class Vec3ScaleSymmetricOuterProduct : public vecops3d_jgTest {
 protected:
}; /*UNQ_TC_TF_01*/

class Vec3ScaleSymmetricOuterProduct_ScalarVec
    : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {
      0.030326499574565464, 0.014906934802982205, 0.13266895883316293,
      0.014906934802982205, 0.007327476244793286, 0.06521318145679289,
      0.13266895883316293,  0.06521318145679289,  0.5803852368321898};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarVec, passive) {
  // Declarations and Initializations:
  T a(a_data);     /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data); /*UNQ_T1F_TFP_01*/
  Mat_t S;         /*UNQ_T1F_TFP_02*/
  // Evaluations:
  A2D::Vec3ScaleSymmetricOuterProduct(a, x, S); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S, S_out); /*UNQ_T1F_TFP_06*/
}

class Vec3ScaleSymmetricOuterProduct_ADScalarVec
    : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {
      0.030326499574565464, 0.014906934802982205, 0.13266895883316293,
      0.014906934802982205, 0.007327476244793286, 0.06521318145679289,
      0.13266895883316293,  0.06521318145679289,  0.5803852368321898};
  const T Sb_out[9] = {
      0.010707229594203083, 0.005263118913176909,  0.046840783545026285,
      0.005263118913176909, 0.0025870763721399562, 0.02302450055963353,
      0.046840783545026285, 0.02302450055963353,   0.20491379061301473};
  const T ab_out = 1.2496403015961977;
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarVec, passive) {
  // Declarations and Initializations:
  T a(a_data);           /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);       /*UNQ_T1F_TFP_01*/
  Mat_t S, Sb;           /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(a, 0); /*UNQ_T1F_TFP_03*/
  ADMat_t S_ad(S, Sb);   /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x, S_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarVec, forward) {
  // Declarations and Initializations:
  T a(a_data), ab(ab_data); /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data);          /*UNQ_T1F_TFF_01*/
  Mat_t S, Sb;              /*UNQ_T1F_TFF_02*/
  ADScalar_t a_ad(a, ab);   /*UNQ_T1F_TFF_03*/
  ADMat_t S_ad(S, Sb);      /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x, S_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                                         /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_MAT_NEAR(3, 3, S_ad.bvalue(), Sb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarVec, reverse) {
  // Declarations and Initializations:
  T a(a_data);           /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data);       /*UNQ_T1F_TFR_01*/
  Mat_t Sb(Sb_data);     /*UNQ_T1F_TFR_01*/
  Mat_t S;               /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(a, 0); /*UNQ_T1F_TFR_03*/
  ADMat_t S_ad(S, Sb);   /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x, S_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                                         /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFR_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out);       /*UNQ_T1F_TFR_07*/
}

class Vec3ScaleSymmetricOuterProduct_ScalarADVec
    : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {
      0.030326499574565464, 0.014906934802982205, 0.13266895883316293,
      0.014906934802982205, 0.007327476244793286, 0.06521318145679289,
      0.13266895883316293,  0.06521318145679289,  0.5803852368321898};
  const T Sb_out[9] = {
      0.2611344213798594,  0.13702366368098343, 0.6746545067516124,
      0.13702366368098343, 0.07161226659930862, 0.36952537529538265,
      0.6746545067516124,  0.36952537529538265, 0.9052432477605453};
  const T xb_out[3] = {1.8166667021956102, 0.714766963705857,
                       1.3831101199754792};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarADVec, passive) {
  // Declarations and Initializations:
  T a(a_data);         /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);     /*UNQ_T1F_TFP_01*/
  Vec_t xb;            /*UNQ_T1F_TFP_02*/
  Mat_t S, Sb;         /*UNQ_T1F_TFP_02*/
  ADVec_t x_ad(x, xb); /*UNQ_T1F_TFP_03*/
  ADMat_t S_ad(S, Sb); /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3ScaleSymmetricOuterProduct(a, x_ad, S_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarADVec, forward) {
  // Declarations and Initializations:
  T a(a_data);                  /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data), xb(xb_data); /*UNQ_T1F_TFF_01*/
  Mat_t S, Sb;                  /*UNQ_T1F_TFF_02*/
  ADVec_t x_ad(x, xb);          /*UNQ_T1F_TFF_03*/
  ADMat_t S_ad(S, Sb);          /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a, x_ad, S_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                                         /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_MAT_NEAR(3, 3, S_ad.bvalue(), Sb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarADVec, reverse) {
  // Declarations and Initializations:
  T a(a_data);         /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data);     /*UNQ_T1F_TFR_01*/
  Mat_t Sb(Sb_data);   /*UNQ_T1F_TFR_01*/
  Vec_t xb;            /*UNQ_T1F_TFR_02*/
  Mat_t S;             /*UNQ_T1F_TFR_02*/
  ADVec_t x_ad(x, xb); /*UNQ_T1F_TFR_03*/
  ADMat_t S_ad(S, Sb); /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a, x_ad, S_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                                         /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFR_06*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
}

class Vec3ScaleSymmetricOuterProduct_ADScalarADVec
    : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {
      0.030326499574565464, 0.014906934802982205, 0.13266895883316293,
      0.014906934802982205, 0.007327476244793286, 0.06521318145679289,
      0.13266895883316293,  0.06521318145679289,  0.5803852368321898};
  const T Sb_out[9] = {
      0.27184165097406254, 0.14228678259416036, 0.7214952902966387,
      0.14228678259416036, 0.07419934297144858, 0.3925498758550162,
      0.7214952902966387,  0.3925498758550162,  1.1101570383735602};
  const T ab_out = 1.2496403015961977;
  const T xb_out[3] = {1.8166667021956102, 0.714766963705857,
                       1.3831101199754792};
}; /*UNQ_T1F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarADVec, passive) {
  // Declarations and Initializations:
  T a(a_data);           /*UNQ_T1F_TFP_01*/
  Vec_t x(x_data);       /*UNQ_T1F_TFP_01*/
  Vec_t xb;              /*UNQ_T1F_TFP_02*/
  Mat_t S, Sb;           /*UNQ_T1F_TFP_02*/
  ADScalar_t a_ad(a, 0); /*UNQ_T1F_TFP_03*/
  ADVec_t x_ad(x, xb);   /*UNQ_T1F_TFP_03*/
  ADMat_t S_ad(S, Sb);   /*UNQ_T1F_TFP_03*/
  // Evaluations:
  A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x_ad, S_ad); /*UNQ_T1F_TFP_04*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarADVec, forward) {
  // Declarations and Initializations:
  T a(a_data), ab(ab_data);     /*UNQ_T1F_TFF_01*/
  Vec_t x(x_data), xb(xb_data); /*UNQ_T1F_TFF_01*/
  Mat_t S, Sb;                  /*UNQ_T1F_TFF_02*/
  ADScalar_t a_ad(a, ab);       /*UNQ_T1F_TFF_03*/
  ADVec_t x_ad(x, xb);          /*UNQ_T1F_TFF_03*/
  ADMat_t S_ad(S, Sb);          /*UNQ_T1F_TFF_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x_ad, S_ad); /*UNQ_T1F_TFF_04*/
  expr.forward();                                            /*UNQ_T1F_TFF_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out);   /*UNQ_T1F_TFF_06*/
  EXPECT_MAT_NEAR(3, 3, S_ad.bvalue(), Sb_out); /*UNQ_T1F_TFF_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ADScalarADVec, reverse) {
  // Declarations and Initializations:
  T a(a_data);           /*UNQ_T1F_TFR_01*/
  Vec_t x(x_data);       /*UNQ_T1F_TFR_01*/
  Mat_t Sb(Sb_data);     /*UNQ_T1F_TFR_01*/
  Vec_t xb;              /*UNQ_T1F_TFR_02*/
  Mat_t S;               /*UNQ_T1F_TFR_02*/
  ADScalar_t a_ad(a, 0); /*UNQ_T1F_TFR_03*/
  ADVec_t x_ad(x, xb);   /*UNQ_T1F_TFR_03*/
  ADMat_t S_ad(S, Sb);   /*UNQ_T1F_TFR_03*/
  // Evaluations:
  auto expr =
      A2D::Vec3ScaleSymmetricOuterProduct(a_ad, x_ad, S_ad); /*UNQ_T1F_TFR_04*/
  expr.reverse();                                            /*UNQ_T1F_TFR_05*/
  // Comparisons:
  EXPECT_MAT_NEAR(3, 3, S_ad.value(), S_out); /*UNQ_T1F_TFR_06*/
  EXPECT_VAL_NEAR(a_ad.bvalue, ab_out);       /*UNQ_T1F_TFR_07*/
  EXPECT_VEC_NEAR(3, x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
}
