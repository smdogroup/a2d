/*
    This is a set of automatically generated unit tests for
    a2d_vecops3d.h A2D classes using Google Test framework.
    These tests were written on 2022-10-19 using the
    A2DTestConstructor package.
*/

#include <gtest/gtest.h>

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2d_vecops3d.h"
#include "test_commons.h"

using T = double;  /*UNQ_TC_TD_01*/
using Vec_t = A2D::Vec<T, 3>;  /*UNQ_TC_TD_01*/
using A2DVec_t = A2D::A2DVec<4, Vec_t>;  /*UNQ_TC_TD_01*/
using Mat_t = A2D::Mat<T, 3, 3>;  /*UNQ_TC_TD_01*/
using A2DMat_t = A2D::A2DMat<4, Mat_t>;  /*UNQ_TC_TD_01*/
using A2DScalar_t = A2D::A2DScalar<4, T>;  /*UNQ_TC_TD_01*/

class vecops3d_a2dTest : public ::testing::Test {
 protected:
  // Input Options:
  const T x_data[3] = {0.16243788741867204, 0.40274570395990394, 0.27523253509698253};
  const T a_data = 0.18405076479370952;
  const T v_data[3] = {0.5810612910693836, 0.33123024458158556, 0.9957906568393305};
  const T y_data[3] = {0.17821294679646638, 0.6163598093229385, 0.8043520507617771};
  const T S_data[9] = {1.863562352437055, 1.3443467625962775, 1.058138824106408,
                       1.3443467625962775, 0.015745248225881303, 0.8996060500091604,
                       1.058138824106408, 0.8996060500091604, 0.9006309693314773};
  const T xb_data[3] = {0.5203786679056438, 0.7385676342426261, 0.7861175322208778};
  const T ab_data = 0.8310608779744129;
  const T vb_data[3] = {0.1902556263764923, 0.4436002353179346, 0.7894772286747886};
  const T yb_data[3] = {0.299853240650116, 0.7001565143143994, 0.8224454419936111};
  const T Sb_data[9] = {0.06670026152713815, 0.8635324135237004, 1.4060804827850246,
                        0.8635324135237004, 0.9965634965451609, 0.5668669985102699,
                        1.4060804827850246, 0.5668669985102699, 1.8380417880847526};
  const T xp0_data[3] = {0.46160113373397027, 0.41045062905653673, 0.4700539443775993};
  const T xp1_data[3] = {0.8911122254787387, 0.949361170617268, 0.5167692781757652};
  const T xp2_data[3] = {0.06502785084463514, 0.3585301972593763, 0.5047695105202312};
  const T xp3_data[3] = {0.9848219195093173, 0.5295612534322616, 0.1479327772454142};
  const T ap0_data = 0.5608739860314739;
  const T ap1_data = 0.8631380642543824;
  const T ap2_data = 0.6416905992000659;
  const T ap3_data = 0.7107021458720669;
  const T yp0_data[3] = {0.2191829210696652, 0.18570803343349163, 0.7174660016775594};
  const T yp1_data[3] = {0.8981869831204218, 0.8606184191151188, 0.7808097295318616};
  const T yp2_data[3] = {0.46013396162004305, 0.2595553895946646, 0.6629536689667503};
  const T yp3_data[3] = {0.12034601412697621, 0.14227338002464518, 0.4556090504915631};
  const T ah0_data = 0.6425352568023172;
  const T ah1_data = 0.14810726064215562;
  const T ah2_data = 0.2852586607737033;
  const T ah3_data = 0.4389438160460044;
  const T vh0_data[3] = {0.6284963202727984, 0.45211785521096726, 0.38502976064339633};
  const T vh1_data[3] = {0.5631258604722853, 0.7356622956967018, 0.04841947244536282};
  const T vh2_data[3] = {0.6926878107409069, 0.9508293682598007, 0.3360382764620453};
  const T vh3_data[3] = {0.7102962012036355, 0.6798507993191757, 0.06371273366041486};
  const T Sh0_data[9] = {1.9841707097671282, 1.5044718748353025, 1.1508619072581627,
                         1.5044718748353025, 0.876045779952817, 0.49868484225370724,
                         1.1508619072581627, 0.49868484225370724, 0.11397213013363072};
  const T Sh1_data[9] = {1.9043267104001387, 1.4523500655861972, 1.3688028431673138,
                         1.4523500655861972, 0.6549482428477225, 0.8855822513223802,
                         1.3688028431673138, 0.8855822513223802, 1.8016647464593099};
  const T Sh2_data[9] = {1.695003537757514, 0.46894972234636934, 0.37128077178030694,
                         0.46894972234636934, 0.47326839606702564, 0.7824026727498676,
                         0.37128077178030694, 0.7824026727498676, 0.31512053477901025};
  const T Sh3_data[9] = {0.6180882418677669, 0.9212798836963006, 0.6667389462458889,
                         0.9212798836963006, 0.15740374169970606, 0.8836150127528192,
                         0.6667389462458889, 0.8836150127528192, 1.6469688448147843};
};  /*UNQ_TC_TIC_01*/

class Vec3Norm : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Norm_A2DVec : public Vec3Norm {
 protected:
  // Results
  const T a_out = 0.5141430906888776;
  const T xb_out[3] = {0.26256459685881517, 0.6509981450589614, 0.44488586233523};
  const T ap0_out = 0.7189883244021776;
  const T ap1_out = 1.301842329394938;
  const T ap2_out = 0.5716083338338123;
  const T ap3_out = 0.8051587830203008;
  const T xh0_out[3] = {0.5819582896210939, 0.25640276295887393, 0.48162139137492566};
  const T xh1_out[3] = {0.8223567336407702, 0.0021969553815232024, -0.2118876676640889};
  const T xh2_out[3] = {-0.09667580004665523, 0.07922143934266154, 0.47400460088537144};
  const T xh3_out[3] = {1.3193639985752674, 0.18034571294580926, -0.22260533784296188};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Norm_A2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
Vec_t xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(0, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Norm(x_a2d, a_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Norm_A2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
T ab(ab_data);  /*UNQ_T2F_TFR_01*/
Vec_t xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Norm(x_a2d, a_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Norm_A2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHF_01*/
T ab(ab_data);  /*UNQ_T2F_TFHF_01*/
Vec_t xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Norm(x_a2d, a_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Norm_A2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHR_01*/
T ab(ab_data), ah0(ah0_data), ah1(ah1_data), ah2(ah2_data), ah3(ah3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.hvalue[0] = ah0;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[1] = ah1;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[2] = ah2;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[3] = ah3;  /*UNQ_T2F_TFHR_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
}
// Evaluations:
auto expr = A2D::Vec3Norm(x_a2d, a_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Scale : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Scale_A2DVecScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.029896817410881075, 0.07412565483120125, 0.050656758580711135};
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T vp0_out[3] = {0.08495804169338061, 0.07554375218791476, 0.08651378795699693};
  const T vp1_out[3] = {0.16400988661638638, 0.1747306495175595, 0.09511178087014281};
  const T vp2_out[3] = {0.011968425680846368, 0.06598775700722775, 0.09290321445579496};
  const T vp3_out[3] = {0.18125722747129888, 0.09746615369932317, 0.027227140790075953};
  const T xh0_out[3] = {0.11567522841624076, 0.08321263702847015, 0.07086502191475601};
  const T xh1_out[3] = {0.10364374529503986, 0.13539920815287404, 0.008911640934476972};
  const T xh2_out[3] = {0.12748972133014425, 0.175000872416536, 0.06184810178279944};
  const T xh3_out[3] = {0.13073055906159567, 0.125127059560309, 0.011726377357297275};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Scale_A2DVecScalar, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Scale(x_a2d, a, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Scale_A2DVecScalar, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Scale_A2DVecScalar, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Scale_A2DVecScalar, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Scale_VecA2DScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.029896817410881075, 0.07412565483120125, 0.050656758580711135};
  const T ab_out = 0.4268526301175164;
  const T vp0_out[3] = {0.0911071853990424, 0.22588958833704328, 0.15437076904539213};
  const T vp1_out[3] = {0.1402063237081239, 0.34762514730272004, 0.23756367756343588};
  const T vp2_out[3] = {0.10423486531048051, 0.2584381320992831, 0.1766141303657359};
  const T vp3_out[3] = {0.11544495515937545, 0.28623223604505993, 0.19560835330723444};
  const T ah0_out = 0.3901528554949718;
  const T ah1_out = 0.4010844184322628;
  const T ah2_out = 0.5879498546032313;
  const T ah3_out = 0.40672182032784615;
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Scale_VecA2DScalar, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb;  /*UNQ_T2F_TFP_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Scale(x, a_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Scale_VecA2DScalar, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v;  /*UNQ_T2F_TFR_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Scale(x, a_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Scale_VecA2DScalar, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v;  /*UNQ_T2F_TFHF_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
// Evaluations:
auto expr = A2D::Vec3Scale(x, a_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Scale_VecA2DScalar, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v;  /*UNQ_T2F_TFHR_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Scale(x, a_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Scale_A2DVecA2DScalar : public Vec3Scale {
 protected:
  // Results
  const T v_out[3] = {0.029896817410881075, 0.07412565483120125, 0.050656758580711135};
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T ab_out = 0.4268526301175164;
  const T vp0_out[3] = {0.176065227092423, 0.301433340524958, 0.24088455700238906};
  const T vp1_out[3] = {0.3042162103245103, 0.5223557968202795, 0.3326754584335787};
  const T vp2_out[3] = {0.11620329099132687, 0.3244258891065109, 0.26951734482153084};
  const T vp3_out[3] = {0.2967021826306743, 0.3836983897443831, 0.2228354940973104};
  const T xh0_out[3] = {0.22238465994704593, 0.33201646921575567, 0.5136622620426942};
  const T ah0_out = 1.031147949300422;
  const T xh1_out[3] = {0.2678606183593153, 0.5182874565680089, 0.6903394878657916};
  const T ah1_out = 1.3997379493485373;
  const T xh2_out[3] = {0.24957496822098196, 0.4596549732230097, 0.568448217705964};
  const T ah2_out = 1.1578698832555172;
  const T xh3_out[3] = {0.2659456409917387, 0.4403946987101389, 0.5728095378936376};
  const T ah3_out = 0.9457927871411109;
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Scale_A2DVecA2DScalar, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Scale(x_a2d, a_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Scale_A2DVecA2DScalar, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Scale_A2DVecA2DScalar, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Scale_A2DVecA2DScalar, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Scale(x_a2d, a_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Axpy_A2DScalarA2DVecVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T ab_out = 0.4268526301175164;
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T vp0_out[3] = {0.176065227092423, 0.301433340524958, 0.24088455700238906};
  const T vp1_out[3] = {0.3042162103245103, 0.5223557968202795, 0.3326754584335787};
  const T vp2_out[3] = {0.11620329099132687, 0.3244258891065109, 0.26951734482153084};
  const T vp3_out[3] = {0.2967021826306743, 0.3836983897443831, 0.2228354940973104};
  const T ah0_out = 1.031147949300422;
  const T xh0_out[3] = {0.22238465994704593, 0.33201646921575567, 0.5136622620426942};
  const T ah1_out = 1.3997379493485373;
  const T xh1_out[3] = {0.2678606183593153, 0.5182874565680089, 0.6903394878657916};
  const T ah2_out = 1.1578698832555172;
  const T xh2_out[3] = {0.24957496822098196, 0.4596549732230097, 0.568448217705964};
  const T ah3_out = 0.9457927871411109;
  const T xh3_out[3] = {0.2659456409917387, 0.4403946987101389, 0.5728095378936376};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_A2DScalarA2DVecVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a_a2d, x_a2d, y, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_ScalarVecA2DVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T yb_out[3] = {0.1902556263764923, 0.4436002353179346, 0.7894772286747886};
  const T vp0_out[3] = {0.2191829210696652, 0.18570803343349163, 0.7174660016775594};
  const T vp1_out[3] = {0.8981869831204218, 0.8606184191151188, 0.7808097295318616};
  const T vp2_out[3] = {0.46013396162004305, 0.2595553895946646, 0.6629536689667503};
  const T vp3_out[3] = {0.12034601412697621, 0.14227338002464518, 0.4556090504915631};
  const T yh0_out[3] = {0.6284963202727984, 0.45211785521096726, 0.38502976064339633};
  const T yh1_out[3] = {0.5631258604722853, 0.7356622956967018, 0.04841947244536282};
  const T yh2_out[3] = {0.6926878107409069, 0.9508293682598007, 0.3360382764620453};
  const T yh3_out[3] = {0.7102962012036355, 0.6798507993191757, 0.06371273366041486};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_ScalarVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a, x, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_ScalarVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_ScalarVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_ScalarVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_ScalarA2DVecA2DVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T yb_out[3] = {0.1902556263764923, 0.4436002353179346, 0.7894772286747886};
  const T vp0_out[3] = {0.3041409627630458, 0.2612517856214064, 0.8039797896345562};
  const T vp1_out[3] = {1.0621968697368083, 1.0353490686326783, 0.8759215104020044};
  const T vp2_out[3] = {0.4721023873008894, 0.3255431466018924, 0.7558568834225453};
  const T vp3_out[3] = {0.3016032415982751, 0.23973953372396836, 0.482836191281639};
  const T xh0_out[3] = {0.11567522841624076, 0.08321263702847015, 0.07086502191475601};
  const T yh0_out[3] = {0.6284963202727984, 0.45211785521096726, 0.38502976064339633};
  const T xh1_out[3] = {0.10364374529503986, 0.13539920815287404, 0.008911640934476972};
  const T yh1_out[3] = {0.5631258604722853, 0.7356622956967018, 0.04841947244536282};
  const T xh2_out[3] = {0.12748972133014425, 0.175000872416536, 0.06184810178279944};
  const T yh2_out[3] = {0.6926878107409069, 0.9508293682598007, 0.3360382764620453};
  const T xh3_out[3] = {0.13073055906159567, 0.125127059560309, 0.011726377357297275};
  const T yh3_out[3] = {0.7102962012036355, 0.6798507993191757, 0.06371273366041486};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_ScalarA2DVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_ScalarA2DVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_ScalarA2DVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_ScalarA2DVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_A2DScalarVecVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T ab_out = 0.4268526301175164;
  const T vp0_out[3] = {0.0911071853990424, 0.22588958833704328, 0.15437076904539213};
  const T vp1_out[3] = {0.1402063237081239, 0.34762514730272004, 0.23756367756343588};
  const T vp2_out[3] = {0.10423486531048051, 0.2584381320992831, 0.1766141303657359};
  const T vp3_out[3] = {0.11544495515937545, 0.28623223604505993, 0.19560835330723444};
  const T ah0_out = 0.3901528554949718;
  const T ah1_out = 0.4010844184322628;
  const T ah2_out = 0.5879498546032313;
  const T ah3_out = 0.40672182032784615;
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_A2DScalarVecVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb;  /*UNQ_T2F_TFP_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a_a2d, x, y, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_A2DScalarVecVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v;  /*UNQ_T2F_TFR_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_A2DScalarVecVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v;  /*UNQ_T2F_TFHF_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_A2DScalarVecVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v;  /*UNQ_T2F_TFHR_02*/
A2DVec_t v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_A2DScalarVecA2DVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T ab_out = 0.4268526301175164;
  const T yb_out[3] = {0.1902556263764923, 0.4436002353179346, 0.7894772286747886};
  const T vp0_out[3] = {0.3102901064687076, 0.41159762177053494, 0.8718367707229515};
  const T vp1_out[3] = {1.0383933068285456, 1.2082435664178388, 1.0183734070952974};
  const T vp2_out[3] = {0.5643688269305235, 0.5179935216939477, 0.8395677993324862};
  const T vp3_out[3] = {0.23579096928635165, 0.4285056160697051, 0.6512174037987974};
  const T ah0_out = 0.3901528554949718;
  const T yh0_out[3] = {0.6284963202727984, 0.45211785521096726, 0.38502976064339633};
  const T ah1_out = 0.4010844184322628;
  const T yh1_out[3] = {0.5631258604722853, 0.7356622956967018, 0.04841947244536282};
  const T ah2_out = 0.5879498546032313;
  const T yh2_out[3] = {0.6926878107409069, 0.9508293682598007, 0.3360382764620453};
  const T ah3_out = 0.40672182032784615;
  const T yh3_out[3] = {0.7102962012036355, 0.6798507993191757, 0.06371273366041486};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_A2DScalarVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a_a2d, x, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_A2DScalarVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_A2DScalarVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_A2DScalarVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_ScalarA2DVecVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T vp0_out[3] = {0.08495804169338061, 0.07554375218791476, 0.08651378795699693};
  const T vp1_out[3] = {0.16400988661638638, 0.1747306495175595, 0.09511178087014281};
  const T vp2_out[3] = {0.011968425680846368, 0.06598775700722775, 0.09290321445579496};
  const T vp3_out[3] = {0.18125722747129888, 0.09746615369932317, 0.027227140790075953};
  const T xh0_out[3] = {0.11567522841624076, 0.08321263702847015, 0.07086502191475601};
  const T xh1_out[3] = {0.10364374529503986, 0.13539920815287404, 0.008911640934476972};
  const T xh2_out[3] = {0.12748972133014425, 0.175000872416536, 0.06184810178279944};
  const T xh3_out[3] = {0.13073055906159567, 0.125127059560309, 0.011726377357297275};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_ScalarA2DVecVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a, x_a2d, y, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_ScalarA2DVecVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_ScalarA2DVecVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_ScalarA2DVecVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a, x_a2d, y, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Axpy_A2DScalarA2DVecA2DVec : public Vec3Axpy {
 protected:
  // Results
  const T v_out[3] = {0.20810976420734745, 0.6904854641541398, 0.8550088093424882};
  const T ab_out = 0.4268526301175164;
  const T xb_out[3] = {0.03501669354089966, 0.08164496257293537, 0.14530388772481315};
  const T yb_out[3] = {0.1902556263764923, 0.4436002353179346, 0.7894772286747886};
  const T vp0_out[3] = {0.39524814816208825, 0.48714137395844964, 0.9583505586799483};
  const T vp1_out[3] = {1.202403193444932, 1.3829742159353984, 1.1134851879654402};
  const T vp2_out[3] = {0.5763372526113699, 0.5839812787011754, 0.9324710137882812};
  const T vp3_out[3] = {0.4170481967576505, 0.5259717697690283, 0.6784445445888735};
  const T ah0_out = 1.031147949300422;
  const T xh0_out[3] = {0.22238465994704593, 0.33201646921575567, 0.5136622620426942};
  const T yh0_out[3] = {0.6284963202727984, 0.45211785521096726, 0.38502976064339633};
  const T ah1_out = 1.3997379493485373;
  const T xh1_out[3] = {0.2678606183593153, 0.5182874565680089, 0.6903394878657916};
  const T yh1_out[3] = {0.5631258604722853, 0.7356622956967018, 0.04841947244536282};
  const T ah2_out = 1.1578698832555172;
  const T xh2_out[3] = {0.24957496822098196, 0.4596549732230097, 0.568448217705964};
  const T yh2_out[3] = {0.6926878107409069, 0.9508293682598007, 0.3360382764620453};
  const T ah3_out = 0.9457927871411109;
  const T xh3_out[3] = {0.2659456409917387, 0.4403946987101389, 0.5728095378936376};
  const T yh3_out[3] = {0.7102962012036355, 0.6798507993191757, 0.06371273366041486};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Axpy_A2DScalarA2DVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Axpy(a_a2d, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Axpy_A2DScalarA2DVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Axpy(a_a2d, x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Dot : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Dot_A2DVecVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.49856865392825395;
  const T xb_out[3] = {0.14810580803107867, 0.512232524284063, 0.668465521506602};
  const T ap0_out = 0.713337423877872;
  const T ap1_out = 1.1596202344005444;
  const T ap2_out = 0.6385847998911346;
  const T ap3_out = 0.620898322288067;
  const T xh0_out[3] = {0.11450810153536521, 0.39603290836594157, 0.5168245514956888};
  const T xh1_out[3] = {0.02639463136099086, 0.0912873629287418, 0.11913037883022691};
  const T xh2_out[3] = {0.05083678653569524, 0.17582197376219658, 0.22944838879088636};
  const T xh3_out[3] = {0.07822547093564451, 0.2705473267615983, 0.3530653586058039};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Dot_A2DVecVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(0, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Dot(x_a2d, y, a_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Dot_A2DVecVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFR_01*/
T ab(ab_data);  /*UNQ_T2F_TFR_01*/
Vec_t xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y, a_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Dot_A2DVecVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHF_01*/
T ab(ab_data);  /*UNQ_T2F_TFHF_01*/
Vec_t xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y, a_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Dot_A2DVecVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHR_01*/
T ab(ab_data), ah0(ah0_data), ah1(ah1_data), ah2(ah2_data), ah3(ah3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.hvalue[0] = ah0;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[1] = ah1;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[2] = ah2;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[3] = ah3;  /*UNQ_T2F_TFHR_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y, a_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Dot_VecA2DVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.49856865392825395;
  const T yb_out[3] = {0.1349957733344704, 0.3347061983333407, 0.22873499226482172};
  const T ap0_out = 0.3078667098006015;
  const T ap1_out = 0.7074142083798535;
  const T ap2_out = 0.3617444257155989;
  const T ap3_out = 0.20224717886674634;
  const T yh0_out[3] = {0.10437206970698233, 0.25877831431990683, 0.17684660761889245};
  const T yh1_out[3] = {0.02405823053007839, 0.059649562948897936, 0.040763936812810035};
  const T yh2_out[3] = {0.04633681422395998, 0.11488670014396457, 0.07851246436311654};
  const T yh3_out[3] = {0.07130110617400316, 0.1767827361922946, 0.12081161925548535};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Dot_VecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(0, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Dot(x, y_a2d, a_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Dot_VecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFR_01*/
T ab(ab_data);  /*UNQ_T2F_TFR_01*/
Vec_t yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Dot(x, y_a2d, a_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Dot_VecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data);  /*UNQ_T2F_TFHF_01*/
T ab(ab_data);  /*UNQ_T2F_TFHF_01*/
Vec_t yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x, y_a2d, a_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Dot_VecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data);  /*UNQ_T2F_TFHR_01*/
T ab(ab_data), ah0(ah0_data), ah1(ah1_data), ah2(ah2_data), ah3(ah3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.hvalue[0] = ah0;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[1] = ah1;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[2] = ah2;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[3] = ah3;  /*UNQ_T2F_TFHR_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x, y_a2d, a_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Dot_A2DVecA2DVec : public Vec3Dot {
 protected:
  // Results
  const T a_out = 0.49856865392825395;
  const T xb_out[3] = {0.14810580803107867, 0.512232524284063, 0.668465521506602};
  const T yb_out[3] = {0.1349957733344704, 0.3347061983333407, 0.22873499226482172};
  const T ap0_out = 1.0212041336784736;
  const T ap1_out = 1.8670344427803978;
  const T ap2_out = 1.0003292256067335;
  const T ap3_out = 0.8231455011548133;
  const T xh0_out[3] = {0.2966624523567004, 0.5503675896780905, 1.1130824767666705};
  const T yh0_out[3] = {0.4879907131826658, 0.5998877744675208, 0.5674900513272312};
  const T xh1_out[3] = {0.7728426941389869, 0.806513661919549, 0.7680307981859807};
  const T yh1_out[3] = {0.7646267390116116, 0.8486264909139831, 0.4702306668421768};
  const T xh2_out[3] = {0.4332361206658768, 0.39152830372174297, 0.780403246978787};
  const T yh2_out[3] = {0.10037891703979668, 0.4128471206575793, 0.4980066569492231};
  const T xh3_out[3] = {0.17824033509683085, 0.3887851668772747, 0.7317042161204349};
  const T yh3_out[3] = {0.889748075251454, 0.6168803764093128, 0.24375276299379783};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Dot_A2DVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t xb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(0, 0);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Dot(x_a2d, y_a2d, a_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Dot_A2DVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFR_01*/
T ab(ab_data);  /*UNQ_T2F_TFR_01*/
Vec_t xb, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y_a2d, a_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Dot_A2DVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data);  /*UNQ_T2F_TFHF_01*/
T ab(ab_data);  /*UNQ_T2F_TFHF_01*/
Vec_t xb, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y_a2d, a_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHF_11*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Dot_A2DVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data);  /*UNQ_T2F_TFHR_01*/
T ab(ab_data), ah0(ah0_data), ah1(ah1_data), ah2(ah2_data), ah3(ah3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t xb, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(0, ab);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.hvalue[0] = ah0;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[1] = ah1;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[2] = ah2;  /*UNQ_T2F_TFHR_05*/
a_a2d.hvalue[3] = ah3;  /*UNQ_T2F_TFHR_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
}
// Evaluations:
auto expr = A2D::Vec3Dot(x_a2d, y_a2d, a_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_val_eq(a_a2d.value, a_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_val_eq(a_a2d.pvalue[0], ap0_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[1], ap1_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[2], ap2_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.pvalue[3], ap3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Cross : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Cross_A2DVecVec : public Vec3Cross {
 protected:
  // Results
  const T v_out[3] = {0.15430706006379918, -0.08160724673272442, 0.028345686603881898};
  const T xb_out[3] = {0.12979127513440708, 0.012337439894056196, -0.038210616460415986};
  const T vp0_out[3] = {0.04042444569001446, -0.287520119972046, 0.2113647706528934};
  const T vp1_out[3] = {0.4451047907392335, -0.6246729701451286, 0.38005730959156864};
  const T vp2_out[3] = {-0.02273513983067892, 0.0376511767393071, -0.02381416922183202};
  const T vp3_out[3] = {0.3347738618266205, -0.7657799944317787, 0.5126299790424196};
  const T xh0_out[3] = {-0.12634505417113728, 0.43691501585905823, -0.30680661684706595};
  const T xh1_out[3] = {-0.5618876594078178, 0.4443224638410546, -0.21598360242235082};
  const T xh2_out[3] = {-0.5576810642789405, 0.4972784896224274, -0.25749482333038753};
  const T xh3_out[3] = {-0.5075694162740358, 0.559973772072363, -0.31663981680807696};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Cross_A2DVecVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Cross(x_a2d, y, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Cross_A2DVecVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Cross_A2DVecVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Cross_A2DVecVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Cross_VecA2DVec : public Vec3Cross {
 protected:
  // Results
  const T v_out[3] = {0.15430706006379918, -0.08160724673272442, 0.028345686603881898};
  const T yb_out[3] = {-0.19586534488676866, 0.0758764748270142, 0.004567151093861771};
  const T vp0_out[3] = {0.23784345708315116, -0.056217390591258375, -0.05810895921457787};
  const T vp1_out[3] = {0.07759757493484024, 0.12037719741424403, -0.22194391092979945};
  const T vp2_out[3] = {0.1955636542245905, 0.018955043297485555, -0.14315534713463596};
  const T vp3_out[3] = {0.14433632470972851, -0.040885033093691646, -0.025358252891220324};
  const T yh0_out[3] = {-0.030631538543512957, -0.11043921461557712, 0.17968312367974684};
  const T yh1_out[3] = {0.18297746410449442, -0.14712540134250543, 0.10729709192739449};
  const T yh2_out[3] = {0.1263612052596281, -0.13606487446069238, 0.12452632608553868};
  const T yh3_out[3] = {0.161457029215053, -0.18514726226956668, 0.17563521597250123};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Cross_VecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Cross(x, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Cross_VecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Cross(x, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Cross_VecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Cross_VecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Cross_A2DVecA2DVec : public Vec3Cross {
 protected:
  // Results
  const T v_out[3] = {0.15430706006379918, -0.08160724673272442, 0.028345686603881898};
  const T xb_out[3] = {0.12979127513440708, 0.012337439894056196, -0.038210616460415986};
  const T yb_out[3] = {-0.19586534488676866, 0.0758764748270142, 0.004567151093861771};
  const T vp0_out[3] = {0.27826790277316565, -0.34373751056330437, 0.15325581143831557};
  const T vp1_out[3] = {0.5227023656740737, -0.5042957727308846, 0.15811339866176916};
  const T vp2_out[3] = {0.17282851439391153, 0.05660622003679268, -0.16696951635646798};
  const T vp3_out[3] = {0.47911018653634896, -0.8066650275254703, 0.4872717261511992};
  const T xh0_out[3] = {-0.29800087777037454, 0.4003770343131332, -0.24490901970692178};
  const T yh0_out[3] = {-0.14615692334108948, 0.16455396157854016, 0.05300729366057741};
  const T xh1_out[3] = {-0.22881639469686843, -0.11622226225171678, 0.01871485824949537};
  const T yh1_out[3] = {-0.33728258848304216, 0.45806914615297645, -0.10737919681556012};
  const T xh2_out[3] = {-0.6468543981727677, 0.260143870331332, -0.10276116290468057};
  const T yh2_out[3] = {0.06722565237904865, -0.18076210638880902, 0.16389234340296338};
  const T xh3_out[3] = {-0.5973561045091121, 0.5516455196414116, -0.29032260765506274};
  const T yh3_out[3] = {-0.1909965067615801, 0.5642021742836232, -0.16048001129037362};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Cross_A2DVecA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb, yb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Cross(x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Cross_A2DVecA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Cross_A2DVecA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Cross_A2DVecA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), y(y_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data), xp3(xp3_data), yp3(yp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb, yb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), y_a2d(y, yb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
y_a2d.pvalue(3)(ii_0) = yp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Cross(x_a2d, y_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(y_a2d.hvalue(3), yh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3Normalize : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3Normalize_A2DVec : public Vec3Normalize {
 protected:
  // Results
  const T v_out[3] = {0.31593906513657255, 0.7833338836087105, 0.5353228314868739};
  const T xb_out[3] = {-0.14012396262140797, -0.40210667629807983, 0.6710986287014249};
  const T vp0_out[3] = {0.45599102452998735, -0.297110454527489, 0.1656408114684889};
  const T vp1_out[3] = {0.9332214818478035, -0.13695805387471324, -0.3503628600824785};
  const T vp2_out[3] = {-0.22477313001766144, -0.17355086628940747, 0.3866132256770997};
  const T vp3_out[3] = {1.4206955603287357, -0.19672909122861396, -0.5505998376799011};
  const T xh0_out[3] = {0.16306451459961013, 0.6346816874496347, -1.3364314569496054};
  const T xh1_out[3] = {-0.43803142967315456, 1.7251015402136651, -1.6851817602560741};
  const T xh2_out[3] = {1.0493281022829306, 0.5517787268521277, -2.1005799050914584};
  const T xh3_out[3] = {-1.0246051088063082, 1.4478517971182696, -0.59957163504673};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3Normalize_A2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
Vec_t v, vb, xb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3Normalize(x_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3Normalize_A2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3Normalize(x_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3Normalize_A2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3Normalize(x_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3Normalize_A2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data), vb(vb_data), vh0(vh0_data), vh1(vh1_data), vh2(vh2_data), vh3(vh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t v, xb;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/
v_a2d.hvalue(3)(ii_0) = vh3(ii_0);  /*UNQ_T2F_TFHR_05*/
}
// Evaluations:
auto expr = A2D::Vec3Normalize(x_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(v_a2d.pvalue(3), vp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3ScaleSymmetricOuterProduct : public vecops3d_a2dTest {
 protected:

};  /*UNQ_TC_TF_01*/

class Vec3ScaleSymmetricOuterProduct_A2DScalarVec : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {0.0048563758607652946, 0.012040814774306012, 0.008228576847328403,
                      0.012040814774306012, 0.029853789036480997, 0.020401791894915408,
                      0.008228576847328403, 0.020401791894915408, 0.013942388083964948};
  const T ab_out = 0.6670299107872728;
  const T Sp0_out[9] = {0.01479925872488173, 0.036693027519342815, 0.02507566160292923,
                        0.036693027519342815, 0.0909760612720154, 0.0621721640500182,
                        0.02507566160292923, 0.0621721640500182, 0.04248785810923407};
  const T Sp1_out[9] = {0.022774819025886118, 0.056467494541458524, 0.038589341910815096,
                        0.056467494541458524, 0.14000453466459928, 0.09567775055558962,
                        0.038589341910815096, 0.09567775055558962, 0.0653852532227466};
  const T Sp2_out[9] = {0.01693169131640428, 0.04198014420663525, 0.02868882622489607,
                        0.04198014420663525, 0.10408484744240842, 0.07113058226341455,
                        0.02868882622489607, 0.07113058226341455, 0.048609954834510444};
  const T Sp3_out[9] = {0.01875263462923227, 0.04649495973428221, 0.031774207672672375,
                        0.04649495973428221, 0.11527880340198506, 0.07878042395315975,
                        0.031774207672672375, 0.07878042395315975, 0.05383778296689637};
  const T ah0_out = 0.6133980634416792;
  const T ah1_out = 0.8017179760897222;
  const T ah2_out = 0.41337558810664504;
  const T ah3_out = 0.5426581694327351;
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Mat_t S, Sb;  /*UNQ_T2F_TFP_02*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x, S_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFR_01*/
Mat_t S;  /*UNQ_T2F_TFR_02*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x, S_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFHF_01*/
Mat_t S;  /*UNQ_T2F_TFHF_02*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x, S_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Mat_t Sb(Sb_data), Sh0(Sh0_data), Sh1(Sh1_data), Sh2(Sh2_data), Sh3(Sh3_data);  /*UNQ_T2F_TFHR_01*/
Mat_t S;  /*UNQ_T2F_TFHR_02*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
S_a2d.hvalue(0)(ii_0, ii_1) = Sh0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(1)(ii_0, ii_1) = Sh1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(2)(ii_0, ii_1) = Sh2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(3)(ii_0, ii_1) = Sh3(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
}
}
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x, S_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3ScaleSymmetricOuterProduct_ScalarA2DVec : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {0.0048563758607652946, 0.012040814774306012, 0.008228576847328403,
                      0.012040814774306012, 0.029853789036480997, 0.020401791894915408,
                      0.008228576847328403, 0.020401791894915408, 0.013942388083964948};
  const T xb_out[3] = {0.27446302144407514, 0.25680687469736846, 0.3543319160917333};
  const T Sp0_out[9] = {0.02760080962380041, 0.04648765382194, 0.03743633414046582,
                        0.04648765382194, 0.060849843309388515, 0.055635154858396596,
                        0.03743633414046582, 0.055635154858396596, 0.047622818360494125};
  const T Sp1_out[9] = {0.05328283899548347, 0.09443715481662533, 0.06059061362757036,
                        0.09443715481662533, 0.1407440368866415, 0.08639742076728603,
                        0.06059061362757036, 0.08639742076728603, 0.052355713132956186};
  const T Sp2_out[9] = {0.003888251566648131, 0.01553914386987501, 0.01838510203186234,
                        0.01553914386987501, 0.053152571297222055, 0.055578348152600036,
                        0.01838510203186234, 0.055578348152600036, 0.05113997446665416};
  const T Sp3_out[9] = {0.05888608221960696, 0.08883276577749036, 0.05431060545196671,
                        0.08883276577749036, 0.07850814936779622, 0.037791470553131416,
                        0.05431060545196671, 0.037791470553131416, 0.014987589966190129};
  const T xh0_out[3] = {0.8433716742130464, 0.6657362512558704, 0.7968859386481212};
  const T xh1_out[3] = {1.0589774467353867, 1.013007216861836, 1.404626552475667};
  const T xh2_out[3] = {0.5853091155632785, 0.4349905939933781, 0.6201081013703943};
  const T xh3_out[3] = {0.5101656369860984, 0.7061178677656276, 1.0580384711890098};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t xb;  /*UNQ_T2F_TFP_02*/
Mat_t S, Sb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3ScaleSymmetricOuterProduct(a, x_a2d, S_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFR_01*/
Vec_t xb;  /*UNQ_T2F_TFR_02*/
Mat_t S;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a, x_a2d, S_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFR_06*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data);  /*UNQ_T2F_TFHF_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t xb;  /*UNQ_T2F_TFHF_02*/
Mat_t S;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a, x_a2d, S_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHF_09*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_ScalarA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data);  /*UNQ_T2F_TFHR_01*/
Mat_t Sb(Sb_data), Sh0(Sh0_data), Sh1(Sh1_data), Sh2(Sh2_data), Sh3(Sh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t xb;  /*UNQ_T2F_TFHR_02*/
Mat_t S;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
S_a2d.hvalue(0)(ii_0, ii_1) = Sh0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(1)(ii_0, ii_1) = Sh1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(2)(ii_0, ii_1) = Sh2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(3)(ii_0, ii_1) = Sh3(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
}
}
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a, x_a2d, S_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHR_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHR_12*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Vec3ScaleSymmetricOuterProduct_A2DScalarA2DVec : public Vec3ScaleSymmetricOuterProduct {
 protected:
  // Results
  const T S_out[9] = {0.0048563758607652946, 0.012040814774306012, 0.008228576847328403,
                      0.012040814774306012, 0.029853789036480997, 0.020401791894915408,
                      0.008228576847328403, 0.020401791894915408, 0.013942388083964948};
  const T ab_out = 0.6670299107872728;
  const T xb_out[3] = {0.27446302144407514, 0.25680687469736846, 0.3543319160917333};
  const T Sp0_out[9] = {0.042400068348682145, 0.08318068134128281, 0.06251199574339505,
                        0.08318068134128281, 0.1518259045814039, 0.11780731890841478,
                        0.06251199574339505, 0.11780731890841478, 0.09011067646972822};
  const T Sp1_out[9] = {0.07605765802136959, 0.15090464935808387, 0.09917995553838545,
                        0.15090464935808387, 0.2807485715512408, 0.18207517132287568,
                        0.09917995553838545, 0.18207517132287568, 0.11774096635570279};
  const T Sp2_out[9] = {0.02081994288305241, 0.05751928807651026, 0.047073928256758414,
                        0.05751928807651026, 0.15723741873963049, 0.12670893041601458,
                        0.047073928256758414, 0.12670893041601458, 0.09974992930116461};
  const T Sp3_out[9] = {0.07763871684883923, 0.13532772551177258, 0.08608481312463909,
                        0.13532772551177258, 0.19378695276978128, 0.11657189450629116,
                        0.08608481312463909, 0.11657189450629116, 0.06882537293308649};
  const T ah0_out = 2.779398842395766;
  const T xh0_out[3] = {1.6797668341482408, 1.4483262916363733, 1.876672563709811};
  const T ah1_out = 4.4501009551138155;
  const T xh1_out[3] = {2.3461195094405145, 2.217347709346151, 3.0663276842361604};
  const T ah2_out = 1.9823813351613593;
  const T xh2_out[3] = {1.5422208724168618, 1.3303444247049852, 1.8554817212928745};
  const T ah3_out = 3.0349568295130216;
  const T xh3_out[3] = {1.5699898573406443, 1.6977638254822685, 2.4262721397625353};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarA2DVec, passive) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
T a(a_data);  /*UNQ_T2F_TFP_01*/
Vec_t xb;  /*UNQ_T2F_TFP_02*/
Mat_t S, Sb;  /*UNQ_T2F_TFP_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFP_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFP_03*/
// Set Derivative Values:
/*None for "passive" tests*/
// Evaluations:
A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x_a2d, S_a2d);  /*UNQ_T2F_TFP_04*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarA2DVec, reverse) {
// Declarations and Initializations:
Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
T a(a_data);  /*UNQ_T2F_TFR_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFR_01*/
Vec_t xb;  /*UNQ_T2F_TFR_02*/
Mat_t S;  /*UNQ_T2F_TFR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFR_03*/
// Set Derivative Values:
/*None for "reverse" tests*/
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x_a2d, S_a2d);  /*UNQ_T2F_TFR_04*/
expr.reverse();  /*UNQ_T2F_TFR_05*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFR_06*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFR_07*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarA2DVec, hforward) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHF_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHF_01*/
Mat_t Sb(Sb_data);  /*UNQ_T2F_TFHF_01*/
Vec_t xb;  /*UNQ_T2F_TFHF_02*/
Mat_t S;  /*UNQ_T2F_TFHF_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHF_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHF_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHF_05*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHF_05*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHF_05*/
}
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x_a2d, S_a2d);  /*UNQ_T2F_TFHF_06*/
expr.reverse();  /*UNQ_T2F_TFHF_07*/
expr.hforward();  /*UNQ_T2F_TFHF_08*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHF_09*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHF_10*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHF_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Vec3ScaleSymmetricOuterProduct_A2DScalarA2DVec, hreverse) {
// Declarations and Initializations:
Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data), xp3(xp3_data);  /*UNQ_T2F_TFHR_01*/
T a(a_data), ap0(ap0_data), ap1(ap1_data), ap2(ap2_data), ap3(ap3_data);  /*UNQ_T2F_TFHR_01*/
Mat_t Sb(Sb_data), Sh0(Sh0_data), Sh1(Sh1_data), Sh2(Sh2_data), Sh3(Sh3_data);  /*UNQ_T2F_TFHR_01*/
Vec_t xb;  /*UNQ_T2F_TFHR_02*/
Mat_t S;  /*UNQ_T2F_TFHR_02*/
A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
A2DScalar_t a_a2d(a, 0);  /*UNQ_T2F_TFHR_03*/
A2DMat_t S_a2d(S, Sb);  /*UNQ_T2F_TFHR_03*/
// Set Derivative Values:
/*UNQ_T2F_CDL_01*/
a_a2d.pvalue[0] = ap0;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[1] = ap1;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[2] = ap2;  /*UNQ_T2F_TFHR_04*/
a_a2d.pvalue[3] = ap3;  /*UNQ_T2F_TFHR_04*/
for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/
x_a2d.pvalue(3)(ii_0) = xp3(ii_0);  /*UNQ_T2F_TFHR_04*/
for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
S_a2d.hvalue(0)(ii_0, ii_1) = Sh0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(1)(ii_0, ii_1) = Sh1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(2)(ii_0, ii_1) = Sh2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
S_a2d.hvalue(3)(ii_0, ii_1) = Sh3(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/
}
}
// Evaluations:
auto expr = A2D::Vec3ScaleSymmetricOuterProduct(a_a2d, x_a2d, S_a2d);  /*UNQ_T2F_TFHR_06*/
expr.reverse();  /*UNQ_T2F_TFHR_07*/
expr.hforward();  /*UNQ_T2F_TFHR_08*/
expr.hreverse();  /*UNQ_T2F_TFHR_09*/
// Comparisons:
expect_mat_eq<3, 3>(S_a2d.value(), S_out);  /*UNQ_T2F_TFHR_10*/
expect_val_eq(a_a2d.bvalue, ab_out);  /*UNQ_T2F_TFHR_11*/
expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
expect_mat_eq<3, 3>(S_a2d.pvalue(0), Sp0_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(1), Sp1_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(2), Sp2_out);  /*UNQ_T2F_TFHR_12*/
expect_mat_eq<3, 3>(S_a2d.pvalue(3), Sp3_out);  /*UNQ_T2F_TFHR_12*/
expect_val_eq(a_a2d.hvalue[0], ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[1], ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[2], ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_val_eq(a_a2d.hvalue[3], ah3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
expect_vec_eq<3>(x_a2d.hvalue(3), xh3_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}
