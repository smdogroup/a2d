/*
    This is a set of automatically generated unit tests for a2dmatvecops3d.h using
    Google Test framework.  These tests were written on 2022-10-23 using the
    A2DTestConstructor package.
*/

#include <gtest/gtest.h>

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2dmatvecops3d.h"
#include "test_commons.h"

using T = double;  /*UNQ_TC_TD_01*/
using Vec_t = A2D::Vec<T, 3>;  /*UNQ_TC_TD_01*/
using A2DVec_t = A2D::A2DVec<3, Vec_t>;  /*UNQ_TC_TD_01*/
using Mat3x2_t = A2D::Mat<T, 3, 2>;  /*UNQ_TC_TD_01*/
using A2DMat3x2_t = A2D::A2DMat<3, Mat3x2_t>;  /*UNQ_TC_TD_01*/
using Mat3x3_t = A2D::Mat<T, 3, 3>;  /*UNQ_TC_TD_01*/
using ADMat3x3_t = A2D::ADMat<Mat3x3_t>;  /*UNQ_TC_TD_01*/
using A2DMat3x3_t = A2D::A2DMat<3, Mat3x3_t>;  /*UNQ_TC_TD_01*/
using ADMat3x2_t = A2D::ADMat<Mat3x2_t>;  /*UNQ_TC_TD_01*/
using ADVec_t = A2D::ADVec<Vec_t>;  /*UNQ_TC_TD_01*/

class matvecops3d_a2dTest : public ::testing::Test {
protected:
    // Input Options:
    const T A_data[6] = {0.4943894693718751, 0.6592561991069723, 
0.9465887250044129, 0.8558301425103879, 
0.6524425854579287, 0.9094094296585536};
    const T x_data[3] = {0.8171547300315154, 0.5671462658111379, 0.6326462008557912};
    const T y_data[3] = {0.9147027276650903, 0.1151828695760968, 0.9238750716313718};
    const T z_data[3] = {0.45879726281241784, 0.6708573392200072, 0.3268990300170723};
    const T Ab_data[6] = {0.3708978579459866, 0.7384936497493327, 
0.3921723853501031, 0.8573007463008531, 
0.5870746787450382, 0.8756208865696513};
    const T ub_data[3] = {0.13637823233780655, 0.46749752767005703, 0.8406096956128962};
    const T vb_data[3] = {0.40175549338424354, 0.5319095796485651, 0.29662884830455416};
    const T xb_data[3] = {0.47956586874092144, 0.4432962905192246, 0.6468568419678706};
    const T yb_data[3] = {0.6051101769012076, 0.6487517388085736, 0.010702071350234266};
    const T zb_data[3] = {0.6119794029525708, 0.8364282494515285, 0.5881527410796049};
    const T Cb_data[9] = {0.6949807472239024, 0.12572152304416717, 0.4188711225520372, 
0.01930685032142776, 0.44132183521619983, 0.2559938992473544, 
0.37224854848033617, 0.8083221322969701, 0.5054872906085297};
    const T Ap0_data[6] = {0.17859684887146865, 0.004938652623771578, 
0.6115790613165083, 0.5626424937353428, 
0.09313402906411528, 0.9435523555925251};
    const T Ap1_data[6] = {0.11449864684271993, 0.37702051702668193, 
0.6875873814395518, 0.8871704439192243, 
0.9565664927098708, 0.6297019578393069};
    const T Ap2_data[6] = {0.5144561070381665, 0.2135760925007788, 
0.6235062784392448, 0.5315930244829316, 
0.1896970725672642, 0.27467004390750815};
    const T xp0_data[3] = {0.5632560445525615, 0.9297696497866957, 0.6639191241480596};
    const T xp1_data[3] = {0.4134336741098551, 0.4417836955915003, 0.2575477044523653};
    const T xp2_data[3] = {0.4841116819654573, 0.41105190129938807, 0.04560950750544501};
    const T yp0_data[3] = {0.308089744949803, 0.8387594032251646, 0.6928243590911211};
    const T yp1_data[3] = {0.4787004182914746, 0.20773010962853278, 0.10660943369769882};
    const T yp2_data[3] = {0.7490976729536114, 0.1850699005947951, 0.20287688570239404};
    const T zp0_data[3] = {0.6338146066065109, 0.8171412122566124, 0.8209797749182447};
    const T zp1_data[3] = {0.22026183805879374, 0.2406212095525121, 0.43882859463385004};
    const T zp2_data[3] = {0.8619937193951389, 0.4250509908063892, 0.9828076665854727};
    const T uh0_data[3] = {0.716513597430403, 0.16067204946487545, 0.9997036542423365};
    const T uh1_data[3] = {0.09431665454471438, 0.6002617692039943, 0.5558037019152263};
    const T uh2_data[3] = {0.6840470174023056, 0.17184068788708617, 0.4491025946840186};
    const T vh0_data[3] = {0.9127791468358138, 0.14747960638043, 0.725439310555214};
    const T vh1_data[3] = {0.4724017927886681, 0.17222000818992, 0.17193289543974233};
    const T vh2_data[3] = {0.39626847163645296, 0.2691361989401383, 0.5302218848378587};
    const T Ch0_data[9] = {0.6622197827410925, 0.12636386010177514, 0.6872364713794958, 
0.8572576025788212, 0.37464379013264515, 0.0582614146530781, 
0.37965962628053496, 0.4516719304600625, 0.14777763366502583};
    const T Ch1_data[9] = {0.709578232730824, 0.7337305800306985, 0.505805922419502, 
0.1530966081596158, 0.09761832055278763, 0.06945786974875079, 
0.8690139240414835, 0.5018004819945504, 0.3773359083397829};
    const T Ch2_data[9] = {0.03899474177603057, 0.6494088076488078, 0.7120318674774812, 
0.7725274677518728, 0.5610718509893741, 0.48218065749992944, 
0.961125836309164, 0.632511767131526, 0.601465326088105};
};  /*UNQ_TC_TIC_01*/

class Mat3x2ToVec3 : public matvecops3d_a2dTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x2ToVec3_A2DMat3x2 : public Mat3x2ToVec3 {
protected:
    // Results
    const T u_out[3] = {0.4943894693718751, 0.9465887250044129, 0.6524425854579287};
    const T v_out[3] = {0.6592561991069723, 0.8558301425103879, 0.9094094296585536};
    const T Ab_out[6] = {0.13637823233780655, 0.40175549338424354, 
0.46749752767005703, 0.5319095796485651, 
0.8406096956128962, 0.29662884830455416};
    const T up0_out[3] = {0.17859684887146862, 0.6115790613165083, 0.09313402906411528};
    const T vp0_out[3] = {0.004938652623771578, 0.5626424937353428, 0.9435523555925251};
    const T up1_out[3] = {0.11449864684271993, 0.6875873814395518, 0.9565664927098708};
    const T vp1_out[3] = {0.3770205170266819, 0.8871704439192243, 0.6297019578393069};
    const T up2_out[3] = {0.5144561070381665, 0.6235062784392448, 0.18969707256726423};
    const T vp2_out[3] = {0.2135760925007788, 0.5315930244829316, 0.27467004390750815};
    const T Ah0_out[6] = {0.716513597430403, 0.9127791468358138, 
0.16067204946487545, 0.14747960638043, 
0.9997036542423365, 0.725439310555214};
    const T Ah1_out[6] = {0.09431665454471438, 0.4724017927886681, 
0.6002617692039943, 0.17222000818992, 
0.5558037019152263, 0.17193289543974233};
    const T Ah2_out[6] = {0.6840470174023056, 0.39626847163645296, 
0.17184068788708617, 0.2691361989401383, 
0.4491025946840186, 0.5302218848378587};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x2ToVec3_A2DMat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFP_01*/
    Vec_t u, v, ub, vb;  /*UNQ_T2F_TFP_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFP_02*/
    A2DVec_t u_a2d(u, ub), v_a2d(v, vb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x2ToVec3(A_a2d, u_a2d, v_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_vec_eq<3>(u_a2d.value(), u_out);  /*UNQ_T2F_TFP_05*/
    expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x2ToVec3_A2DMat3x2, reverse) {
    // Declarations and Initializations:
    Vec_t ub(ub_data), vb(vb_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFR_01*/
    Vec_t u, v;  /*UNQ_T2F_TFR_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFR_02*/
    A2DVec_t u_a2d(u, ub), v_a2d(v, vb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x2ToVec3(A_a2d, u_a2d, v_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_vec_eq<3>(u_a2d.value(), u_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFR_06*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x2ToVec3_A2DMat3x2, hforward) {
    // Declarations and Initializations:
    Vec_t ub(ub_data), vb(vb_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t u, v;  /*UNQ_T2F_TFHF_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t u_a2d(u, ub), v_a2d(v, vb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        for (int ii_1 = 0; ii_1 < 2; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x2ToVec3(A_a2d, u_a2d, v_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_vec_eq<3>(u_a2d.value(), u_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHF_09*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(u_a2d.pvalue(0), up0_out);  /*UNQ_T2F_TFHF_11*/
    expect_vec_eq<3>(u_a2d.pvalue(1), up1_out);  /*UNQ_T2F_TFHF_11*/
    expect_vec_eq<3>(u_a2d.pvalue(2), up2_out);  /*UNQ_T2F_TFHF_11*/
    expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x2ToVec3_A2DMat3x2, hreverse) {
    // Declarations and Initializations:
    Vec_t ub(ub_data), vb(vb_data), uh0(uh0_data), vh0(vh0_data), uh1(uh1_data), vh1(vh1_data), uh2(uh2_data), vh2(vh2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t u, v;  /*UNQ_T2F_TFHR_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t u_a2d(u, ub), v_a2d(v, vb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {u_a2d.hvalue(0)(ii_0) = uh0(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        if (ii_0 < 3) {v_a2d.hvalue(0)(ii_0) = vh0(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        if (ii_0 < 3) {u_a2d.hvalue(1)(ii_0) = uh1(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        if (ii_0 < 3) {v_a2d.hvalue(1)(ii_0) = vh1(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        if (ii_0 < 3) {u_a2d.hvalue(2)(ii_0) = uh2(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        if (ii_0 < 3) {v_a2d.hvalue(2)(ii_0) = vh2(ii_0);  /*UNQ_T2F_TFHR_05*/ }
        for (int ii_1 = 0; ii_1 < 2; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x2ToVec3(A_a2d, u_a2d, v_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_vec_eq<3>(u_a2d.value(), u_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(v_a2d.value(), v_out);  /*UNQ_T2F_TFHR_10*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(u_a2d.pvalue(0), up0_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(u_a2d.pvalue(1), up1_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(u_a2d.pvalue(2), up2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(v_a2d.pvalue(0), vp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(v_a2d.pvalue(1), vp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(v_a2d.pvalue(2), vp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(0), Ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(1), Ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(2), Ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3 : public matvecops3d_a2dTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromThreeVec3_VecA2DVecA2DVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T yb_out[3] = {0.12572152304416717, 0.44132183521619983, 0.8083221322969701};
    const T zb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.0, 0.308089744949803, 0.6338146066065109, 
0.0, 0.8387594032251646, 0.8171412122566124, 
0.0, 0.6928243590911211, 0.8209797749182446};
    const T Cp1_out[9] = {0.0, 0.4787004182914746, 0.22026183805879374, 
0.0, 0.20773010962853278, 0.2406212095525121, 
0.0, 0.10660943369769882, 0.4388285946338501};
    const T Cp2_out[9] = {0.0, 0.7490976729536114, 0.8619937193951389, 
0.0, 0.1850699005947951, 0.4250509908063892, 
0.0, 0.20287688570239404, 0.9828076665854727};
    const T yh0_out[3] = {0.12636386010177514, 0.37464379013264515, 0.4516719304600625};
    const T zh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T yh1_out[3] = {0.7337305800306985, 0.09761832055278763, 0.5018004819945504};
    const T zh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T yh2_out[3] = {0.6494088076488078, 0.5610718509893741, 0.632511767131526};
    const T zh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecA2DVecA2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t yb, zb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecA2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t yb, zb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecA2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yp0(yp0_data), zp0(zp0_data), yp1(yp1_data), zp1(zp1_data), yp2(yp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t yb, zb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecA2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yp0(yp0_data), zp0(zp0_data), yp1(yp1_data), zp1(zp1_data), yp2(yp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t yb, zb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(0), zh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(1), zh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(2), zh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_A2DVecVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T xb_out[3] = {0.6949807472239024, 0.01930685032142776, 0.37224854848033617};
    const T Cp0_out[9] = {0.5632560445525615, 0.0, 0.0, 
0.9297696497866957, 0.0, 0.0, 
0.6639191241480596, 0.0, 0.0};
    const T Cp1_out[9] = {0.4134336741098551, 0.0, 0.0, 
0.4417836955915003, 0.0, 0.0, 
0.2575477044523653, 0.0, 0.0};
    const T Cp2_out[9] = {0.48411168196545734, 0.0, 0.0, 
0.41105190129938807, 0.0, 0.0, 
0.04560950750544501, 0.0, 0.0};
    const T xh0_out[3] = {0.6622197827410925, 0.8572576025788212, 0.37965962628053496};
    const T xh1_out[3] = {0.709578232730824, 0.1530966081596158, 0.8690139240414835};
    const T xh2_out[3] = {0.03899474177603057, 0.7725274677518728, 0.961125836309164};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_A2DVecVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_a2d, y, z, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_VecVecA2DVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T zb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.0, 0.0, 0.6338146066065109, 
0.0, 0.0, 0.8171412122566124, 
0.0, 0.0, 0.8209797749182446};
    const T Cp1_out[9] = {0.0, 0.0, 0.22026183805879374, 
0.0, 0.0, 0.2406212095525121, 
0.0, 0.0, 0.4388285946338501};
    const T Cp2_out[9] = {0.0, 0.0, 0.8619937193951389, 
0.0, 0.0, 0.4250509908063892, 
0.0, 0.0, 0.9828076665854727};
    const T zh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T zh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T zh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecVecA2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t zb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t z_a2d(z, zb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y, z_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecVecA2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t zb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t z_a2d(z, zb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y, z_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecVecA2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), zp0(zp0_data), zp1(zp1_data), zp2(zp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t zb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t z_a2d(z, zb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y, z_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_VecVecA2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), zp0(zp0_data), zp1(zp1_data), zp2(zp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t zb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t z_a2d(z, zb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y, z_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(z_a2d.hvalue(0), zh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(1), zh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(2), zh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_A2DVecVecA2DVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T xb_out[3] = {0.6949807472239024, 0.01930685032142776, 0.37224854848033617};
    const T zb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.5632560445525615, 0.0, 0.6338146066065109, 
0.9297696497866957, 0.0, 0.8171412122566124, 
0.6639191241480596, 0.0, 0.8209797749182446};
    const T Cp1_out[9] = {0.4134336741098551, 0.0, 0.22026183805879374, 
0.4417836955915003, 0.0, 0.2406212095525121, 
0.2575477044523653, 0.0, 0.4388285946338501};
    const T Cp2_out[9] = {0.48411168196545734, 0.0, 0.8619937193951389, 
0.41105190129938807, 0.0, 0.4250509908063892, 
0.04560950750544501, 0.0, 0.9828076665854727};
    const T xh0_out[3] = {0.6622197827410925, 0.8572576025788212, 0.37965962628053496};
    const T zh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T xh1_out[3] = {0.709578232730824, 0.1530966081596158, 0.8690139240414835};
    const T zh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T xh2_out[3] = {0.03899474177603057, 0.7725274677518728, 0.961125836309164};
    const T zh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_A2DVecVecA2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb, zb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb), z_a2d(z, zb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_a2d, y, z_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecA2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb, zb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb), z_a2d(z, zb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecA2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), zp0(zp0_data), xp1(xp1_data), zp1(zp1_data), xp2(xp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb, zb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb), z_a2d(z, zb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecVecA2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), zp0(zp0_data), xp1(xp1_data), zp1(zp1_data), xp2(xp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb, zb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb), z_a2d(z, zb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y, z_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(0), zh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(1), zh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(2), zh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_VecA2DVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T yb_out[3] = {0.12572152304416717, 0.44132183521619983, 0.8083221322969701};
    const T Cp0_out[9] = {0.0, 0.308089744949803, 0.0, 
0.0, 0.8387594032251646, 0.0, 
0.0, 0.6928243590911211, 0.0};
    const T Cp1_out[9] = {0.0, 0.4787004182914746, 0.0, 
0.0, 0.20773010962853278, 0.0, 
0.0, 0.10660943369769882, 0.0};
    const T Cp2_out[9] = {0.0, 0.7490976729536114, 0.0, 
0.0, 0.1850699005947951, 0.0, 
0.0, 0.20287688570239404, 0.0};
    const T yh0_out[3] = {0.12636386010177514, 0.37464379013264515, 0.4516719304600625};
    const T yh1_out[3] = {0.7337305800306985, 0.09761832055278763, 0.5018004819945504};
    const T yh2_out[3] = {0.6494088076488078, 0.5610718509893741, 0.632511767131526};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecA2DVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t yb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y_a2d, z, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t yb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t yb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_VecA2DVecVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yp0(yp0_data), yp1(yp1_data), yp2(yp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t yb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t y_a2d(y, yb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_a2d, z, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_A2DVecA2DVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T xb_out[3] = {0.6949807472239024, 0.01930685032142776, 0.37224854848033617};
    const T yb_out[3] = {0.12572152304416717, 0.44132183521619983, 0.8083221322969701};
    const T Cp0_out[9] = {0.5632560445525615, 0.308089744949803, 0.0, 
0.9297696497866957, 0.8387594032251646, 0.0, 
0.6639191241480596, 0.6928243590911211, 0.0};
    const T Cp1_out[9] = {0.4134336741098551, 0.4787004182914746, 0.0, 
0.4417836955915003, 0.20773010962853278, 0.0, 
0.2575477044523653, 0.10660943369769882, 0.0};
    const T Cp2_out[9] = {0.48411168196545734, 0.7490976729536114, 0.0, 
0.41105190129938807, 0.1850699005947951, 0.0, 
0.04560950750544501, 0.20287688570239404, 0.0};
    const T xh0_out[3] = {0.6622197827410925, 0.8572576025788212, 0.37965962628053496};
    const T yh0_out[3] = {0.12636386010177514, 0.37464379013264515, 0.4516719304600625};
    const T xh1_out[3] = {0.709578232730824, 0.1530966081596158, 0.8690139240414835};
    const T yh1_out[3] = {0.7337305800306985, 0.09761832055278763, 0.5018004819945504};
    const T xh2_out[3] = {0.03899474177603057, 0.7725274677518728, 0.961125836309164};
    const T yh2_out[3] = {0.6494088076488078, 0.5610718509893741, 0.632511767131526};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb, yb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb, yb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb, yb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), yp0(yp0_data), xp1(xp1_data), yp1(yp1_data), xp2(xp2_data), yp2(yp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb, yb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromThreeVec3_A2DVecA2DVecA2DVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.8171547300315154, 0.9147027276650903, 0.45879726281241784, 
0.5671462658111379, 0.1151828695760968, 0.6708573392200072, 
0.6326462008557912, 0.9238750716313718, 0.3268990300170723};
    const T xb_out[3] = {0.6949807472239024, 0.01930685032142776, 0.37224854848033617};
    const T yb_out[3] = {0.12572152304416717, 0.44132183521619983, 0.8083221322969701};
    const T zb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.5632560445525615, 0.308089744949803, 0.6338146066065109, 
0.9297696497866957, 0.8387594032251646, 0.8171412122566124, 
0.6639191241480596, 0.6928243590911211, 0.8209797749182446};
    const T Cp1_out[9] = {0.4134336741098551, 0.4787004182914746, 0.22026183805879374, 
0.4417836955915003, 0.20773010962853278, 0.2406212095525121, 
0.2575477044523653, 0.10660943369769882, 0.4388285946338501};
    const T Cp2_out[9] = {0.48411168196545734, 0.7490976729536114, 0.8619937193951389, 
0.41105190129938807, 0.1850699005947951, 0.4250509908063892, 
0.04560950750544501, 0.20287688570239404, 0.9828076665854727};
    const T xh0_out[3] = {0.6622197827410925, 0.8572576025788212, 0.37965962628053496};
    const T yh0_out[3] = {0.12636386010177514, 0.37464379013264515, 0.4516719304600625};
    const T zh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T xh1_out[3] = {0.709578232730824, 0.1530966081596158, 0.8690139240414835};
    const T yh1_out[3] = {0.7337305800306985, 0.09761832055278763, 0.5018004819945504};
    const T zh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T xh2_out[3] = {0.03899474177603057, 0.7725274677518728, 0.961125836309164};
    const T yh2_out[3] = {0.6494088076488078, 0.5610718509893741, 0.632511767131526};
    const T zh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecA2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb, yb, zb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecA2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb, yb, zb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecA2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), yp0(yp0_data), zp0(zp0_data), xp1(xp1_data), yp1(yp1_data), zp1(zp1_data), xp2(xp2_data), yp2(yp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb, yb, zb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromThreeVec3_A2DVecA2DVecA2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xp0(xp0_data), yp0(yp0_data), zp0(zp0_data), xp1(xp1_data), yp1(yp1_data), zp1(zp1_data), xp2(xp2_data), yp2(yp2_data), zp2(zp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb, yb, zb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb), y_a2d(y, yb), z_a2d(z, zb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(0)(ii_0) = yp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(0)(ii_0) = zp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(1)(ii_0) = yp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(1)(ii_0) = zp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {y_a2d.pvalue(2)(ii_0) = yp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {z_a2d.pvalue(2)(ii_0) = zp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_a2d, y_a2d, z_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(y_a2d.bvalue(), yb_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(z_a2d.bvalue(), zb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(0), yh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(1), yh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(y_a2d.hvalue(2), yh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(0), zh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(1), zh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(z_a2d.hvalue(2), zh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromMat3x2AndVec3 : public matvecops3d_a2dTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromMat3x2AndVec3_A2DMat3x2Vec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.4943894693718751, 0.6592561991069723, 0.8171547300315154, 
0.9465887250044129, 0.8558301425103879, 0.5671462658111379, 
0.6524425854579287, 0.9094094296585536, 0.6326462008557912};
    const T Ab_out[6] = {0.6949807472239024, 0.12572152304416717, 
0.01930685032142776, 0.44132183521619983, 
0.37224854848033617, 0.8083221322969701};
    const T Cp0_out[9] = {0.17859684887146862, 0.004938652623771578, 0.0, 
0.6115790613165083, 0.5626424937353428, 0.0, 
0.09313402906411528, 0.9435523555925251, 0.0};
    const T Cp1_out[9] = {0.11449864684271993, 0.3770205170266819, 0.0, 
0.6875873814395518, 0.8871704439192243, 0.0, 
0.9565664927098708, 0.6297019578393069, 0.0};
    const T Cp2_out[9] = {0.5144561070381665, 0.2135760925007788, 0.0, 
0.6235062784392448, 0.5315930244829316, 0.0, 
0.18969707256726423, 0.27467004390750815, 0.0};
    const T Ah0_out[6] = {0.6622197827410925, 0.12636386010177514, 
0.8572576025788212, 0.37464379013264515, 
0.37965962628053496, 0.4516719304600625};
    const T Ah1_out[6] = {0.709578232730824, 0.7337305800306985, 
0.1530966081596158, 0.09761832055278763, 
0.8690139240414835, 0.5018004819945504};
    const T Ah2_out[6] = {0.03899474177603057, 0.6494088076488078, 
0.7725274677518728, 0.5610718509893741, 
0.961125836309164, 0.632511767131526};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2Vec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFP_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2Vec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2Vec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        for (int ii_1 = 0; ii_1 < 2; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2Vec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(0), Ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(1), Ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(2), Ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromMat3x2AndVec3_Mat3x2A2DVec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.4943894693718751, 0.6592561991069723, 0.8171547300315154, 
0.9465887250044129, 0.8558301425103879, 0.5671462658111379, 
0.6524425854579287, 0.9094094296585536, 0.6326462008557912};
    const T xb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.0, 0.0, 0.5632560445525615, 
0.0, 0.0, 0.9297696497866957, 
0.0, 0.0, 0.6639191241480596};
    const T Cp1_out[9] = {0.0, 0.0, 0.4134336741098551, 
0.0, 0.0, 0.4417836955915003, 
0.0, 0.0, 0.2575477044523653};
    const T Cp2_out[9] = {0.0, 0.0, 0.48411168196545734, 
0.0, 0.0, 0.41105190129938807, 
0.0, 0.0, 0.04560950750544501};
    const T xh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T xh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T xh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2A2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A, x_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2A2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A, x_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2A2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A, x_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2A2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A, x_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromMat3x2AndVec3_A2DMat3x2A2DVec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.4943894693718751, 0.6592561991069723, 0.8171547300315154, 
0.9465887250044129, 0.8558301425103879, 0.5671462658111379, 
0.6524425854579287, 0.9094094296585536, 0.6326462008557912};
    const T Ab_out[6] = {0.6949807472239024, 0.12572152304416717, 
0.01930685032142776, 0.44132183521619983, 
0.37224854848033617, 0.8083221322969701};
    const T xb_out[3] = {0.4188711225520372, 0.2559938992473544, 0.5054872906085297};
    const T Cp0_out[9] = {0.17859684887146862, 0.004938652623771578, 0.5632560445525615, 
0.6115790613165083, 0.5626424937353428, 0.9297696497866957, 
0.09313402906411528, 0.9435523555925251, 0.6639191241480596};
    const T Cp1_out[9] = {0.11449864684271993, 0.3770205170266819, 0.4134336741098551, 
0.6875873814395518, 0.8871704439192243, 0.4417836955915003, 
0.9565664927098708, 0.6297019578393069, 0.2575477044523653};
    const T Cp2_out[9] = {0.5144561070381665, 0.2135760925007788, 0.48411168196545734, 
0.6235062784392448, 0.5315930244829316, 0.41105190129938807, 
0.18969707256726423, 0.27467004390750815, 0.04560950750544501};
    const T Ah0_out[6] = {0.6622197827410925, 0.12636386010177514, 
0.8572576025788212, 0.37464379013264515, 
0.37965962628053496, 0.4516719304600625};
    const T xh0_out[3] = {0.6872364713794958, 0.0582614146530781, 0.14777763366502583};
    const T Ah1_out[6] = {0.709578232730824, 0.7337305800306985, 
0.1530966081596158, 0.09761832055278763, 
0.8690139240414835, 0.5018004819945504};
    const T xh1_out[3] = {0.505805922419502, 0.06945786974875079, 0.3773359083397829};
    const T Ah2_out[6] = {0.03899474177603057, 0.6494088076488078, 
0.7725274677518728, 0.5610718509893741, 
0.961125836309164, 0.632511767131526};
    const T xh2_out[3] = {0.7120318674774812, 0.48218065749992944, 0.601465326088105};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2A2DVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFP_01*/
    Vec_t xb;  /*UNQ_T2F_TFP_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFP_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2A2DVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Vec_t xb;  /*UNQ_T2F_TFR_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFR_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFR_07*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2A2DVec, hforward) {
    // Declarations and Initializations:
    Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Vec_t xb;  /*UNQ_T2F_TFHF_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHF_05*/ }
        for (int ii_1 = 0; ii_1 < 2; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHF_10*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_A2DMat3x2A2DVec, hreverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), xp0(xp0_data), xp1(xp1_data), xp2(xp2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Vec_t xb;  /*UNQ_T2F_TFHR_02*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DVec_t x_a2d(x, xb);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        if (ii_0 < 3) {x_a2d.pvalue(0)(ii_0) = xp0(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(1)(ii_0) = xp1(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        if (ii_0 < 3) {x_a2d.pvalue(2)(ii_0) = xp2(ii_0);  /*UNQ_T2F_TFHR_04*/ }
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_a2d, x_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHR_11*/
    expect_vec_eq<3>(x_a2d.bvalue(), xb_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(0), Ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(1), Ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(2), Ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(0), xh0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(1), xh1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_vec_eq<3>(x_a2d.hvalue(2), xh2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}

class Mat3x3FromMat3x2 : public matvecops3d_a2dTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromMat3x2_A2DMat3x2 : public Mat3x3FromMat3x2 {
protected:
    // Results
    const T C_out[9] = {0.4943894693718751, 0.6592561991069723, 0.0, 
0.9465887250044129, 0.8558301425103879, 0.0, 
0.6524425854579287, 0.9094094296585536, 0.0};
    const T Ab_out[6] = {0.6949807472239024, 0.12572152304416717, 
0.01930685032142776, 0.44132183521619983, 
0.37224854848033617, 0.8083221322969701};
    const T Cp0_out[9] = {0.17859684887146862, 0.004938652623771578, 0.0, 
0.6115790613165083, 0.5626424937353428, 0.0, 
0.09313402906411528, 0.9435523555925251, 0.0};
    const T Cp1_out[9] = {0.11449864684271993, 0.3770205170266819, 0.0, 
0.6875873814395518, 0.8871704439192243, 0.0, 
0.9565664927098708, 0.6297019578393069, 0.0};
    const T Cp2_out[9] = {0.5144561070381665, 0.2135760925007788, 0.0, 
0.6235062784392448, 0.5315930244829316, 0.0, 
0.18969707256726423, 0.27467004390750815, 0.0};
    const T Ah0_out[6] = {0.6622197827410925, 0.12636386010177514, 
0.8572576025788212, 0.37464379013264515, 
0.37965962628053496, 0.4516719304600625};
    const T Ah1_out[6] = {0.709578232730824, 0.7337305800306985, 
0.1530966081596158, 0.09761832055278763, 
0.8690139240414835, 0.5018004819945504};
    const T Ah2_out[6] = {0.03899474177603057, 0.6494088076488078, 
0.7725274677518728, 0.5610718509893741, 
0.961125836309164, 0.632511767131526};
};  /*UNQ_T2F_FTV_01*/

TEST_F(Mat3x3FromMat3x2_A2DMat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFP_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T2F_TFP_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFP_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFP_03*/
    // Set Derivative Values:
        /*None for "passive" tests*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2(A_a2d, C_a2d);  /*UNQ_T2F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2_A2DMat3x2, reverse) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T2F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFR_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFR_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFR_03*/
    // Set Derivative Values:
        /*None for "reverse" tests*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2(A_a2d, C_a2d);  /*UNQ_T2F_TFR_04*/
    expr.reverse();  /*UNQ_T2F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFR_06*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFR_07*/
}

TEST_F(Mat3x3FromMat3x2_A2DMat3x2, hforward) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T2F_TFHF_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHF_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHF_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHF_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHF_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        for (int ii_1 = 0; ii_1 < 2; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHF_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2(A_a2d, C_a2d);  /*UNQ_T2F_TFHF_06*/
    expr.reverse();  /*UNQ_T2F_TFHF_07*/
    expr.hforward();  /*UNQ_T2F_TFHF_08*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHF_09*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHF_10*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHF_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHF_11*/
}

TEST_F(Mat3x3FromMat3x2_A2DMat3x2, hreverse) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data), Ap0(Ap0_data), Ap1(Ap1_data), Ap2(Ap2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x3_t Cb(Cb_data), Ch0(Ch0_data), Ch1(Ch1_data), Ch2(Ch2_data);  /*UNQ_T2F_TFHR_01*/
    Mat3x2_t Ab;  /*UNQ_T2F_TFHR_02*/
    Mat3x3_t C;  /*UNQ_T2F_TFHR_02*/
    A2DMat3x2_t A_a2d(A, Ab);  /*UNQ_T2F_TFHR_03*/
    A2DMat3x3_t C_a2d(C, Cb);  /*UNQ_T2F_TFHR_03*/
    // Set Derivative Values:
    /*UNQ_T2F_CDL_01*/
    for (int ii_0 = 0; ii_0 < 3; ii_0++) {  /*UNQ_T2F_CDL_02*/
        for (int ii_1 = 0; ii_1 < 3; ii_1++) {  /*UNQ_T2F_CDL_02*/
            if (ii_1 < 2) {A_a2d.pvalue(0)(ii_0, ii_1) = Ap0(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(1)(ii_0, ii_1) = Ap1(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 2) {A_a2d.pvalue(2)(ii_0, ii_1) = Ap2(ii_0, ii_1);  /*UNQ_T2F_TFHR_04*/ }
            if (ii_1 < 3) {C_a2d.hvalue(0)(ii_0, ii_1) = Ch0(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(1)(ii_0, ii_1) = Ch1(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
            if (ii_1 < 3) {C_a2d.hvalue(2)(ii_0, ii_1) = Ch2(ii_0, ii_1);  /*UNQ_T2F_TFHR_05*/ }
        }
    }
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2(A_a2d, C_a2d);  /*UNQ_T2F_TFHR_06*/
    expr.reverse();  /*UNQ_T2F_TFHR_07*/
    expr.hforward();  /*UNQ_T2F_TFHR_08*/
    expr.hreverse();  /*UNQ_T2F_TFHR_09*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_a2d.value(), C_out);  /*UNQ_T2F_TFHR_10*/
    expect_mat_eq<3, 2>(A_a2d.bvalue(), Ab_out);  /*UNQ_T2F_TFHR_11*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(0), Cp0_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(1), Cp1_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 3>(C_a2d.pvalue(2), Cp2_out);  /*UNQ_T2F_TFHR_12*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(0), Ah0_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(1), Ah1_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
    expect_mat_eq<3, 2>(A_a2d.hvalue(2), Ah2_out, 1e-8);  /*UNQ_T2F_TFHR_13*/
}
