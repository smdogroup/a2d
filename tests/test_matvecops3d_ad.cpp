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
using Mat3x3_t = A2D::Mat<T, 3, 3>;  /*UNQ_TC_TD_01*/
using ADMat3x3_t = A2D::ADMat<Mat3x3_t>;  /*UNQ_TC_TD_01*/
using Mat3x2_t = A2D::Mat<T, 3, 2>;  /*UNQ_TC_TD_01*/
using ADMat3x2_t = A2D::ADMat<Mat3x2_t>;  /*UNQ_TC_TD_01*/
using ADVec_t = A2D::ADVec<Vec_t>;  /*UNQ_TC_TD_01*/

class matvecops3d_adTest : public ::testing::Test {
protected:
    // Input Options:
    const T A_data[6] = {0.8696716749452119, 0.3123934153539396, 
0.539756796471562, 0.22324745061187035, 
0.4750063052837804, 0.8867632661556809};
    const T x_data[3] = {0.5091865469526378, 0.25444157059373695, 0.6223444159078971};
    const T y_data[3] = {0.8447415869436454, 0.05660059371215864, 0.5300425913528337};
    const T z_data[3] = {0.5428871900961088, 0.663815030391408, 0.4622604491401766};
    const T Ab_data[6] = {0.4371569663083663, 0.17663070533524516, 
0.0905213166025215, 0.48408509902693586, 
0.09676444216081104, 0.30652444030178005};
    const T ub_data[3] = {0.25853859707099647, 0.03636698777385605, 0.8749643717152243};
    const T vb_data[3] = {0.2986821714640696, 0.8340399047875401, 0.4880218904037843};
    const T xb_data[3] = {0.7107834004392464, 0.34530490358404275, 0.7580821768228523};
    const T yb_data[3] = {0.9027412776983588, 0.7566224901549622, 0.28291098882869137};
    const T zb_data[3] = {0.8009402504883824, 0.7020964963497749, 0.9934665217962456};
    const T Cb_data[9] = {0.49709266900490845, 0.3931450974708073, 0.4250556626414279, 
0.3529705062922772, 0.0831448824605504, 0.633714170902854, 
0.1365698081639385, 0.38416420211810287, 0.33488110986345954};
};  /*UNQ_TC_TIC_01*/

class Mat3x2ToVec3 : public matvecops3d_adTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x2ToVec3_Mat3x2 : public Mat3x2ToVec3 {
protected:
    // Results
    const T u_out[3] = {0.8696716749452119, 0.539756796471562, 0.4750063052837804};
    const T v_out[3] = {0.3123934153539396, 0.22324745061187035, 0.8867632661556809};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x2ToVec3_Mat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Vec_t u, v;  /*UNQ_T1F_TFP_02*/
    // Evaluations:
    A2D::Mat3x2ToVec3(A, u, v);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_vec_eq<3>(u, u_out);  /*UNQ_T1F_TFP_06*/
    expect_vec_eq<3>(v, v_out);  /*UNQ_T1F_TFP_06*/
}

class Mat3x2ToVec3_ADMat3x2 : public Mat3x2ToVec3 {
protected:
    // Results
    const T u_out[3] = {0.8696716749452119, 0.539756796471562, 0.4750063052837804};
    const T v_out[3] = {0.3123934153539396, 0.22324745061187035, 0.8867632661556809};
    const T ub_out[3] = {0.4371569663083663, 0.0905213166025215, 0.09676444216081102};
    const T vb_out[3] = {0.17663070533524516, 0.4840850990269358, 0.30652444030178005};
    const T Ab_out[6] = {0.25853859707099647, 0.2986821714640696, 
0.03636698777385605, 0.8340399047875401, 
0.8749643717152243, 0.4880218904037843};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x2ToVec3_ADMat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Vec_t u, v, ub, vb;  /*UNQ_T1F_TFP_02*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFP_02*/
    ADVec_t u_ad(u, ub), v_ad(v, vb);  /*UNQ_T1F_TFP_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x2ToVec3(A_ad, u_ad, v_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_vec_eq<3>(u_ad.value(), u_out);  /*UNQ_T1F_TFP_05*/
    expect_vec_eq<3>(v_ad.value(), v_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x2ToVec3_ADMat3x2, forward) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data), Ab(Ab_data);  /*UNQ_T1F_TFF_01*/
    Vec_t u, v, ub, vb;  /*UNQ_T1F_TFF_02*/
    ADVec_t u_ad(u, ub), v_ad(v, vb);  /*UNQ_T1F_TFF_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x2ToVec3(A_ad, u_ad, v_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_vec_eq<3>(u_ad.value(), u_out);  /*UNQ_T1F_TFF_06*/
    expect_vec_eq<3>(v_ad.value(), v_out);  /*UNQ_T1F_TFF_06*/
    expect_vec_eq<3>(u_ad.bvalue(), ub_out);  /*UNQ_T1F_TFF_07*/
    expect_vec_eq<3>(v_ad.bvalue(), vb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x2ToVec3_ADMat3x2, reverse) {
    // Declarations and Initializations:
    Vec_t ub(ub_data), vb(vb_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFR_01*/
    Vec_t u, v;  /*UNQ_T1F_TFR_02*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFR_02*/
    ADVec_t u_ad(u, ub), v_ad(v, vb);  /*UNQ_T1F_TFR_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x2ToVec3(A_ad, u_ad, v_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_vec_eq<3>(u_ad.value(), u_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(v_ad.value(), v_out);  /*UNQ_T1F_TFR_06*/
    expect_mat_eq<3, 2>(A_ad.bvalue(), Ab_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3 : public matvecops3d_adTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromThreeVec3_VecVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Mat3x3_t C;  /*UNQ_T1F_TFP_02*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y, z, C);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C, C_out);  /*UNQ_T1F_TFP_06*/
}

class Mat3x3FromThreeVec3_VecADVecADVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.0, 0.9027412776983588, 0.8009402504883825, 
0.0, 0.7566224901549622, 0.7020964963497749, 
0.0, 0.28291098882869137, 0.9934665217962456};
    const T yb_out[3] = {0.3931450974708073, 0.0831448824605504, 0.38416420211810287};
    const T zb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecADVecADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t yb, zb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecADVecADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yb(yb_data), zb(zb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecADVecADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t yb, zb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(y_ad.bvalue(), yb_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(z_ad.bvalue(), zb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_ADVecVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.7107834004392463, 0.0, 0.0, 
0.3453049035840427, 0.0, 0.0, 
0.7580821768228522, 0.0, 0.0};
    const T xb_out[3] = {0.49709266900490845, 0.3529705062922772, 0.1365698081639385};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_ADVecVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_ad, y, z, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecVecVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xb(xb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y, z, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y, z, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_VecVecADVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.0, 0.0, 0.8009402504883825, 
0.0, 0.0, 0.7020964963497749, 
0.0, 0.0, 0.9934665217962456};
    const T zb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecVecADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t zb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t z_ad(z, zb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y, z_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecVecADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), zb(zb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t z_ad(z, zb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y, z_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecVecADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t zb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t z_ad(z, zb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y, z_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(z_ad.bvalue(), zb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_ADVecVecADVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.7107834004392463, 0.0, 0.8009402504883825, 
0.3453049035840427, 0.0, 0.7020964963497749, 
0.7580821768228522, 0.0, 0.9934665217962456};
    const T xb_out[3] = {0.49709266900490845, 0.3529705062922772, 0.1365698081639385};
    const T zb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_ADVecVecADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb, zb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb), z_ad(z, zb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_ad, y, z_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecVecADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xb(xb_data), zb(zb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb), z_ad(z, zb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y, z_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecVecADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb, zb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb), z_ad(z, zb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y, z_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(z_ad.bvalue(), zb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_VecADVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.0, 0.9027412776983588, 0.0, 
0.0, 0.7566224901549622, 0.0, 
0.0, 0.28291098882869137, 0.0};
    const T yb_out[3] = {0.3931450974708073, 0.0831448824605504, 0.38416420211810287};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_VecADVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t yb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t y_ad(y, yb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x, y_ad, z, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_VecADVecVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), yb(yb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t y_ad(y, yb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_ad, z, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_VecADVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t yb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t y_ad(y, yb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x, y_ad, z, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(y_ad.bvalue(), yb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_ADVecADVecVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.7107834004392463, 0.9027412776983588, 0.0, 
0.3453049035840427, 0.7566224901549622, 0.0, 
0.7580821768228522, 0.28291098882869137, 0.0};
    const T xb_out[3] = {0.49709266900490845, 0.3529705062922772, 0.1365698081639385};
    const T yb_out[3] = {0.3931450974708073, 0.0831448824605504, 0.38416420211810287};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_ADVecADVecVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb, yb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecADVecVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xb(xb_data), yb(yb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecADVecVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb, yb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(y_ad.bvalue(), yb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromThreeVec3_ADVecADVecADVec : public Mat3x3FromThreeVec3 {
protected:
    // Results
    const T C_out[9] = {0.5091865469526378, 0.8447415869436454, 0.5428871900961088, 
0.25444157059373695, 0.05660059371215864, 0.663815030391408, 
0.6223444159078971, 0.5300425913528337, 0.4622604491401766};
    const T Cb_out[9] = {0.7107834004392463, 0.9027412776983588, 0.8009402504883825, 
0.3453049035840427, 0.7566224901549622, 0.7020964963497749, 
0.7580821768228522, 0.28291098882869137, 0.9934665217962456};
    const T xb_out[3] = {0.49709266900490845, 0.3529705062922772, 0.1365698081639385};
    const T yb_out[3] = {0.3931450974708073, 0.0831448824605504, 0.38416420211810287};
    const T zb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromThreeVec3_ADVecADVecADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb, yb, zb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecADVecADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data), xb(xb_data), yb(yb_data), zb(zb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromThreeVec3_ADVecADVecADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data), y(y_data), z(z_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb, yb, zb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb), y_ad(y, yb), z_ad(z, zb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromThreeVec3(x_ad, y_ad, z_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(y_ad.bvalue(), yb_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(z_ad.bvalue(), zb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromMat3x2AndVec3 : public matvecops3d_adTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromMat3x2AndVec3_Mat3x2Vec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.5091865469526378, 
0.539756796471562, 0.22324745061187035, 0.25444157059373695, 
0.4750063052837804, 0.8867632661556809, 0.6223444159078971};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2Vec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Mat3x3_t C;  /*UNQ_T1F_TFP_02*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A, x, C);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C, C_out);  /*UNQ_T1F_TFP_06*/
}

class Mat3x3FromMat3x2AndVec3_Mat3x2ADVec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.5091865469526378, 
0.539756796471562, 0.22324745061187035, 0.25444157059373695, 
0.4750063052837804, 0.8867632661556809, 0.6223444159078971};
    const T Cb_out[9] = {0.0, 0.0, 0.7107834004392463, 
0.0, 0.0, 0.3453049035840427, 
0.0, 0.0, 0.7580821768228522};
    const T xb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2ADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A, x_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2ADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), xb(xb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A, x_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_Mat3x2ADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A, x_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromMat3x2AndVec3_ADMat3x2Vec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.5091865469526378, 
0.539756796471562, 0.22324745061187035, 0.25444157059373695, 
0.4750063052837804, 0.8867632661556809, 0.6223444159078971};
    const T Cb_out[9] = {0.4371569663083663, 0.17663070533524516, 0.0, 
0.0905213166025215, 0.4840850990269358, 0.0, 
0.09676444216081102, 0.30652444030178005, 0.0};
    const T Ab_out[6] = {0.49709266900490845, 0.3931450974708073, 
0.3529705062922772, 0.0831448824605504, 
0.1365698081639385, 0.38416420211810287};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2Vec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A_ad, x, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2Vec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFF_01*/
    Mat3x2_t A(A_data), Ab(Ab_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_ad, x, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2Vec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_ad, x, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_mat_eq<3, 2>(A_ad.bvalue(), Ab_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromMat3x2AndVec3_ADMat3x2ADVec : public Mat3x3FromMat3x2AndVec3 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.5091865469526378, 
0.539756796471562, 0.22324745061187035, 0.25444157059373695, 
0.4750063052837804, 0.8867632661556809, 0.6223444159078971};
    const T Cb_out[9] = {0.4371569663083663, 0.17663070533524516, 0.7107834004392463, 
0.0905213166025215, 0.4840850990269358, 0.3453049035840427, 
0.09676444216081102, 0.30652444030178005, 0.7580821768228522};
    const T Ab_out[6] = {0.49709266900490845, 0.3931450974708073, 
0.3529705062922772, 0.0831448824605504, 
0.1365698081639385, 0.38416420211810287};
    const T xb_out[3] = {0.4250556626414279, 0.633714170902854, 0.33488110986345954};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2ADVec, passive) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Vec_t xb;  /*UNQ_T1F_TFP_02*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFP_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2AndVec3(A_ad, x_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2ADVec, forward) {
    // Declarations and Initializations:
    Vec_t x(x_data), xb(xb_data);  /*UNQ_T1F_TFF_01*/
    Mat3x2_t A(A_data), Ab(Ab_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFF_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_ad, x_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromMat3x2AndVec3_ADMat3x2ADVec, reverse) {
    // Declarations and Initializations:
    Vec_t x(x_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Vec_t xb;  /*UNQ_T1F_TFR_02*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADVec_t x_ad(x, xb);  /*UNQ_T1F_TFR_03*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2AndVec3(A_ad, x_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_mat_eq<3, 2>(A_ad.bvalue(), Ab_out);  /*UNQ_T1F_TFR_07*/
    expect_vec_eq<3>(x_ad.bvalue(), xb_out);  /*UNQ_T1F_TFR_07*/
}

class Mat3x3FromMat3x2 : public matvecops3d_adTest {
protected:
    
};  /*UNQ_TC_TF_01*/

class Mat3x3FromMat3x2_Mat3x2 : public Mat3x3FromMat3x2 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.0, 
0.539756796471562, 0.22324745061187035, 0.0, 
0.4750063052837804, 0.8867632661556809, 0.0};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2_Mat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Mat3x3_t C;  /*UNQ_T1F_TFP_02*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2(A, C);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C, C_out);  /*UNQ_T1F_TFP_06*/
}

class Mat3x3FromMat3x2_ADMat3x2 : public Mat3x3FromMat3x2 {
protected:
    // Results
    const T C_out[9] = {0.8696716749452119, 0.3123934153539396, 0.0, 
0.539756796471562, 0.22324745061187035, 0.0, 
0.4750063052837804, 0.8867632661556809, 0.0};
    const T Cb_out[9] = {0.4371569663083663, 0.17663070533524516, 0.0, 
0.0905213166025215, 0.4840850990269358, 0.0, 
0.09676444216081102, 0.30652444030178005, 0.0};
    const T Ab_out[6] = {0.49709266900490845, 0.3931450974708073, 
0.3529705062922772, 0.0831448824605504, 
0.1365698081639385, 0.38416420211810287};
};  /*UNQ_T1F_FTV_01*/

TEST_F(Mat3x3FromMat3x2_ADMat3x2, passive) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFP_01*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFP_02*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFP_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFP_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFP_03*/
    // Evaluations:
    A2D::Mat3x3FromMat3x2(A_ad, C_ad);  /*UNQ_T1F_TFP_04*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFP_05*/
}

TEST_F(Mat3x3FromMat3x2_ADMat3x2, forward) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data), Ab(Ab_data);  /*UNQ_T1F_TFF_01*/
    Mat3x3_t C, Cb;  /*UNQ_T1F_TFF_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFF_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFF_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2(A_ad, C_ad);  /*UNQ_T1F_TFF_04*/
    expr.forward();  /*UNQ_T1F_TFF_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFF_06*/
    expect_mat_eq<3, 3>(C_ad.bvalue(), Cb_out);  /*UNQ_T1F_TFF_07*/
}

TEST_F(Mat3x3FromMat3x2_ADMat3x2, reverse) {
    // Declarations and Initializations:
    Mat3x2_t A(A_data);  /*UNQ_T1F_TFR_01*/
    Mat3x3_t Cb(Cb_data);  /*UNQ_T1F_TFR_01*/
    Mat3x2_t Ab;  /*UNQ_T1F_TFR_02*/
    Mat3x3_t C;  /*UNQ_T1F_TFR_02*/
    ADMat3x2_t A_ad(A, Ab);  /*UNQ_T1F_TFR_03*/
    ADMat3x3_t C_ad(C, Cb);  /*UNQ_T1F_TFR_03*/
    // Evaluations:
    auto expr = A2D::Mat3x3FromMat3x2(A_ad, C_ad);  /*UNQ_T1F_TFR_04*/
    expr.reverse();  /*UNQ_T1F_TFR_05*/
    // Comparisons:
    expect_mat_eq<3, 3>(C_ad.value(), C_out);  /*UNQ_T1F_TFR_06*/
    expect_mat_eq<3, 2>(A_ad.bvalue(), Ab_out);  /*UNQ_T1F_TFR_07*/
}
