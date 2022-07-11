/*
  This is a set of unit tests for a2dtmp2d.h using Google Test framework.
 */
#include "a2dtmp2d.h"
#include "test_commons.h"

// Global typenames
typedef A2D::Mat<T, 2, 2> Mat;
typedef A2D::SymmMat<T, 2> SMat;
typedef A2D::ADMat<Mat> ADMat;
typedef A2D::ADMat<SMat> ADSMat;
typedef A2D::A2DMat<1, Mat> A2DMat;
typedef A2D::A2DMat<1, SMat> A2DSMat;
typedef A2D::ADScalar<T> ADScalar;
typedef A2D::A2DScalar<1, T> A2DScalar;

class ADExpressionTest : public ::testing::Test {
 protected:
  // Inputs
  const T A_data[4] = {0.96155402, 0.02855176, 0.95787560, 0.45439794};
  const T B_data[4] = {0.80766462, 0.60212270, 0.86418474, 0.65304149};
  const T dA_data[4] = {0.69626791, 0.53987667, 0.65374594, 0.45510150};
  const T dB_data[4] = {0.14453098, 0.50253050, 0.44495442, 0.81206852};
  const T Cb_data[4] = {0.95793759, 0.25352624, 0.24823463, 0.26332010};
  const T Ch_data[4] = {0.66773557, 0.09141460, 0.95390316, 0.44180188};

  // Symmetric matrices
  const T S_data[3] = {0.16665530, 0.09586054, 0.29580571};
  const T E_data[3] = {0.28147380, 0.47227010, 0.39857781};
  const T dS_data[3] = {0.14503781, 0.66755052, 0.83809867};
  const T dE_data[3] = {0.87502649, 0.41741302, 0.66047610};
  const T Sh_data[3] = {0.52028303, 0.73053919, 0.44468893};

  const T sb_data = 0.33622324;
  const T sh_data = 0.73157930;

  const T mu = 0.84010356;
  const T lam = 0.76764533;

  const T mub = 0.45061325;
  const T lamb = 0.30576147;

  // Outputs
  const T AB_data[4] = {0.8012871574649149, 0.5976189866107764,
                        1.1663259981167076, 0.8734993503266507};

  const T ABT_data[4] = {0.7938048249937244, 0.8496057946621772,
                         1.0472455469885100, 1.1245221841288746};

  const T ATB_data[4] = {1.6043946385111165, 1.2045060117768980,
                         0.4157440120261668, 0.3139323706114826};

  const T ATBT_data[4] = {1.3533718047088925, 1.4564928198282989,
                          0.2966635608979692, 0.3214147030826730};

  const T detA_data = 0.4095791316456628;
  const T invA_data[4] = {1.1094264939091454, -0.0697099969065339,
                          -2.3386826280707149, 2.3476636032127356};

  const T trSE = 0.2553548262821431;
};

// Test suite: pure mat-mat multiplication, no AD
class MatxMat : public ADExpressionTest {
 protected:
  void SetUp() {
    Mat A(A_data), B(B_data);
    this->A = A;
    this->B = B;
  };

  Mat A, B;
  Mat AB, ABT, ATB, ATBT;
};

TEST_F(MatxMat, AB) {
  A2D::MatMatMult<T, false, false>(A, B, AB);
  expect_mat_eq<2, 2, Mat>(AB, AB_data);
}

TEST_F(MatxMat, ATB) {
  A2D::MatMatMult<T, true, false>(A, B, ATB);
  expect_mat_eq<2, 2, Mat>(ATB, ATB_data);
}

TEST_F(MatxMat, ABT) {
  A2D::MatMatMult<T, false, true>(A, B, ABT);
  expect_mat_eq<2, 2, Mat>(ABT, ABT_data);
}

TEST_F(MatxMat, ATBT) {
  A2D::MatMatMult<T, true, true>(A, B, ATBT);
  expect_mat_eq<2, 2, Mat>(ATBT, ATBT_data);
}

// Test suite: C = AB, where A and B are both AD-able
class ADMatxADMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[4] = {1.1805827132886788, 1.2781967868198005,
                          1.2619283086857356, 1.5411993990007218};
  const T AB_dA_out[4] = {1.7231461889181325, 1.8549570993098863,
                          1.9472059912665993, 2.0970063392507581};
  const T AB_dB_out[4] = {2.3439643897445710, 2.7053325577551179,
                          0.6071253381844496, 0.7368125999252162};

  // ATB
  const T ATB_dC_out[4] = {1.6924935489969912, 2.1072327799153299,
                           1.0356440426689451, 1.0056225527241729};
  const T ATB_dA_out[4] = {2.6357798500742304, 1.4419612188046487,
                           2.8387375319643913, 1.5517010279550085};
  const T ATB_dB_out[4] = {1.6569934360138374, 2.0549304443793286,
                           2.0917927933036635, 2.4754196797764703};

  // ABT
  const T ABT_dC_out[4] = {1.0407454302179728, 1.4052996645210951,
                           1.1688259333339508, 1.6573706715031915};
  const T ABT_dA_out[4] = {2.0550117876199856, 1.5443754352708634,
                           2.3762937961289219, 1.7861084396098033};
  const T ABT_dB_out[4] = {2.1203227944105398, 0.5608272100702049,
                           2.9388264681134331, 0.7932295976969538};

  // ATBT
  const T ATBT_dC_out[4] = {1.5763222764945217, 2.2343356576166244,
                            0.9425416673171600, 1.1454598357948786};
  const T ATBT_dA_out[4] = {3.2040185117926345, 1.7511464679447204,
                            2.4082533122031182, 1.3155585314901530};
  const T ATBT_dB_out[4] = {1.5426302452540983, 1.9382096383836389,
                            2.1811393279318634, 2.6607101983788501};
};

TEST_F(ADMatxADMat, AB) {
  auto dC_out = AB_dC_out;
  auto dA_out = AB_dA_out;
  auto dB_out = AB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, false>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dA_out = ATB_dA_out;
  auto dB_out = ATB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, false>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dA_out = ABT_dA_out;
  auto dB_out = ABT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, true>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dA_out = ATBT_dA_out;
  auto dB_out = ATBT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, true>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

// Test suite: C = AB, where A and B are both A2D-able
class A2DMatxA2DMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[4] = {1.1805827132886788, 1.2781967868198005,
                          1.2619283086857356, 1.5411993990007218};
  const T AB_Ab_out[4] = {0.9263462037607140, 0.9933982006740742,
                          0.3590413376860606, 0.3864795296364952};
  const T AB_Bb_out[4] = {1.1588866357256398, 0.4960070740270449,
                          0.1401481086733206, 0.1268907313587764};
  const T AB_Ah_out[4] = {0.8602055292507194, 1.2688636627963492,
                          1.2046567507447885, 1.4371505735590697};
  const T AB_Bh_out[4] = {2.3850479684377248, 0.8597579484903980,
                          1.0826567652694772, 0.4600741865747527};

  // ATB
  const T ATB_Cp_out[4] = {1.6924935489969912, 2.1072327799153299,
                           1.0356440426689451, 1.0056225527241729};
  const T ATB_Ab_out[4] = {0.9263462037607140, 0.3590413376860606,
                           0.9933982006740742, 0.3864795296364952};
  const T ATB_Bb_out[4] = {0.9281962761530607, 0.2512974275458608,
                           1.0303823482924663, 0.3624987102563381};
  const T ATB_Ah_out[4] = {0.8602055292507194, 1.2046567507447885,
                           1.2688636627963492, 1.4371505735590697};
  const T ATB_Bh_out[4] = {1.4702967247379561, 0.4191968613354567,
                           1.8122790031863218, 0.5738968015060868};

  // ABT
  const T ABT_Cp_out[4] = {1.0407454302179728, 1.4052996645210951,
                           1.1688259333339508, 1.6573706715031915};
  const T ABT_Ab_out[4] = {0.9927858074086435, 0.7423591216459906,
                           0.4280475402650646, 0.3214266561000500};
  const T ABT_Bb_out[4] = {1.1588866357256398, 0.1401481086733206,
                           0.4960070740270449, 0.1268907313587764};
  const T ABT_Ah_out[4] = {0.8695647786668232, 1.1490298087204351,
                           1.3052753134781960, 1.2014611424705672};
  const T ABT_Bh_out[4] = {2.3850479684377248, 1.0826567652694772,
                           0.8597579484903980, 0.4600741865747527};

  // ATBT
  const T ATBT_Cp_out[4] = {1.5763222764945217, 2.2343356576166244,
                            0.9425416673171600, 1.1454598357948786};
  const T ATBT_Ab_out[4] = {0.9927858074086435, 0.4280475402650646,
                            0.7423591216459906, 0.3214266561000500};
  const T ATBT_Bb_out[4] = {0.9281962761530607, 1.0303823482924663,
                            0.2512974275458608, 0.3624987102563381};
  const T ATBT_Ah_out[4] = {0.8695647786668232, 1.3052753134781960,
                            1.1490298087204351, 1.2014611424705672};
  const T ATBT_Bh_out[4] = {1.4702967247379561, 1.8122790031863218,
                            0.4191968613354567, 0.5738968015060868};
};

TEST_F(A2DMatxA2DMat, AB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;
  const T *Bp_in = dB_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = AB_Cp_out;
  auto Ab_out = AB_Ab_out;
  auto Bb_out = AB_Bb_out;
  auto Ah_out = AB_Ah_out;
  auto Bh_out = AB_Bh_out;

  Mat Ab, Ah;
  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, false>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out, 1e-8);
}

TEST_F(A2DMatxA2DMat, ATB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;
  const T *Bp_in = dB_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATB_Cp_out;
  auto Ab_out = ATB_Ab_out;
  auto Bb_out = ATB_Bb_out;
  auto Ah_out = ATB_Ah_out;
  auto Bh_out = ATB_Bh_out;

  Mat Ab, Ah;
  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, false>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out, 1e-8);
}

TEST_F(A2DMatxA2DMat, ABT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;
  const T *Bp_in = dB_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ABT_Cp_out;
  auto Ab_out = ABT_Ab_out;
  auto Bb_out = ABT_Bb_out;
  auto Ah_out = ABT_Ah_out;
  auto Bh_out = ABT_Bh_out;

  Mat Ab, Ah;
  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, true>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out, 1e-8);
}

TEST_F(A2DMatxA2DMat, ATBT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;
  const T *Bp_in = dB_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATBT_Cp_out;
  auto Ab_out = ATBT_Ab_out;
  auto Bb_out = ATBT_Bb_out;
  auto Ah_out = ATBT_Ah_out;
  auto Bh_out = ATBT_Bh_out;

  Mat Ab, Ah;
  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, true>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out, 1e-8);
}

// Test suite: C = AB, where only A is AD-able
class ADMatxMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[4] = {1.0289041366443599, 0.7718005788855953,
                          0.9212992376577528, 0.6908354321680730};
  const T AB_dA_out[4] = {1.2957281169594528, 1.3931810538292422,
                          1.1600684943618456, 1.2473169421252956};

  // ATB
  const T ATB_dC_out[4] = {1.1273082221332997, 0.8461619366316077,
                           0.8293310569735255, 0.6222721598686439};
  const T ATB_dA_out[4] = {1.4199802767740197, 1.0445055460196604,
                           1.5267814147233185, 1.1230647823107287};

  // ABT
  const T ABT_dC_out[4] = {0.8874229551557533, 0.9542659677667317,
                           0.8020344101606928, 0.8621574268461906};
  const T ABT_dA_out[4] = {1.5414022111004901, 1.1575127752470395,
                           1.3928381088675044, 1.0459476951810660};

  // ATBT
  const T ATBT_dC_out[4] = {0.9559862274551822, 1.0286273255127438,
                            0.7100662294764654, 0.7637533413572508};
  const T ATBT_dA_out[4] = {1.6610402909779491, 1.2335193541298892,
                            1.2473573297456857, 0.9263096153036066};
};

TEST_F(ADMatxMat, AB) {
  auto dC_out = AB_dC_out;
  auto dA_out = AB_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, false>(A_, B, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dA_out = ATB_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, false>(A_, B, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dA_out = ABT_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, true>(A_, B, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dA_out = ATBT_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, true>(A_, B, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, dA_out);
}

// Test suite: C = AB, where only B is AD-able
class MatxADMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[4] = {0.1516785766443188, 0.5063962079342051,
                          0.3406290710279828, 0.8503639668326488};
  const T AB_dB_out[4] = {0.4721274209085946, 1.3014702044000945,
                          0.1591118384967193, 0.4008621377728315};

  // ATB
  const T ATB_dC_out[4] = {0.5651853268636916, 1.2610708432837221,
                           0.2063129856954196, 0.3833503928555288};
  const T ATB_dB_out[4] = {0.5493468219432557, 1.2235330672769698,
                           0.6351254297760029, 1.3821426194646445};

  // ABT
  const T ABT_dC_out[4] = {0.1533224750622196, 0.4510336967543636,
                           0.3667915231732580, 0.7952132446570008};
  const T ABT_dB_out[4] = {0.4987684925869255, 0.1710469390499732,
                           1.1954086281233909, 0.3742210660945006};

  // ATBT
  const T ATBT_dC_out[4] = {0.6203360490393396, 1.2057083321038804,
                            0.2324754378406948, 0.3817064944376280};
  const T ATBT_dB_out[4] = {0.6031242046118166, 0.6998411252305967,
                            1.1702520859016057, 1.3283652367960836};
};

TEST_F(MatxADMat, AB) {
  auto dC_out = AB_dC_out;
  auto dB_out = AB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, false>(A, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dB_out = ATB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, false>(A, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dB_out = ABT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, false, true>(A, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dB_out = ATBT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = A2D::MatMatMult<T, true, true>(A, B_, C_);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<2, 2, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dB, dB_out);
}

// Test suite: C = AB, where only A is A2D-able
class A2DMatxMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[4] = {1.0289041366443599, 0.7718005788855953,
                          0.9212992376577528, 0.6908354321680730};
  const T AB_Ab_out[4] = {0.9263462037607140, 0.9933982006740742,
                          0.3590413376860606, 0.3864795296364952};
  const T AB_Ah_out[4] = {0.5943492011759534, 0.6367444165409557,
                          1.0364527740888752, 1.1128635123097796};

  // ATB
  const T ATB_Cp_out[4] = {1.1273082221332997, 0.8461619366316077,
                           0.8293310569735255, 0.6222721598686439};
  const T ATB_Ab_out[4] = {0.9263462037607140, 0.3590413376860606,
                           0.9933982006740742, 0.3864795296364952};
  const T ATB_Ah_out[4] = {0.5943492011759534, 1.0364527740888752,
                           0.6367444165409557, 1.1128635123097796};

  // ABT
  const T ABT_Cp_out[4] = {0.8874229551557533, 0.9542659677667317,
                           0.8020344101606928, 0.8621574268461906};
  const T ABT_Ab_out[4] = {0.9927858074086435, 0.7423591216459906,
                           0.4280475402650646, 0.3214266561000500};
  const T ABT_Ah_out[4] = {0.6183054977377374, 0.4617562708861930,
                           1.1522322760375103, 0.8628817042377331};

  // ATBT
  const T ATBT_Cp_out[4] = {0.9559862274551822, 1.0286273255127438,
                            0.7100662294764654, 0.7637533413572508};
  const T ATBT_Ab_out[4] = {0.9927858074086435, 0.4280475402650646,
                            0.7423591216459906, 0.3214266561000500};
  const T ATBT_Ah_out[4] = {0.6183054977377374, 1.1522322760375103,
                            0.4617562708861930, 0.8628817042377331};
};

TEST_F(A2DMatxMat, AB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = AB_Cp_out;
  auto Ab_out = AB_Ab_out;
  auto Ah_out = AB_Ah_out;

  Mat Ab, Ah;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, false>(A__, B, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out);
}

TEST_F(A2DMatxMat, ATB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATB_Cp_out;
  auto Ab_out = ATB_Ab_out;
  auto Ah_out = ATB_Ah_out;

  Mat Ab, Ah;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, false>(A__, B, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out);
}

TEST_F(A2DMatxMat, ABT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ABT_Cp_out;
  auto Ab_out = ABT_Ab_out;
  auto Ah_out = ABT_Ah_out;

  Mat Ab, Ah;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, true>(A__, B, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out);
}

TEST_F(A2DMatxMat, ATBT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Ap_in = dA_data;

  Mat A(A_data), Ap(Ap_in);
  Mat B(B_data);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATBT_Cp_out;
  auto Ab_out = ATBT_Ab_out;
  auto Ah_out = ATBT_Ah_out;

  Mat Ab, Ah;
  Mat C, Cp;

  // A2D data types
  A2DMat A__(A, Ab), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, true>(A__, B, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out);
}

// Test suite: C = AB, where only B is A2D-able
class MatxA2DMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[4] = {0.1516785766443188, 0.5063962079342051,
                          0.3406290710279828, 0.8503639668326488};
  const T AB_Bb_out[4] = {1.1588866357256398, 0.4960070740270449,
                          0.1401481086733206, 0.1268907313587764};
  const T AB_Bh_out[4] = {1.5557843833573874, 0.5110913170028201,
                          0.4525166566015936, 0.2033639118798232};

  // ATB
  const T ATB_Cp_out[4] = {0.5651853268636916, 1.2610708432837221,
                           0.2063129856954196, 0.3833503928555288};
  const T ATB_Bb_out[4] = {0.9281962761530607, 0.2512974275458608,
                           1.0303823482924663, 0.3624987102563381};
  const T ATB_Bh_out[4] = {0.6692994357180529, 0.1005142973620008,
                           1.0730592406185824, 0.2883176789838872};

  // ABT
  const T ABT_Cp_out[4] = {0.1533224750622196, 0.4510336967543636,
                           0.3667915231732580, 0.7952132446570008};
  const T ABT_Bb_out[4] = {1.1588866357256398, 0.1401481086733206,
                           0.4960070740270449, 0.1268907313587764};
  const T ABT_Bh_out[4] = {1.5557843833573874, 0.4525166566015936,
                           0.5110913170028201, 0.2033639118798232};

  // ATBT
  const T ATBT_Cp_out[4] = {0.6203360490393396, 1.2057083321038804,
                            0.2324754378406948, 0.3817064944376280};
  const T ATBT_Bb_out[4] = {0.9281962761530607, 1.0303823482924663,
                            0.2512974275458608, 0.3624987102563381};
  const T ATBT_Bh_out[4] = {0.6692994357180529, 1.0730592406185824,
                            0.1005142973620008, 0.2883176789838872};
};

TEST_F(MatxA2DMat, AB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Bp_in = dB_data;

  Mat A(A_data);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = AB_Cp_out;
  auto Bb_out = AB_Bb_out;
  auto Bh_out = AB_Bh_out;

  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, false>(A, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out);
}

TEST_F(MatxA2DMat, ATB) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Bp_in = dB_data;

  Mat A(A_data);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATB_Cp_out;
  auto Bb_out = ATB_Bb_out;
  auto Bh_out = ATB_Bh_out;

  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, false>(A, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out);
}

TEST_F(MatxA2DMat, ABT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Bp_in = dB_data;

  Mat A(A_data);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ABT_Cp_out;
  auto Bb_out = ABT_Bb_out;
  auto Bh_out = ABT_Bh_out;

  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, false, true>(A, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out);
}

TEST_F(MatxA2DMat, ATBT) {
  // Set inputs
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;
  const T *Bp_in = dB_data;

  Mat A(A_data);
  Mat B(B_data), Bp(Bp_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Cp_out = ATBT_Cp_out;
  auto Bb_out = ATBT_Bb_out;
  auto Bh_out = ATBT_Bh_out;

  Mat Bb, Bh;
  Mat C, Cp;

  // A2D data types
  A2DMat B__(B, Bb), C__(C, Cb);

  auto expr = A2D::MatMatMult<1, T, true, true>(A, B__, C__);

  // Check expression result
  expect_mat_eq<2, 2, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(B__.hvalue(0), Bh_out);
}

// Test suite: C = det(A)
class DetA : public ADExpressionTest {
 protected:
  // AD
  const T AD_sb_out = 0.2181870944430290;
  const T AD_Ab_out[4] = {0.0991437662494978, -0.2089960940018731,
                          -0.0062296255556347, 0.2097986777738142};

  // A2D
  const T A2D_sp_out = 0.2181870944430290;
  const T A2D_Ab_out[4] = {0.1527791476361256, -0.3220600377489440,
                           -0.0095997652549024, 0.3232968080394248};
  const T A2D_Ah_out[4] = {0.4854438271673693, -0.9205665401420510,
                           -0.2024069598158007, 0.9375544706684048};
};

TEST_F(DetA, A) {
  Mat A(A_data);
  T det;
  A2D::MatDet(A, det);
  expect_val_eq(det, detA_data);
}

TEST_F(DetA, AD) {
  // Set inputs
  Mat A(A_data), dA(dA_data);
  ADMat A_(A, dA);

  // Set outputs
  auto sb_out = AD_sb_out;
  auto Ab_out = AD_Ab_out;
  ADScalar det;

  auto expr = A2D::MatDet(A_, det);

  // Check expression result
  expect_val_eq(det.value, detA_data);

  // Check forward AD result
  expr.forward();
  expect_val_eq(det.bvalue, sb_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(dA, Ab_out);
}

TEST_F(DetA, A2D) {
  // Set inputs
  const T *A_in = A_data;
  const T *Ap_in = dA_data;
  const T sb = sb_data;
  const T sh = sh_data;

  Mat A(A_in), Ap(Ap_in);

  // Set outputs
  auto sp_out = A2D_sp_out;
  auto Ab_out = A2D_Ab_out;
  auto Ah_out = A2D_Ah_out;

  Mat Ab;

  // A2D data types
  A2DMat A__(A, Ab);
  A2DScalar s__(0.0, sb);

  auto expr = A2D::MatDet(A__, s__);

  // Check expression result
  expect_val_eq(s__.value, detA_data);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(s__.pvalue[0], sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out);

  // Check hreverse
  s__.hvalue[0] = sh;
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
}

// Test suite: C = inv(A)
class InvA : public ADExpressionTest {
 protected:
  // AD
  const T AD_invAb_out[4] = {0.5201411408336892, -1.2809901867149651,
                             -0.3503000066931826, 0.4493344398331727};
  const T AD_Ab_out[4] = {-1.7214151608155448, 9.0689598641803251,
                          1.0323652505658460, -4.6942585519301279};

  // A2D
  const T A2D_Ab_out[4] = {-0.5583086126097820, 1.9131768361965673,
                           -0.5305951394897133, -0.2030589394010456};
  const T A2D_Ah_out[4] = {0.4552365695998450, -0.1090881922485685,
                           -0.6368853232331602, 0.4950724611175659};
  const T A2D_Cp_out[4] = {0.5201411408336892, -1.2809901867149651,
                           -0.3503000066931826, 0.4493344398331727};
};

TEST_F(InvA, A) {
  Mat A(A_data), invA;
  A2D::MatInverse(A, invA);
  expect_mat_eq<2, 2, Mat>(invA, invA_data, 1e-14);
}

TEST_F(InvA, AD) {
  // Set inputs
  const T *A_in = A_data;
  const T *Ab_in = dA_data;
  Mat A(A_in), Ab(Ab_in);

  // Set outputs
  auto invAb_out = AD_invAb_out;
  auto Ab_out = AD_Ab_out;
  Mat invA, invAb;

  // AD data type
  ADMat A_(A, Ab), invA_(invA, invAb);

  // Check expression result
  auto expr = A2D::MatInverse(A_, invA_);
  expect_mat_eq<2, 2, Mat>(invA, invA_data, 1e-14);

  // Check forward
  expr.forward();
  expect_mat_eq<2, 2, Mat>(invAb, invAb_out, 1e-14);

  // Check reverse
  Ab.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out, 1e-12);
}

TEST_F(InvA, A2D) {
  // Set inputs
  const T *A_in = A_data;
  const T *Ap_in = dA_data;
  const T *Cb_in = Cb_data;
  const T *Ch_in = Ch_data;

  Mat A(A_in), Ap(Ap_in);
  Mat Cb(Cb_in), Ch(Ch_in);

  // Set outputs
  auto Ab_out = A2D_Ab_out;
  auto Ah_out = A2D_Ah_out;
  auto Cp_out = A2D_Cp_out;

  Mat Ab, C;

  // A2D data types
  A2DMat A__(A, Ab), C__(C, Cb);

  auto expr = A2D::MatInverse(A__, C__);

  // Check exprssion
  expect_mat_eq<2, 2, Mat>(C, invA_data, 1e-14);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, Mat>(C__.pvalue(0), Cp_out, 1e-14);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Ab, Ab_out, 1e-13);

  // Check hreverse
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(A__.hvalue(0), Ah_out, 1e-8);
}

// Test suite: trace = tr(S)
class SymmTracesSuite : public ADExpressionTest {
 protected:
  const T tr_out = 0.4624610100000000;

  // AD
  const T AD_trb_out = 0.9831364800000001;
  const T AD_Sb_out[4] = {0.9831364800000001, 0.0000000000000000,
                          0.0000000000000000, 0.9831364800000001};

  // A2D
  const T A2D_Sb_out[4] = {0.3362232400000000, 0.0000000000000000,
                           0.0000000000000000, 0.3362232400000000};
  const T A2D_Sh_out[4] = {0.7315793000000000, 0.0000000000000000,
                           0.0000000000000000, 0.7315793000000000};
  const T A2D_trp_out = 0.9831364800000001;
};

TEST_F(SymmTracesSuite, A) {
  // Set inputs
  const T *S_in = S_data;
  SMat S(S_in);

  // Set outputs
  T trace;

  // Compute
  A2D::SymmTrace(S, trace);
  expect_val_eq(trace, tr_out);
}

TEST_F(SymmTracesSuite, AD) {
  // Set inputs
  const T *S_in = S_data;
  const T *Sb_in = dS_data;
  SMat S(S_in), Sb(Sb_in);

  // Set outputs
  auto trb_out = AD_trb_out;
  auto Sb_out = AD_Sb_out;

  // AD types
  ADSMat S_(S, Sb);
  ADScalar trace_;

  auto expr = A2D::SymmTrace(S_, trace_);

  // Check expression result
  expect_val_eq(trace_.value, tr_out);

  // Check forward
  expr.forward();
  expect_val_eq(trace_.bvalue, trb_out);

  // Check reverse
  Sb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Sb, Sb_out);
}

TEST_F(SymmTracesSuite, A2D) {
  // Set inputs
  const T *S_in = S_data;
  const T *Sp_in = dS_data;
  const T trb = sb_data;
  const T trh = sh_data;

  SMat S(S_in), Sp(Sp_in);

  // Set outputs
  auto Sb_out = A2D_Sb_out;
  auto Sh_out = A2D_Sh_out;
  auto trp_out = A2D_trp_out;

  SMat Sb, Eb;

  // A2D types
  A2DSMat S__(S, Sb);
  A2DScalar trace__(0.0, trb);

  auto expr = A2D::SymmTrace(S__, trace__);

  // Check expression
  expect_val_eq(trace__.value, tr_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = i; j < 2; j++) {
      S__.pvalue(0)(i, j) = Sp(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(trace__.pvalue[0], trp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Sb, Sb_out);

  // Check hreverse
  trace__.hvalue[0] = trh;
  expr.hreverse();
  expect_mat_eq<2, 2, SMat>(S__.hvalue(0), Sh_out);
}

// Test suite: trace = tr(SE)
class SymmMultTrace : public ADExpressionTest {
 protected:
  // AD
  const T AD_Cb_out = 1.4266274565456845;
  const T AD_Sb_out[4] = {0.4015582513782487, 1.3475069831311521,
                          1.3475069831311521, 0.5686220473158491};
  const T AD_Eb_out[4] = {0.2377550267588580, 0.2735145567265917,
                          0.2735145567265917, 0.4220045476889903};

  // A2D
  const T A2D_Sb_out[4] = {0.0946380330111120, 0.3175763663542480,
                           0.3175763663542480, 0.1340111226703044};
  const T A2D_Sh_out[4] = {0.5001246470729986, 0.9716939742924943,
                           0.9716939742924943, 0.5136586894799157};
  const T A2D_Eb_out[4] = {0.0560333849291720, 0.0644610826938992,
                           0.0644610826938992, 0.0994567542267004};
  const T A2D_Eh_out[4] = {0.1706866501072147, 0.5891511705403434,
                           0.5891511705403434, 0.4981935844741603};
  const T A2D_sp_out = 1.4266274565456845;
};

TEST_F(SymmMultTrace, A) {
  // Set inputs
  const T *S_in = S_data;
  const T *E_in = E_data;
  SMat S(S_in), E(E_in);

  // Set outputs
  T trace;

  // Compute
  A2D::SymmSymmMultTrace(S, E, trace);
  expect_val_eq(trace, trSE);
}

TEST_F(SymmMultTrace, AD) {
  // Set inputs
  const T *S_in = S_data;
  const T *E_in = E_data;
  const T *Sb_in = dS_data;
  const T *Eb_in = dE_data;
  SMat S(S_in), E(E_in), Sb(Sb_in), Eb(Eb_in);

  // Set outputs
  auto Cb_out = AD_Cb_out;
  auto Sb_out = AD_Sb_out;
  auto Eb_out = AD_Eb_out;

  // AD types
  ADSMat S_(S, Sb), E_(E, Eb);
  ADScalar trace_;

  auto expr = A2D::SymmSymmMultTrace(S_, E_, trace_);

  // Check expression result
  expect_val_eq(trace_.value, trSE);

  // Check forward
  expr.forward();
  expect_val_eq(trace_.bvalue, Cb_out);

  // Check reverse
  Sb.zero();
  Eb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Sb, Sb_out);
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);
}

TEST_F(SymmMultTrace, A2D) {
  // Set inputs
  const T *S_in = S_data;
  const T *E_in = E_data;
  const T *Sp_in = dS_data;
  const T *Ep_in = dE_data;
  const T sb = sb_data;
  const T sh = sh_data;

  SMat S(S_in), E(E_in), Sp(Sp_in), Ep(Ep_in);

  // Set outputs
  auto Sb_out = A2D_Sb_out;
  auto Sh_out = A2D_Sh_out;
  auto Eb_out = A2D_Eb_out;
  auto Eh_out = A2D_Eh_out;
  auto sp_out = A2D_sp_out;

  SMat Sb, Eb;

  // A2D types
  A2DSMat S__(S, Sb);
  A2DSMat E__(E, Eb);
  A2DScalar trace__(0.0, sb);

  auto expr = A2D::SymmSymmMultTrace(S__, E__, trace__);

  // Check expression
  expect_val_eq(trace__.value, trSE);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = i; j < 2; j++) {
      S__.pvalue(0)(i, j) = Sp(i, j);
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(trace__.pvalue[0], sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Sb, Sb_out);
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);

  // Check hreverse
  trace__.hvalue[0] = sh;
  expr.hreverse();
  expect_mat_eq<2, 2, SMat>(S__.hvalue(0), Sh_out, 1e-9);
  expect_mat_eq<2, 2, SMat>(E__.hvalue(0), Eh_out, 1e-9);
}

// Test suite: symmetric isotropic constitutive, S = 2*mu*E + lambda * tr(E) * I
class SymmIsoConstitutive : public ADExpressionTest {
 protected:
  const T S_out[4] = {0.9949727254289373, 0.7935115845831120,
                      0.7935115845831120, 1.1917317168114885};

  // AD
  const T AD_Sb_out[4] = {2.6489471311030135, 0.7013403281847024,
                          0.7013403281847024, 2.2884580382262367};
  const T AD_Eb_out[4] = {8.2409558507363148, 1.1783970129590735,
                          1.1783970129590735, 7.6352595102024132};

  // A2D
  const T A2D_Sp_out[4] = {2.6489471311030135, 0.7013403281847024,
                           0.7013403281847024, 2.2884580382262367};
  const T A2D_Eb_out[4] = {0.9983936886558455, 1.1216231366637024,
                           1.1216231366637024, 2.1628794802211684};
  const T A2D_Eh_out[4] = {1.6149394700961202, 1.2274571484770327,
                           1.2274571484770327, 1.4879257250461282};
};

TEST_F(SymmIsoConstitutive, A) {
  // Set inputs
  const T *E_in = E_data;
  SMat E(E_in);

  // Set outputs
  SMat S;

  // Compute
  A2D::SymmIsotropicConstitutive(mu, lam, E, S);
  expect_mat_eq<2, 2, SMat>(S, S_out);
}

TEST_F(SymmIsoConstitutive, AD) {
  // Set inputs
  const T *E_in = E_data;
  const T *Eb_in = dE_data;

  SMat E(E_in), Eb(Eb_in);

  // Set outputs
  auto Sb_out = AD_Sb_out;
  auto Eb_out = AD_Eb_out;
  SMat S, Sb;

  // AD datatypes
  ADSMat E_(E, Eb), S_(S, Sb);

  auto expr = A2D::SymmIsotropicConstitutive(mu, lam, E_, S_);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(S, S_out);

  // Check forward AD
  expr.forward();
  expect_mat_eq<2, 2, SMat>(Sb, Sb_out);

  // Check reverse AD
  Eb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out, 1e-14);
}

TEST_F(SymmIsoConstitutive, A2D) {
  // Set inputs
  const T *E_in = E_data;
  const T *Ep_in = dE_data;
  const T *Sb_in = dS_data;
  const T *Sh_in = Sh_data;
  SMat E(E_in), Ep(Ep_in), Sb(Sb_in), Sh(Sh_in);

  // Set outputs
  auto Sp_out = A2D_Sp_out;
  auto Eb_out = A2D_Eb_out;
  auto Eh_out = A2D_Eh_out;

  SMat S, Sp, Eb, Eh;

  // A2D types
  A2DSMat E__(E, Eb), S__(S, Sb);

  auto expr = A2D::SymmIsotropicConstitutive(mu, lam, E__, S__);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(S, S_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, SMat>(S__.pvalue(0), Sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);

  // Check hreverse
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      S__.hvalue(0)(i, j) = Sh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(E__.hvalue(0), Eh_out);
}

// Test suite: symmetric isotropic energy w = f(E)
// S = 2*mu*E + lambda * tr(E) * I
// W = 1/2 * tr(ES)
class SymmIsoEnergy : public ADExpressionTest {
 protected:
  const T w_out = 0.7522800812607762;

  // AD_E_only
  const T AD_E_only_wb_out = 2.3201819419954175;
  const T AD_E_only_Eb_out[4] = {2.3085177503181851, 3.6821824986278116,
                                 3.6821824986278116, 2.7650344090492123};

  // AD_E_mu_lam
  const T AD_E_mu_lam_wb_out = 2.6991807315127199;
  const T AD_E_mu_lam_mub_out = 1.8466981178202271;
  const T AD_E_mu_lam_lamb_out = 0.6241453159284354;
  const T AD_E_mu_lam_Eb_out[4] = {2.6856112088584831, 4.2836623586777236,
                                   4.2836623586777236, 3.2166992871501430};

  // A2D_E_only
  const T A2D_E_only_wp_out = 2.3201819419954175;
  const T A2D_E_only_Eb_out[4] = {0.3345329534553476, 0.5335940718921359,
                                  0.5335940718921359, 0.4006878990371211};
  const T A2D_E_only_Eh_out[4] = {1.6185390415909939, 1.6326471345872511,
                                  1.6326471345872511, 1.6412790358718623};

  // A2D_E_mu_lam
  const T A2D_E_mu_lam_wp_out = 2.6991807315127199;
  const T A2D_E_mu_lam_Eb_out[4] = {0.3345329534553476, 0.5335940718921359,
                                    0.5335940718921359, 0.4006878990371211};
  const T A2D_E_mu_lam_Eh_out[4] = {1.7737414473551882, 1.9188553725435744,
                                    1.9188553725435744, 1.8319655105427990};
  const T A2D_E_mu_lam_mub_out = 0.2300338088614917;
  const T A2D_E_mu_lam_muh_out = 1.1082893898437427;
  const T A2D_E_mu_lam_lamb_out = 0.0777466132231439;
  const T A2D_E_mu_lam_lamh_out = 0.5202581814022481;
};

TEST_F(SymmIsoEnergy, A) {
  // Set inputs
  const T *E_in = E_data;
  SMat E(E_in);

  // Set outputs
  T w;

  // Compute
  A2D::SymmIsotropicEnergy(mu, lam, E, w);
  expect_val_eq(w, w_out);
}

TEST_F(SymmIsoEnergy, AD_E_only) {
  // Set inputs
  const T *E_in = E_data;
  const T *Eb_in = dE_data;
  SMat E(E_in), Eb(Eb_in);

  // Set outputs
  auto wb_out = AD_E_only_wb_out;
  auto Eb_out = AD_E_only_Eb_out;

  // AD data types
  ADSMat E_(E, Eb);
  ADScalar w_;

  // Compute
  auto expr = A2D::SymmIsotropicEnergy(mu, lam, E_, w_);

  // Check expression
  expect_val_eq(w_.value, w_out);

  // Check forward
  expr.forward();
  expect_val_eq(w_.bvalue, wb_out);

  // Check reverse
  Eb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);
}

TEST_F(SymmIsoEnergy, AD_E_mu_lam) {
  // Set inputs
  const T *E_in = E_data;
  const T *Eb_in = dE_data;
  SMat E(E_in), Eb(Eb_in);

  // Set outputs
  auto wb_out = AD_E_mu_lam_wb_out;
  auto Eb_out = AD_E_mu_lam_Eb_out;
  auto mub_out = AD_E_mu_lam_mub_out;
  auto lamb_out = AD_E_mu_lam_lamb_out;

  // AD data types
  ADScalar mu_(mu, mub);
  ADScalar lam_(lam, lamb);
  ADSMat E_(E, Eb);
  ADScalar w_;

  // Compute
  auto expr = A2D::SymmIsotropicEnergy(mu_, lam_, E_, w_);

  // Check expression
  expect_val_eq(w_.value, w_out);

  // Check forward
  expr.forward();
  expect_val_eq(w_.bvalue, wb_out);

  // Check reverse
  Eb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);
  expect_val_eq(mu_.bvalue, mub_out, 1e-14);
  expect_val_eq(lam_.bvalue, lamb_out, 1e-14);
}

TEST_F(SymmIsoEnergy, A2D_E_only) {
  // Set inputs
  const T *E_in = E_data;
  const T *Ep_in = dE_data;
  const T wb = sb_data;
  const T wh = sh_data;
  SMat E(E_in), Ep(Ep_in);

  // Set outputs
  auto wp_out = A2D_E_only_wp_out;
  auto Eb_out = A2D_E_only_Eb_out;
  auto Eh_out = A2D_E_only_Eh_out;
  SMat Eb, Eh;

  // A2D data types
  A2DSMat E__(E, Eb);
  A2DScalar w__(0.0, wb);

  // Compute
  auto expr = A2D::SymmIsotropicEnergy(mu, lam, E__, w__);

  // Check expression
  expect_val_eq(w__.value, w_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(w__.pvalue[0], wp_out);

  // Check hreverse
  w__.hvalue[0] = wh;
  expr.hreverse();
  expect_mat_eq<2, 2, SMat>(E__.hvalue(0), Eh_out, 1e-8);
}

TEST_F(SymmIsoEnergy, A2D_E_mu_lam) {
  // Set inputs
  const T *E_in = E_data;
  const T *Ep_in = dE_data;
  const T mu_in = mu;
  const T mup_in = mub;
  const T lam_in = lam;
  const T lamp_in = lamb;
  const T wb_in = sb_data;
  const T wh_in = sh_data;

  SMat E(E_in), Ep(Ep_in);

  // Set outputs
  auto wp_out = A2D_E_mu_lam_wp_out;
  auto Eb_out = A2D_E_mu_lam_Eb_out;
  auto Eh_out = A2D_E_mu_lam_Eh_out;
  auto mub_out = A2D_E_mu_lam_mub_out;
  auto muh_out = A2D_E_mu_lam_muh_out;
  auto lamb_out = A2D_E_mu_lam_lamb_out;
  auto lamh_out = A2D_E_mu_lam_lamh_out;
  SMat Eb;

  // A2D data types
  A2DScalar mu__(mu_in, 0.0);
  A2DScalar lam__(lam_in, 0.0);
  A2DSMat E__(E, Eb);
  A2DScalar w__(0.0, wb_in);

  // Compute
  auto expr = A2D::SymmIsotropicEnergy(mu__, lam__, E__, w__);

  // Check expression
  expect_val_eq(w__.value, w_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);
  expect_val_eq(mu__.bvalue, mub_out);
  expect_val_eq(lam__.bvalue, lamb_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  mu__.pvalue[0] = mup_in;
  lam__.pvalue[0] = lamp_in;
  expr.hforward();
  expect_val_eq(w__.pvalue[0], wp_out);

  // Check hreverse
  w__.hvalue[0] = wh_in;
  expr.hreverse();
  expect_mat_eq<2, 2, SMat>(E__.hvalue(0), Eh_out, 1e-8);
  expect_val_eq(mu__.hvalue[0], muh_out, 1e-8);
  expect_val_eq(lam__.hvalue[0], lamh_out, 1e-8);
}

class GreenStrain : public ADExpressionTest {
 protected:
  const T E_out[4] = {1.8826099192267602, 0.7245690595111697,
                      0.7245690595111697, 0.5580442854376706};

  // AD
  const T AD_Eb_out[4] = {1.9919744023825618, 1.2328071486987491,
                          1.2328071486987491, 0.6773131132023492};
  const T AD_Uxb_out[4] = {3.9249648036485771, 1.2284473906603905,
                           2.8045597646092055, 1.5755207401985354};

  // A2D
  const T A2D_Uxb_out[4] = {0.2940293703749538, 0.6786473951117044,
                            0.6243702298434004, 1.5386441566023961};
  const T A2D_Uxh_out[4] = {1.3121750081321486, 1.4140596171119642,
                            1.2763332655667154, 1.5962616755510131};
  const T A2D_Ep_out[4] = {1.9919744023825618, 1.2328071486987491,
                           1.2328071486987491, 0.6773131132023492};
};

TEST_F(GreenStrain, A) {
  // Set inputs
  const T *Ux_in = A_data;
  Mat Ux(Ux_in);

  // Set outputs
  SMat E;

  // Compute
  A2D::MatGreenStrain(Ux, E);
  expect_mat_eq<2, 2, SMat>(E, E_out);
}

TEST_F(GreenStrain, AD) {
  // Set inputs
  const T *Ux_in = A_data;
  const T *Uxb_in = dA_data;
  Mat Ux(Ux_in), Uxb(Uxb_in);

  // Set outputs
  SMat E, Eb;
  auto Eb_out = AD_Eb_out;
  auto Uxb_out = AD_Uxb_out;

  // AD data type
  ADMat Ux_(Ux, Uxb);
  ADSMat E_(E, Eb);

  // Compute
  auto expr = A2D::MatGreenStrain(Ux_, E_);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(E, E_out);

  // Check forward
  expr.forward();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);

  // Check reverse
  Uxb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Uxb, Uxb_out, 1e-14);
}

TEST_F(GreenStrain, A2D) {
  // Set inputs
  const T *Ux_in = A_data;
  const T *Uxp_in = dA_data;
  const T *E_in = E_data;
  const T *Eb_in = dS_data;
  const T *Eh_in = Sh_data;

  Mat Ux(Ux_in), Uxp(Uxp_in);
  SMat E(E_in), Eb(Eb_in), Eh(Eh_in);

  // Set outputs
  auto Ep_out = A2D_Ep_out;
  auto Uxb_out = A2D_Uxb_out;
  auto Uxh_out = A2D_Uxh_out;
  Mat Uxb;

  // A2D data type
  A2DMat Ux__(Ux, Uxb);
  A2DSMat E__(E, Eb);

  // Compute
  auto expr = A2D::MatGreenStrain(Ux__, E__);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(E, E_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      Ux__.pvalue(0)(i, j) = Uxp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, SMat>(E__.pvalue(0), Ep_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Uxb, Uxb_out);

  // Check hreverse
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      E__.hvalue(0)(i, j) = Eh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(Ux__.hvalue(0), Uxh_out, 1e-8);
}

class LinearGreenStrain : public ADExpressionTest {
 protected:
  const T E_out[4] = {0.9615540200000000, 0.4932136800000000,
                      0.4932136800000000, 0.4543979400000000};

  // AD
  const T AD_Eb_out[4] = {0.6962679100000000, 0.5968113049999999,
                          0.5968113049999999, 0.4551015000000000};
  const T AD_Uxb_out[4] = {0.6962679100000000, 0.2984056525000000,
                           0.2984056525000000, 0.4551015000000000};

  // A2D
  const T A2D_Uxb_out[4] = {0.1450378100000000, 0.3337752600000000,
                            0.3337752600000000, 0.8380986700000000};
  const T A2D_Uxh_out[4] = {0.5202830300000000, 0.3652695950000000,
                            0.3652695950000000, 0.4446889300000000};
  const T A2D_Ep_out[4] = {0.6962679100000000, 0.5968113049999999,
                           0.5968113049999999, 0.4551015000000000};
};

TEST_F(LinearGreenStrain, A) {
  // Set inputs
  const T *Ux_in = A_data;
  Mat Ux(Ux_in);

  // Set outputs
  SMat E;

  // Compute
  A2D::MatLinearGreenStrain(Ux, E);
  expect_mat_eq<2, 2, SMat>(E, E_out);
}

TEST_F(LinearGreenStrain, AD) {
  // Set inputs
  const T *Ux_in = A_data;
  const T *Uxb_in = dA_data;
  Mat Ux(Ux_in), Uxb(Uxb_in);

  // Set outputs
  SMat E, Eb;
  auto Eb_out = AD_Eb_out;
  auto Uxb_out = AD_Uxb_out;

  // AD data type
  ADMat Ux_(Ux, Uxb);
  ADSMat E_(E, Eb);

  // Compute
  auto expr = A2D::MatLinearGreenStrain(Ux_, E_);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(E, E_out);

  // Check forward
  expr.forward();
  expect_mat_eq<2, 2, SMat>(Eb, Eb_out);

  // Check reverse
  Uxb.zero();
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Uxb, Uxb_out, 1e-14);
}

TEST_F(LinearGreenStrain, A2D) {
  // Set inputs
  const T *Ux_in = A_data;
  const T *Uxp_in = dA_data;
  const T *E_in = E_data;
  const T *Eb_in = dS_data;
  const T *Eh_in = Sh_data;

  Mat Ux(Ux_in), Uxp(Uxp_in);
  SMat E(E_in), Eb(Eb_in), Eh(Eh_in);

  // Set outputs
  auto Ep_out = A2D_Ep_out;
  auto Uxb_out = A2D_Uxb_out;
  auto Uxh_out = A2D_Uxh_out;
  Mat Uxb;

  // A2D data type
  A2DMat Ux__(Ux, Uxb);
  A2DSMat E__(E, Eb);

  // Compute
  auto expr = A2D::MatLinearGreenStrain(Ux__, E__);

  // Check expression result
  expect_mat_eq<2, 2, SMat>(E, E_out);

  // Check hforward
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      Ux__.pvalue(0)(i, j) = Uxp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<2, 2, SMat>(E__.pvalue(0), Ep_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<2, 2, Mat>(Uxb, Uxb_out);

  // Check hreverse
  for (I i = 0; i < 2; i++) {
    for (I j = 0; j < 2; j++) {
      E__.hvalue(0)(i, j) = Eh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<2, 2, Mat>(Ux__.hvalue(0), Uxh_out, 1e-8);
}