/*
  This is a set of unit tests for a2dtmp3d.h using Google Test framework.
 */
#include "a2dtmp3d.h"
#include "test_commons.h"

// Global typenames
typedef A2D::Mat<T, 3, 3> Mat;
typedef A2D::SymmMat<T, 3> SMat;
typedef A2D::ADMat<Mat> ADMat;
typedef A2D::ADMat<SMat> ADSMat;
typedef A2D::A2DMat<1, Mat> A2DMat;
typedef A2D::A2DMat<1, SMat> A2DSMat;
typedef A2D::ADScalar<T> ADScalar;
typedef A2D::A2DScalar<1, T> A2DScalar;

class ADExpressionTest : public ::testing::Test {
 protected:
  const T A_data[9] = {0.54881350, 0.71518937, 0.60276338,
                       0.54488318, 0.42365480, 0.64589411,
                       0.43758721, 0.89177300, 0.96366276};

  const T B_data[9] = {0.63263687, 0.19519233, 0.13896359,
                       0.83668126, 0.50298017, 0.20476354,
                       0.57524785, 0.83518588, 0.58878908};

  const T dA_data[9] = {0.69626791, 0.53987667, 0.65374594,
                        0.45510150, 0.98931708, 0.74223721,
                        0.40945653, 0.74347344, 0.88483858};

  const T dB_data[9] = {0.14453098, 0.50253050, 0.44495442,
                        0.81206852, 0.59942161, 0.72829853,
                        0.99976353, 0.17018835, 0.75282306};

  const T Cb_data[9] = {0.95793759, 0.25352624, 0.24823463,
                        0.26332010, 0.32297026, 0.40684991,
                        0.86658516, 0.82340961, 0.60516882};

  const T Ch_data[9] = {0.66773557, 0.09141460, 0.95390316,
                        0.44180188, 0.53497833, 0.96218661,
                        0.45728781, 0.18113042, 0.98689272};

  // Symmetric matrices
  const T S_data[6] = {0.16665530, 0.09586054, 0.29580571,
                       0.16117561, 0.23525668, 0.79122117};

  const T E_data[6] = {0.28147380, 0.47227010, 0.39857781,
                       0.92539657, 0.00320169, 0.15774955};

  const T dS_data[6] = {0.14503781, 0.66755052, 0.83809867,
                        0.61659135, 0.99207311, 0.10221423};

  const T dE_data[6] = {0.87502649, 0.41741302, 0.66047610,
                        0.84442668, 0.65391596, 0.25608671};

  const T Sh_data[6] = {0.52028303, 0.73053919, 0.44468893,
                        0.68831122, 0.54873182, 0.19163654};

  const T sb_data = 0.33622324;
  const T sh_data = 0.73157930;

  const T mu = 0.84010356;
  const T lam = 0.76764533;

  const T mub = 0.45061325;
  const T lamb = 0.30576147;

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

  const T detA_data = -0.0843033050078113;
  const T invA_data[9] = {
      1.9896085246745707,  1.7991377032023534,  -2.4503547327977944,
      2.8759089284053432,  -3.1447116618538922, 0.3088820759058467,
      -3.5648209297984907, 2.0931485459459962,  1.8645435624787421};

  const T trSE = 0.6799787610430104;
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
  MatMatMult<T, false, false>(A, B, AB);
  expect_mat_eq<3, 3, Mat>(AB, AB_data);
}

TEST_F(MatxMat, ATB) {
  MatMatMult<T, true, false>(A, B, ATB);
  expect_mat_eq<3, 3, Mat>(ATB, ATB_data);
}

TEST_F(MatxMat, ABT) {
  MatMatMult<T, false, true>(A, B, ABT);
  expect_mat_eq<3, 3, Mat>(ABT, ABT_data);
}

TEST_F(MatxMat, ATBT) {
  MatMatMult<T, true, true>(A, B, ATBT);
  expect_mat_eq<3, 3, Mat>(ATBT, ATBT_data);
}

// Test suite: C = AB, where A and B are both AD-able
class ADMatxADMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[9] = {
      2.5309795609466685, 1.7605315843727500, 1.8110639486651825,
      2.6111580106063448, 1.8440371446992301, 1.7400790974735909,
      3.1409496540879167, 2.1113329613916307, 2.2997699723047305};
  const T AB_dA_out[9] = {
      2.1965051974896728, 3.3739755089803416, 3.9926463475464047,
      2.2536643760547048, 3.4685258467790412, 4.0667163879885972,
      2.7187808496552606, 4.1608413667410273, 4.9242594590101207};
  const T AB_dB_out[9] = {
      4.1862512274560597, 2.8948806240783922, 2.9484261022963034,
      5.7173733983030974, 3.9231783914018710, 4.0833193140120247,
      6.2389295868764929, 4.2868496476293565, 4.4317525459621798};

  // ATB
  const T ATB_dC_out[9] = {
      2.0160854005628668, 1.3836680186544550, 1.4014879702188572,
      2.9359363150181164, 1.9890468669008348, 2.0134694285170829,
      3.1186643787516175, 2.0940155245376690, 2.2278868036613981};
  const T ATB_dA_out[9] = {
      1.7402871416558603, 2.5254271934239512, 2.6913129887950107,
      2.6697520865063527, 3.8691791544457992, 4.1187663150826097,
      3.1405495964941834, 4.5356237232696088, 4.8546525988150435};
  const T ATB_dB_out[9] = {
      5.0860220105016625, 3.4441167391378968, 3.5520560302714399,
      4.3566814898644850, 2.9491189762036796, 3.0556421742207389,
      6.5057426437481567, 4.3971784993448066, 4.5557624295975181};

  // ABT
  const T ABT_dC_out[9] = {
      1.3436386740022619, 2.3013307183401648, 2.3605172536202854,
      1.1632104393382681, 2.1971989388738971, 2.6281804026654543,
      1.4672920005589434, 2.4894526595805879, 2.6921814582489407};
  const T ABT_dA_out[9] = {
      4.1333981252622696, 3.3912623591136111, 2.0477922604834835,
      4.0861001141873992, 3.5272164141322744, 2.1589940527434077,
      4.5598130022054537, 3.7870214066295822, 2.2987763472095413};
  const T ABT_dB_out[9] = {
      2.0132890593902797, 2.7622471720175295, 2.9751816190839739,
      3.5495707548289639, 4.7967678102123426, 5.2053085562222154,
      3.9055992043166552, 5.2024728259130075, 5.7147146149165238};

  // ATBT
  const T ATBT_dC_out[9] = {
      1.1340639558972678, 1.9862861521408786, 1.9925467324846766,
      1.3510335479306113, 2.5857553995624358, 3.0330442317201616,
      1.5219086178529759, 2.6799724795708406, 2.9549673224909254};
  const T ATBT_dA_out[9] = {
      3.5255472958187832, 4.7629188940374441, 4.9049368544803498,
      2.8845700286044087, 4.0974507923272014, 4.1129848861414917,
      1.7375023400430898, 2.5030362240132145, 2.5001030284477861};
  const T ATBT_dB_out[9] = {
      2.5059852234013720, 2.1732960343839944, 3.1676737817093574,
      4.5547947004382259, 3.9087400411359310, 5.7576699418982340,
      5.0438836305287351, 4.2792649361147781, 6.4242918848177846};
};

TEST_F(ADMatxADMat, AB) {
  auto dC_out = AB_dC_out;
  auto dA_out = AB_dA_out;
  auto dB_out = AB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, false, false>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dA_out = ATB_dA_out;
  auto dB_out = ATB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, true, false>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dA_out = ABT_dA_out;
  auto dB_out = ABT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, false, true>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(ADMatxADMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dA_out = ATBT_dA_out;
  auto dB_out = ATBT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dB(dB_data), dC;
  ADMat A_(A, dA), B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, true, true>(A_, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

// Test suite: C = AB, where A and B are both A2D-able
class A2DMatxA2DMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[9] = {
      2.5309795609466685, 1.7605315843727500, 1.8110639486651825,
      2.6111580106063448, 1.8440371446992301, 1.7400790974735909,
      3.1409496540879167, 2.1113329613916307, 2.2997699723047305};
  const T AB_Ab_out[9] = {
      0.6900085914418043, 0.9798365026866145, 0.9089509143610131,
      0.2861646455269697, 0.4660706571513516, 0.6607633064056966,
      0.7930533953144042, 1.2631307790643582, 1.5425181226326985};
  const T AB_Bb_out[9] = {
      1.0484143573871867, 0.6754331993575550, 0.6227343243778110,
      1.5694608536815782, 1.0524416309088469, 0.8895719000108111,
      1.5825824482727067, 1.1549100993633834, 0.9955873604820227};
  const T AB_Ah_out[9] = {
      0.9491448004736809, 1.9106551453921896, 2.2098444365044738,
      0.8990227148307055, 1.5394890917929669, 1.8919865033261383,
      1.2701039395508038, 2.3138284075115210, 2.4575030667653941};
  const T AB_Bh_out[9] = {
      1.9489435412107077, 1.0815874944545052, 2.0854320025701192,
      2.4944840717793531, 1.5221281781350136, 2.9563872527080322,
      2.7169967870828575, 1.7092362527036455, 3.1472185644883757};

  // ATB
  const T ATB_Cp_out[9] = {
      2.0160854005628668, 1.3836680186544550, 1.4014879702188572,
      2.9359363150181164, 1.9890468669008348, 2.0134694285170829,
      3.1186643787516175, 2.0940155245376690, 2.2278868036613981};
  const T ATB_Ab_out[9] = {
      0.6900085914418043, 0.2861646455269697, 0.7930533953144042,
      0.9798365026866145, 0.4660706571513516, 1.2631307790643582,
      0.9089509143610131, 0.6607633064056966, 1.5425181226326985};
  const T ATB_Bb_out[9] = {
      1.2363986180762427, 0.8664446795424580, 0.7919828503427733,
      1.1932431552396239, 0.8068055019872884, 0.6984977682252416,
      1.4890988399601655, 1.1924451750774940, 1.0546207193046553};
  const T ATB_Ah_out[9] = {
      0.9491448004736809, 0.8990227148307055, 1.2701039395508038,
      1.9106551453921896, 1.5394890917929669, 2.3138284075115210,
      2.2098444365044738, 1.8919865033261383, 2.4575030667653941};
  const T ATB_Bh_out[9] = {
      2.3337387717260243, 1.4311461538750090, 2.5946366042394642,
      2.1860465339635566, 1.4395129084375375, 2.5294832946630956,
      2.4816440504319157, 1.7641421050954142, 3.1660997360901537};

  // ABT
  const T ABT_Cp_out[9] = {
      1.3436386740022619, 2.3013307183401648, 2.3605172536202854,
      1.1632104393382681, 2.1971989388738971, 2.6281804026654543,
      1.4672920005589434, 2.4894526595805879, 2.6921814582489407};
  const T ABT_Ab_out[9] = {
      0.9609437297222513, 0.5218227993843698, 0.3311892163094781,
      0.6708487039516081, 0.5536410002458479, 0.3422732243744622,
      1.5852871757937947, 1.0887379356215181, 0.6453448442604295};
  const T ABT_Bb_out[9] = {
      1.0484143573871867, 1.5694608536815782, 1.5825824482727067,
      0.6754331993575550, 1.0524416309088469, 1.1549100993633834,
      0.6227343243778110, 0.8895719000108111, 0.9955873604820227};
  const T ABT_Ah_out[9] = {
      1.6401580355322301, 1.6486116482647004, 1.4709151802624600,
      1.9876878791512917, 1.5540873510089139, 1.3961336044378592,
      2.4074926203474400, 2.0366506204060517, 2.1225708690842846};
  const T ABT_Bh_out[9] = {
      1.9489435412107077, 2.4944840717793531, 2.7169967870828575,
      1.0815874944545052, 1.5221281781350136, 1.7092362527036455,
      2.0854320025701192, 2.9563872527080322, 3.1472185644883757};

  // ATBT
  const T ATBT_Cp_out[9] = {
      1.1340639558972678, 1.9862861521408786, 1.9925467324846766,
      1.3510335479306113, 2.5857553995624358, 3.0330442317201616,
      1.5219086178529759, 2.6799724795708406, 2.9549673224909254};
  const T ATBT_Ab_out[9] = {
      0.9609437297222513, 0.6708487039516081, 1.5852871757937947,
      0.5218227993843698, 0.5536410002458479, 1.0887379356215181,
      0.3311892163094781, 0.3422732243744622, 0.6453448442604295};
  const T ATBT_Bb_out[9] = {
      1.2363986180762427, 1.1932431552396239, 1.4890988399601655,
      0.8664446795424580, 0.8068055019872884, 1.1924451750774940,
      0.7919828503427733, 0.6984977682252416, 1.0546207193046553};
  const T ATBT_Ah_out[9] = {
      1.6401580355322301, 1.9876878791512917, 2.4074926203474400,
      1.6486116482647004, 1.5540873510089139, 2.0366506204060517,
      1.4709151802624600, 1.3961336044378592, 2.1225708690842846};
  const T ATBT_Bh_out[9] = {
      2.3337387717260243, 2.1860465339635566, 2.4816440504319157,
      1.4311461538750090, 1.4395129084375375, 1.7641421050954142,
      2.5946366042394642, 2.5294832946630956, 3.1660997360901537};
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

  auto expr = MatMatMult<1, T, false, false>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out, 1e-8);
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

  auto expr = MatMatMult<1, T, true, false>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out, 1e-8);
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

  auto expr = MatMatMult<1, T, false, true>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out, 1e-8);
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

  auto expr = MatMatMult<1, T, true, true>(A__, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out, 1e-8);
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out, 1e-8);
}

// Test suite: C = AB, where only A is AD-able
class ADMatxMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[9] = {
      1.2682553901952749, 0.9534527931080914, 0.5922214170543438,
      1.5426274087687242, 1.2063452326563935, 0.7028397697553150,
      1.3900890828380486, 1.1928798594613501, 0.7301190962948266};
  const T AB_dA_out[9] = {
      1.0707490068343024, 1.6619587196547210, 1.8745649998149054,
      1.3090614498062219, 2.0413711335872371, 2.3087399863868585,
      1.2137225763148520, 1.9125581706215586, 2.2258083223674259};

  // ATB
  const T ATB_dC_out[9] = {
      1.0567986362626920, 0.7067854978241817, 0.4310276161733993,
      1.5969704456208476, 1.2239251774552717, 0.7153483104457437,
      1.5436012400805454, 1.2399414794180161, 0.7638132948763544};
  const T ATB_dA_out[9] = {
      0.8664260345286076, 1.3486105606532113, 1.3247083610296713,
      1.3279614449291188, 2.0982425908284346, 2.0715693207538246,
      1.4520027649459017, 2.3620481153309756, 2.3732598374319229};

  // ABT
  const T ABT_dC_out[9] = {
      0.6367114191641073, 0.9879649044270280, 1.2363424605432483,
      0.5841650417794853, 1.0303648881695169, 1.5250809794832723,
      0.5271179562342785, 0.8977196826466857, 1.3774608012606941};
  const T ABT_dA_out[9] = {
      1.9406221826762378, 1.6537837069171670, 1.0187238357354746,
      2.1089508910118697, 1.9060037421440181, 1.1901098601869264,
      1.8769608535409716, 1.7048603921186007, 1.0681043416085370};

  // ATBT
  const T ATBT_dC_out[9] = {
      0.5862166227930794, 0.8953031106325379, 1.0217044986170061,
      0.6379687309408689, 1.1015478190528851, 1.5745755924835247,
      0.6814231413367092, 1.1014902548345833, 1.5169552773005301};
  const T ATBT_dA_out[9] = {
      1.7076789002364965, 2.2310181825343474, 2.2253149193868276,
      1.4180578700525981, 1.9936466141847144, 1.9539759544963609,
      0.8663566524997005, 1.2413041705566648, 1.2134047520075819};
};

TEST_F(ADMatxMat, AB) {
  auto dC_out = AB_dC_out;
  auto dA_out = AB_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = MatMatMult<T, false, false>(A_, B, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dA_out = ATB_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = MatMatMult<T, true, false>(A_, B, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dA_out = ABT_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = MatMatMult<T, false, true>(A_, B, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
}

TEST_F(ADMatxMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dA_out = ATBT_dA_out;

  Mat A(A_data), B(B_data), C;
  Mat dA(dA_data), dC;
  ADMat A_(A, dA), C_(C, dC);

  auto expr = MatMatMult<T, true, true>(A_, B, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, dA_out);
}

// Test suite: C = AB, where only B is AD-able
class MatxADMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_dC_out[9] = {
      1.2627241707513939, 0.8070787912646588, 1.2188425316108389,
      1.0685306018376208, 0.6376919120428365, 1.0372393277182761,
      1.7508605712498686, 0.9184531019302811, 1.5696508760099037};
  const T AB_dB_out[9] = {
      2.0413786164135028, 1.1923066634934252, 1.9209506465376314,
      2.9171452067840677, 1.6664270897907307, 2.7109069129005614,
      3.1385206417173066, 1.7834380412688020, 2.9172345118872216};

  // ATB
  const T ATB_dC_out[9] = {
      0.9592867643001749, 0.6768825208302733, 0.9704603540454581,
      1.3389658693972686, 0.7651216894455629, 1.2981211180713395,
      1.5750631386710723, 0.8540740451196530, 1.4640735087850434};
  const T ATB_dB_out[9] = {
      2.4334740643837729, 1.4334937226001876, 2.3434940648557960,
      2.1072825444158747, 1.2446107720007960, 2.0243792223597974,
      3.1316549204160795, 1.8015493497141932, 2.9931635208868759};

  // ABT
  const T ABT_dC_out[9] = {
      0.7069272548381547, 1.3133658139131372, 1.1241747930770374,
      0.5790453975587826, 1.1668340507043797, 1.1030994231821820,
      0.9401740443246649, 1.5917329769339021, 1.3147206569882468};
  const T ABT_dB_out[9] = {
      1.1148914555297598, 1.5893240481467552, 1.7061225876804622,
      2.0531031296356019, 2.8531046074557045, 3.0791938516032742,
      1.7933275685809114, 2.4437636119833894, 2.6570441551059920};

  // ATBT
  const T ATBT_dC_out[9] = {
      0.5478473331041881, 1.0909830415083410, 0.9708422338676705,
      0.7130648169897426, 1.4842075805095507, 1.4584686392366362,
      0.8404854765162666, 1.5784822247362573, 1.4380120451903955};
  const T ATBT_dB_out[9] = {
      1.3172562562444901, 1.1434707482675552, 1.6855674910799765,
      2.6116869869566339, 2.2427823463592396, 3.3221010092849750,
      2.4426735924168170, 2.0756863534277197, 3.1112097546677155};
};

TEST_F(MatxADMat, AB) {
  auto dC_out = AB_dC_out;
  auto dB_out = AB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, false, false>(A, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ATB) {
  auto dC_out = ATB_dC_out;
  auto dB_out = ATB_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, true, false>(A, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ABT) {
  auto dC_out = ABT_dC_out;
  auto dB_out = ABT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, false, true>(A, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

TEST_F(MatxADMat, ATBT) {
  auto dC_out = ATBT_dC_out;
  auto dB_out = ATBT_dB_out;

  Mat A(A_data), B(B_data), C;
  Mat dB(dB_data), dC;
  ADMat B_(B, dB), C_(C, dC);

  auto expr = MatMatMult<T, true, true>(A, B_, C_);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  expr.forward();
  expect_mat_eq<3, 3, Mat>(dC, dC_out);

  // Check reverse AD results
  dB.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dB, dB_out);
}

// Test suite: C = AB, where only A is A2D-able
class A2DMatxMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[9] = {
      1.2682553901952749, 0.9534527931080914, 0.5922214170543438,
      1.5426274087687242, 1.2063452326563935, 0.7028397697553150,
      1.3900890828380486, 1.1928798594613501, 0.7301190962948266};
  const T AB_Ab_out[9] = {
      0.6900085914418043, 0.9798365026866145, 0.9089509143610131,
      0.2861646455269697, 0.4660706571513516, 0.6607633064056966,
      0.7930533953144042, 1.2631307790643582, 1.5425181226326985};
  const T AB_Ah_out[9] = {
      0.5728353773884283, 0.7999861569616865, 1.0221093981423652,
      0.5176327308310544, 0.8357515814026844, 1.2674768978081572,
      0.4617945528372981, 0.6757887974446408, 0.9954030554236757};

  // ATB
  const T ATB_Cp_out[9] = {
      1.0567986362626920, 0.7067854978241817, 0.4310276161733993,
      1.5969704456208476, 1.2239251774552717, 0.7153483104457437,
      1.5436012400805454, 1.2399414794180161, 0.7638132948763544};
  const T ATB_Ab_out[9] = {
      0.6900085914418043, 0.2861646455269697, 0.7930533953144042,
      0.9798365026866145, 0.4660706571513516, 1.2631307790643582,
      0.9089509143610131, 0.6607633064056966, 1.5425181226326985};
  const T ATB_Ah_out[9] = {
      0.5728353773884283, 0.5176327308310544, 0.4617945528372981,
      0.7999861569616865, 0.8357515814026844, 0.6757887974446408,
      1.0221093981423652, 1.2674768978081572, 0.9954030554236757};

  // ABT
  const T ABT_Cp_out[9] = {
      0.6367114191641073, 0.9879649044270280, 1.2363424605432483,
      0.5841650417794853, 1.0303648881695169, 1.5250809794832723,
      0.5271179562342785, 0.8977196826466857, 1.3774608012606941};
  const T ABT_Ab_out[9] = {
      0.9609437297222513, 0.5218227993843698, 0.3311892163094781,
      0.6708487039516081, 0.5536410002458479, 0.3422732243744622,
      1.5852871757937947, 1.0887379356215181, 0.6453448442604295};
  const T ABT_Ah_out[9] = {
      1.0476497656010679, 0.9730030429000408, 0.6731570730670731,
      1.2806022804417001, 1.1589245003223632, 0.7374634008778562,
      1.0085534721981360, 1.0046029473770621, 0.6817069184092223};

  // ATBT
  const T ATBT_Cp_out[9] = {
      0.5862166227930794, 0.8953031106325379, 1.0217044986170061,
      0.6379687309408689, 1.1015478190528851, 1.5745755924835247,
      0.6814231413367092, 1.1014902548345833, 1.5169552773005301};
  const T ATBT_Ab_out[9] = {
      0.9609437297222513, 0.6708487039516081, 1.5852871757937947,
      0.5218227993843698, 0.5536410002458479, 1.0887379356215181,
      0.3311892163094781, 0.3422732243744622, 0.6453448442604295};
  const T ATBT_Ah_out[9] = {
      1.0476497656010679, 1.2806022804417001, 1.0085534721981360,
      0.9730030429000408, 1.1589245003223632, 1.0046029473770621,
      0.6731570730670731, 0.7374634008778562, 0.6817069184092223};
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

  auto expr = MatMatMult<1, T, false, false>(A__, B, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out);
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

  auto expr = MatMatMult<1, T, true, false>(A__, B, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out);
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

  auto expr = MatMatMult<1, T, false, true>(A__, B, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out);
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

  auto expr = MatMatMult<1, T, true, true>(A__, B, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out);
}

// Test suite: C = AB, where only B is A2D-able
class MatxA2DMat : public ADExpressionTest {
 protected:
  // AB
  const T AB_Cp_out[9] = {
      1.2627241707513939, 0.8070787912646588, 1.2188425316108389,
      1.0685306018376208, 0.6376919120428365, 1.0372393277182761,
      1.7508605712498686, 0.9184531019302811, 1.5696508760099037};
  const T AB_Bb_out[9] = {
      1.0484143573871867, 0.6754331993575550, 0.6227343243778110,
      1.5694608536815782, 1.0524416309088469, 0.8895719000108111,
      1.5825824482727067, 1.1549100993633834, 0.9955873604820227};
  const T AB_Bh_out[9] = {
      0.8072960054954835, 0.4209306153925176, 1.4796458636249910,
      1.0725257909330448, 0.4535521056179460, 1.9699406574561971,
      1.1285150142973088, 0.5751893660591434, 2.1474803194132548};

  // ATB
  const T ATB_Cp_out[9] = {
      0.9592867643001749, 0.6768825208302733, 0.9704603540454581,
      1.3389658693972686, 0.7651216894455629, 1.2981211180713395,
      1.5750631386710723, 0.8540740451196530, 1.4640735087850434};
  const T ATB_Bb_out[9] = {
      1.2363986180762427, 0.8664446795424580, 0.7919828503427733,
      1.1932431552396239, 0.8068055019872884, 0.6984977682252416,
      1.4890988399601655, 1.1924451750774940, 1.0546207193046553};
  const T ATB_Bh_out[9] = {
      0.9580706494566084, 0.5419591655534717, 1.8065233589335892,
      0.8463688709455357, 0.3934474867667382, 1.5648289581049559,
      1.1268507661262552, 0.6916297305035152, 2.2264996245332207};

  // ABT
  const T ABT_Cp_out[9] = {
      0.7069272548381547, 1.3133658139131372, 1.1241747930770374,
      0.5790453975587826, 1.1668340507043797, 1.1030994231821820,
      0.9401740443246649, 1.5917329769339021, 1.3147206569882468};
  const T ABT_Bb_out[9] = {
      1.0484143573871867, 1.5694608536815782, 1.5825824482727067,
      0.6754331993575550, 1.0524416309088469, 1.1549100993633834,
      0.6227343243778110, 0.8895719000108111, 0.9955873604820227};
  const T ABT_Bh_out[9] = {
      0.8072960054954835, 1.0725257909330448, 1.1285150142973088,
      0.4209306153925176, 0.4535521056179460, 0.5751893660591434,
      1.4796458636249910, 1.9699406574561971, 2.1474803194132548};

  // ATBT
  const T ATBT_Cp_out[9] = {
      0.5478473331041881, 1.0909830415083410, 0.9708422338676705,
      0.7130648169897426, 1.4842075805095507, 1.4584686392366362,
      0.8404854765162666, 1.5784822247362573, 1.4380120451903955};
  const T ATBT_Bb_out[9] = {
      1.2363986180762427, 1.1932431552396239, 1.4890988399601655,
      0.8664446795424580, 0.8068055019872884, 1.1924451750774940,
      0.7919828503427733, 0.6984977682252416, 1.0546207193046553};
  const T ATBT_Bh_out[9] = {
      0.9580706494566084, 0.8463688709455357, 1.1268507661262552,
      0.5419591655534717, 0.3934474867667382, 0.6916297305035152,
      1.8065233589335892, 1.5648289581049559, 2.2264996245332207};
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

  auto expr = MatMatMult<1, T, false, false>(A, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, AB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out);
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

  auto expr = MatMatMult<1, T, true, false>(A, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATB_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out);
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

  auto expr = MatMatMult<1, T, false, true>(A, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ABT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out);
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

  auto expr = MatMatMult<1, T, true, true>(A, B__, C__);

  // Check expression result
  expect_mat_eq<3, 3, Mat>(C, ATBT_data);

  // Check forward AD result
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      B__.pvalue(0)(i, j) = Bp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out);

  // Check reverse AD results
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Bb, Bb_out);

  // Check reverse A2D results
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(B__.hvalue(0), Bh_out);
}

// Test suite: C = det(A)
class DetA : public ADExpressionTest {
 protected:
  // AD
  const T AD_sb_out = -0.0627964774418258;
  const T AD_Ab_out[9] = {
      0.0105328892254463,  0.0152249197717527,  -0.0188719858688089,
      0.0095245461075082,  -0.0166479481614073, 0.0110810249822685,
      -0.0129720568863309, 0.0016352064483511,  0.0098708014948920};

  // A2D
  const T A2D_sp_out = -0.0627964774418258;
  const T A2D_Ab_out[9] = {
      -0.0563949171388059, -0.0815168630938085, 0.1010438880098533,
      -0.0509960730641692, 0.0891360040914050,  -0.0593297311224298,
      0.0694546441689654,  -0.0087551791522515, -0.0528499845088277};
  const T A2D_Ah_out[9] = {
      -0.0601297024768627, -0.2888078775111398, 0.2886408254311685,
      -0.1519649990997223, 0.4036519783269401,  -0.2971582469942756,
      0.1532281527518056,  -0.0952153063131011, -0.0416059848892236};
};

TEST_F(DetA, A) {
  Mat A(A_data);
  T det;
  MatDet(A, det);
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

  auto expr = MatDet(A_, det);

  // Check expression result
  expect_val_eq(det.value, detA_data);

  // Check forward AD result
  expr.forward();
  expect_val_eq(det.bvalue, sb_out);

  // Check reverse AD results
  dA.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(dA, Ab_out);
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

  auto expr = MatDet(A__, s__);

  // Check expression result
  expect_val_eq(s__.value, detA_data);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(s__.pvalue[0], sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out);

  // Check hreverse
  s__.hvalue[0] = sh;
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out, 1e-8);
}

// Test suite: C = inv(A)
class InvA : public ADExpressionTest {
 protected:
  // AD
  const T AD_invAb_out[9] = {
      -3.6897999072316887, 0.1064624509411591,  1.7510230980042913,
      1.7892774874694848,  -5.0558603245912472, 2.4570181642106101,
      0.2287635114258432,  4.3701340615693320,  -3.9780303671725092};
  const T AD_Ab_out[9] = {
      120.3973246131449883, -93.0289820006620118, 5.7545319340518617,
      -53.5731920528893255, 117.2485960258310485, -70.7854458380928122,
      -58.1465838256700067, -5.5434985887259032,  42.8953522403000278};

  // A2D
  const T A2D_Ab_out[9] = {
      2.3410906920289722,  -3.3459619930684923, 2.5452604094640368,
      -6.4214289213302544, -4.2652437638722427, 6.4126931000537191,
      1.0519501300050025,  4.8585038067463699,  -5.6438538826213520};
  const T A2D_Ah_out[9] = {
      4.2846057796144192, 2.3423048586980575,  -6.0023168756884964,
      4.3152713577318380, -7.3351514870213510, 3.5374811782671349,
      0.1741911137413386, 4.6233386521636808,  -2.1532599075282208};
  const T A2D_Cp_out[9] = {
      -3.6897999072316887, 0.1064624509411591,  1.7510230980042913,
      1.7892774874694848,  -5.0558603245912472, 2.4570181642106101,
      0.2287635114258432,  4.3701340615693320,  -3.9780303671725092};
};

TEST_F(InvA, A) {
  Mat A(A_data), invA;
  MatInverse(A, invA);
  expect_mat_eq<3, 3, Mat>(invA, invA_data, 1e-14);
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
  auto expr = MatInverse(A_, invA_);
  expect_mat_eq<3, 3, Mat>(invA, invA_data, 1e-14);

  // Check forward
  expr.forward();
  expect_mat_eq<3, 3, Mat>(invAb, invAb_out, 1e-14);

  // Check reverse
  Ab.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out, 1e-12);
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

  auto expr = MatInverse(A__, C__);

  // Check exprssion
  expect_mat_eq<3, 3, Mat>(C, invA_data, 1e-14);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      A__.pvalue(0)(i, j) = Ap(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, Mat>(C__.pvalue(0), Cp_out, 1e-14);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Ab, Ab_out, 1e-13);

  // Check hreverse
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      C__.hvalue(0)(i, j) = Ch(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(A__.hvalue(0), Ah_out,
                           1e-5);  // TODO: double-check this tolerance
}

// Test suite: trace = tr(S)
class SymmTraceSuite : public ADExpressionTest {
 protected:
  const T tr_out = 1.2536821800000000;

  // AD
  const T AD_trb_out = 1.0853507099999999;
  const T AD_Sb_out[9] = {
      1.0853507099999999, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 1.0853507099999999, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 1.0853507099999999};

  // A2D
  const T A2D_Sb_out[9] = {
      0.3362232400000000, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 0.3362232400000000, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 0.3362232400000000};
  const T A2D_Sh_out[9] = {
      0.7315793000000000, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 0.7315793000000000, 0.0000000000000000,
      0.0000000000000000, 0.0000000000000000, 0.7315793000000000};
  const T A2D_trp_out = 1.0853507099999999;
};

TEST_F(SymmTraceSuite, A) {
  // Set inputs
  const T *S_in = S_data;
  SMat S(S_in);

  // Set outputs
  T trace;

  // Compute
  SymmTrace(S, trace);
  expect_val_eq(trace, tr_out);
}

TEST_F(SymmTraceSuite, AD) {
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

  auto expr = SymmTrace(S_, trace_);

  // Check expression result
  expect_val_eq(trace_.value, tr_out);

  // Check forward
  expr.forward();
  expect_val_eq(trace_.bvalue, trb_out);

  // print_mat<3, 3, SMat>(Sb);
  // print_mat<3, 3, SMat>(Eb);

  // Check reverse
  Sb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Sb, Sb_out);
}

TEST_F(SymmTraceSuite, A2D) {
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

  auto expr = SymmTrace(S__, trace__);

  // Check expression
  expect_val_eq(trace__.value, tr_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = i; j < 3; j++) {
      S__.pvalue(0)(i, j) = Sp(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(trace__.pvalue[0], trp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Sb, Sb_out);

  // Check hreverse
  trace__.hvalue[0] = trh;
  expr.hreverse();
  expect_mat_eq<3, 3, SMat>(S__.hvalue(0), Sh_out);
}

// Test suite: trace = tr(SE)
class SymmMultTrace : public ADExpressionTest {
 protected:
  // AD
  const T AD_Cb_out = 3.3727867595096575;
  const T AD_Sb_out[9] = {
      0.9493511057888694, 3.1857326803846040, 6.2423305971833027,
      3.1857326803846040, 1.3443179602023558, 0.0215972352801089,
      6.2423305971833027, 0.0215972352801089, 0.5320555935586067};
  const T AD_Eb_out[9] = {
      0.5620927892421098, 0.6466343201428918, 1.0872219267277847,
      0.6466343201428918, 0.9976895820753535, 1.5869412307804009,
      1.0872219267277847, 1.5869412307804009, 2.6686202860197401};

  // A2D
  const T A2D_Sb_out[9] = {
      0.0946380330111120, 0.3175763663542480, 0.6222796661005735,
      0.3175763663542480, 0.1340111226703044, 0.0021529651705512,
      0.6222796661005735, 0.0021529651705512, 0.0530390648095420};
  const T A2D_Sh_out[9] = {
      0.5001246470729986, 0.9716939742924943, 1.9218336982878552,
      0.9716939742924943, 0.5136586894799157, 0.4444080657813740,
      1.9218336982878552, 0.4444080657813740, 0.2015086087059534};
  const T A2D_Eb_out[9] = {
      0.0560333849291720, 0.0644610826938992, 0.1083819716063528,
      0.0644610826938992, 0.0994567542267004, 0.1581975263624864,
      0.1083819716063528, 0.1581975263624864, 0.2660269453339908};
  const T A2D_Eh_out[9] = {
      0.1706866501072147, 0.5891511705403434, 0.6504501627130445,
      0.5891511705403434, 0.4981935844741603, 1.0113339051514922,
      0.6504501627130445, 1.0113339051514922, 0.6132078294417402};
  const T A2D_sp_out = 3.3727867595096575;
};

TEST_F(SymmMultTrace, A) {
  // Set inputs
  const T *S_in = S_data;
  const T *E_in = E_data;
  SMat S(S_in), E(E_in);

  // Set outputs
  T trace;

  // Compute
  SymmSymmMultTrace(S, E, trace);
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

  auto expr = SymmSymmMultTrace(S_, E_, trace_);

  // Check expression result
  expect_val_eq(trace_.value, trSE);

  // Check forward
  expr.forward();
  expect_val_eq(trace_.bvalue, Cb_out);

  // print_mat<3, 3, SMat>(Sb);
  // print_mat<3, 3, SMat>(Eb);

  // Check reverse
  Sb.zero();
  Eb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Sb, Sb_out);
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);
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

  auto expr = SymmSymmMultTrace(S__, E__, trace__);

  // Check expression
  expect_val_eq(trace__.value, trSE);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = i; j < 3; j++) {
      S__.pvalue(0)(i, j) = Sp(i, j);
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(trace__.pvalue[0], sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Sb, Sb_out);
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);

  // Check hreverse
  trace__.hvalue[0] = sh;
  expr.hreverse();
  expect_mat_eq<3, 3, SMat>(S__.hvalue(0), Sh_out, 1e-9);
  expect_mat_eq<3, 3, SMat>(E__.hvalue(0), Eh_out, 1e-9);
}

// Test suite: symmetric isotropic constitutive, S = 2*mu*E + lambda * tr(E) * I
class SymmIsoConstitutive : public ADExpressionTest {
 protected:
  const T S_out[9] = {
      1.1160684307960387, 0.7935115845831120, 1.5548579057375782,
      0.7935115845831120, 1.3128274221785898, 0.0053795023340328,
      1.5548579057375782, 0.0053795023340328, 0.9081860650293788};

  // AD
  const T AD_Sb_out[9] = {
      2.8455308981095779, 0.7013403281847024, 1.4188117200539616,
      0.7013403281847024, 2.4850418052328012, 1.0987142518736353,
      1.4188117200539616, 1.0987142518736353, 1.8055838729023441};
  const T AD_Eb_out[9] = {
      10.2591185450867570, 1.1783970129590735, 2.3838975539741130,
      1.1783970129590735,  9.6534222045528555, 1.8460675088435552,
      2.3838975539741130,  1.8460675088435552, 8.5117921489107431};

  // A2D
  const T A2D_Sp_out[9] = {
      2.8455308981095779, 0.7013403281847024, 1.4188117200539616,
      0.7013403281847024, 2.4850418052328012, 1.0987142518736353,
      1.4188117200539616, 1.0987142518736353, 1.8055838729023441};
  const T A2D_Eb_out[9] = {
      1.0768579649748915, 1.1216231366637024, 1.0360011764004120,
      1.1216231366637024, 2.2413437565402141, 1.6668883029825432,
      1.0360011764004120, 1.6668883029825432, 1.0049054809550020};
  const T A2D_Eh_out[9] = {
      1.7620483650844783, 1.2274571484770327, 1.1565054126198864,
      1.2274571484770327, 1.6350346200344863, 0.9219831109345583,
      1.1565054126198864, 0.9219831109345583, 1.2098541926234696};
};

TEST_F(SymmIsoConstitutive, A) {
  // Set inputs
  const T *E_in = E_data;
  SMat E(E_in);

  // Set outputs
  SMat S;

  // Compute
  SymmIsotropicConstitutive(mu, lam, E, S);
  expect_mat_eq<3, 3, SMat>(S, S_out);
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

  auto expr = SymmIsotropicConstitutive(mu, lam, E_, S_);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(S, S_out);

  // Check forward AD
  expr.forward();
  expect_mat_eq<3, 3, SMat>(Sb, Sb_out);

  // Check reverse AD
  Eb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out, 1e-14);
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

  auto expr = SymmIsotropicConstitutive(mu, lam, E__, S__);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(S, S_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, SMat>(S__.pvalue(0), Sp_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);

  // Check hreverse
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      S__.hvalue(0)(i, j) = Sh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(E__.hvalue(0), Eh_out);
}

// Test suite: symmetric isotropic energy w = f(E)
// S = 2*mu*E + lambda * tr(E) * I
// W = 1/2 * tr(ES)
class SymmIsoEnergy : public ADExpressionTest {
 protected:
  const T w_out = 2.3039661138033614;

  // AD_E_only
  const T AD_E_only_wb_out = 5.3716615759793331;
  const T AD_E_only_Eb_out[9] = {
      5.9951419058706303,  8.5249513779991553, 16.7043409367164912,
      8.5249513779991553,  7.0520646196087293, 0.0577937319712303,
      16.7043409367164912, 0.0577937319712303, 4.8784681893581814};

  // AD_E_mu_lam
  const T AD_E_mu_lam_wb_out = 6.5702621870009370;
  const T AD_E_mu_lam_mub_out = 15.9118143886609182;
  const T AD_E_mu_lam_lamb_out = 2.3058689403874317;
  const T AD_E_mu_lam_Eb_out[9] = {
      7.3328622089646851,  10.4271583182672334, 20.4316482084541562,
      10.4271583182672334, 8.6256203699979039,  0.0706894815403580,
      20.4316482084541562, 0.0706894815403580,  5.9670205618237011};

  // A2D_E_only
  const T A2D_E_only_wp_out = 5.3716615759793331;
  const T A2D_E_only_Eb_out[9] = {
      0.3752481438639599, 0.5335940718921359, 1.0455587256134062,
      0.5335940718921359, 0.4414030894457333, 0.0036174274086721,
      1.0455587256134062, 0.0036174274086721, 0.3053532613070284};
  const T A2D_E_only_Eh_out[9] = {
      1.7732261794453628, 1.6326471345872511, 3.2290786693440907,
      1.6326471345872511, 1.7959661710272448, 0.7466975962914767,
      3.2290786693440907, 0.7466975962914767, 1.2714893852494025};

  // A2D_E_mu_lam
  const T A2D_E_mu_lam_wp_out = 6.5702621870009370;
  const T A2D_E_mu_lam_Eb_out[9] = {
      0.3752481438639599, 0.5335940718921359, 1.0455587256134062,
      0.5335940718921359, 0.4414030894457333, 0.0036174274086721,
      1.0455587256134062, 0.0036174274086721, 0.3053532613070284};
  const T A2D_E_mu_lam_Eh_out[9] = {
      1.9446458871820334, 1.9188553725435744, 3.7898935972354293,
      1.9188553725435744, 2.0028699512248536, 0.7486379055661978,
      3.7898935972354293, 0.7486379055661978, 1.4054190018701651};
  const T A2D_E_mu_lam_mub_out = 0.8142630591848906;
  const T A2D_E_mu_lam_muh_out = 3.4604184768609523;
  const T A2D_E_mu_lam_lamb_out = 0.1179993589428304;
  const T A2D_E_mu_lam_lamh_out = 0.7614213002495509;
};

TEST_F(SymmIsoEnergy, A) {
  // Set inputs
  const T *E_in = E_data;
  SMat E(E_in);

  // Set outputs
  T w;

  // Compute
  SymmIsotropicEnergy(mu, lam, E, w);
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
  auto expr = SymmIsotropicEnergy(mu, lam, E_, w_);

  // Check expression
  expect_val_eq(w_.value, w_out);

  // Check forward
  expr.forward();
  expect_val_eq(w_.bvalue, wb_out);

  // Check reverse
  Eb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);
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
  auto expr = SymmIsotropicEnergy(mu_, lam_, E_, w_);

  // Check expression
  expect_val_eq(w_.value, w_out);

  // Check forward
  expr.forward();
  expect_val_eq(w_.bvalue, wb_out);

  // Check reverse
  Eb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);
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
  auto expr = SymmIsotropicEnergy(mu, lam, E__, w__);

  // Check expression
  expect_val_eq(w__.value, w_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      E__.pvalue(0)(i, j) = Ep(i, j);
    }
  }
  expr.hforward();
  expect_val_eq(w__.pvalue[0], wp_out);

  // Check hreverse
  w__.hvalue[0] = wh;
  expr.hreverse();
  expect_mat_eq<3, 3, SMat>(E__.hvalue(0), Eh_out, 1e-8);
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
  auto expr = SymmIsotropicEnergy(mu__, lam__, E__, w__);

  // Check expression
  expect_val_eq(w__.value, w_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);
  expect_val_eq(mu__.bvalue, mub_out);
  expect_val_eq(lam__.bvalue, lamb_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
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
  expect_mat_eq<3, 3, SMat>(E__.hvalue(0), Eh_out, 1e-7);
  expect_val_eq(mu__.hvalue[0], muh_out, 1e-7);
  expect_val_eq(lam__.hvalue[0], lamh_out, 1e-7);
}

class GreenStrain : public ADExpressionTest {
 protected:
  const T E_out[9] = {
      0.9436017519923732, 1.1368244824910445, 1.0723893026894997,
      1.1368244824910445, 1.1667739540265185, 1.5508808212339893,
      1.0723893026894997, 1.5508808212339893, 1.8182371643062669};

  // AD
  const T AD_Eb_out[9] = {
      1.5055392317465364, 1.6057890503612342, 1.6609119355915014,
      1.6057890503612342, 2.4575696051681013, 2.3688299685949716,
      1.6609119355915014, 2.3688299685949716, 2.6109853227747908};
  const T AD_Uxb_out[9] = {
      3.4065895626387768, 3.7150835165849525, 3.7071087588020268,
      2.4997742670895411, 4.7012211489701663, 3.8251196074462523,
      3.0055398291496793, 4.8687221917585912, 6.5467208590069497};

  // A2D
  const T A2D_Uxb_out[9] = {
      0.6491783791618027, 1.4153475589445037, 0.8938635694579302,
      0.7533357747742581, 1.6954188087124340, 0.9401895193906458,
      0.9665071916249530, 1.8674980606845262, 0.7779725280334064};
  const T A2D_Uxh_out[9] = {
      1.7572324071129239, 2.0582978076623135, 1.3940450133788342,
      1.6508443808732531, 2.3585478594177447, 1.4088143291783681,
      1.8095439459254061, 2.2938438511143504, 1.3570465181240248};
  const T A2D_Ep_out[9] = {
      1.5055392317465364, 1.6057890503612342, 1.6609119355915014,
      1.6057890503612342, 2.4575696051681013, 2.3688299685949716,
      1.6609119355915014, 2.3688299685949716, 2.6109853227747908};
};

TEST_F(GreenStrain, A) {
  // Set inputs
  const T *Ux_in = A_data;
  Mat Ux(Ux_in);

  // Set outputs
  SMat E;

  // Compute
  MatGreenStrain(Ux, E);
  expect_mat_eq<3, 3, SMat>(E, E_out);
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
  auto expr = MatGreenStrain(Ux_, E_);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(E, E_out);

  // Check forward
  expr.forward();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);

  // Check reverse
  Uxb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Uxb, Uxb_out, 1e-14);
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
  auto expr = MatGreenStrain(Ux__, E__);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(E, E_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      Ux__.pvalue(0)(i, j) = Uxp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, SMat>(E__.pvalue(0), Ep_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Uxb, Uxb_out);

  // Check hreverse
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      E__.hvalue(0)(i, j) = Eh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(Ux__.hvalue(0), Uxh_out, 1e-8);
}

class LinearGreenStrain : public ADExpressionTest {
 protected:
  const T E_out[9] = {
      0.5488135000000000, 0.6300362750000000, 0.5201752950000000,
      0.6300362750000000, 0.4236548000000000, 0.7688335550000001,
      0.5201752950000000, 0.7688335550000001, 0.9636627600000000};

  // AD
  const T AD_Eb_out[9] = {
      0.6962679100000000, 0.4974890850000000, 0.5316012349999999,
      0.4974890850000000, 0.9893170799999998, 0.7428553250000000,
      0.5316012349999999, 0.7428553250000000, 0.8848385800000000};
  const T AD_Uxb_out[9] = {
      0.6962679100000000, 0.2487445425000000, 0.2658006174999999,
      0.2487445425000000, 0.9893170799999998, 0.3714276625000000,
      0.2658006174999999, 0.3714276625000000, 0.8848385800000000};

  // A2D
  const T A2D_Uxb_out[9] = {
      0.1450378100000000, 0.3337752600000000, 0.3082956750000000,
      0.3337752600000000, 0.8380986700000000, 0.4960365550000000,
      0.3082956750000000, 0.4960365550000000, 0.1022142300000000};
  const T A2D_Uxh_out[9] = {
      0.5202830300000000, 0.3652695950000000, 0.3441556100000000,
      0.3652695950000000, 0.4446889300000000, 0.2743659100000000,
      0.3441556100000000, 0.2743659100000000, 0.1916365400000000};
  const T A2D_Ep_out[9] = {
      0.6962679100000000, 0.4974890850000000, 0.5316012349999999,
      0.4974890850000000, 0.9893170799999998, 0.7428553250000000,
      0.5316012349999999, 0.7428553250000000, 0.8848385800000000};
};

TEST_F(LinearGreenStrain, A) {
  // Set inputs
  const T *Ux_in = A_data;
  Mat Ux(Ux_in);

  // Set outputs
  SMat E;

  // Compute
  MatLinearGreenStrain(Ux, E);
  expect_mat_eq<3, 3, SMat>(E, E_out);
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
  auto expr = MatLinearGreenStrain(Ux_, E_);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(E, E_out);

  // Check forward
  expr.forward();
  expect_mat_eq<3, 3, SMat>(Eb, Eb_out);

  // Check reverse
  Uxb.zero();
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Uxb, Uxb_out, 1e-14);
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
  auto expr = MatLinearGreenStrain(Ux__, E__);

  // Check expression result
  expect_mat_eq<3, 3, SMat>(E, E_out);

  // Check hforward
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      Ux__.pvalue(0)(i, j) = Uxp(i, j);
    }
  }
  expr.hforward();
  expect_mat_eq<3, 3, SMat>(E__.pvalue(0), Ep_out);

  // Check reverse
  expr.reverse();
  expect_mat_eq<3, 3, Mat>(Uxb, Uxb_out);

  // Check hreverse
  for (I i = 0; i < 3; i++) {
    for (I j = 0; j < 3; j++) {
      E__.hvalue(0)(i, j) = Eh(i, j);
    }
  }
  expr.hreverse();
  expect_mat_eq<3, 3, Mat>(Ux__.hvalue(0), Uxh_out, 1e-8);
}