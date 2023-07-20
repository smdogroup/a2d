/*
  This is a set of unit tests for a2dvecops.h using Google Test framework.
 */

#include <gtest/gtest.h>

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2dvecops3d.h"
#include "test_commons.h"

using T = double;
using I = int;
using Vec_t = A2D::Vec<T, 3>;
using ADVec_t = A2D::ADVec<A2D::Vec<T, 3>>;
using ADScalar_t = A2D::ADScalar<T>;

class VecOpsTest : public ::testing::Test {
 protected:
  const T x_data[3] = {0.96155402, 0.02855176, 0.95787560};
  const T y_data[3] = {0.80766462, 0.60212270, 0.86418474};
  const T alpha = 0.33622324;
  const T scale = 0.45061325;

  const T axpy_out[3] = {1.1309614280394249, 0.6117224652549024,
                         1.1862447777489440};
  const T saxpy_out[3] = {0.9533464453852714, 0.6064484814207487,
                          1.0093092603051743};
};

class Axpy : public VecOpsTest {
 protected:
  // Forward
  const T ab_in = 0.73157930;
  const T xb_in[3] = {0.16665530, 0.09586054, 0.29580571};
  const T yb_in[3] = {0.28147380, 0.47227010, 0.39857781};

  const T AXX_axpyb_out[3] = {0.7034530168637860, 0.0208878765945680,
                              0.7007619609350800};
  const T AXX_saxpyb_out[3] = {0.3169852501512955, 0.0094123539578772,
                               0.3157726246933295};

  const T AAX_axpyb_out[3] = {0.7594864017929580, 0.0531184179415176,
                              0.8002187151617804};
  const T AAX_saxpyb_out[3] = {0.3422346358427307, 0.0239358629434856,
                               0.3605891559498742};

  const T XAA_axpyb_out[3] = {0.3375071849291720, 0.5045006413469496,
                              0.4980345642267003};
  const T XAA_saxpyb_out[3] = {0.3067231856914352, 0.4867936089856084,
                               0.4433943412565447};

  const T AAA_axpyb_out[3] = {1.0409602017929580, 0.5253885179415176,
                              1.1987965251617805};
  const T AAA_saxpyb_out[3] = {0.6237084358427307, 0.4962059629434856,
                               0.7591669659498742};
  // Reverse
  const T axpyb_in[3] = {0.14503781, 0.66755052, 0.83809867};

  const T ab_out = 0.9613156978778634;
  const T xb_out[3] = {0.0487650824007044, 0.2244459986980848,
                       0.2817882502670908};
  const T yb_out[3] = {0.1450378100000000, 0.6675505200000000,
                       0.8380986700000000};

  const T sab_out = 0.4331815908967622;
  const T sxb_out[3] = {0.0219741922670992, 0.1011383409228398,
                        0.1269775192646672};
  const T syb_out[3] = {0.1450378100000000, 0.6675505200000000,
                        0.8380986700000000};
};

TEST_F(Axpy, all_tests) {
  /*
    Allocate
  */
  Vec_t x, y, axpy, saxpy;
  Vec_t xb, yb, axpy_b, saxpy_b;
  ADVec_t xobj(x, xb), yobj(y, yb), axpy_obj(axpy, axpy_b),
      saxpy_obj(saxpy, saxpy_b);

  ADScalar_t aobj;

  /*
    Initialize
  */
  aobj.value = alpha;
  for (int i = 0; i != 3; i++) {
    x(i) = x_data[i];
    y(i) = y_data[i];
  }

  /*
    XXX (passive)
  */
  A2D::Vec3Axpy(alpha, x, y, axpy);
  A2D::Vec3Axpy(scale, alpha, x, y, saxpy);
  EXPECT_VEC_NEAR(3, axpy, axpy_out);
  EXPECT_VEC_NEAR(3, saxpy, saxpy_out);

  /*
    Expression (axpy) and forward (axpyb)
  */
  aobj.bvalue = ab_in;
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = xb_in[i];
    yobj.bvalue()(i) = yb_in[i];
  }

  // AXX
  auto axx = A2D::Vec3Axpy(aobj, x, y, axpy_obj);
  axx.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), axpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AXX_axpyb_out);

  auto saxx = A2D::Vec3Axpy(scale, aobj, x, y, axpy_obj);
  saxx.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), saxpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AXX_saxpyb_out);

  // AAX
  auto aax = A2D::Vec3Axpy(aobj, xobj, y, axpy_obj);
  aax.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), axpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AAX_axpyb_out);

  auto saax = A2D::Vec3Axpy(scale, aobj, xobj, y, axpy_obj);
  saax.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), saxpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AAX_saxpyb_out);

  // XAA
  auto xaa = A2D::Vec3Axpy(alpha, xobj, yobj, axpy_obj);
  xaa.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), axpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), XAA_axpyb_out);

  auto sxaa = A2D::Vec3Axpy(scale, alpha, xobj, yobj, axpy_obj);
  sxaa.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), saxpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), XAA_saxpyb_out);

  // AAA
  auto aaa = A2D::Vec3Axpy(aobj, xobj, yobj, axpy_obj);
  aaa.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), axpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AAA_axpyb_out);

  auto saaa = A2D::Vec3Axpy(scale, aobj, xobj, yobj, axpy_obj);
  saaa.forward();
  EXPECT_VEC_NEAR(3, axpy_obj.value(), saxpy_out);
  EXPECT_VEC_NEAR(3, axpy_obj.bvalue(), AAA_saxpyb_out);

  /*
    Reverse
  */

  // Set reverse seed
  for (int i = 0; i != 3; i++) {
    axpy_obj.bvalue()(i) = axpyb_in[i];
  }

  // AXX
  aobj.bvalue = 0.0;
  axx.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, ab_out);

  aobj.bvalue = 0.0;
  saxx.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, sab_out);

  // AAX
  aobj.bvalue = 0.0;
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
  }
  aax.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, ab_out);
  EXPECT_VEC_NEAR(3, xobj.bvalue(), xb_out);

  aobj.bvalue = 0.0;
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
  }
  saax.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, sab_out);
  EXPECT_VEC_NEAR(3, xobj.bvalue(), sxb_out);

  // XAA
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
    yobj.bvalue()(i) = 0.0;
  }
  xaa.reverse();
  EXPECT_VEC_NEAR(3, xobj.bvalue(), xb_out);
  EXPECT_VEC_NEAR(3, yobj.bvalue(), yb_out);

  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
    yobj.bvalue()(i) = 0.0;
  }
  sxaa.reverse();
  EXPECT_VEC_NEAR(3, xobj.bvalue(), sxb_out);
  EXPECT_VEC_NEAR(3, yobj.bvalue(), syb_out);

  // AAA
  aobj.bvalue = 0.0;
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
    yobj.bvalue()(i) = 0.0;
  }
  aaa.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, ab_out);
  EXPECT_VEC_NEAR(3, xobj.bvalue(), xb_out);
  EXPECT_VEC_NEAR(3, yobj.bvalue(), yb_out);

  aobj.bvalue = 0.0;
  for (int i = 0; i != 3; i++) {
    xobj.bvalue()(i) = 0.0;
    yobj.bvalue()(i) = 0.0;
  }
  saaa.reverse();
  EXPECT_VAL_NEAR(aobj.bvalue, sab_out);
  EXPECT_VEC_NEAR(3, xobj.bvalue(), sxb_out);
  EXPECT_VEC_NEAR(3, yobj.bvalue(), syb_out);
}

class Cross : public VecOpsTest {
 protected:
  // x cross y
  const T cross_out[3] = {-0.5520846472439777, -0.0573180782883828,
                          0.5559132563275228};
};

TEST_F(Cross, passive) {
  Vec_t x(x_data), y(y_data), cross;
  A2D::Vec3Cross(x, y, cross);
  EXPECT_VEC_NEAR(3, cross, cross_out);
}

TEST_F(Cross, forward) {
  const T xb_data[3] = {0.69626791, 0.53987667, 0.65374594};
  const T yb_data[3] = {0.14453098, 0.50253050, 0.44495442};
  const T crossb_out[3] = {-0.3957395632058430, -0.3631016486487311,
                           0.4622830369642267};
  Vec_t x(x_data), xb(xb_data), y(y_data), yb(yb_data);
  Vec_t cross, crossb;
  ADVec_t x_ad(x, xb), y_ad(y, yb), cross_ad(cross, crossb);
  auto expr = A2D::Vec3Cross(x_ad, y_ad, cross_ad);
  expr.forward();
  EXPECT_VEC_NEAR(3, cross, cross_out);
  EXPECT_VEC_NEAR(3, crossb, crossb_out);
}

TEST_F(Cross, reverse) {
  const T vb_data[3] = {0.95793759, 0.25352624, 0.24823463};
  const T xb_out[3] = {-0.0696258021484766, 0.6273447190405860,
                       -0.3720317938326643};
  const T yb_out[3] = {0.2357590636762953, -0.6788940374040915,
                       -0.2164283710828264};
  Vec_t x(x_data), y(y_data), cross;
  Vec_t xb, yb, crossb(vb_data);
  ADVec_t x_ad(x, xb), y_ad(y, yb), cross_ad(cross, crossb);
  auto expr = A2D::Vec3Cross(x_ad, y_ad, cross_ad);
  expr.reverse();
  EXPECT_VEC_NEAR(3, cross, cross_out);
  EXPECT_VEC_NEAR(3, xb, xb_out);
  EXPECT_VEC_NEAR(3, yb, yb_out);
}
