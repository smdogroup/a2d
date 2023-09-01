#include "multiphysics/feelementvector.h"
#include "test_commons.h"

using namespace A2D;

struct Basis {
  static constexpr index_t ndof = 42;
};
struct VecType;

using EVempty = ElementVector_Empty;
using EVserial = ElementVector_Serial<T, Basis, VecType>;
using EVparallel = ElementVector_Parallel<T, Basis, VecType>;

TEST(HaveSameEvtype_Value, ZeroArg) {
  EXPECT_TRUE((have_same_evtype<>::value));
}

TEST(HaveSameEvtype_Value, OneArg) {
  EXPECT_TRUE((have_same_evtype<EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel>::value));
}

TEST(HaveSameEvtype_Value, TwoArgs) {
  EXPECT_TRUE((have_same_evtype<EVserial, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel>::value));
  EXPECT_TRUE((have_same_evtype<EVserial, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVparallel>::value));
  EXPECT_FALSE((have_same_evtype<EVserial, EVparallel>::value));
}

TEST(HaveSameEvtype_Value, ThreeArgs) {
  EXPECT_TRUE((have_same_evtype<EVserial, EVserial, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVserial, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVserial, EVempty, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVserial, EVserial, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVserial, EVempty, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVserial, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVserial>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVempty>::value));

  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel, EVparallel>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVparallel, EVparallel>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVempty, EVparallel>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVempty, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVparallel, EVempty>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVparallel>::value));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVempty>::value));

  EXPECT_FALSE((have_same_evtype<EVserial, EVserial, EVparallel>::value));
  EXPECT_FALSE((have_same_evtype<EVserial, EVparallel, EVserial>::value));
  EXPECT_FALSE((have_same_evtype<EVparallel, EVserial, EVserial>::value));

  EXPECT_FALSE((have_same_evtype<EVserial, EVempty, EVparallel>::value));
  EXPECT_FALSE((have_same_evtype<EVempty, EVserial, EVparallel>::value));
  EXPECT_FALSE((have_same_evtype<EVempty, EVparallel, EVserial>::value));
  EXPECT_FALSE((have_same_evtype<EVserial, EVparallel, EVempty>::value));
  EXPECT_FALSE((have_same_evtype<EVparallel, EVempty, EVserial>::value));
  EXPECT_FALSE((have_same_evtype<EVparallel, EVserial, EVempty>::value));
}

TEST(HaveSameEvtype_Evtype, ZeroArg) {
  EXPECT_TRUE((have_same_evtype<>::evtype == ElemVecType::Empty));
}

TEST(HaveSameEvtype_Evtype, OneArg) {
  EXPECT_TRUE((have_same_evtype<EVserial>::evtype == ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVparallel>::evtype == ElemVecType::Parallel));
}

TEST(HaveSameEvtype_Evtype, TwoArg) {
  EXPECT_TRUE(
      (have_same_evtype<EVserial, EVserial>::evtype == ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE(
      (have_same_evtype<EVserial, EVempty>::evtype == ElemVecType::Serial));
  EXPECT_TRUE(
      (have_same_evtype<EVparallel, EVempty>::evtype == ElemVecType::Parallel));
  EXPECT_TRUE(
      (have_same_evtype<EVempty, EVserial>::evtype == ElemVecType::Serial));
  EXPECT_TRUE(
      (have_same_evtype<EVempty, EVparallel>::evtype == ElemVecType::Parallel));
}

TEST(HaveSameEvtype_Evtype, ThreeArgs) {
  EXPECT_TRUE((have_same_evtype<EVserial, EVserial, EVserial>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVempty, EVserial, EVserial>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVserial, EVempty, EVserial>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVserial, EVserial, EVempty>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVserial, EVempty, EVempty>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVempty, EVserial, EVempty>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVserial>::evtype ==
               ElemVecType::Serial));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVempty>::evtype ==
               ElemVecType::Empty));

  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel, EVparallel>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVempty, EVparallel, EVparallel>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVempty, EVparallel>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVparallel, EVempty>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVparallel, EVempty, EVempty>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVempty, EVparallel, EVempty>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVparallel>::evtype ==
               ElemVecType::Parallel));
  EXPECT_TRUE((have_same_evtype<EVempty, EVempty, EVempty>::evtype ==
               ElemVecType::Empty));
}
