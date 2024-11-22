#include <string>

#include "a2dtuple.h"
#include "ad/a2dvartuple.h"
#include "test_commons.h"

TEST(test_a2dtuple, tuple_of_objs_default_construction) {
  A2D::a2d_tuple<int, std::string, double> t;

  EXPECT_EQ(A2D::a2d_get<0>(t), 0);
  EXPECT_EQ(A2D::a2d_get<1>(t), "");
  EXPECT_DOUBLE_EQ(A2D::a2d_get<2>(t), 0.0);
}

TEST(test_a2dtuple, tuple_of_objs) {
  int x = 5;
  std::string s = "hello";
  double y = 42.1;

  auto t = A2D::a2d_tuple<int, std::string, double>{x, s, y};

  EXPECT_EQ(A2D::a2d_get<0>(t), 5);
  EXPECT_EQ(A2D::a2d_get<1>(t), "hello");
  EXPECT_DOUBLE_EQ(A2D::a2d_get<2>(t), 42.1);
}

TEST(test_a2dtuple, tuple_of_refs) {
  int x = 5;
  std::string s = "hello";
  double y = 42.1;

  auto t = A2D::a2d_tuple<int&, std::string&, double&>{x, s, y};

  EXPECT_EQ(A2D::a2d_get<0>(t), 5);
  EXPECT_EQ(A2D::a2d_get<1>(t), "hello");
  EXPECT_DOUBLE_EQ(A2D::a2d_get<2>(t), 42.1);

  A2D::a2d_get<0>(t) += 2;

  EXPECT_EQ(A2D::a2d_get<0>(t), 7);
  EXPECT_EQ(x, 7);
}

TEST(test_a2dtuple, TieTuple) {
  using T = double;

  A2D::Mat<T, 3, 3> Ux;
  A2D::A2DObj<A2D::Mat<T, 3, 3>> Uxb;

  auto g = A2D::MakeTieTuple<T, A2D::ADseed::b>(Uxb);
  g.zero();
}
