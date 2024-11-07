#include <gtest/gtest.h>

#include <any>

#include "a2ddefs.h"
#include "ad/a2dgemm.h"
#include "ad/a2dmat.h"
#include "test_commons.h"

using namespace A2D;

template <bool additive, MatOp opA, MatOp opB, class Atype, class Btype,
          class Ctype>
void test_gemm_core(Atype A, Btype B, Ctype C, bool scale) {
  using T = typename get_object_numeric_type<Atype>::type;

  T alpha = scale ? 1.234 : 1.0;

  Ctype C_expect;

  constexpr ADObjType At = get_a2d_object_type<Atype>::value;
  constexpr ADObjType Bt = get_a2d_object_type<Btype>::value;
  constexpr ADObjType Ct = get_a2d_object_type<Ctype>::value;

  static_assert(
      std::is_same_v<typename get_object_numeric_type<Atype>::type, double>,
      "only support double type");
  static_assert(
      std::is_same_v<typename get_object_numeric_type<Btype>::type, double>,
      "only support double type");
  static_assert(
      std::is_same_v<typename get_object_numeric_type<Ctype>::type, double>,
      "only support double type");

  // Check if shapes are consistent
  if constexpr (opA == MatOp::TRANSPOSE && opB == MatOp::TRANSPOSE) {
    static_assert(Atype::nrows == Btype::ncols &&
                      Atype::ncols == Ctype::nrows &&
                      Btype::nrows == Ctype::ncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opA == MatOp::TRANSPOSE) {
    static_assert(Atype::nrows == Btype::nrows &&
                      Atype::ncols == Ctype::nrows &&
                      Btype::ncols == Ctype::ncols,
                  "Matrix dimensions must agree.");
  } else if constexpr (opB == MatOp::TRANSPOSE) {
    static_assert(Atype::ncols == Btype::ncols &&
                      Atype::nrows == Ctype::nrows &&
                      Btype::nrows == Ctype::ncols,
                  "Matrix dimensions must agree.");
  } else {
    static_assert(Atype::ncols == Btype::nrows &&
                      Atype::nrows == Ctype::nrows &&
                      Btype::ncols == Ctype::ncols,
                  "Matrix dimensions must agree.");
  }

  for (int i = 0; i < A.ncomp; i++) {
    A[i] = static_cast<T>(rand()) / RAND_MAX;
  }
  for (int i = 0; i < B.ncomp; i++) {
    B[i] = static_cast<T>(rand()) / RAND_MAX;
  }
  for (int i = 0; i < C.ncomp; i++) {
    C[i] = static_cast<T>(rand()) / RAND_MAX;
    C_expect[i] = C[i];
  }

  // Compute ground truth
  int idim = opA == MatOp::NORMAL ? A.nrows : A.ncols;
  int jdim = opA == MatOp::NORMAL ? A.ncols : A.nrows;
  int kdim = opB == MatOp::NORMAL ? B.ncols : B.nrows;
  for (int i = 0; i < idim; i++) {
    for (int k = 0; k < kdim; k++) {
      T value = 0.0;
      for (int j = 0; j < jdim; j++) {
        value += A(opA == MatOp::NORMAL ? i : j, opA == MatOp::NORMAL ? j : i) *
                 B(opB == MatOp::NORMAL ? j : k, opB == MatOp::NORMAL ? k : j);
      }
      if constexpr (additive) {
        C_expect(i, k) += alpha * value;
      } else {
        C_expect(i, k) = alpha * value;
      }
    }
  }

  // Compute test values
  if constexpr (At == ADObjType::MATRIX and Bt == ADObjType::MATRIX) {
    if (scale) {
      MatMatMultScaleCore<T, Atype::nrows, Atype::ncols, Btype::nrows,
                          Btype::ncols, Ctype::nrows, Ctype::ncols, opA, opB,
                          additive>(alpha, get_data(A), get_data(B),
                                    get_data(C));
    } else {
      MatMatMultCore<T, Atype::nrows, Atype::ncols, Btype::nrows, Btype::ncols,
                     Ctype::nrows, Ctype::ncols, opA, opB, additive>(
          get_data(A), get_data(B), get_data(C));
    }
  } else if constexpr (At == ADObjType::SYMMAT and Bt == ADObjType::SYMMAT) {
    if (scale) {
      SMatSMatMultScaleCore<T, Atype::nrows, Btype::nrows, Ctype::nrows,
                            Ctype::ncols, additive>(alpha, get_data(A),
                                                    get_data(B), get_data(C));
    } else {
      SMatSMatMultCore<T, Atype::nrows, Btype::nrows, Ctype::nrows,
                       Ctype::ncols, additive>(get_data(A), get_data(B),
                                               get_data(C));
    }
  } else if constexpr (At == ADObjType::SYMMAT and Bt == ADObjType::MATRIX) {
    if (scale) {
      SMatMatMultScaleCore<T, Atype::nrows, Btype::nrows, Btype::ncols,
                           Ctype::nrows, Ctype::ncols, opB, additive>(
          alpha, get_data(A), get_data(B), get_data(C));
    } else {
      SMatMatMultCore<T, Atype::nrows, Btype::nrows, Btype::ncols, Ctype::nrows,
                      Ctype::ncols, opB, additive>(get_data(A), get_data(B),
                                                   get_data(C));
    }
  } else if constexpr (At == ADObjType::MATRIX and Bt == ADObjType::SYMMAT) {
    if (scale) {
      MatSMatMultScaleCore<T, Atype::nrows, Atype::ncols, Btype::nrows,
                           Ctype::nrows, Ctype::ncols, opA, additive>(
          alpha, get_data(A), get_data(B), get_data(C));
    } else {
      MatSMatMultCore<T, Atype::nrows, Atype::ncols, Btype::nrows, Ctype::nrows,
                      Ctype::ncols, opA, additive>(get_data(A), get_data(B),
                                                   get_data(C));
    }
  } else {
    throw std::runtime_error("Not implemented");
  }

  // Check value match
  for (int i = 0; i < C.nrows; i++) {
    for (int j = 0; j < C.ncols; j++) {
      EXPECT_DOUBLE_EQ(C_expect(i, j), C(i, j))
          << "i: " + std::to_string(i) + ", j: " + std::to_string(j);
    }
  }

  // // Debug
  // std::printf("A:\n");
  // print_mat<Atype::nrows, Atype::ncols>(A);
  //
  // std::printf("B:\n");
  // print_mat<Btype::nrows, Btype::ncols>(B);
  //
  // std::printf("C:\n");
  // print_mat<Ctype::nrows, Ctype::ncols>(C);
  //
  // std::printf("C_expect:\n");
  // print_mat<Ctype::nrows, Ctype::ncols>(C_expect);
}

template <bool additive, MatOp opA, MatOp opB, class case_t>
void run_single_test(case_t c) {
  using T = double;

  for (auto scale : {true, false}) {
    auto A = std::get<0>(c);
    auto B = std::get<1>(c);
    auto C = std::get<2>(c);
    test_gemm_core<additive, opA, opB>(A, B, C, scale);
  }
}

template <bool additive, MatOp opA, MatOp opB, class cases_t>
void run_tests(cases_t cases) {
  using T = double;

  for (auto& c_v : cases) {
    for (auto scale : {true, false}) {
      std::visit(
          [&](const auto& c) {
            auto A = std::get<0>(c);
            auto B = std::get<1>(c);
            auto C = std::get<2>(c);
            test_gemm_core<additive, opA, opB>(A, B, C, scale);
          },
          c_v);
    }
  }
}

TEST(test_a2dgemmcore, square_matrices) {
  using case1 = std::tuple<Mat<T, 2, 2>, Mat<T, 2, 2>, Mat<T, 2, 2>>;
  using case2 = std::tuple<Mat<T, 3, 3>, Mat<T, 3, 3>, Mat<T, 3, 3>>;
  using case3 = std::tuple<Mat<T, 4, 4>, Mat<T, 4, 4>, Mat<T, 4, 4>>;

  using case4 = std::tuple<SymMat<T, 2>, SymMat<T, 2>, Mat<T, 2, 2>>;
  using case5 = std::tuple<SymMat<T, 3>, SymMat<T, 3>, Mat<T, 3, 3>>;
  using case6 = std::tuple<SymMat<T, 4>, SymMat<T, 4>, Mat<T, 4, 4>>;

  using case7 = std::tuple<SymMat<T, 2>, Mat<T, 2, 2>, Mat<T, 2, 2>>;
  using case8 = std::tuple<SymMat<T, 3>, Mat<T, 3, 3>, Mat<T, 3, 3>>;
  using case9 = std::tuple<SymMat<T, 4>, Mat<T, 4, 4>, Mat<T, 4, 4>>;

  using case10 = std::tuple<Mat<T, 2, 2>, SymMat<T, 2>, Mat<T, 2, 2>>;
  using case11 = std::tuple<Mat<T, 3, 3>, SymMat<T, 3>, Mat<T, 3, 3>>;
  using case12 = std::tuple<Mat<T, 4, 4>, SymMat<T, 4>, Mat<T, 4, 4>>;

  std::vector<std::variant<case1, case2, case3, case4, case5, case6, case7,
                           case8, case9, case10, case11, case12>>
      cases = {case1(), case2(), case3(), case4(),  case5(),  case6(),
               case7(), case8(), case9(), case10(), case11(), case12()};

  // Regular and Scale
  run_tests<false, MatOp::NORMAL, MatOp::NORMAL>(cases);
  run_tests<false, MatOp::NORMAL, MatOp::TRANSPOSE>(cases);
  run_tests<false, MatOp::TRANSPOSE, MatOp::NORMAL>(cases);
  run_tests<false, MatOp::TRANSPOSE, MatOp::TRANSPOSE>(cases);

  // Add and ScaleAdd
  run_tests<true, MatOp::NORMAL, MatOp::NORMAL>(cases);
  run_tests<true, MatOp::NORMAL, MatOp::TRANSPOSE>(cases);
  run_tests<true, MatOp::TRANSPOSE, MatOp::NORMAL>(cases);
  run_tests<true, MatOp::TRANSPOSE, MatOp::TRANSPOSE>(cases);
}

TEST(test_a2dgemmcore, rectangular_matrices) {
  using case1 = std::tuple<Mat<T, 2, 3>, Mat<T, 3, 4>, Mat<T, 2, 4>>;
  using case2 = std::tuple<Mat<T, 3, 2>, Mat<T, 3, 4>, Mat<T, 2, 4>>;
  using case3 = std::tuple<Mat<T, 2, 3>, Mat<T, 4, 3>, Mat<T, 2, 4>>;
  using case4 = std::tuple<Mat<T, 3, 2>, Mat<T, 4, 3>, Mat<T, 2, 4>>;

  using case5 = std::tuple<SymMat<T, 4>, Mat<T, 4, 5>, Mat<T, 4, 5>>;
  using case6 = std::tuple<SymMat<T, 2>, Mat<T, 3, 2>, Mat<T, 2, 3>>;

  using case7 = std::tuple<Mat<T, 4, 5>, SymMat<T, 5>, Mat<T, 4, 5>>;
  using case8 = std::tuple<Mat<T, 5, 4>, SymMat<T, 5>, Mat<T, 4, 5>>;

  // Regular and Scale
  run_single_test<false, MatOp::NORMAL, MatOp::NORMAL>(case1());
  run_single_test<false, MatOp::TRANSPOSE, MatOp::NORMAL>(case2());
  run_single_test<false, MatOp::NORMAL, MatOp::TRANSPOSE>(case3());
  run_single_test<false, MatOp::TRANSPOSE, MatOp::TRANSPOSE>(case4());
  run_single_test<false, MatOp::NORMAL, MatOp::NORMAL>(case5());
  run_single_test<false, MatOp::NORMAL, MatOp::TRANSPOSE>(case6());
  run_single_test<false, MatOp::NORMAL, MatOp::NORMAL>(case7());
  run_single_test<false, MatOp::TRANSPOSE, MatOp::NORMAL>(case8());

  // Add and ScaleAdd
  run_single_test<true, MatOp::NORMAL, MatOp::NORMAL>(case1());
  run_single_test<true, MatOp::TRANSPOSE, MatOp::NORMAL>(case2());
  run_single_test<true, MatOp::NORMAL, MatOp::TRANSPOSE>(case3());
  run_single_test<true, MatOp::TRANSPOSE, MatOp::TRANSPOSE>(case4());
  run_single_test<true, MatOp::NORMAL, MatOp::NORMAL>(case5());
  run_single_test<true, MatOp::NORMAL, MatOp::TRANSPOSE>(case6());
  run_single_test<true, MatOp::NORMAL, MatOp::NORMAL>(case7());
  run_single_test<true, MatOp::TRANSPOSE, MatOp::NORMAL>(case8());
}
