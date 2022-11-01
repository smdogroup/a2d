//
// Created by James on 10/25/2022.
//

#include <stdio.h>
#include <iostream>
#include "a2dtypes.h"
#include "shell_development.h"

/*template <int N, typename T, int D>
A2D::Vec<T, D> VecScalarAssembly(A2D::A2DVec<N, A2D::Vec<T, D>>& v, A2D::A2DScalar<N, T> x[D]) {
  return;
};*/

//template <int N, typename T>
//class A2DInstanceTest {
// public:
//  A2DInstanceTest(A2D::A2DVec<N, A2D::Vec<T, 3>>& x, A2D::A2DVec<N, A2D::Vec<T, 3>>& y,
//                  A2D::A2DVec<N, A2D::Vec<T, 3>>& v)
//      : x_cross_y(x_cross_y_v, x_cross_y_bv),
//        x_cross_y_expression(x, y, x_cross_y),
//        x_dot_y_expression(x, y, x_dot_y),
//        v_expression(x_cross_y, x_dot_y, v) {
//  };
//
// private:
//  /* Instantiating objects: */
//  A2D::Vec<T, 3> x_cross_y_v, x_cross_y_bv;
//
//  /* A2D Objects: */
//  A2D::A2DVec<N, A2D::Vec<T, 3>> x_cross_y;
//  A2D::A2DScalar<N, T> x_dot_y;
//
//  /* Expressions: */
//  A2D::A2DVec3CrossA2DVecExpr<N, T> x_cross_y_expression;
//  A2D::A2DVec3DotA2DVecExpr<N, T> x_dot_y_expression;
//  A2D::A2DVec3A2DScaleExpr<N, T> v_expression;
//};
//
//template <typename T, int P>
//void print_vector(const A2D::Vec<T, P>& x) {
//  for (int i = 0; i < P; ++i) {
//    std::cout << x(i) << ", ";
//  }
//  std::cout << std::endl;
//};

int main() {
  /*/// Test to see if Vec3Scale (and other operators) are cast to non-a2d types if the output is constant
  const double x_data[3] = {0.20970639762168464, 0.10308079207850707, 0.9173999579053238};
  const double xb_data[3] = {0.4242524109372384, 0.28692392009681067, 0.9492874123274755};
  A2D::Vec<double, 3> x(x_data), xb(xb_data);
  A2D::ADVec<A2D::Vec<double, 3>> x_ad(x, xb);

//  A2D::Vec<double, 3> y;
//  A2D::Vec3Scale(x_ad, 5.0, y);
  // ^ error because no matching signature*/

  /*/// Test to make sure LinearIsotropicMaterial class works as intended.
  A2D::LinearIsotropicMaterial<double> test_mat1(5.4, 0.5);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      std::cout << test_mat1.D(i, j) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/

  /*/// Test to see if derivatives act as expected
  const double u_data[3] = {0.5, -1.2, 3.7};
  const double v_data[3] = {1.3, 0.82, -1.01};
//  const double norm_ub_data[3] = {1, 1, 1};

  A2D::Vec<double, 3> u(u_data), v(v_data);
//  norm_ub(norm_ub_data)
  A2D::Vec<double, 3> ub, norm_u, norm_ub;

  A2D::A2DVec<3, A2D::Vec<double, 3>> u_a2d(u, ub);
  A2D::A2DVec<3, A2D::Vec<double, 3>> norm_u_a2d(norm_u, norm_ub);
  A2D::A2DScalar<3, double> W(0, 1);

  auto norm_expr = A2D::Vec3Normalize(u_a2d, norm_u_a2d);
  auto dot_expr = A2D::Vec3Dot(norm_u_a2d, v, W);

  dot_expr.reverse();
  norm_expr.reverse();
  for (int j = 0; j < 3; ++j) {
    std::cout << u_a2d.bvalue()(j) << ", ";
  }
  std::cout << std::endl << std::endl;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == j) { u_a2d.pvalue(i)(j) = 1; } else { u_a2d.pvalue(i)(j) = 0; }
    }
  }
  norm_expr.hforward();
  dot_expr.hforward();
  dot_expr.hreverse();
  norm_expr.hreverse();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      std::cout << u_a2d.hvalue(i)(j) << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/

  /*/// Test to see if an array of matrices behaves as expected
  A2D::Mat<double, 2, 2> test[3]{{1, 2, 3, 4}};*/

  /*/// Test to check that Ni_rj_sk returns correct values
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        std::cout << A2D::ShellElementMITC4<24, double>::Ni_rj_sk(i, j, k) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }*/

  /*/// Test of branch-less program
  int r_ind;
  int s_ind;
  for (int alpha = 0; alpha < 2; ++alpha) {
    for (int n_alpha_val_ind = 0; n_alpha_val_ind < 2; ++n_alpha_val_ind) {
      /// Branch-less logic
      r_ind = (4 * n_alpha_val_ind * alpha) + (2 * (1 - alpha));
      s_ind = (2 * alpha) + (4 * n_alpha_val_ind * (1 - alpha));
      /// Branching logic:
      if (alpha == 0) {
        std::cout << (r_ind == 2) << std::endl;
        std::cout << (s_ind == 4 * n_alpha_val_ind) << std::endl;
      } else {
        std::cout << (s_ind == 2) << std::endl;
        std::cout << (r_ind == 4 * n_alpha_val_ind) << std::endl;
      }
    }
  }*/

  /*/// Instantiation test
  const double x_data[3]{-2.3, 1.2, 5};
  const double y_data[3]{8, -32, 0.6};
  A2D::Vec<double, 3> xv{x_data}, xbv;
  A2D::Vec<double, 3> yv{y_data}, ybv;
  A2D::Vec<double, 3> vv, vbv;
  A2D::A2DVec<2, A2D::Vec<double, 3>> x{xv, xbv}, y{yv, ybv}, v{vv, vbv};
  A2DInstanceTest<2, double> cross_dot_expr(x, y, v);
  for (int i = 0; i < 3; ++i) {
    std::cout << v.value()(i) << ", ";
  }
  std::cout << std::endl;*/

  /*/// Test the behavior of += with two comma separated outputs
  double x{0}, y{0};
  x, y += 4.5;
  std::cout << x << std::endl;
  std::cout << y << std::endl;*/

  /*/// Test to see resetting behavior of loops with declared vectors
  for (int i = 0; i < 3; ++i) {
    A2D::Vec<double, 3> x;
    for (int j = 0; j < 3; ++j) {
      std::cout << x(j) << ", ";
    }
    std::cout << std::endl;
    x(0) += 2;
    x(1) *= 2;
    x(2) -= 1;
    for (int j = 0; j < 3; ++j) {
      std::cout << x(j) << ", ";
    }
    std::cout << std::endl;
  }*/

  /// Test shell node constructor
  const double initial_position1_data[3] = {0, 0, 0};
  const double initial_position2_data[3] = {0, 1, 0};
  const double initial_position3_data[3] = {1, 0, 0};
  const double initial_position4_data[3] = {1, 1, 0};
  A2D::Vec<double, 3> initial_position1(initial_position1_data);
  A2D::Vec<double, 3> initial_position2(initial_position2_data);
  A2D::Vec<double, 3> initial_position3(initial_position3_data);
  A2D::Vec<double, 3> initial_position4(initial_position4_data);
  const double initial_shell_director_data[3] = {0, 0, 1};
  A2D::Vec<double, 3> initial_shell_director(initial_shell_director_data);

  A2D::ShellNodeMITC<24, double> n1_1{initial_position1, initial_shell_director};
  A2D::ShellNodeMITC<24, double> n1_2{initial_position2, initial_shell_director};
  A2D::ShellNodeMITC<24, double> n1_3{initial_position3, initial_shell_director};
  A2D::ShellNodeMITC<24, double> n1_4{initial_position4, initial_shell_director};

  n1_1.displacement.value()(0) += 0.01;
  n1_1.displacement.value()(1) += 0.01;
  n1_1.displacement.value()(2) += 0.01;
  n1_3.displacement.value()(0) += 0.01;
  n1_3.displacement.value()(1) += 0.01;
  n1_3.displacement.value()(2) += 0.01;

  A2D::ShellElementMITC4 x1(n1_1,
                            n1_2,
                            n1_3,
                            n1_4);
  /*A2D::Vec<double, 3> position_vec;
  x.position(0, 0, 0, position_vec);
  print_vector(position_vec);*/
  std::cout << x1.strain_energy.value << std::endl << std::endl;

//  const double initial_position2_1_data[3] = {0, 0, 0};
//  const double initial_position2_2_data[3] = {0, 1, 0};
//  const double initial_position2_3_data[3] = {1, 0, 0};
//  const double initial_position2_4_data[3] = {1, 1, 0};
//  A2D::Vec<double, 3> initial_position2_1(initial_position2_1_data);
//  A2D::Vec<double, 3> initial_position2_2(initial_position2_2_data);
//  A2D::Vec<double, 3> initial_position2_3(initial_position2_3_data);
//  A2D::Vec<double, 3> initial_position2_4(initial_position2_4_data);
//  const double initial_shell_director2_data[3] = {0, 0, 1};
//  A2D::Vec<double, 3> initial_shell_director2(initial_shell_director2_data);
//
//  A2D::ShellNodeMITC<24, double> n2_1{initial_position2_1, initial_shell_director2};
//  A2D::ShellNodeMITC<24, double> n2_2{initial_position2_2, initial_shell_director2};
//  A2D::ShellNodeMITC<24, double> n2_3{initial_position2_3, initial_shell_director2};
//  A2D::ShellNodeMITC<24, double> n2_4{initial_position2_4, initial_shell_director2};
//  double move_amount[3] = {2.3, 1.2, 3.9};
//  for (int i = 0; i < 3; ++i) {
//    n2_1.position(i) += move_amount[i];
//    n2_2.position(i) += move_amount[i];
//    n2_3.position(i) += move_amount[i];
//    n2_4.position(i) += move_amount[i];
//  }
//  A2D::ShellElementMITC4 x2(n2_1, n2_2, n2_3, n2_4);
//  std::cout << x2.strain_energy.value << std::endl << std::endl;
////  print_vector(n2_1.displacement.value());
////  print_vector(n2_2.displacement.value());
////  print_vector(n2_3.displacement.value());
////  print_vector(n2_4.displacement.value());
////  print_vector(n2_1.position);
////  print_vector(n2_2.position);
////  print_vector(n2_3.position);
////  print_vector(n2_4.position);
////  print_vector(n1_1.position);
////  print_vector(n1_2.position);
////  print_vector(n1_3.position);
////  print_vector(n1_4.position);
//
//  const double initial_position3_1_data[3] = {0, 0, 0};
//  const double initial_position3_2_data[3] = {0, 1, 0};
//  const double initial_position3_3_data[3] = {1, 0, 0};
//  const double initial_position3_4_data[3] = {1, 1, 0};
//  A2D::Vec<double, 3> initial_position3_1(initial_position3_1_data);
//  A2D::Vec<double, 3> initial_position3_2(initial_position3_2_data);
//  A2D::Vec<double, 3> initial_position3_3(initial_position3_3_data);
//  A2D::Vec<double, 3> initial_position3_4(initial_position3_4_data);
//  const double initial_shell_director3_data[3] = {0, 0, 1};
//  A2D::Vec<double, 3> initial_shell_director3(initial_shell_director3_data);
//
//  A2D::ShellNodeMITC<24, double> n3_1{initial_position3_1, initial_shell_director3};
//  A2D::ShellNodeMITC<24, double> n3_2{initial_position3_2, initial_shell_director3};
//  A2D::ShellNodeMITC<24, double> n3_3{initial_position3_3, initial_shell_director3};
//  A2D::ShellNodeMITC<24, double> n3_4{initial_position3_4, initial_shell_director3};
////  A2D::ShellNodeMITC<24, double> n3_1{n1_1};
////  A2D::ShellNodeMITC<24, double> n3_2{n1_2};
////  A2D::ShellNodeMITC<24, double> n3_3{n1_3};
////  A2D::ShellNodeMITC<24, double> n3_4{n1_4};
//  double move_amount3[3] = {-2.3, 51.2, -0.9};
//  for (int i = 0; i < 3; ++i) {
//    n3_1.position(i) += move_amount3[i];
//    n3_2.position(i) += move_amount3[i];
//    n3_3.position(i) += move_amount3[i];
//    n3_4.position(i) += move_amount3[i];
//  }
//  A2D::ShellElementMITC4 x3(n3_1,
//                            n3_2,
//                            n3_3,
//                            n3_4);
//  std::cout << x3.strain_energy.value << std::endl << std::endl;
//
//  const double initial_position4_1_data[3] = {0, 0, 0};
//  const double initial_position4_2_data[3] = {0, 1, 0};
//  const double initial_position4_3_data[3] = {1, 0, 0};
//  const double initial_position4_4_data[3] = {1, 1, 0};
//  A2D::Vec<double, 3> initial_position4_1(initial_position4_1_data);
//  A2D::Vec<double, 3> initial_position4_2(initial_position4_2_data);
//  A2D::Vec<double, 3> initial_position4_3(initial_position4_3_data);
//  A2D::Vec<double, 3> initial_position4_4(initial_position4_4_data);
//  const double initial_shell_director4_data[3] = {0, 0, 1};
//  A2D::Vec<double, 3> initial_shell_director4(initial_shell_director4_data);
//
//  A2D::ShellNodeMITC<24, double> n4_1{initial_position4_1, initial_shell_director4};
//  A2D::ShellNodeMITC<24, double> n4_2{initial_position4_2, initial_shell_director4};
//  A2D::ShellNodeMITC<24, double> n4_3{initial_position4_3, initial_shell_director4};
//  A2D::ShellNodeMITC<24, double> n4_4{initial_position4_4, initial_shell_director4};
////  A2D::ShellNodeMITC<24, double> n4_1{n1_1};
////  A2D::ShellNodeMITC<24, double> n4_2{n1_2};
////  A2D::ShellNodeMITC<24, double> n4_3{n1_3};
////  A2D::ShellNodeMITC<24, double> n4_4{n1_4};
//  double move_amount4[3] = {-2.3, 51.2, -0.9};
//  for (int i = 0; i < 3; ++i) {
//    n4_1.position(i) += move_amount4[i];
//    n4_2.position(i) += move_amount4[i];
//    n4_3.position(i) += move_amount4[i];
//    n4_4.position(i) += move_amount4[i];
//  }
//  A2D::ShellElementMITC4 x4(n4_1,
//                            n4_2,
//                            n4_3,
//                            n4_4);
//  std::cout << x4.strain_energy.value << std::endl << std::endl;

//  x.strain_energy.bvalue = 1;
//  x.reverse();
//  std::cout << std::endl;
//  print_vector(x.node1.displacement.bvalue());
//  print_vector(x.node1.rotation.bvalue());
//  print_vector(x.node2.displacement.bvalue());
//  print_vector(x.node2.rotation.bvalue());
//  print_vector(x.node3.displacement.bvalue());
//  print_vector(x.node3.rotation.bvalue());
//  print_vector(x.node4.displacement.bvalue());
//  print_vector(x.node4.rotation.bvalue());
//  std::cout << std::endl;
//
//  /* Zero out p values: */
//  for (int i = 0; i < 24; ++i) {
//    x.node1.displacement.pvalue(i).zero();
//    x.node1.rotation.pvalue(i).zero();
//    x.node2.displacement.pvalue(i).zero();
//    x.node2.rotation.pvalue(i).zero();
//    x.node3.displacement.pvalue(i).zero();
//    x.node3.rotation.pvalue(i).zero();
//    x.node4.displacement.pvalue(i).zero();
//    x.node4.rotation.pvalue(i).zero();
//  }
//  /* Assign unit p values to degrees of freedom: */
//  x.node1.displacement.pvalue(0)(0) = 1;
//  x.node1.displacement.pvalue(1)(1) = 1;
//  x.node1.displacement.pvalue(2)(2) = 1;
//  x.node1.rotation.pvalue(3)(0) = 1;
//  x.node1.rotation.pvalue(4)(1) = 1;
//  x.node1.rotation.pvalue(5)(2) = 1;
//  x.node2.displacement.pvalue(6)(0) = 1;
//  x.node2.displacement.pvalue(7)(1) = 1;
//  x.node2.displacement.pvalue(8)(2) = 1;
//  x.node2.rotation.pvalue(9)(0) = 1;
//  x.node2.rotation.pvalue(10)(1) = 1;
//  x.node2.rotation.pvalue(11)(2) = 1;
//  x.node3.displacement.pvalue(12)(0) = 1;
//  x.node3.displacement.pvalue(13)(1) = 1;
//  x.node3.displacement.pvalue(14)(2) = 1;
//  x.node3.rotation.pvalue(15)(0) = 1;
//  x.node3.rotation.pvalue(16)(1) = 1;
//  x.node3.rotation.pvalue(17)(2) = 1;
//  x.node4.displacement.pvalue(18)(0) = 1;
//  x.node4.displacement.pvalue(19)(1) = 1;
//  x.node4.displacement.pvalue(20)(2) = 1;
//  x.node4.rotation.pvalue(21)(0) = 1;
//  x.node4.rotation.pvalue(22)(1) = 1;
//  x.node4.rotation.pvalue(23)(2) = 1;
//  /* forward and reverse passes */
//  x.hforward();
//  for (int i = 0; i < 24; ++i) {
//    std::cout << x.strain_energy.pvalue[i] << std::endl;
//  }
//  std::cout << std::endl;
//  x.hreverse();
//
//  for (int i = 0; i < 24; ++i) {
//    printf(""
//           "% 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, % 5.2f, "
//           //           "% 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, % 7.2f, "
//           "\n",
//           x.node1.displacement.hvalue(i)(0),
//           x.node1.displacement.hvalue(i)(1),
//           x.node1.displacement.hvalue(i)(2),
//           x.node1.rotation.hvalue(i)(0),
//           x.node1.rotation.hvalue(i)(1),
//           x.node1.rotation.hvalue(i)(2),
//           x.node2.displacement.hvalue(i)(0),
//           x.node2.displacement.hvalue(i)(1),
//           x.node2.displacement.hvalue(i)(2),
//           x.node2.rotation.hvalue(i)(0),
//           x.node2.rotation.hvalue(i)(1),
//           x.node2.rotation.hvalue(i)(2)
//        /*x.node3.displacement.hvalue(i)(0),
//        x.node3.displacement.hvalue(i)(1),
//        x.node3.displacement.hvalue(i)(2),
//        x.node3.rotation.hvalue(i)(0),
//        x.node3.rotation.hvalue(i)(1),
//        x.node3.rotation.hvalue(i)(2),
//        x.node4.displacement.hvalue(i)(0),
//        x.node4.displacement.hvalue(i)(1),
//        x.node4.displacement.hvalue(i)(2),
//        x.node4.rotation.hvalue(i)(0),
//        x.node4.rotation.hvalue(i)(1),
//        x.node4.rotation.hvalue(i)(2)*/
//    );
//  }


  /*print_vector(n1_1.rotation.value());
  print_vector(n1_2.rotation.value());
  print_vector(n1_3.rotation.value());
  print_vector(n1_4.rotation.value());

  std::cout << std::endl;
  n1_1.rotation.value()(0) += 1;
  print_vector(n1_1.rotation.value());
  print_vector(n1_2.rotation.value());
  print_vector(n1_3.rotation.value());
  print_vector(n1_4.rotation.value());*/

};
