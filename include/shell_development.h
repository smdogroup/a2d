//
// Created by James on 10/23/2022.
//

#ifndef A2D_SHELL_DEVELOPMENT_H_
#define A2D_SHELL_DEVELOPMENT_H_

#include "a2dtypes.h"
#include "a2dvecops3d.h"
#include "a2dmatops3d.h"
#include "shell_dev_ops.h"

namespace A2D {

template <typename T>
class LinearIsotropicMaterial {
 private:
  const T temp1;
  const T temp2;
  const T temp3;
  const T D_init[25];

 public:
  LinearIsotropicMaterial(T E, T nu)
      : temp1(E / (1 - (nu * nu))), temp2(temp1 * nu), temp3(E / (2 * (1 + nu))),
        D_init{temp1, temp2, 0, 0, 0,
               temp2, temp1, 0, 0, 0,
               0, 0, temp3, 0, 0,
               0, 0, 0, temp3 * 5 / 6, 0,
               0, 0, 0, 0, temp3 * 5 / 6}, D(D_init) {};

  const Mat<T, 5, 5> D;
};

/**********************************************************************************************************************
 * @brief A class to contain all necessary properties of a MITC Shell element node.
 *
 * A MITC shell node has 13 values, six of which are degrees of freedom, associated with the node.  The six degrees of
 * freedom come, three each, from the displacement and rotation vectors.  The other seven values are design parameters,
 * and they are the position of the node (3), thickness of the shell at the node (1), and the shell director vector (3)
 * describing how the cross section is angled in relation to the shell's mid-surface.
 **********************************************************************************************************************/
template <int N, typename T>
class ShellNodeMITC {
 public:
  ShellNodeMITC(Vec<T, 3> disp, Vec<T, 3> dispb, Vec<T, 3> rot, Vec<T, 3> rotb)
      : displacement(disp, dispb), rotation(rot, rotb) {
    // TODO: delete me
    position.zero();
    shell_director.zero();
  }

//  Vec<T, 3>& position;  /**< <u> X </u> <sup> I </sup>: The position vector at the node.*/
//  A2DScalar<N, T>& thickness;  /**< h <sup> I </sup>: The thickness of the shell at the node. */
//  Vec<T, 3>& shell_director;  /**< <u> v </u> <sup> I </sup>: The shell director vector at the node. */
//  A2DVec<N, Vec<T, 3>>& displacement;  /**< <u> u </u> <sup> I </sup>: The displacement vector at the node.*/
//  A2DVec<N, Vec<T, 3>>& rotation;  /**< <u> &theta </u> <sup> I </sup>: The rotation vector at the node. */

//  A2DVec<N, Vec<T, 3>>& position;  /**< <u> X </u> <sup> I </sup>: The position vector at the node.*/
  Vec<T, 3> position;  /**< <u> X </u> <sup> I </sup>: The position vector at the node.*/
  A2DScalar<N, T> thickness;  /**< h <sup> I </sup>: The thickness of the shell at the node. */
//  A2DVec<N, Vec<T, 3>>& shell_director;  /**< <u> v </u> <sup> I </sup>: The shell director vector at the node. */
  Vec<T, 3> shell_director;  /**< <u> v </u> <sup> I </sup>: The shell director vector at the node. */
  A2DVec<N, Vec<T, 3>> displacement;  /**< <u> u </u> <sup> I </sup>: The displacement vector at the node.*/
  A2DVec<N, Vec<T, 3>> rotation;  /**< <u> &theta </u> <sup> I </sup>: The rotation vector at the node. */
};

template <int N, typename T>
class g_alpha_expr;

template <int N, typename T>
class u_alpha_expr;

template <int N, typename T>
class g_t_expr;

template <int N, typename T>
class u_t_expr;

template <int N, typename T>
class tying_shear_expr;

template <int N, typename T>
class contravariant_basis_expr;

template <int N, typename T>
class cartesian_local_basis_expr;

template <int N, typename T>
class in_plane_strain_components_expr;

template <int N, typename T>
class local_strain_expr;

template <int N, typename T>
class local_strains_expr;

template <int N, typename T>
class strain_energy_expr;

/**********************************************************************************************************************
 * @brief This is the basic MITC4 element as given by the original paper <it> A continuum mechanics based four-node
 * shell element for general non-linear analysis </it> Dvorkin, E.N. and Bathe, K.J. (1983).
 *
 * The paper from which this element is based can be found at
 * <a href="https://doi.org/10.1108/eb023562"> https://doi.org/10.1108/eb023562 </a>.  Furthermore, a significant
 * amount of useful information regarding this element can be found
 * <a href="https://www.sesamx.io/blog/shell_finite_element/"> here </a>.
 *
 * <h1>Assumptions:</h1>
 *
 * <h2>Stress assumption:</h2>
 *
 * <p>
 * <center> &sigma <sub>3,3</sub> = 0 </center>
 *
 * The 3 direction is collinear to the director vector at each point of the shell mid-surface.
 * </p>
 *
 * <h2>Kinematic assumption:</h2>
 *
 * <p>
 * The position of every point in the space can be described by the parametric coordinates r, s, and t having values
 * from -1 to 1.
 *
 * <center>
 * <u>X</u>(r,s,t) = <u>X</u><sub>0</sub>(r,s) + <sup>t</sup>&frasl;<sub>2</sub> h(r,s) <u>v</u>(r,s)
 * </center>
 *
 * Which becomes the following after the shape functions are applied, which is summed over I=1,2,3,4 for the
 * contribution of each of the four nodes.
 *
 * <center>
 * <u>X</u>(r,s,t) = N<sup>I</sup>(r,s) <u>X</u><sup>I</sup> + <sup>t</sup>&frasl;<sub>2</sub> N<sup>I</sup>(r,s)
 * h<sup>I</sup> <u>v</u><sup>I</sup>
 * </center>
 * </p>
 *
 * <p>
 * Similarly, the displacements are interpolated from the values at the nodes using the same parametric coordinates: r,
 * s, and t, and is summed over I=1,2,3,4 for the contribution of the four nodes.
 *
 * <center>
 * <u>u</u>(r,s,t) = N<sup>I</sup>(r,s) <u>u</u><sup>I</sup> + <sup>t</sup>&frasl;<sub>2</sub> N<sup>I</sup>(r,s)
 * h<sup>I</sup> (<u>&theta</u><sup>I</sup> &times <u>v</u><sup>I</sup>)
 * </center>
 * </p>
 **********************************************************************************************************************/
template <int N, typename T>  // N = 6 * 4;
class ShellElementMITC4 {
 public:
  ShellElementMITC4(ShellNodeMITC<N, T>& node1,
                    ShellNodeMITC<N, T>& node2,
                    ShellNodeMITC<N, T>& node3,
                    ShellNodeMITC<N, T>& node4)
      : node1(node1), node2(node2), node3(node3), node4(node4),
      /* Private Initializations: */
        gr_rAs0t0(gr_rAs0t0_v, gr_rAs0t0_bv),
        gr_rAs0t1(gr_rAs0t1_v, gr_rAs0t1_bv),
        gr_rAs1t0(gr_rAs1t0_v, gr_rAs1t0_bv),
        gr_rAs1t1(gr_rAs1t1_v, gr_rAs1t1_bv),
        gs_r0sAt0(gs_r0sAt0_v, gs_r0sAt0_bv),
        gs_r0sAt1(gs_r0sAt1_v, gs_r0sAt1_bv),
        gs_r1sAt0(gs_r1sAt0_v, gs_r1sAt0_bv),
        gs_r1sAt1(gs_r1sAt1_v, gs_r1sAt1_bv),
        gt_r0s0tA(gt_r0s0tA_v, gt_r0s0tA_bv),
        gt_r0s1tA(gt_r0s1tA_v, gt_r0s1tA_bv),
        gt_r1s0tA(gt_r1s0tA_v, gt_r1s0tA_bv),
        gt_r1s1tA(gt_r1s1tA_v, gt_r1s1tA_bv),

        ur_rAs0t0(ur_rAs0t0_v, ur_rAs0t0_bv),
        ur_rAs0t1(ur_rAs0t1_v, ur_rAs0t1_bv),
        ur_rAs1t0(ur_rAs1t0_v, ur_rAs1t0_bv),
        ur_rAs1t1(ur_rAs1t1_v, ur_rAs1t1_bv),
        us_r0sAt0(us_r0sAt0_v, us_r0sAt0_bv),
        us_r0sAt1(us_r0sAt1_v, us_r0sAt1_bv),
        us_r1sAt0(us_r1sAt0_v, us_r1sAt0_bv),
        us_r1sAt1(us_r1sAt1_v, us_r1sAt1_bv),
        ut_r0s0tA(ut_r0s0tA_v, ut_r0s0tA_bv),
        ut_r0s1tA(ut_r0s1tA_v, ut_r0s1tA_bv),
        ut_r1s0tA(ut_r1s0tA_v, ut_r1s0tA_bv),
        ut_r1s1tA(ut_r1s1tA_v, ut_r1s1tA_bv),
      /* Energy Operations: */
      /* Tying points: */
        e_rt_A_expression(0, 1, node1, node2, node3, node4, e_rt_A),
        e_rt_B_expression(0, 0, node1, node2, node3, node4, e_rt_B),
        e_st_C_expression(1, 1, node1, node2, node3, node4, e_st_C),
        e_st_D_expression(1, 0, node1, node2, node3, node4, e_st_D),

      /* e_rt_r... and e_st_r... values */
        e_rt_rAs0tA_expression((1 + quad_0) * 0.5, e_rt_A, (1 - quad_0) * 0.5, e_rt_B, e_rt_rAs0tA),
        e_st_r0sAtA_expression((1 + quad_0) * 0.5, e_st_C, (1 - quad_0) * 0.5, e_st_D, e_st_r0sAtA),
        e_rt_rAs1tA_expression((1 + quad_1) * 0.5, e_rt_A, (1 - quad_1) * 0.5, e_rt_B, e_rt_rAs1tA),
        e_st_r1sAtA_expression((1 + quad_1) * 0.5, e_st_C, (1 - quad_1) * 0.5, e_st_D, e_st_r1sAtA),

      /* gr calculations */
        gr_rAs0t0_expression(0, 1, quad_0, node1, node2, node3, node4, gr_rAs0t0),
        gr_rAs0t1_expression(0, 1, quad_1, node1, node2, node3, node4, gr_rAs0t1),
        gr_rAs1t0_expression(0, 2, quad_0, node1, node2, node3, node4, gr_rAs1t0),
        gr_rAs1t1_expression(0, 2, quad_1, node1, node2, node3, node4, gr_rAs1t1),
      /* gs calculations */
        gs_r0sAt0_expression(1, 1, quad_0, node1, node2, node3, node4, gs_r0sAt0),
        gs_r0sAt1_expression(1, 1, quad_1, node1, node2, node3, node4, gs_r0sAt1),
        gs_r1sAt0_expression(1, 2, quad_0, node1, node2, node3, node4, gs_r1sAt0),
        gs_r1sAt1_expression(1, 2, quad_1, node1, node2, node3, node4, gs_r1sAt1),

      /* ur calculations */
        ur_rAs0t0_expression(0, 1, quad_0, node1, node2, node3, node4, ur_rAs0t0),
        ur_rAs0t1_expression(0, 1, quad_1, node1, node2, node3, node4, ur_rAs0t1),
        ur_rAs1t0_expression(0, 2, quad_0, node1, node2, node3, node4, ur_rAs1t0),
        ur_rAs1t1_expression(0, 2, quad_1, node1, node2, node3, node4, ur_rAs1t1),
      /* us calculations */
        us_r0sAt0_expression(1, 1, quad_0, node1, node2, node3, node4, us_r0sAt0),
        us_r0sAt1_expression(1, 1, quad_1, node1, node2, node3, node4, us_r0sAt1),
        us_r1sAt0_expression(1, 2, quad_0, node1, node2, node3, node4, us_r1sAt0),
        us_r1sAt1_expression(1, 2, quad_1, node1, node2, node3, node4, us_r1sAt1),

      /* Energy at the quadrature points: */
        strain_energy_r0s0t0_expression(gr_rAs0t0, gs_r0sAt0, gt_r0s0tA, ur_rAs0t0, us_r0sAt0, ut_r0s0tA,
                                        e_rt_rAs0tA, e_st_r0sAtA, material, energy_r0s0t0),
        strain_energy_r0s0t1_expression(gr_rAs0t1, gs_r0sAt1, gt_r0s0tA, ur_rAs0t1, us_r0sAt1, ut_r0s0tA,
                                        e_rt_rAs0tA, e_st_r0sAtA, material, energy_r0s0t1),
        strain_energy_r0s1t0_expression(gr_rAs1t0, gs_r0sAt0, gt_r0s1tA, ur_rAs1t0, us_r0sAt0, ut_r0s1tA,
                                        e_rt_rAs1tA, e_st_r0sAtA, material, energy_r0s1t0),
        strain_energy_r0s1t1_expression(gr_rAs1t1, gs_r0sAt1, gt_r0s1tA, ur_rAs1t1, us_r0sAt1, ut_r0s1tA,
                                        e_rt_rAs1tA, e_st_r0sAtA, material, energy_r0s1t1),
        strain_energy_r1s0t0_expression(gr_rAs0t0, gs_r1sAt0, gt_r1s0tA, ur_rAs0t0, us_r1sAt0, ut_r1s0tA,
                                        e_rt_rAs0tA, e_st_r1sAtA, material, energy_r1s0t0),
        strain_energy_r1s0t1_expression(gr_rAs0t1, gs_r1sAt1, gt_r1s0tA, ur_rAs0t1, us_r1sAt1, ut_r1s0tA,
                                        e_rt_rAs0tA, e_st_r1sAtA, material, energy_r1s0t1),
        strain_energy_r1s1t0_expression(gr_rAs1t0, gs_r1sAt0, gt_r1s1tA, ur_rAs1t0, us_r1sAt0, ut_r1s1tA,
                                        e_rt_rAs1tA, e_st_r1sAtA, material, energy_r1s1t0),
        strain_energy_r1s1t1_expression(gr_rAs1t1, gs_r1sAt1, gt_r1s1tA, ur_rAs1t1, us_r1sAt1, ut_r1s1tA,
                                        e_rt_rAs1tA, e_st_r1sAtA, material, energy_r1s1t1),

      /* Summation of the strain energy: */
        sum_12_expression(1, energy_r0s0t0, energy_r0s0t1, sum_12),
        sum_123_expression(1, sum_12, energy_r0s1t0, sum_123),
        sum_1234_expression(1, sum_123, energy_r0s1t1, sum_1234),
        sum_12345_expression(1, sum_1234, energy_r1s0t0, sum_12345),
        sum_123456_expression(1, sum_12345, energy_r1s0t1, sum_123456),
        sum_1234567_expression(1, sum_123456, energy_r1s1t0, sum_1234567),
        sum_12345678_expression(1, sum_1234567, energy_r1s1t1, strain_energy) {};

  /**
   * @brief Position of a point on the undeformed shell element defined by parameters r, s, and t.
   *
   * @param r &isin [-1, 1]: one of two mid-surface curvilinear coordinates.  A value of 1 corresponds to the edge
   * between the second and third nodes, a value of -1 corresponds to the edge between the first and fourth nodes.
   * @param s &isin [-1, 1]: one of two mid-surface curvilinear coordinates.  A value of 1 corresponds to the edge
   * between the third and fourth nodes, a value of -1 corresponds to the edge between the first and second nodes.
   * @param t &isin [-1, 1]: local coordinate describing distance along the director vector from the mid-surface. A
   * value of 1 corresponds to the exterior surface, a value of -1 corresponds to the interior surface.
   * @param X: a three vector output indicating the position of the point in the undeformed shell element as defined by
   * the given values of the parameters.
   * */
  template <typename S>
  // might want to replace S with T or double
  void position(const S r, const S s, const S t, Vec<T, 3>& X) {
    S N1, N2, N3, N4;
    S tNh1, tNh2, tNh3, tNh4;
    S s1{1 - s}, s2{1 + s};

    N1 = (1 / 4);
    N2 = N1 * (1 + r);
    N1 *= (1 - r);
    N3 = N2 * s2;
    N4 = N1 * s2;
    N1 *= s1;
    N2 *= s1;

    tNh1 = t / 2;
    tNh2 = tNh1 * N2 * node2.thickness.value;
    tNh3 = tNh1 * N3 * node3.thickness.value;
    tNh4 = tNh1 * N4 * node2.thickness.value;
    tNh1 *= N1 * node1.thickness.value;

    /** N <sup> I </sup> (r, s) <u>X</u> <sup> I </sup> + ...*/
    Vec3Scale(node1.position, N1, X);
    Vec3Axpy(N2, node2.position, X, X);
    Vec3Axpy(N3, node3.position, X, X);
    Vec3Axpy(N4, node4.position, X, X);
    /** ... + (t/2) N <sup> I </sup> (r,s) h <sup> I </sup> <u>v</u> <sup> I </sup>*/
    Vec3Axpy(tNh1, node1.shell_director, X, X);
    Vec3Axpy(tNh2, node2.shell_director, X, X);
    Vec3Axpy(tNh3, node3.shell_director, X, X);
    Vec3Axpy(tNh4, node4.shell_director, X, X);
  }

  void generate_stiffness_matrix() {
    /* Generate the internal energy of the shell element, that will then be used to calculate the global stiffness
     * matrix
     * */

    // code for tying points
    /*tying_shear_expr e_rt_A_expression(0, 1, this, e_rt_A);
    tying_shear_expr e_rt_B_expression(0, 0, this, e_rt_B);
    tying_shear_expr e_st_C_expression(1, 1, this, e_st_C);
    tying_shear_expr e_st_D_expression(1, 0, this, e_st_D);*/

    /* e_rt_r... and e_st_r... values */
    /*ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T>
        e_rt_rAs0tA_expression = ScalarAxpby((1 + quad_0) * 0.5, e_rt_A, (1 - quad_0) * 0.5, e_rt_B, e_rt_rAs0tA);
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T>
        e_st_r0sAtA_expression = ScalarAxpby((1 + quad_0) * 0.5, e_st_C, (1 - quad_0) * 0.5, e_st_D, e_st_r0sAtA);
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T>
        e_rt_rAs1tA_expression = ScalarAxpby((1 + quad_1) * 0.5, e_rt_A, (1 - quad_1) * 0.5, e_rt_B, e_rt_rAs1tA);
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T>
        e_st_r1sAtA_expression = ScalarAxpby((1 + quad_1) * 0.5, e_st_C, (1 - quad_1) * 0.5, e_st_D, e_st_r1sAtA);*/

    /* gr and gs calculations */
    /*g_alpha_expr gr_rAs0t0_expression(0, 1, quad_0, this, gr_rAs0t0);
    g_alpha_expr gr_rAs0t1_expression(0, 1, quad_1, this, gr_rAs0t1);
    g_alpha_expr gr_rAs1t0_expression(0, 2, quad_0, this, gr_rAs1t0);
    g_alpha_expr gr_rAs1t1_expression(0, 2, quad_1, this, gr_rAs1t1);

    g_alpha_expr gs_r0sAt0_expression(1, 1, quad_0, this, gs_r0sAt0);
    g_alpha_expr gs_r0sAt1_expression(1, 1, quad_1, this, gs_r0sAt1);
    g_alpha_expr gs_r1sAt0_expression(1, 2, quad_0, this, gs_r1sAt0);
    g_alpha_expr gs_r1sAt1_expression(1, 2, quad_1, this, gs_r1sAt1);*/

    /* ur and us calculations */
    /*u_alpha_expr ur_rAs0t0_expression(0, 1, quad_0, this, ur_rAs0t0);
    u_alpha_expr ur_rAs0t1_expression(0, 1, quad_1, this, ur_rAs0t1);
    u_alpha_expr ur_rAs1t0_expression(0, 2, quad_0, this, ur_rAs1t0);
    u_alpha_expr ur_rAs1t1_expression(0, 2, quad_1, this, ur_rAs1t1);

    u_alpha_expr us_r0sAt0_expression(1, 1, quad_0, this, us_r0sAt0);
    u_alpha_expr us_r0sAt1_expression(1, 1, quad_1, this, us_r0sAt1);
    u_alpha_expr us_r1sAt0_expression(1, 2, quad_0, this, us_r1sAt0);
    u_alpha_expr us_r1sAt1_expression(1, 2, quad_1, this, us_r1sAt1);*/

    /* Energy at the quadrature points: */
    /*strain_energy_expr strain_energy_r0s0t0_expression(gr_rAs0t0, gs_r0sAt0, gt_r0s0tA, ur_rAs0t0, us_r0sAt0, ut_r0s0tA, e_rt_rAs0tA, e_st_r0sAtA, this, energy_r0s0t0);
    strain_energy_expr strain_energy_r0s0t1_expression(gr_rAs0t1, gs_r0sAt1, gt_r0s0tA, ur_rAs0t1, us_r0sAt1, ut_r0s0tA, e_rt_rAs0tA, e_st_r0sAtA, this, energy_r0s0t1);
    strain_energy_expr strain_energy_r0s1t0_expression(gr_rAs1t0, gs_r0sAt0, gt_r0s1tA, ur_rAs1t0, us_r0sAt0, ut_r0s1tA, e_rt_rAs1tA, e_st_r0sAtA, this, energy_r0s1t0);
    strain_energy_expr strain_energy_r0s1t1_expression(gr_rAs1t1, gs_r0sAt1, gt_r0s1tA, ur_rAs1t1, us_r0sAt1, ut_r0s1tA, e_rt_rAs1tA, e_st_r0sAtA, this, energy_r0s1t1);
    strain_energy_expr strain_energy_r1s0t0_expression(gr_rAs0t0, gs_r1sAt0, gt_r1s0tA, ur_rAs0t0, us_r1sAt0, ut_r1s0tA, e_rt_rAs0tA, e_st_r1sAtA, this, energy_r1s0t0);
    strain_energy_expr strain_energy_r1s0t1_expression(gr_rAs0t1, gs_r1sAt1, gt_r1s0tA, ur_rAs0t1, us_r1sAt1, ut_r1s0tA, e_rt_rAs0tA, e_st_r1sAtA, this, energy_r1s0t1);
    strain_energy_expr strain_energy_r1s1t0_expression(gr_rAs1t0, gs_r1sAt0, gt_r1s1tA, ur_rAs1t0, us_r1sAt0, ut_r1s1tA, e_rt_rAs1tA, e_st_r1sAtA, this, energy_r1s1t0);
    strain_energy_expr strain_energy_r1s1t1_expression(gr_rAs1t1, gs_r1sAt1, gt_r1s1tA, ur_rAs1t1, us_r1sAt1, ut_r1s1tA, e_rt_rAs1tA, e_st_r1sAtA, this, energy_r1s1t1);*/

    /* Summation of the strain energy: */
    /*ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12_expression(1, energy_r0s0t0, energy_r0s0t1, sum_12);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_123_expression(1, sum_12, energy_r0s1t0, sum_123);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_1234_expression(1, sum_123, energy_r0s1t1, sum_1234);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12345_expression(1, sum_1234, energy_r1s0t0, sum_12345);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_123456_expression(1, sum_12345, energy_r1s0t1, sum_123456);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_1234567_expression(1, sum_123456, energy_r1s1t0, sum_1234567);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12345678_expression(1, sum_1234567, energy_r1s1t1, strain_energy);*/
  };

  ShellNodeMITC<N, T>& node1;  /**< The top left node. */
  ShellNodeMITC<N, T>& node2;  /**< The top right node.*/
  ShellNodeMITC<N, T>& node3;  /**< The bottom left node. */
  ShellNodeMITC<N, T>& node4;  /**< The bottom right node.*/
  // TODO: make sure this is correct ^

  const LinearIsotropicMaterial<T> material{5.2, 0.5};  /**< I'm using a linear isotropic material assumption here */

 private:
  /** <h1>Instantiating objects:</h1> */
  Vec<T, 3>
      gr_rAs0t0_v, gr_rAs0t0_bv,
      gr_rAs0t1_v, gr_rAs0t1_bv,
      gr_rAs1t0_v, gr_rAs1t0_bv,
      gr_rAs1t1_v, gr_rAs1t1_bv,
      gs_r0sAt0_v, gs_r0sAt0_bv,
      gs_r0sAt1_v, gs_r0sAt1_bv,
      gs_r1sAt0_v, gs_r1sAt0_bv,
      gs_r1sAt1_v, gs_r1sAt1_bv,
      gt_r0s0tA_v, gt_r0s0tA_bv,
      gt_r0s1tA_v, gt_r0s1tA_bv,
      gt_r1s0tA_v, gt_r1s0tA_bv,
      gt_r1s1tA_v, gt_r1s1tA_bv;
  Vec<T, 3>
      ur_rAs0t0_v, ur_rAs0t0_bv,
      ur_rAs0t1_v, ur_rAs0t1_bv,
      ur_rAs1t0_v, ur_rAs1t0_bv,
      ur_rAs1t1_v, ur_rAs1t1_bv,
      us_r0sAt0_v, us_r0sAt0_bv,
      us_r0sAt1_v, us_r0sAt1_bv,
      us_r1sAt0_v, us_r1sAt0_bv,
      us_r1sAt1_v, us_r1sAt1_bv,
      ut_r0s0tA_v, ut_r0s0tA_bv,
      ut_r0s1tA_v, ut_r0s1tA_bv,
      ut_r1s0tA_v, ut_r1s0tA_bv,
      ut_r1s1tA_v, ut_r1s1tA_bv;

 private:
  /** <h1>A2D Objects:</h1> */
  /** gr vector evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      gr_rAs0t0,
      gr_rAs0t1,
      gr_rAs1t0,
      gr_rAs1t1;
  /** gs vector evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      gs_r0sAt0,
      gs_r0sAt1,
      gs_r1sAt0,
      gs_r1sAt1;
  /** gt vector evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      gt_r0s0tA,
      gt_r0s1tA,
      gt_r1s0tA,
      gt_r1s1tA;

  /** derivatives of u with respect to r evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      ur_rAs0t0,
      ur_rAs0t1,
      ur_rAs1t0,
      ur_rAs1t1;
  /** derivatives of u with respect to s evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      us_r0sAt0,
      us_r0sAt1,
      us_r1sAt0,
      us_r1sAt1;
  /** derivatives of u with respect to t evaluated at the various quadrature points */
  A2DVec<N, Vec<T, 3>>
      ut_r0s0tA,
      ut_r0s1tA,
      ut_r1s0tA,
      ut_r1s1tA;

  /** Coefficients for the tying scheme. */
  A2DScalar<N, T>
      e_rt_A, e_rt_B, e_st_C, e_st_D;
  /** Covariant transverse shear strains evaluated (using the tying points) at the various quadrature points. NOTE: the
   * values are the same for multiple quadrature points because we're assuming constant covariant transverse shear
   * strain conditions along the edges.*/
  A2DScalar<N, T>
      e_rt_rAs0tA, e_rt_rAs1tA, e_st_r0sAtA, e_st_r1sAtA;

  /** Energy contribution from the various quadrature points */
  A2DScalar<N, T>
      energy_r0s0t0,
      energy_r0s0t1,
      energy_r0s1t0,
      energy_r0s1t1,
      energy_r1s0t0,
      energy_r1s0t1,
      energy_r1s1t0,
      energy_r1s1t1;

  /** Strain energy summation objects */
  A2DScalar<N, T> sum_12, sum_123, sum_1234, sum_12345, sum_123456, sum_1234567;
  /** Strain energy object */
  A2DScalar<N, T> strain_energy;

 private:
  /** <h1>Expressions:</h1> */
  tying_shear_expr<N, T> e_rt_A_expression; /**< e_rt_A calculations: */
  tying_shear_expr<N, T> e_rt_B_expression; /**< e_rt_B calculations: */
  tying_shear_expr<N, T> e_st_C_expression; /**< e_st_C calculations: */
  tying_shear_expr<N, T> e_st_D_expression; /**< e_st_D calculations: */

  ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_rt_rAs0tA_expression; /**< Evaluate e_rt at s = quad_0; */
  ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_st_r0sAtA_expression; /**< Evaluate e_st at r = quad_0; */
  ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_rt_rAs1tA_expression; /**< Evaluate e_rt at s = quad_1; */
  ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_st_r1sAtA_expression; /**< Evaluate e_st at r = quad_1; */

  /* gr and gs calculations */
  g_alpha_expr<N, T> gr_rAs0t0_expression;
  g_alpha_expr<N, T> gr_rAs0t1_expression;
  g_alpha_expr<N, T> gr_rAs1t0_expression;
  g_alpha_expr<N, T> gr_rAs1t1_expression;

  g_alpha_expr<N, T> gs_r0sAt0_expression;
  g_alpha_expr<N, T> gs_r0sAt1_expression;
  g_alpha_expr<N, T> gs_r1sAt0_expression;
  g_alpha_expr<N, T> gs_r1sAt1_expression;

  /* ur and us calculations */
  u_alpha_expr<N, T> ur_rAs0t0_expression;
  u_alpha_expr<N, T> ur_rAs0t1_expression;
  u_alpha_expr<N, T> ur_rAs1t0_expression;
  u_alpha_expr<N, T> ur_rAs1t1_expression;

  u_alpha_expr<N, T> us_r0sAt0_expression;
  u_alpha_expr<N, T> us_r0sAt1_expression;
  u_alpha_expr<N, T> us_r1sAt0_expression;
  u_alpha_expr<N, T> us_r1sAt1_expression;

  /* Energy at the quadrature points: */
  strain_energy_expr<N, T> strain_energy_r0s0t0_expression;
  strain_energy_expr<N, T> strain_energy_r0s0t1_expression;
  strain_energy_expr<N, T> strain_energy_r0s1t0_expression;
  strain_energy_expr<N, T> strain_energy_r0s1t1_expression;
  strain_energy_expr<N, T> strain_energy_r1s0t0_expression;
  strain_energy_expr<N, T> strain_energy_r1s0t1_expression;
  strain_energy_expr<N, T> strain_energy_r1s1t0_expression;
  strain_energy_expr<N, T> strain_energy_r1s1t1_expression;

  /* Summation of the strain energy: */
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_123_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_1234_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12345_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_123456_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_1234567_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> sum_12345678_expression;

 private:
  /** <h1> Static variables: </h1>*/

  /** quadrature point values */
  constexpr static const T quad_0{-0.5773502691896257645091488}, quad_1{0.5773502691896257645091488};

  /** shape functions evaluated at the various points of interest */
  constexpr static const T Ni_rj_sk_data[100]{
      (1 - (-1)) * (1 - (-1)) * 0.25,      /**< N1 evaluated at r=-1,     s=-1 */
      (1 - (-1)) * (1 - quad_0) * 0.25,    /**< N1 evaluated at r=-1,     s=quad_0 */
      (1 - (-1)) * (1 - (0)) * 0.25,       /**< N1 evaluated at r=-1,     s=0 */
      (1 - (-1)) * (1 - quad_1) * 0.25,    /**< N1 evaluated at r=-1,     s=quad_1 */
      (1 - (-1)) * (1 - (1)) * 0.25,       /**< N1 evaluated at r=-1,     s=1 */

      (1 - quad_0) * (1 - (-1)) * 0.25,    /**< N1 evaluated at r=quad_0, s=-1 */
      (1 - quad_0) * (1 - quad_0) * 0.25,  /**< N1 evaluated at r=quad_0, s=quad_0 */
      (1 - quad_0) * (1 - (0)) * 0.25,     /**< N1 evaluated at r=quad_0, s=0 */
      (1 - quad_0) * (1 - quad_1) * 0.25,  /**< N1 evaluated at r=quad_0, s=quad_1 */
      (1 - quad_0) * (1 - (1)) * 0.25,     /**< N1 evaluated at r=quad_0, s=1 */

      (1 - (0)) * (1 - (-1)) * 0.25,       /**< N1 evaluated at r=0,      s=-1 */
      (1 - (0)) * (1 - quad_0) * 0.25,     /**< N1 evaluated at r=0,      s=quad_0 */
      (1 - (0)) * (1 - (0)) * 0.25,        /**< N1 evaluated at r=0,      s=0 */
      (1 - (0)) * (1 - quad_1) * 0.25,     /**< N1 evaluated at r=0,      s=quad_1 */
      (1 - (0)) * (1 - (1)) * 0.25,        /**< N1 evaluated at r=0,      s=1 */

      (1 - quad_1) * (1 - (-1)) * 0.25,    /**< N1 evaluated at r=quad_1, s=-1 */
      (1 - quad_1) * (1 - quad_0) * 0.25,  /**< N1 evaluated at r=quad_1, s=quad_0 */
      (1 - quad_1) * (1 - (0)) * 0.25,     /**< N1 evaluated at r=quad_1, s=0 */
      (1 - quad_1) * (1 - quad_1) * 0.25,  /**< N1 evaluated at r=quad_1, s=quad_1 */
      (1 - quad_1) * (1 - (1)) * 0.25,     /**< N1 evaluated at r=quad_1, s=1 */

      (1 - (1)) * (1 - (-1)) * 0.25,       /**< N1 evaluated at r=1,      s=-1 */
      (1 - (1)) * (1 - quad_0) * 0.25,     /**< N1 evaluated at r=1,      s=quad_0 */
      (1 - (1)) * (1 - (0)) * 0.25,        /**< N1 evaluated at r=1,      s=0 */
      (1 - (1)) * (1 - quad_1) * 0.25,     /**< N1 evaluated at r=1,      s=quad_1 */
      (1 - (1)) * (1 - (1)) * 0.25,        /**< N1 evaluated at r=1,      s=1 */



      (1 + (-1)) * (1 - (-1)) * 0.25,      /**< N2 evaluated at r=-1,     s=-1 */
      (1 + (-1)) * (1 - quad_0) * 0.25,    /**< N2 evaluated at r=-1,     s=quad_0 */
      (1 + (-1)) * (1 - (0)) * 0.25,       /**< N2 evaluated at r=-1,     s=0 */
      (1 + (-1)) * (1 - quad_1) * 0.25,    /**< N2 evaluated at r=-1,     s=quad_1 */
      (1 + (-1)) * (1 - (1)) * 0.25,       /**< N2 evaluated at r=-1,     s=1 */

      (1 + quad_0) * (1 - (-1)) * 0.25,    /**< N2 evaluated at r=quad_0, s=-1 */
      (1 + quad_0) * (1 - quad_0) * 0.25,  /**< N2 evaluated at r=quad_0, s=quad_0 */
      (1 + quad_0) * (1 - (0)) * 0.25,     /**< N2 evaluated at r=quad_0, s=0 */
      (1 + quad_0) * (1 - quad_1) * 0.25,  /**< N2 evaluated at r=quad_0, s=quad_1 */
      (1 + quad_0) * (1 - (1)) * 0.25,     /**< N2 evaluated at r=quad_0, s=1 */

      (1 + (0)) * (1 - (-1)) * 0.25,       /**< N2 evaluated at r=0,      s=-1 */
      (1 + (0)) * (1 - quad_0) * 0.25,     /**< N2 evaluated at r=0,      s=quad_0 */
      (1 + (0)) * (1 - (0)) * 0.25,        /**< N2 evaluated at r=0,      s=0 */
      (1 + (0)) * (1 - quad_1) * 0.25,     /**< N2 evaluated at r=0,      s=quad_1 */
      (1 + (0)) * (1 - (1)) * 0.25,        /**< N2 evaluated at r=0,      s=1 */

      (1 + quad_1) * (1 - (-1)) * 0.25,    /**< N2 evaluated at r=quad_1, s=-1 */
      (1 + quad_1) * (1 - quad_0) * 0.25,  /**< N2 evaluated at r=quad_1, s=quad_0 */
      (1 + quad_1) * (1 - (0)) * 0.25,     /**< N2 evaluated at r=quad_1, s=0 */
      (1 + quad_1) * (1 - quad_1) * 0.25,  /**< N2 evaluated at r=quad_1, s=quad_1 */
      (1 + quad_1) * (1 - (1)) * 0.25,     /**< N2 evaluated at r=quad_1, s=1 */

      (1 + (1)) * (1 - (-1)) * 0.25,       /**< N2 evaluated at r=1,      s=-1 */
      (1 + (1)) * (1 - quad_0) * 0.25,     /**< N2 evaluated at r=1,      s=quad_0 */
      (1 + (1)) * (1 - (0)) * 0.25,        /**< N2 evaluated at r=1,      s=0 */
      (1 + (1)) * (1 - quad_1) * 0.25,     /**< N2 evaluated at r=1,      s=quad_1 */
      (1 + (1)) * (1 - (1)) * 0.25,        /**< N2 evaluated at r=1,      s=1 */



      (1 + (-1)) * (1 + (-1)) * 0.25,      /**< N3 evaluated at r=-1,     s=-1 */
      (1 + (-1)) * (1 + quad_0) * 0.25,    /**< N3 evaluated at r=-1,     s=quad_0 */
      (1 + (-1)) * (1 + (0)) * 0.25,       /**< N3 evaluated at r=-1,     s=0 */
      (1 + (-1)) * (1 + quad_1) * 0.25,    /**< N3 evaluated at r=-1,     s=quad_1 */
      (1 + (-1)) * (1 + (1)) * 0.25,       /**< N3 evaluated at r=-1,     s=1 */

      (1 + quad_0) * (1 + (-1)) * 0.25,    /**< N3 evaluated at r=quad_0, s=-1 */
      (1 + quad_0) * (1 + quad_0) * 0.25,  /**< N3 evaluated at r=quad_0, s=quad_0 */
      (1 + quad_0) * (1 + (0)) * 0.25,     /**< N3 evaluated at r=quad_0, s=0 */
      (1 + quad_0) * (1 + quad_1) * 0.25,  /**< N3 evaluated at r=quad_0, s=quad_1 */
      (1 + quad_0) * (1 + (1)) * 0.25,     /**< N3 evaluated at r=quad_0, s=1 */

      (1 + (0)) * (1 + (-1)) * 0.25,       /**< N3 evaluated at r=0,      s=-1 */
      (1 + (0)) * (1 + quad_0) * 0.25,     /**< N3 evaluated at r=0,      s=quad_0 */
      (1 + (0)) * (1 + (0)) * 0.25,        /**< N3 evaluated at r=0,      s=0 */
      (1 + (0)) * (1 + quad_1) * 0.25,     /**< N3 evaluated at r=0,      s=quad_1 */
      (1 + (0)) * (1 + (1)) * 0.25,        /**< N3 evaluated at r=0,      s=1 */

      (1 + quad_1) * (1 + (-1)) * 0.25,    /**< N3 evaluated at r=quad_1, s=-1 */
      (1 + quad_1) * (1 + quad_0) * 0.25,  /**< N3 evaluated at r=quad_1, s=quad_0 */
      (1 + quad_1) * (1 + (0)) * 0.25,     /**< N3 evaluated at r=quad_1, s=0 */
      (1 + quad_1) * (1 + quad_1) * 0.25,  /**< N3 evaluated at r=quad_1, s=quad_1 */
      (1 + quad_1) * (1 + (1)) * 0.25,     /**< N3 evaluated at r=quad_1, s=1 */

      (1 + (1)) * (1 + (-1)) * 0.25,       /**< N3 evaluated at r=1,      s=-1 */
      (1 + (1)) * (1 + quad_0) * 0.25,     /**< N3 evaluated at r=1,      s=quad_0 */
      (1 + (1)) * (1 + (0)) * 0.25,        /**< N3 evaluated at r=1,      s=0 */
      (1 + (1)) * (1 + quad_1) * 0.25,     /**< N3 evaluated at r=1,      s=quad_1 */
      (1 + (1)) * (1 + (1)) * 0.25,        /**< N3 evaluated at r=1,      s=1 */



      (1 - (-1)) * (1 + (-1)) * 0.25,      /**< N4 evaluated at r=-1,     s=-1 */
      (1 - (-1)) * (1 + quad_0) * 0.25,    /**< N4 evaluated at r=-1,     s=quad_0 */
      (1 - (-1)) * (1 + (0)) * 0.25,       /**< N4 evaluated at r=-1,     s=0 */
      (1 - (-1)) * (1 + quad_1) * 0.25,    /**< N4 evaluated at r=-1,     s=quad_1 */
      (1 - (-1)) * (1 + (1)) * 0.25,       /**< N4 evaluated at r=-1,     s=1 */

      (1 - quad_0) * (1 + (-1)) * 0.25,    /**< N4 evaluated at r=quad_0, s=-1 */
      (1 - quad_0) * (1 + quad_0) * 0.25,  /**< N4 evaluated at r=quad_0, s=quad_0 */
      (1 - quad_0) * (1 + (0)) * 0.25,     /**< N4 evaluated at r=quad_0, s=0 */
      (1 - quad_0) * (1 + quad_1) * 0.25,  /**< N4 evaluated at r=quad_0, s=quad_1 */
      (1 - quad_0) * (1 + (1)) * 0.25,     /**< N4 evaluated at r=quad_0, s=1 */

      (1 - (0)) * (1 + (-1)) * 0.25,       /**< N4 evaluated at r=0,      s=-1 */
      (1 - (0)) * (1 + quad_0) * 0.25,     /**< N4 evaluated at r=0,      s=quad_0 */
      (1 - (0)) * (1 + (0)) * 0.25,        /**< N4 evaluated at r=0,      s=0 */
      (1 - (0)) * (1 + quad_1) * 0.25,     /**< N4 evaluated at r=0,      s=quad_1 */
      (1 - (0)) * (1 + (1)) * 0.25,        /**< N4 evaluated at r=0,      s=1 */

      (1 - quad_1) * (1 + (-1)) * 0.25,    /**< N4 evaluated at r=quad_1, s=-1 */
      (1 - quad_1) * (1 + quad_0) * 0.25,  /**< N4 evaluated at r=quad_1, s=quad_0 */
      (1 - quad_1) * (1 + (0)) * 0.25,     /**< N4 evaluated at r=quad_1, s=0 */
      (1 - quad_1) * (1 + quad_1) * 0.25,  /**< N4 evaluated at r=quad_1, s=quad_1 */
      (1 - quad_1) * (1 + (1)) * 0.25,     /**< N4 evaluated at r=quad_1, s=1 */

      (1 - (1)) * (1 + (-1)) * 0.25,       /**< N4 evaluated at r=1,      s=-1 */
      (1 - (1)) * (1 + quad_0) * 0.25,     /**< N4 evaluated at r=1,      s=quad_0 */
      (1 - (1)) * (1 + (0)) * 0.25,        /**< N4 evaluated at r=1,      s=0 */
      (1 - (1)) * (1 + quad_1) * 0.25,     /**< N4 evaluated at r=1,      s=quad_1 */
      (1 - (1)) * (1 + (1)) * 0.25,        /**< N4 evaluated at r=1,      s=1 */
  };
 public:
  /**
   * @brief getter for shape functions evaluated at the various points of interest
   *
   * @param i:  denotes which shape function is being queried
   * @param j:  denotes the value for r (0 for r=-1; 1 for r=quad_0; 2 for r=0; 3 for r=quad_1; and 4 for r=1)
   * @param k:  denotes the value for s (0 for s=-1; 1 for s=quad_0; 2 for s=0; 3 for s=quad_1; and 4 for s=1)
   * */
  constexpr static T Ni_rj_sk(const int i, const int j, const int k) {
    return Ni_rj_sk_data[i * 25 + j * 5 + k];
  }

 private:
  /** shape function derivatives with respect to r and s, evaluated at the various points of interest */
  constexpr static const T dNi_d_alpha_alpha_j_data[32]{
      /** shape function derivative with respect to r, evaluated at s=-1, quad_0, quad_1, and 1 */
      -0.25 * (1 - (-1)),   /**< dN1/dr evaluated at s=-1 */
      -0.25 * (1 - quad_0), /**< dN1/dr evaluated at s=quad_0 */
      -0.25 * (1 - quad_1), /**< dN1/dr evaluated at s=quad_1 */
      -0.25 * (1 - (1)),    /**< dN1/dr evaluated at s=1 */

      0.25 * (1 - (-1)),    /**< dN2/dr evaluated at s=-1 */
      0.25 * (1 - quad_0),  /**< dN2/dr evaluated at s=quad_0 */
      0.25 * (1 - quad_1),  /**< dN2/dr evaluated at s=quad_1 */
      0.25 * (1 - (1)),     /**< dN2/dr evaluated at s=1 */

      0.25 * (1 + (-1)),    /**< dN3/dr evaluated at s=-1 */
      0.25 * (1 + quad_0),  /**< dN3/dr evaluated at s=quad_0 */
      0.25 * (1 + quad_1),  /**< dN3/dr evaluated at s=quad_1 */
      0.25 * (1 + (1)),     /**< dN3/dr evaluated at s=1 */

      -0.25 * (1 + (-1)),   /**< dN4/dr evaluated at s=-1 */
      -0.25 * (1 + quad_0), /**< dN4/dr evaluated at s=quad_0 */
      -0.25 * (1 + quad_1), /**< dN4/dr evaluated at s=quad_1 */
      -0.25 * (1 + (1)),    /**< dN4/dr evaluated at s=1 */

      /** shape function derivative with respect to s, evaluated at r=-1, quad_0, quad_1, and 1 */
      -0.25 * (1 - (-1)),   /**< dN1/ds evaluated at r=-1 */
      -0.25 * (1 - quad_0), /**< dN1/ds evaluated at r=quad_0 */
      -0.25 * (1 - quad_1), /**< dN1/ds evaluated at r=quad_1 */
      -0.25 * (1 - (1)),    /**< dN1/ds evaluated at r=1 */

      -0.25 * (1 + (-1)),   /**< dN2/ds evaluated at r=-1 */
      -0.25 * (1 + quad_0), /**< dN2/ds evaluated at r=quad_0 */
      -0.25 * (1 + quad_1), /**< dN2/ds evaluated at r=quad_1 */
      -0.25 * (1 + (1)),    /**< dN2/ds evaluated at r=1 */

      0.25 * (1 + (-1)),    /**< dN3/ds evaluated at r=-1 */
      0.25 * (1 + quad_0),  /**< dN3/ds evaluated at r=quad_0 */
      0.25 * (1 + quad_1),  /**< dN3/ds evaluated at r=quad_1 */
      0.25 * (1 + (1)),     /**< dN3/ds evaluated at r=1 */

      0.25 * (1 - (-1)),    /**< dN4/ds evaluated at r=-1 */
      0.25 * (1 - quad_0),  /**< dN4/ds evaluated at r=quad_0 */
      0.25 * (1 - quad_1),  /**< dN4/ds evaluated at r=quad_1 */
      0.25 * (1 - (1))      /**< dN4/ds evaluated at r=1 */
  };
 public:
  /** @brief getter for shape function derivatives evaluated at the various points of interest
   *
   * @param alpha:  denotes the derivative of the shape function, a value of 0 corresponds to  dNi/dr while a value of
   *                1 corresponds to dNi/ds.
   * @param i:      denotes which shape function to use.
   * @param j:      denotes the index of the value of the parameter <u>not</u> associated with the alpha parameter.  A
   *                value of 0 corresponds with -1; 1 with quad_0; 2 with quad_1; and 3 with 1.  For example,
   *                dNi_d_alpha_alpha_j(0, i, 1) corresponds to dNi/dr evaluated at s=quad_0.
   * */
  constexpr static T dNi_d_alpha_alpha_j(const int alpha, const int i, const int j) {
    // TODO: see if template <int alpha, int i, int j> improves performance (or is viable)
    return dNi_d_alpha_alpha_j_data[alpha * 16 + i * 4 + j];
  }
};

/* Helper Classes */

template <int N, typename T>
class g_alpha_expr {
 public:
  /**
   * @brief Computes the g<sub>&alpha</sub> vector for the given situation and element.
   *
   *
   * @param alpha:            denotes the variant of the g<sub>&alpha</sub> vector, a value of 0 corresponds to the
   *                          g<sub>r</sub> vector while a value of 1 corresponds to the g<sub>s</sub> vector.
   * @param n_alpha_var_ind:  denotes which value to use (0 for ~&alpha =-1; 1 for ~&alpha =quad_0; 2 for
   *                          ~&alpha =quad_1; or 3 for ~&alpha=1).  This corresponds to the value of the index
   *                          <u>not</u> represented by the alpha parameter.  For example, alpha=0, n_alpha_var_ind=0
   *                          corresponds to evaluating g<sub>r</sub>(s,t) with s=-1; and alpha=1, n_alpha_var_ind=2
   *                          corresponds to evaluating g<sub>s</sub>(r,t) with r=quad_1.
   * @param t:                the value of the t parametric coordinate.
   * @param element:          the MITC4 element object for which the g<sub>&alpha</sub> vector is being computed.
   * @param result:           an A2DVec where the resulting g<sub>&alpha</sub> vector should be stored.
   * */
  g_alpha_expr(const int alpha,
               const int n_alpha_var_ind,
               const T& t,
//               const ShellElementMITC4<N, T>& element,
               ShellNodeMITC<N, T>& node1, ShellNodeMITC<N, T>& node2,
               ShellNodeMITC<N, T>& node3, ShellNodeMITC<N, T>& node4,
               A2DVec<N, Vec<T, 3>>& result)
      : /* Instantiations: */
      node1_contribution_unscaled(node1_contribution_unscaled_v, node1_contribution_unscaled_bv),
      node2_contribution_unscaled(node2_contribution_unscaled_v, node2_contribution_unscaled_bv),
      node3_contribution_unscaled(node3_contribution_unscaled_v, node3_contribution_unscaled_bv),
      node4_contribution_unscaled(node4_contribution_unscaled_v, node4_contribution_unscaled_bv),
      node1_contribution(node1_contribution_v, node1_contribution_bv),
      node2_contribution(node2_contribution_v, node2_contribution_bv),
      node3_contribution(node3_contribution_v, node3_contribution_bv),
      node4_contribution(node4_contribution_v, node4_contribution_bv),
      n1c_n2c(n1c_n2c_v, n1c_n2c_bv),
      n1c_n2c_n3c(n1c_n2c_n3c_v, n1c_n2c_n3c_bv),
      /* Operations: */
      /* Calculations for first node. */
      node1_thickness_scale_expression(node1.thickness, t * 0.5, node1_scaled_thickness),
      node1_addition_expression(node1_scaled_thickness,
                                node1.shell_director,
                                node1.position,
                                node1_contribution_unscaled),
      node1_scale_expression(node1_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 0, n_alpha_var_ind),
                             node1_contribution),
      /* Calculations for second node. */
      node2_thickness_scale_expression(node2.thickness, t * 0.5, node2_scaled_thickness),
      node2_addition_expression(node2_scaled_thickness,
                                node2.shell_director,
                                node2.position,
                                node2_contribution_unscaled),
      node2_scale_expression(node2_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 1, n_alpha_var_ind),
                             node2_contribution),
      /* Calculations for third node. */
      node3_thickness_scale_expression(node3.thickness, t * 0.5, node3_scaled_thickness),
      node3_addition_expression(node3_scaled_thickness,
                                node3.shell_director,
                                node3.position,
                                node3_contribution_unscaled),
      node3_scale_expression(node3_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 2, n_alpha_var_ind),
                             node3_contribution),
      /* Calculations for fourth node. */
      node4_thickness_scale_expression(node4.thickness, t * 0.5, node4_scaled_thickness),
      node4_addition_expression(node4_scaled_thickness,
                                node4.shell_director,
                                node4.position,
                                node4_contribution_unscaled),
      node4_scale_expression(node4_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 3, n_alpha_var_ind),
                             node4_contribution),

      /* Sum together the components. */
      sum_12_expression(1.0, node1_contribution, node2_contribution, n1c_n2c),
      sum_123_expression(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c),
      sum_1234_expression(1.0, n1c_n2c_n3c, node4_contribution, result) {};

  void reverse() {
    /* Sum component calls: */
    sum_1234_expression.reverse();
    sum_123_expression.reverse();
    sum_12_expression.reverse();

    /* Node expression calls: */
    node4_scale_expression.reverse();
    node4_addition_expression.reverse();
    node4_thickness_scale_expression.reverse();
    node3_scale_expression.reverse();
    node3_addition_expression.reverse();
    node3_thickness_scale_expression.reverse();
    node2_scale_expression.reverse();
    node2_addition_expression.reverse();
    node2_thickness_scale_expression.reverse();
    node1_scale_expression.reverse();
    node1_addition_expression.reverse();
    node1_thickness_scale_expression.reverse();
  };

  void hforward() {
    /* Node expression calls: */
    node1_thickness_scale_expression.hforward();
    node1_addition_expression.hforward();
    node1_scale_expression.hforward();
    node2_thickness_scale_expression.hforward();
    node2_addition_expression.hforward();
    node2_scale_expression.hforward();
    node3_thickness_scale_expression.hforward();
    node3_addition_expression.hforward();
    node3_scale_expression.hforward();
    node4_thickness_scale_expression.hforward();
    node4_addition_expression.hforward();
    node4_scale_expression.hforward();

    /* Sum component calls: */
    sum_12_expression.hforward();
    sum_123_expression.hforward();
    sum_1234_expression.hforward();
  };

  void hreverse() {
    /* Sum component calls: */
    sum_1234_expression.hreverse();
    sum_123_expression.hreverse();
    sum_12_expression.hreverse();

    /* Node expression calls: */
    node4_scale_expression.hhreverse();
    node4_addition_expression.hreverse();
    node4_thickness_scale_expression.hreverse();
    node3_scale_expression.hreverse();
    node3_addition_expression.hreverse();
    node3_thickness_scale_expression.hreverse();
    node2_scale_expression.hreverse();
    node2_addition_expression.hreverse();
    node2_thickness_scale_expression.hreverse();
    node1_scale_expression.hreverse();
    node1_addition_expression.hreverse();
    node1_thickness_scale_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      node1_contribution_unscaled_v, node1_contribution_unscaled_bv,
      node2_contribution_unscaled_v, node2_contribution_unscaled_bv,
      node3_contribution_unscaled_v, node3_contribution_unscaled_bv,
      node4_contribution_unscaled_v, node4_contribution_unscaled_bv,
      node1_contribution_v, node1_contribution_bv,
      node2_contribution_v, node2_contribution_bv,
      node3_contribution_v, node3_contribution_bv,
      node4_contribution_v, node4_contribution_bv,
      n1c_n2c_v, n1c_n2c_bv,
      n1c_n2c_n3c_v, n1c_n2c_n3c_bv;

  /* A2D Objects: */
  A2DVec<N, Vec<T, 3>>
      node1_contribution_unscaled,
      node2_contribution_unscaled,
      node3_contribution_unscaled,
      node4_contribution_unscaled;
  A2DVec<N, Vec<T, 3>>
      node1_contribution,
      node2_contribution,
      node3_contribution,
      node4_contribution;
  A2DScalar<N, T>
      node1_scaled_thickness,
      node2_scaled_thickness,
      node3_scaled_thickness,
      node4_scaled_thickness;
  A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

  /* Expressions: */
  A2DScalarScalarMultExpr<N, T> node1_thickness_scale_expression;
  Vec3VecA2DScalarAxpyExpr<N, T> node1_addition_expression;
  A2DVec3ScaleExpr<N, T> node1_scale_expression;

  A2DScalarScalarMultExpr<N, T> node2_thickness_scale_expression;
  Vec3VecA2DScalarAxpyExpr<N, T> node2_addition_expression;
  A2DVec3ScaleExpr<N, T> node2_scale_expression;

  A2DScalarScalarMultExpr<N, T> node3_thickness_scale_expression;
  Vec3VecA2DScalarAxpyExpr<N, T> node3_addition_expression;
  A2DVec3ScaleExpr<N, T> node3_scale_expression;

  A2DScalarScalarMultExpr<N, T> node4_thickness_scale_expression;
  Vec3VecA2DScalarAxpyExpr<N, T> node4_addition_expression;
  A2DVec3ScaleExpr<N, T> node4_scale_expression;

  A2DVec3A2DVecScalarAxpyExpr<N, T>
      sum_12_expression,
      sum_123_expression,
      sum_1234_expression;
};

template <int N, typename T>
class u_alpha_expr {
 public:
  /**
   * @brief Computes the du/d&alpha vector for the given situation and element.
   *
   *
   * @param alpha:            denotes the variant of the du/d&alpha vector, a value of 0 corresponds to the du/dr
   *                          vector while a value of 1 corresponds to the du/ds vector.
   * @param n_alpha_var_ind:  denotes which value to use (0 for ~&alpha =-1; 1 for ~&alpha =quad_0; 2 for
   *                          ~&alpha =quad_1; or 3 for ~&alpha=1).  This corresponds to the value of the index
   *                          <u>not</u> represented by the alpha parameter.  For example, alpha=0, n_alpha_var_ind=0
   *                          corresponds to evaluating du/dr(s,t) with s=-1; and alpha=1, n_alpha_var_ind=2
   *                          corresponds to evaluating du/ds(r,t) with r=quad_1.
   * @param t:                the value of the t parametric coordinate.
   * @param element:          the MITC4 element object for which the du/d&alpha vector is being computed.
   * @param result:           an A2DVec where the resulting du/d&alpha vector should be stored.
   * */
  u_alpha_expr(const int alpha,
               const int n_alpha_var_ind,
               const T& t,
//               const ShellElementMITC4<N, T>& element,
               ShellNodeMITC<N, T>& node1, ShellNodeMITC<N, T>& node2,
               ShellNodeMITC<N, T>& node3, ShellNodeMITC<N, T>& node4,
               A2DVec<N, Vec<T, 3>>& result)
      : /* Instantiations: */
      node1_phi(node1_phi_v, node1_phi_bv),
      node2_phi(node2_phi_v, node2_phi_bv),
      node3_phi(node3_phi_v, node3_phi_bv),
      node4_phi(node4_phi_v, node4_phi_bv),
      node1_contribution_unscaled(node1_contribution_unscaled_v, node1_contribution_unscaled_bv),
      node2_contribution_unscaled(node2_contribution_unscaled_v, node2_contribution_unscaled_bv),
      node3_contribution_unscaled(node3_contribution_unscaled_v, node3_contribution_unscaled_bv),
      node4_contribution_unscaled(node4_contribution_unscaled_v, node4_contribution_unscaled_bv),
      node1_contribution(node1_contribution_v, node1_contribution_bv),
      node2_contribution(node2_contribution_v, node2_contribution_bv),
      node3_contribution(node3_contribution_v, node3_contribution_bv),
      node4_contribution(node4_contribution_v, node4_contribution_bv),
      n1c_n2c(n1c_n2c_v, n1c_n2c_bv),
      n1c_n2c_n3c(n1c_n2c_n3c_v, n1c_n2c_n3c_bv),
      /* Operations: */
      /* Calculations for first node. */
      node1_phi_expression(node1.rotation, node1.shell_director, node1_phi),
      node1_thickness_scale_expression(node1.thickness, t * 0.5, node1_scaled_thickness),
      node1_addition_expression(node1_scaled_thickness,
                                node1_phi,
                                node1.displacement,
                                node1_contribution_unscaled),
      node1_scale_expression(node1_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 0, n_alpha_var_ind),
                             node1_contribution),
      /* Calculations for second node. */
      node2_phi_expression(node2.rotation, node2.shell_director, node2_phi),
      node2_thickness_scale_expression(node2.thickness, t * 0.5, node2_scaled_thickness),
      node2_addition_expression(node2_scaled_thickness,
                                node2_phi,
                                node2.displacement,
                                node2_contribution_unscaled),
      node2_scale_expression(node2_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 1, n_alpha_var_ind),
                             node2_contribution),
      /* Calculations for third node. */
      node3_phi_expression(node3.rotation, node3.shell_director, node3_phi),
      node3_thickness_scale_expression(node3.thickness, t * 0.5, node3_scaled_thickness),
      node3_addition_expression(node3_scaled_thickness,
                                node3_phi,
                                node3.displacement,
                                node3_contribution_unscaled),
      node3_scale_expression(node3_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 2, n_alpha_var_ind),
                             node3_contribution),
      /* Calculations for fourth node. */
      node4_phi_expression(node4.rotation, node4.shell_director, node4_phi),
      node4_thickness_scale_expression(node4.thickness, t * 0.5, node4_scaled_thickness),
      node4_addition_expression(node4_scaled_thickness,
                                node4_phi,
                                node4.displacement,
                                node4_contribution_unscaled),
      node4_scale_expression(node4_contribution_unscaled,
                             ShellElementMITC4<N, T>::dNi_d_alpha_alpha_j(alpha, 3, n_alpha_var_ind),
                             node4_contribution),

      /* Sum together the components. */
      sum_12_expression(1.0, node1_contribution, node2_contribution, n1c_n2c),
      sum_123_expression(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c),
      sum_1234_expression(1.0, n1c_n2c_n3c, node4_contribution, result) {};

  void reverse() {
    /* Sum component calls: */
    sum_1234_expression.reverse();
    sum_123_expression.reverse();
    sum_12_expression.reverse();

    /* Node expression calls: */
    node4_scale_expression.reverse();
    node4_addition_expression.reverse();
    node4_thickness_scale_expression.reverse();
    node4_phi_expression.reverse();
    node3_scale_expression.reverse();
    node3_addition_expression.reverse();
    node3_thickness_scale_expression.reverse();
    node3_phi_expression.reverse();
    node2_scale_expression.reverse();
    node2_addition_expression.reverse();
    node2_thickness_scale_expression.reverse();
    node2_phi_expression.reverse();
    node1_scale_expression.reverse();
    node1_addition_expression.reverse();
    node1_thickness_scale_expression.reverse();
    node1_phi_expression.reverse();
  };

  void hforward() {
    /* Node expression calls: */
    node1_phi_expression.hforward();
    node1_thickness_scale_expression.hforward();
    node1_addition_expression.hforward();
    node1_scale_expression.hforward();
    node2_phi_expression.hforward();
    node2_thickness_scale_expression.hforward();
    node2_addition_expression.hforward();
    node2_scale_expression.hforward();
    node3_phi_expression.hforward();
    node3_thickness_scale_expression.hforward();
    node3_addition_expression.hforward();
    node3_scale_expression.hforward();
    node4_phi_expression.hforward();
    node4_thickness_scale_expression.hforward();
    node4_addition_expression.hforward();
    node4_scale_expression.hforward();

    /* Sum component calls: */
    sum_12_expression.hforward();
    sum_123_expression.hforward();
    sum_1234_expression.hforward();
  };

  void hreverse() {
    /* Sum component calls: */
    sum_1234_expression.hreverse();
    sum_123_expression.hreverse();
    sum_12_expression.hreverse();

    /* Node expression calls: */
    node4_scale_expression.hhreverse();
    node4_addition_expression.hreverse();
    node4_thickness_scale_expression.hreverse();
    node4_phi_expression.hreverse();
    node3_scale_expression.hreverse();
    node3_addition_expression.hreverse();
    node3_thickness_scale_expression.hreverse();
    node3_phi_expression.hreverse();
    node2_scale_expression.hreverse();
    node2_addition_expression.hreverse();
    node2_thickness_scale_expression.hreverse();
    node2_phi_expression.hreverse();
    node1_scale_expression.hreverse();
    node1_addition_expression.hreverse();
    node1_thickness_scale_expression.hreverse();
    node1_phi_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      node1_phi_v, node1_phi_bv,
      node2_phi_v, node2_phi_bv,
      node3_phi_v, node3_phi_bv,
      node4_phi_v, node4_phi_bv;
  Vec<T, 3>
      node1_contribution_unscaled_v, node1_contribution_unscaled_bv,
      node2_contribution_unscaled_v, node2_contribution_unscaled_bv,
      node3_contribution_unscaled_v, node3_contribution_unscaled_bv,
      node4_contribution_unscaled_v, node4_contribution_unscaled_bv,
      node1_contribution_v, node1_contribution_bv,
      node2_contribution_v, node2_contribution_bv,
      node3_contribution_v, node3_contribution_bv,
      node4_contribution_v, node4_contribution_bv,
      n1c_n2c_v, n1c_n2c_bv,
      n1c_n2c_n3c_v, n1c_n2c_n3c_bv;

  /* A2D Objects: */
  A2DVec<N, Vec<T, 3>> node1_phi;
  A2DVec<N, Vec<T, 3>> node2_phi;
  A2DVec<N, Vec<T, 3>> node3_phi;
  A2DVec<N, Vec<T, 3>> node4_phi;
  A2DVec<N, Vec<T, 3>>
      node1_contribution_unscaled,
      node2_contribution_unscaled,
      node3_contribution_unscaled,
      node4_contribution_unscaled;
  A2DVec<N, Vec<T, 3>>
      node1_contribution,
      node2_contribution,
      node3_contribution,
      node4_contribution;
  A2DScalar<N, T>
      node1_scaled_thickness,
      node2_scaled_thickness,
      node3_scaled_thickness,
      node4_scaled_thickness;
  A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

  /* Expressions: */
  A2DVec3CrossVecExpr<N, T> node1_phi_expression;
  A2DScalarScalarMultExpr<N, T> node1_thickness_scale_expression;
  A2DVec3A2DVecA2DScalarAxpyExpr<N, T> node1_addition_expression;
  A2DVec3ScaleExpr<N, T> node1_scale_expression;

  A2DVec3CrossVecExpr<N, T> node2_phi_expression;
  A2DScalarScalarMultExpr<N, T> node2_thickness_scale_expression;
  A2DVec3A2DVecA2DScalarAxpyExpr<N, T> node2_addition_expression;
  A2DVec3ScaleExpr<N, T> node2_scale_expression;

  A2DVec3CrossVecExpr<N, T> node3_phi_expression;
  A2DScalarScalarMultExpr<N, T> node3_thickness_scale_expression;
  A2DVec3A2DVecA2DScalarAxpyExpr<N, T> node3_addition_expression;
  A2DVec3ScaleExpr<N, T> node3_scale_expression;

  A2DVec3CrossVecExpr<N, T> node4_phi_expression;
  A2DScalarScalarMultExpr<N, T> node4_thickness_scale_expression;
  A2DVec3A2DVecA2DScalarAxpyExpr<N, T> node4_addition_expression;
  A2DVec3ScaleExpr<N, T> node4_scale_expression;

  A2DVec3A2DVecScalarAxpyExpr<N, T>
      sum_12_expression,
      sum_123_expression,
      sum_1234_expression;
};

template <int N, typename T>
class g_t_expr {
 public:
  /**
   * @brief Computes the g<sub>t</sub> vector for the given situation and element.
   *
   *
   * @param r_ind:      denotes which value to use for r (0 for r=-1; 1 for r=quad_0; 2 for r=0; 3 for r=quad_1; and 4
   *                    for r=1)
   * @param s_ind:      denotes which value to use for s (0 for s=-1; 1 for s=quad_0; 2 for s=0; 3 for s=quad_1; and 4
   *                    for s=1)
   * @param element:    the MITC4 element object for which the g<sub>t</sub> vector is being computed.
   * @param result:     an A2DVec where the resulting g<sub>t</sub> vector should be stored.
   * */
  g_t_expr(const int r_ind,
           const int s_ind,
//           const ShellElementMITC4<N, T>& element,
           ShellNodeMITC<N, T>& node1, ShellNodeMITC<N, T>& node2,
           ShellNodeMITC<N, T>& node3, ShellNodeMITC<N, T>& node4,
           A2DVec<N, Vec<T, 3>>& result)
      : /* Instantiations: */
      node1_contribution(node1_contribution_v, node1_contribution_bv),
      node2_contribution(node2_contribution_v, node2_contribution_bv),
      node3_contribution(node3_contribution_v, node3_contribution_bv),
      node4_contribution(node4_contribution_v, node4_contribution_bv),
      n1c_n2c(n1c_n2c_v, n1c_n2c_bv),
      n1c_n2c_n3c(n1c_n2c_n3c_v, n1c_n2c_n3c_bv),
      /* Operations: */
      /* Calculations for first node contribution */
      node1_thickness_scale_expression(node1.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(0, r_ind, s_ind) * 0.5,
                                       node1_scaled_thickness),
      node1_scale_expression(node1.shell_director, node1_scaled_thickness, node1_contribution),
      /* Calculations for second node contribution */
      node2_thickness_scale_expression(node2.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(1, r_ind, s_ind) * 0.5,
                                       node2_scaled_thickness),
      node2_scale_expression(node2.shell_director, node2_scaled_thickness, node2_contribution),
      /* Calculations for third node contribution */
      node3_thickness_scale_expression(node3.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(2, r_ind, s_ind) * 0.5,
                                       node3_scaled_thickness),
      node3_scale_expression(node3.shell_director, node3_scaled_thickness, node3_contribution),
      /* Calculations for fourth node contribution */
      node4_thickness_scale_expression(node4.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(3, r_ind, s_ind) * 0.5,
                                       node4_scaled_thickness),
      node4_scale_expression(node4.shell_director, node4_scaled_thickness, node4_contribution),
      /* Sum together the components. */
      sum_12_expression(1.0, node1_contribution, node2_contribution, n1c_n2c),
      sum_123_expression(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c),
      sum_1234_expression(1.0, n1c_n2c_n3c, node4_contribution, result) {};

  void reverse() {
    /* Sum component calls: */
    sum_1234_expression.reverse();
    sum_123_expression.reverse();
    sum_12_expression.reverse();

    /* Node expression calls: */
    node4_scale_expression.reverse();
    node4_thickness_scale_expression.reverse();
    node3_scale_expression.reverse();
    node3_thickness_scale_expression.reverse();
    node2_scale_expression.reverse();
    node2_thickness_scale_expression.reverse();
    node1_scale_expression.reverse();
    node1_thickness_scale_expression.reverse();
  };

  void hforward() {
    /* Node expression calls: */
    node1_thickness_scale_expression.hforward();
    node1_scale_expression.hforward();
    node2_thickness_scale_expression.hforward();
    node2_scale_expression.hforward();
    node3_thickness_scale_expression.hforward();
    node3_scale_expression.hforward();
    node4_thickness_scale_expression.hforward();
    node4_scale_expression.hforward();

    /* Sum component calls: */
    sum_12_expression.hforward();
    sum_123_expression.hforward();
    sum_1234_expression.hforward();
  };

  void hreverse() {
    /* Sum component calls: */
    sum_1234_expression.hreverse();
    sum_123_expression.hreverse();
    sum_12_expression.hreverse();

    /* Node expression calls: */
    node4_scale_expression.hhreverse();
    node4_thickness_scale_expression.hreverse();
    node3_scale_expression.hreverse();
    node3_thickness_scale_expression.hreverse();
    node2_scale_expression.hreverse();
    node2_thickness_scale_expression.hreverse();
    node1_scale_expression.hreverse();
    node1_thickness_scale_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      node1_contribution_v, node1_contribution_bv,
      node2_contribution_v, node2_contribution_bv,
      node3_contribution_v, node3_contribution_bv,
      node4_contribution_v, node4_contribution_bv,
      n1c_n2c_v, n1c_n2c_bv,
      n1c_n2c_n3c_v, n1c_n2c_n3c_bv;

  /* A2D Objects: */
  A2DScalar<N, T>
      node1_scaled_thickness,
      node2_scaled_thickness,
      node3_scaled_thickness,
      node4_scaled_thickness;
  A2DVec<N, Vec<T, 3>>
      node1_contribution,
      node2_contribution,
      node3_contribution,
      node4_contribution;
  A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

  /* Expressions: */
  A2DScalarScalarMultExpr<N, T> node1_thickness_scale_expression;
  Vec3A2DScaleExpr<N, T> node1_scale_expression;

  A2DScalarScalarMultExpr<N, T> node2_thickness_scale_expression;
  Vec3A2DScaleExpr<N, T> node2_scale_expression;

  A2DScalarScalarMultExpr<N, T> node3_thickness_scale_expression;
  Vec3A2DScaleExpr<N, T> node3_scale_expression;

  A2DScalarScalarMultExpr<N, T> node4_thickness_scale_expression;
  Vec3A2DScaleExpr<N, T> node4_scale_expression;

  A2DVec3A2DVecScalarAxpyExpr<N, T>
      sum_12_expression,
      sum_123_expression,
      sum_1234_expression;
};

template <int N, typename T>
class u_t_expr {
 public:
  /**
   * @brief Computes the du/dt vector for the given situation and element.
   *
   *
   * @param r_ind:      denotes which value to use for r (0 for r=-1; 1 for r=quad_0; 2 for r=0; 3 for r=quad_1; and 4
   *                    for r=1)
   * @param s_ind:      denotes which value to use for s (0 for s=-1; 1 for s=quad_0; 2 for s=0; 3 for s=quad_1; and 4
   *                    for s=1)
   * @param element:    the MITC4 element object for which the du/dt vector is being computed.
   * @param result:     an A2DVec where the resulting du/dt vector should be stored.
   * */
  u_t_expr(const int r_ind,
           const int s_ind,
//           const ShellElementMITC4<N, T>& element,
           ShellNodeMITC<N, T>& node1, ShellNodeMITC<N, T>& node2,
           ShellNodeMITC<N, T>& node3, ShellNodeMITC<N, T>& node4,
           A2DVec<N, Vec<T, 3>>& result)
      : /* Instantiations: */
      node1_phi(node1_phi_v, node1_phi_bv),
      node2_phi(node2_phi_v, node2_phi_bv),
      node3_phi(node3_phi_v, node3_phi_bv),
      node4_phi(node4_phi_v, node4_phi_bv),
      node1_contribution(node1_contribution_v, node1_contribution_bv),
      node2_contribution(node2_contribution_v, node2_contribution_bv),
      node3_contribution(node3_contribution_v, node3_contribution_bv),
      node4_contribution(node4_contribution_v, node4_contribution_bv),
      n1c_n2c(n1c_n2c_v, n1c_n2c_bv),
      n1c_n2c_n3c(n1c_n2c_n3c_v, n1c_n2c_n3c_bv),
      /* Operations: */
      /* Calculations for first node. */
      node1_phi_expression(node1.rotation, node1.shell_director, node1_phi),
      node1_thickness_scale_expression(node1.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(0, r_ind, s_ind) * 0.5,
                                       node1_scaled_thickness),
      node1_scale_expression(node1_phi, node1_scaled_thickness, node1_contribution),
      /* Calculations for second node. */
      node2_phi_expression(node2.rotation, node2.shell_director, node2_phi),
      node2_thickness_scale_expression(node2.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(1, r_ind, s_ind) * 0.5,
                                       node2_scaled_thickness),
      node2_scale_expression(node2_phi, node2_scaled_thickness, node2_contribution),
      /* Calculations for third node. */
      node3_phi_expression(node3.rotation, node3.shell_director, node3_phi),
      node3_thickness_scale_expression(node3.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(2, r_ind, s_ind) * 0.5,
                                       node3_scaled_thickness),
      node3_scale_expression(node3_phi, node3_scaled_thickness, node3_contribution),
      /* Calculations for fourth node. */
      node4_phi_expression(node4.rotation, node4.shell_director, node4_phi),
      node4_thickness_scale_expression(node4.thickness,
                                       ShellElementMITC4<N, T>::Ni_rj_sk(3, r_ind, s_ind) * 0.5,
                                       node4_scaled_thickness),
      node4_scale_expression(node4_phi, node4_scaled_thickness, node4_contribution),
      /* Sum together the components. */
      sum_12_expression(1.0, node1_contribution, node2_contribution, n1c_n2c),
      sum_123_expression(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c),
      sum_1234_expression(1.0, n1c_n2c_n3c, node4_contribution, result) {};

  void reverse() {
    /* Sum component calls: */
    sum_1234_expression.reverse();
    sum_123_expression.reverse();
    sum_12_expression.reverse();

    /* Node expression calls: */
    node4_scale_expression.reverse();
    node4_thickness_scale_expression.reverse();
    node4_phi_expression.reverse();
    node3_scale_expression.reverse();
    node3_thickness_scale_expression.reverse();
    node3_phi_expression.reverse();
    node2_scale_expression.reverse();
    node2_thickness_scale_expression.reverse();
    node2_phi_expression.reverse();
    node1_scale_expression.reverse();
    node1_thickness_scale_expression.reverse();
    node1_phi_expression.reverse();
  };

  void hforward() {
    /* Node expression calls: */
    node1_phi_expression.hforward();
    node1_thickness_scale_expression.hforward();
    node1_scale_expression.hforward();
    node2_phi_expression.hforward();
    node2_thickness_scale_expression.hforward();
    node2_scale_expression.hforward();
    node3_phi_expression.hforward();
    node3_thickness_scale_expression.hforward();
    node3_scale_expression.hforward();
    node4_phi_expression.hforward();
    node4_thickness_scale_expression.hforward();
    node4_scale_expression.hforward();

    /* Sum component calls: */
    sum_12_expression.hforward();
    sum_123_expression.hforward();
    sum_1234_expression.hforward();
  };

  void hreverse() {
    /* Sum component calls: */
    sum_1234_expression.hreverse();
    sum_123_expression.hreverse();
    sum_12_expression.hreverse();

    /* Node expression calls: */
    node4_scale_expression.hhreverse();
    node4_thickness_scale_expression.hreverse();
    node4_phi_expression.hreverse();
    node3_scale_expression.hreverse();
    node3_thickness_scale_expression.hreverse();
    node3_phi_expression.hreverse();
    node2_scale_expression.hreverse();
    node2_thickness_scale_expression.hreverse();
    node2_phi_expression.hreverse();
    node1_scale_expression.hreverse();
    node1_thickness_scale_expression.hreverse();
    node1_phi_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      node1_phi_v, node1_phi_bv,
      node2_phi_v, node2_phi_bv,
      node3_phi_v, node3_phi_bv,
      node4_phi_v, node4_phi_bv,
      node1_contribution_v, node1_contribution_bv,
      node2_contribution_v, node2_contribution_bv,
      node3_contribution_v, node3_contribution_bv,
      node4_contribution_v, node4_contribution_bv,
      n1c_n2c_v, n1c_n2c_bv,
      n1c_n2c_n3c_v, n1c_n2c_n3c_bv;

  /* A2D Objects: */
  A2DVec<N, Vec<T, 3>>
      node1_phi,
      node2_phi,
      node3_phi,
      node4_phi;
  A2DVec<N, Vec<T, 3>>
      node1_contribution,
      node2_contribution,
      node3_contribution,
      node4_contribution;
  A2DScalar<N, T>
      node1_scaled_thickness,
      node2_scaled_thickness,
      node3_scaled_thickness,
      node4_scaled_thickness;
  A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

  /* Expressions: */

  A2DVec3CrossVecExpr<N, T> node1_phi_expression;
  A2DScalarScalarMultExpr<N, T> node1_thickness_scale_expression;
  A2DVec3A2DScaleExpr<N, T> node1_scale_expression;

  A2DVec3CrossVecExpr<N, T> node2_phi_expression;
  A2DScalarScalarMultExpr<N, T> node2_thickness_scale_expression;
  A2DVec3A2DScaleExpr<N, T> node2_scale_expression;

  A2DVec3CrossVecExpr<N, T> node3_phi_expression;
  A2DScalarScalarMultExpr<N, T> node3_thickness_scale_expression;
  A2DVec3A2DScaleExpr<N, T> node3_scale_expression;

  A2DVec3CrossVecExpr<N, T> node4_phi_expression;
  A2DScalarScalarMultExpr<N, T> node4_thickness_scale_expression;
  A2DVec3A2DScaleExpr<N, T> node4_scale_expression;

  A2DVec3A2DVecScalarAxpyExpr<N, T>
      sum_12_expression,
      sum_123_expression,
      sum_1234_expression;
};

template <int N, typename T>
class tying_shear_expr {
 public:
  /**
   * @brief Computes the e<sub>&alpha t</sub> covariant strain at a tying point.
   *
   * @param TODO
   * */
  tying_shear_expr(const int alpha, /* values are 0 (corresponding to r) or 1 (corresponding to s)*/
                   const int n_alpha_val_ind, /* values are 0 (corresponding to -1) or 1 (corresponding to 1)*/
//                   const ShellElementMITC4<N, T>& element,
                   ShellNodeMITC<N, T>& node1, ShellNodeMITC<N, T>& node2,
                   ShellNodeMITC<N, T>& node3, ShellNodeMITC<N, T>& node4,
                   A2DScalar<N, T>& e_alpha_t)
      : /* Initializations: */
      r_ind((4 * n_alpha_val_ind * alpha) + (2 * (1 - alpha))),
      s_ind((2 * alpha) + (4 * n_alpha_val_ind * (1 - alpha))),
      n_alpha_var_ind(3 * n_alpha_val_ind),
      g_alpha(g_alpha_v, g_alpha_bv),
      u_alpha(u_alpha_v, u_alpha_bv),
      gt(gt_v, gt_bv),
      ut(ut_v, ut_bv),
      /* Operations: */
      g_alpha_expression(alpha, n_alpha_var_ind, 0, node1, node2, node3, node4, g_alpha),
      gt_expression(r_ind, s_ind, node1, node2, node3, node4, gt),
      u_alpha_expression(alpha, n_alpha_var_ind, 0, node1, node2, node3, node4, u_alpha),
      ut_expression(r_ind, s_ind, node1, node2, node3, node4, ut),
      g_alpha_ut_expression(g_alpha, ut, g_alpha_ut),
      gt_u_alpha_expression(gt, u_alpha, gt_u_alpha),
      e_alpha_t_expression(0.5, g_alpha_ut, gt_u_alpha, e_alpha_t) {};

  void reverse() {
    e_alpha_t_expression.reverse();
    gt_u_alpha_expression.reverse();
    g_alpha_ut_expression.reverse();
    ut_expression.reverse();
    u_alpha_expression.reverse();
    gt_expression.reverse();
    g_alpha_expression.reverse();
  };

  void hforward() {
    g_alpha_expression.hforward();
    gt_expression.hforward();
    u_alpha_expression.hforward();
    ut_expression.hforward();
    g_alpha_ut_expression.hforward();
    gt_u_alpha_expression.hforward();
    e_alpha_t_expression.hforward();
  };

  void hreverse() {
    e_alpha_t_expression.reverse();
    gt_u_alpha_expression.reverse();
    g_alpha_ut_expression.reverse();
    ut_expression.reverse();
    u_alpha_expression.reverse();
    gt_expression.reverse();
    g_alpha_expression.reverse();
  };

 private:
  /* Instantiating objects: */
  int r_ind, s_ind, n_alpha_var_ind;
  Vec<T, 3>
      g_alpha_v, g_alpha_bv,
      u_alpha_v, u_alpha_bv,
      gt_v, gt_bv,
      ut_v, ut_bv;

  /* A2D Objects: */
  A2DScalar<N, T>
      g_alpha_ut, gt_u_alpha;
  A2DVec<N, Vec<T, 3>>
      g_alpha,
      u_alpha,
      gt,
      ut;

  /* Expressions: */
  g_alpha_expr<N, T> g_alpha_expression;
  g_t_expr<N, T> gt_expression;
  u_alpha_expr<N, T> u_alpha_expression;
  u_t_expr<N, T> ut_expression;
  A2DVec3DotA2DVecExpr<N, T> g_alpha_ut_expression;
  A2DVec3DotA2DVecExpr<N, T> gt_u_alpha_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_alpha_t_expression;
};

template <int N, typename T>
class contravariant_basis_expr {
 public:
  /**
   * @brief Constructs the contravariant basis vectors (g<sup>r</sup>, g<sup>s</sup>, and g<sup>t</sup>) from the
   * covariant basis vectors (g<sub>r</sub>, g<sub>s</sub>, and g<sub>t</sub>).
   *
   * @note The covariant basis vectors must be evaluated at the desired point of interest.
   *
   * @param gr: the g<sub>r</sub> covariant basis vector
   * @param gs: the g<sub>s</sub> covariant basis vector
   * @param gt: the g<sub>t</sub> covariant basis vector
   * @param Gr: the g<sup>r</sup> contravariant basis vector (output)
   * @param Gs: the g<sup>s</sup> contravariant basis vector (output)
   * @param Gt: the g<sup>t</sup> contravariant basis vector (output)
   * */
  contravariant_basis_expr(A2DVec<N, Vec<T, 3>>& gr, A2DVec<N, Vec<T, 3>>& gs, A2DVec<N, Vec<T, 3>>& gt,
                           A2DVec<N, Vec<T, 3>>& Gr, A2DVec<N, Vec<T, 3>>& Gs, A2DVec<N, Vec<T, 3>>& Gt)
      : /* Initializations: */
      gs_cross_gt(gs_cross_gt_v, gs_cross_gt_bv),
      gt_cross_gr(gt_cross_gr_v, gt_cross_gr_bv),
      gr_cross_gs(gr_cross_gs_v, gr_cross_gs_bv),
      /* Operations: */
      gs_cross_gt_expression(gs, gt, gs_cross_gt),
      gt_cross_gr_expression(gt, gr, gt_cross_gr),
      gr_cross_gs_expression(gr, gs, gr_cross_gs),

      gr_dot_gs_cross_gt_expression(gr, gs_cross_gt, gr_dot_gs_cross_gt),

      Gr_expression(gs_cross_gt, gr_dot_gs_cross_gt, Gr),
      Gs_expression(gt_cross_gr, gr_dot_gs_cross_gt, Gs),
      Gt_expression(gr_cross_gs, gr_dot_gs_cross_gt, Gt) {};

  void reverse() {
    Gt_expression.reverse();
    Gs_expression.reverse();
    Gr_expression.reverse();

    gr_dot_gs_cross_gt_expression.reverse();

    gr_cross_gs_expression.reverse();
    gt_cross_gr_expression.reverse();
    gs_cross_gt_expression.reverse();
  };

  void hforward() {
    gs_cross_gt_expression.reverse();
    gt_cross_gr_expression.reverse();
    gr_cross_gs_expression.reverse();

    gr_dot_gs_cross_gt_expression.reverse();

    Gr_expression.reverse();
    Gs_expression.reverse();
    Gt_expression.reverse();
  };

  void hreverse() {
    Gt_expression.hreverse();
    Gs_expression.hreverse();
    Gr_expression.hreverse();

    gr_dot_gs_cross_gt_expression.hreverse();

    gr_cross_gs_expression.hreverse();
    gt_cross_gr_expression.hreverse();
    gs_cross_gt_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      gs_cross_gt_v, gs_cross_gt_bv,
      gt_cross_gr_v, gt_cross_gr_bv,
      gr_cross_gs_v, gr_cross_gs_bv;

  /* A2D Objects: */
  A2DScalar<N, T>
      gr_dot_gs_cross_gt;
  A2DVec<N, Vec<T, 3>>
      gs_cross_gt,
      gt_cross_gr,
      gr_cross_gs;

  /* Expressions: */
  A2DVec3CrossA2DVecExpr<N, T>
      gs_cross_gt_expression,
      gt_cross_gr_expression,
      gr_cross_gs_expression;
  A2DVec3DotA2DVecExpr<N, T>
      gr_dot_gs_cross_gt_expression;
  A2DVec3A2DScaleDivExpr<N, T>
      Gr_expression,
      Gs_expression,
      Gt_expression;
};

template <int N, typename T>
class cartesian_local_basis_expr {
 public:
  /**
   * @brief Constructs the cartesian local basis vectors (e<sub>1</sub>, e<sub>2</sub>, and e<sub>3</sub>) from two of
   * the covariant basis vectors (g<sub>s</sub> and g<sub>t</sub>).
   *
   * @note The covariant basis vectors, g<sub>s</sub> and g<sub>t</sub>, must be evaluated at the desired point of
   * interest.
   *
   * @param gs: the g<sub>s</sub> covariant basis vector
   * @param gt: the g<sub>t</sub> covariant basis vector
   * @param e1: the e<sub>1</sub> cartesian local basis vector (output)
   * @param e2: the e<sub>2</sub> cartesian local basis vector (output)
   * @param e3: the e<sub>3</sub> cartesian local basis vector (output)
   * */
  cartesian_local_basis_expr(A2DVec<N, Vec<T, 3>>& gs, A2DVec<N, Vec<T, 3>>& gt,
                             A2DVec<N, Vec<T, 3>>& e1, A2DVec<N, Vec<T, 3>>& e2, A2DVec<N, Vec<T, 3>>& e3)
      : /* Initializations: */
      gs_cross_e3(gs_cross_e3_v, gs_cross_e3_bv),
      /* Operations: */
      e3_expression(gt, e3),
      gs_cross_e3_expression(gs, e3, gs_cross_e3),
      e1_expression(gs_cross_e3, e1),
      e2_expression(e3, e1, e2) {};

  void reverse() {
    e2_expression.reverse();
    e1_expression.reverse();
    gs_cross_e3_expression.reverse();
    e3_expression.reverse();
  };

  void hforward() {
    e3_expression.hforward();
    gs_cross_e3_expression.hforward();
    e1_expression.hforward();
    e2_expression.hforward();
  };

  void hreverse() {
    e2_expression.hreverse();
    e1_expression.hreverse();
    gs_cross_e3_expression.hreverse();
    e3_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      gs_cross_e3_v, gs_cross_e3_bv;

  /* A2D Objects: */
  A2DVec<N, Vec<T, 3>> gs_cross_e3;

  /* Expressions: */
  A2DVec3NormalizeExpr<N, T> e3_expression;
  A2DVec3CrossA2DVecExpr<N, T> gs_cross_e3_expression;
  A2DVec3NormalizeExpr<N, T> e1_expression;
  A2DVec3CrossA2DVecExpr<N, T> e2_expression;
};

template <int N, typename T>
class in_plane_strain_components_expr {
 public:
  /**
   * @brief Calculates the in-plane stain components (e<sub>rr</sub>, e<sub>ss</sub>, and e<sub>rs</sub>) based on the
   * values of the covariant basis vectors g<sub>r</sub> and g<sub>s</sub>, and the displacement derivative vectors
   * du/dr and du/ds.
   *
   * @note The covariant basis vectors, g<sub>r</sub> and g<sub>s</sub>, and displacement derivative vectors, du/dr and
   * du/ds, must be evaluated at the desired point of interest (i.e. a quadrature point).
   *
   * @param gr: the g<sub>r</sub> covariant basis vector
   * @param gs: the g<sub>s</sub> covariant basis vector
   * @param ur: the derivative of the displacement vector with respect to r (i.e. du/dr)
   * @param us: the derivative of the displacement vector with respect to s (i.e. du/ds)
   * */
  in_plane_strain_components_expr(A2DVec<N, Vec<T, 3>>& gr,
                                  A2DVec<N, Vec<T, 3>>& gs,
                                  A2DVec<N, Vec<T, 3>>& ur,
                                  A2DVec<N, Vec<T, 3>>& us,
                                  A2DScalar<N, T>& e_rr,
                                  A2DScalar<N, T>& e_ss,
                                  A2DScalar<N, T>& e_rs)
      : e_rr_expression(gr, ur, e_rr),
        e_ss_expression(gs, us, e_ss),
        gr_us_expression(gr, us, gr_us),
        gs_ur_expression(gs, ur, gs_ur),
        e_rs_expression(0.5, gr_us, gs_ur, e_rs) {};

  void reverse() {
    e_rs_expression.reverse();
    gs_ur_expression.reverse();
    gr_us_expression.reverse();
    e_ss_expression.reverse();
    e_rr_expression.reverse();
  };

  void hforward() {
    e_rr_expression.hforward();
    e_ss_expression.hforward();
    gr_us_expression.hforward();
    gs_ur_expression.hforward();
    e_rs_expression.hforward();
  };

  void hreverse() {
    e_rs_expression.hreverse();
    gs_ur_expression.hreverse();
    gr_us_expression.hreverse();
    e_ss_expression.hreverse();
    e_rr_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  /* None necessary */

  /* A2D Objects: */
  A2DScalar<N, T>
      gr_us, gs_ur;

  /* Expressions: */
  A2DVec3DotA2DVecExpr<N, T>
      e_rr_expression, e_ss_expression,
      gr_us_expression, gs_ur_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T>
      e_rs_expression;
};

template <int N, typename T>
class local_strain_expr {
 public:
  /**
   * @brief Computes the i,j<sup>th</sup> local strain (&epsilon <sub>ij</sub>) from the covariant strains and dot
   * products of the contravariant basis vectors and cartesian basis vectors.
   *
   * @param TODO
   * */
  local_strain_expr(A2DScalar<N, T>& Gr_ei, A2DScalar<N, T>& Gs_ei, A2DScalar<N, T>& Gt_ei,
                    A2DScalar<N, T>& Gr_ej, A2DScalar<N, T>& Gs_ej, A2DScalar<N, T>& Gt_ej,
                    A2DScalar<N, T>& e_rr, A2DScalar<N, T>& e_ss, A2DScalar<N, T>& e_rs,
                    A2DScalar<N, T>& e_rt, A2DScalar<N, T>& e_st,
                    A2DScalar<N, T>& e_ij)
      : Gr_ej_Gs_ei_expression(Gr_ej, Gs_ei, Gr_ej_Gs_ei),
        Gr_ei_Gs_ej_expression(Gr_ei, Gs_ej, Gr_ei_Gs_ej),

        Gr_ej_Gt_ei_expression(Gr_ej, Gt_ei, Gr_ej_Gt_ei),
        Gr_ei_Gt_ej_expression(Gr_ei, Gt_ej, Gr_ei_Gt_ej),

        Gs_ej_Gt_ei_expression(Gs_ej, Gt_ei, Gs_ej_Gt_ei),
        Gs_ei_Gt_ej_expression(Gs_ei, Gt_ej, Gs_ei_Gt_ej),

        e_rr_multiplier_expression(Gr_ei, Gr_ej, e_rr_multiplier),
        e_ss_multiplier_expression(Gs_ei, Gs_ej, e_ss_multiplier),
        e_rs_multiplier_expression(Gr_ej, Gs_ei, Gr_ei, Gs_ej, e_rs_multiplier),
        e_rt_multiplier_expression(Gr_ej, Gt_ei, Gr_ei, Gt_ej, e_rt_multiplier),
        e_st_multiplier_expression(Gs_ej, Gt_ei, Gs_ei, Gt_ej, e_st_multiplier),

        e_rr_expression(e_rr_multiplier, e_rr, scaled_e_rr),
        e_ss_expression(e_ss_multiplier, e_ss, scaled_e_ss),
        e_rs_expression(e_rs_multiplier, e_rs, scaled_e_rs),
        e_rt_expression(e_rt_multiplier, e_rt, scaled_e_rt),
        e_st_expression(e_st_multiplier, e_st, scaled_e_st),

        sum_12_expression(1, scaled_e_rr, scaled_e_ss, sum_12),
        sum_123_expression(1, sum_12, scaled_e_rs, sum_123),
        sum_1234_expression(1, sum_123, scaled_e_rt, sum_1234),
        e_ij_expression(1, sum_1234, scaled_e_st, e_ij) {};

  void reverse() {
    e_ij_expression.reverse();
    sum_1234_expression.reverse();
    sum_123_expression.reverse();
    sum_12_expression.reverse();

    e_st_expression.reverse();
    e_rt_expression.reverse();
    e_rs_expression.reverse();
    e_ss_expression.reverse();
    e_rr_expression.reverse();

    e_st_multiplier_expression.reverse();
    e_rt_multiplier_expression.reverse();
    e_rs_multiplier_expression.reverse();
    e_ss_multiplier_expression.reverse();
    e_rr_multiplier_expression.reverse();

    Gs_ei_Gt_ej_expression.reverse();
    Gs_ej_Gt_ei_expression.reverse();

    Gr_ei_Gt_ej_expression.reverse();
    Gr_ej_Gt_ei_expression.reverse();

    Gr_ei_Gs_ej_expression.reverse();
    Gr_ej_Gs_ei_expression.reverse();
  };

  void hforward() {
    Gr_ej_Gs_ei_expression.hforward();
    Gr_ei_Gs_ej_expression.hforward();

    Gr_ej_Gt_ei_expression.hforward();
    Gr_ei_Gt_ej_expression.hforward();

    Gs_ej_Gt_ei_expression.hforward();
    Gs_ei_Gt_ej_expression.hforward();

    e_rr_multiplier_expression.hforward();
    e_ss_multiplier_expression.hforward();
    e_rs_multiplier_expression.hforward();
    e_rt_multiplier_expression.hforward();
    e_st_multiplier_expression.hforward();

    e_rr_expression.hforward();
    e_ss_expression.hforward();
    e_rs_expression.hforward();
    e_rt_expression.hforward();
    e_st_expression.hforward();

    sum_12_expression.hforward();
    sum_123_expression.hforward();
    sum_1234_expression.hforward();
    e_ij_expression.hforward();
  };

  void hreverse() {
    e_ij_expression.hreverse();
    sum_1234_expression.hreverse();
    sum_123_expression.hreverse();
    sum_12_expression.hreverse();

    e_st_expression.hreverse();
    e_rt_expression.hreverse();
    e_rs_expression.hreverse();
    e_ss_expression.hreverse();
    e_rr_expression.hreverse();

    e_st_multiplier_expression.hreverse();
    e_rt_multiplier_expression.hreverse();
    e_rs_multiplier_expression.hreverse();
    e_ss_multiplier_expression.hreverse();
    e_rr_multiplier_expression.hreverse();

    Gs_ei_Gt_ej_expression.hreverse();
    Gs_ej_Gt_ei_expression.hreverse();

    Gr_ei_Gt_ej_expression.hreverse();
    Gr_ej_Gt_ei_expression.hreverse();

    Gr_ei_Gs_ej_expression.hreverse();
    Gr_ej_Gs_ei_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  /* None necessary */

  /* A2D Objects: */
  A2DScalar<N, T>
      e_rr_multiplier, e_ss_multiplier, e_rs_multiplier, e_rt_multiplier, e_st_multiplier,
      scaled_e_rr, scaled_e_ss, scaled_e_rs, scaled_e_rt, scaled_e_st;
  A2DScalar<N, T>
      Gr_ej_Gs_ei, Gr_ei_Gs_ej,
      Gr_ej_Gt_ei, Gr_ei_Gt_ej,
      Gs_ej_Gt_ei, Gs_ei_Gt_ej;
  A2DScalar<N, T>
      sum_12, sum_123, sum_1234;

  /* Expressions: */
  A2DScalarA2DScalarMultExpr<N, T>
      Gr_ej_Gs_ei_expression,
      Gr_ei_Gs_ej_expression,
      Gr_ej_Gt_ei_expression,
      Gr_ei_Gt_ej_expression,
      Gs_ej_Gt_ei_expression,
      Gs_ei_Gt_ej_expression;
  A2DScalarA2DScalarMultExpr<N, T>
      e_rr_multiplier_expression,
      e_ss_multiplier_expression;
  A2DScalarA2DScalarA2DScalarA2DScalarAxpbyExpr<N, T>
      e_rs_multiplier_expression,
      e_rt_multiplier_expression,
      e_st_multiplier_expression;
  A2DScalarA2DScalarMultExpr<N, T>
      e_rr_expression,
      e_ss_expression,
      e_rs_expression,
      e_rt_expression,
      e_st_expression;
  ScalarA2DScalarA2DScalarAxpayExpr<N, T>
      sum_12_expression,
      sum_123_expression,
      sum_1234_expression,
      e_ij_expression;
};

template <int N, typename T>
class local_strains_expr {
 public:
  /**
   * @brief Computes the local strains from the contravariant basis vectors, cartesian basis vectors, and covariant
   * strains.
   *
   * @param: TODO
   * */
  local_strains_expr(A2DVec<N, Vec<T, 3>>& Gr, A2DVec<N, Vec<T, 3>>& Gs, A2DVec<N, Vec<T, 3>>& Gt,
                     A2DVec<N, Vec<T, 3>>& e1, A2DVec<N, Vec<T, 3>>& e2, A2DVec<N, Vec<T, 3>>& e3,
                     A2DScalar<N, T>& e_rr, A2DScalar<N, T>& e_ss, A2DScalar<N, T>& e_rs,
                     A2DScalar<N, T>& e_rt, A2DScalar<N, T>& e_st,
                     A2DScalar<N, T>& e_11, A2DScalar<N, T>& e_22, A2DScalar<N, T>& e_12,
                     A2DScalar<N, T>& e_13, A2DScalar<N, T>& e_23)
      : Gr_e1_expression(Gr, e1, Gr_e1),
        Gr_e2_expression(Gr, e2, Gr_e2),
        Gr_e3_expression(Gr, e3, Gr_e3),
        Gs_e1_expression(Gs, e1, Gs_e1),
        Gs_e2_expression(Gs, e2, Gs_e2),
        Gs_e3_expression(Gs, e3, Gs_e3),
        Gt_e1_expression(Gt, e1, Gt_e1),
        Gt_e2_expression(Gt, e2, Gt_e2),
        Gt_e3_expression(Gt, e3, Gt_e3),

        e_11_expression(Gr_e1, Gs_e1, Gt_e1, Gr_e1, Gs_e1, Gt_e1, e_rr, e_ss, e_rs, e_rt, e_st, e_11),
        e_22_expression(Gr_e2, Gs_e2, Gt_e2, Gr_e2, Gs_e2, Gt_e2, e_rr, e_ss, e_rs, e_rt, e_st, e_22),
        e_12_expression(Gr_e1, Gs_e1, Gt_e1, Gr_e2, Gs_e2, Gt_e2, e_rr, e_ss, e_rs, e_rt, e_st, e_12),
        e_13_expression(Gr_e1, Gs_e1, Gt_e1, Gr_e3, Gs_e3, Gt_e3, e_rr, e_ss, e_rs, e_rt, e_st, e_13),
        e_23_expression(Gr_e2, Gs_e2, Gt_e2, Gr_e3, Gs_e3, Gt_e3, e_rr, e_ss, e_rs, e_rt, e_st, e_23) {};

  void reverse() {
    e_23_expression.reverse();
    e_13_expression.reverse();
    e_12_expression.reverse();
    e_22_expression.reverse();
    e_11_expression.reverse();

    Gt_e3_expression.reverse();
    Gt_e2_expression.reverse();
    Gt_e1_expression.reverse();
    Gs_e3_expression.reverse();
    Gs_e2_expression.reverse();
    Gs_e1_expression.reverse();
    Gr_e3_expression.reverse();
    Gr_e2_expression.reverse();
    Gr_e1_expression.reverse();
  };

  void hforward() {
    Gr_e1_expression.hforward();
    Gr_e2_expression.hforward();
    Gr_e3_expression.hforward();
    Gs_e1_expression.hforward();
    Gs_e2_expression.hforward();
    Gs_e3_expression.hforward();
    Gt_e1_expression.hforward();
    Gt_e2_expression.hforward();
    Gt_e3_expression.hforward();

    e_11_expression.hforward();
    e_22_expression.hforward();
    e_12_expression.hforward();
    e_13_expression.hforward();
    e_23_expression.hforward();
  };

  void hreverse() {
    e_23_expression.hreverse();
    e_13_expression.hreverse();
    e_12_expression.hreverse();
    e_22_expression.hreverse();
    e_11_expression.hreverse();

    Gt_e3_expression.hreverse();
    Gt_e2_expression.hreverse();
    Gt_e1_expression.hreverse();
    Gs_e3_expression.hreverse();
    Gs_e2_expression.hreverse();
    Gs_e1_expression.hreverse();
    Gr_e3_expression.hreverse();
    Gr_e2_expression.hreverse();
    Gr_e1_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  /* None necessary */

  /* A2D Objects: */
  A2DScalar<N, T>
      Gr_e1, Gr_e2, Gr_e3,
      Gs_e1, Gs_e2, Gs_e3,
      Gt_e1, Gt_e2, Gt_e3;

  /* Expressions: */
  A2DVec3DotA2DVecExpr<N, T>
      Gr_e1_expression,
      Gr_e2_expression,
      Gr_e3_expression,
      Gs_e1_expression,
      Gs_e2_expression,
      Gs_e3_expression,
      Gt_e1_expression,
      Gt_e2_expression,
      Gt_e3_expression;
  local_strain_expr<N, T>
      e_11_expression,
      e_22_expression,
      e_12_expression,
      e_13_expression,
      e_23_expression;
};

template <int N, typename T>
class strain_energy_expr {
 public:
  /**
   * @brief Compute the strain energy of the element for some point within the element.
   *
   * @note All inputs (gr, gs, gt, ur, us, ut, e_rt, and e_st) must be evaluated at the point of interest in order for
   * the calculations to be correct.
   *
   * @param gr:         the g<sub>r</sub> covariant basis vector.
   * @param gs:         the g<sub>s</sub> covariant basis vector.
   * @param gt:         the g<sub>t</sub> covariant basis vector.
   * @param ur:         the derivative of the displacement vector with respect to r (i.e. du/dr).
   * @param us:         the derivative of the displacement vector with respect to s (i.e. du/ds).
   * @param ut:         the derivative of the displacement vector with respect to t (i.e. du/dt).
   * @param e_rt:       the covariant shear strain component between the r and t directions (e<sub>rt</sub>).
   * @param e_st:       the covariant shear strain component between the s and t directions (e<sub>st</sub>).
   * @param element:    the MITC4 element object for which the strain energy is being computed.
   * @param energy:     the strain energy as an A2DScalar object to store the values (output).
   * */
  strain_energy_expr(A2DVec<N, Vec<T, 3>>& gr, A2DVec<N, Vec<T, 3>>& gs, A2DVec<N, Vec<T, 3>>& gt,
                     A2DVec<N, Vec<T, 3>>& ur, A2DVec<N, Vec<T, 3>>& us, A2DVec<N, Vec<T, 3>>& ut,
                     A2DScalar<N, T>& e_rt, A2DScalar<N, T>& e_st,
                     const LinearIsotropicMaterial<T>& material,
                     A2DScalar<N, T>& energy)
      : /* Initializations: */
      Gr(Gr_v, Gr_bv), Gs(Gs_v, Gs_bv), Gt(Gt_v, Gt_bv),
      e1(e1_v, e1_bv), e2(e2_v, e2_bv), e3(e3_v, e3_bv),
      local_strains_vec(lsv_v, lsv_bv),
      /* In plane strain calculations: */
      strain_expression(gr, gs, ur, us, e_rr, e_ss, e_rs),
      /* Contravariant basis: */
      contravariant_basis_expression(gr, gs, gt, Gr, Gs, Gt),
      /* Cartesian local basis: */
      local_basis_expression(gs, gt, e1, e2, e3),
      /* Calculate local strains */
      local_strains_expression(Gr, Gs, Gt, e1, e2, e3, e_rr, e_ss, e_rs, e_rt, e_st, e_11, e_22, e_12, e_13, e_23),
      /* Assemble local strain vector */
      local_strains_vec_expression(local_strains_vec, e_11, e_22, e_12, e_13, e_23),
      /* Calculate strain energy */
      strain_energy_expression(material.D, local_strains_vec, local_strains_vec, energy) {};

  void reverse() {
    strain_energy_expression.reverse();
    local_strains_vec_expression.reverse();
    local_strains_expression.reverse();
    local_basis_expression.reverse();
    contravariant_basis_expression.reverse();
    strain_expression.reverse();
  };

  void hforward() {
    strain_expression.hforward();
    contravariant_basis_expression.hforward();
    local_basis_expression.hforward();
    local_strains_expression.hforward();
    local_strains_vec_expression.hforward();
    strain_energy_expression.hforward();
  };

  void hreverse() {
    strain_energy_expression.hreverse();
    local_strains_vec_expression.hreverse();
    local_strains_expression.hreverse();
    local_basis_expression.hreverse();
    contravariant_basis_expression.hreverse();
    strain_expression.hreverse();
  };

 private:
  /* Instantiating objects: */
  Vec<T, 3>
      Gr_v, Gr_bv,
      Gs_v, Gs_bv,
      Gt_v, Gt_bv,
      e1_v, e1_bv,
      e2_v, e2_bv,
      e3_v, e3_bv;
  Vec<T, 5>
      lsv_v, lsv_bv;

  /* A2D Objects: */
  A2DScalar<N, T>
      e_rr, e_ss, e_rs;
  A2DVec<N, Vec<T, 3>>
      Gr, Gs, Gt,
      e1, e2, e3;
  A2DScalar<N, T>
      e_11, e_22, e_12, e_13, e_23;
  A2DVec<N, Vec<T, 5>>
      local_strains_vec;

  /* Expressions: */
  in_plane_strain_components_expr<N, T> strain_expression;
  contravariant_basis_expr<N, T> contravariant_basis_expression;
  cartesian_local_basis_expr<N, T> local_basis_expression;
  local_strains_expr<N, T> local_strains_expression;
  A2DScalar5VecAssemblyExpr<N, T> local_strains_vec_expression;
  MatA2DVecA2DVecInnerProductExpr<N, T, 5, 5> strain_energy_expression;
};

}  // namespace A2D

namespace A2D::TEST {

int main() {
  Vec<double, 3>
      d1, db1, d2, db2, d3, db3, d4, db4,
      r1, rb1, r2, rb2, r3, rb3, r4, rb4;
  ShellNodeMITC<24, double>
      n1{d1, db1, r1, rb1},
      n2{d2, db2, r2, rb2},
      n3{d3, db3, r3, rb3},
      n4{d4, db4, r4, rb4};
  ShellElementMITC4<24, double> x(n1, n2, n3, n4);
//  x.generate_stiffness_matrix();

//  x.test;
//  x.Ni_rj_sk(0, 0, 0);

  return 0;
}

}

#endif //A2D_SHELL_DEVELOPMENT_H_
