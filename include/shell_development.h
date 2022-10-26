//
// Created by James on 10/23/2022.
//

#ifndef A2D_SHELL_DEVELOPMENT_H_
#define A2D_SHELL_DEVELOPMENT_H_

#include "a2dtypes.h"
#include "a2dvecops3d.h"
#include "a2dmatops3d.h"

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
               0, 0, 0, temp3, 0,
               0, 0, 0, 0, temp3}, D(D_init) {
  };

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
      : node1(node1), node2(node2), node3(node3), node4(node4) {
  };

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

  /*class InternalEnergy {
    InternalEnergy() {

    };

    void forward();
    void reverse();
    void hforward();
    void hreverse();

    T value;
  };

  InternalEnergy internal_energy();*/

  class g_alpha {
    g_alpha(const int alpha,  /* alpha can be 0 (corresponding to r) or 1 (corresponding to s) */
            const T t,
            const ShellElementMITC4<N, T>& element,
            A2DVec<N, Vec<T, 3>>& result) {
      A2DVec<N, Vec<T, 3>>
          node1_contribution, node1_contribution_unscaled,
          node2_contribution, node2_contribution_unscaled,
          node3_contribution, node3_contribution_unscaled,
          node4_contribution, node4_contribution_unscaled;

      /* Calculations for first node. */

    }

    void reverse();
    void hforward();
    void hreverse();
  };

  void generate_energy() {
    /* Generate the internal energy of the shell element, that will then be used to calculate the global stiffness
     * matrix
     * */

    A2DVec<N, Vec<T, 3>>
    /** gr vector evaluated at the various quadrature points */
    gr_r0s0t0, gr_r0s0t1, gr_r0s1t0, gr_r0s1t1, gr_r1s0t0, gr_r1s0t1, gr_r1s1t0, gr_r1s1t1,
    /** gs vector evaluated at the various quadrature points */
    gs_r0s0t0, gs_r0s0t1, gs_r0s1t0, gs_r0s1t1, gs_r1s0t0, gs_r1s0t1, gs_r1s1t0, gs_r1s1t1,
    /** gt vector evaluated at the tying points (s={1, -1} with r=t=0, r={1, -1} with s=t=0)*/
    gt_r0_sp1_t0, gt_r0_sn1_t0, gt_rp1_s0_t0, gt_rn1_s0_t0,
    /** derivatives of u with respect to r evaluated at the various quadrature points */
    ur_r0s0t0, ur_r0s0t1, ur_r0s1t0, ur_r0s1t1, ur_r1s0t0, ur_r1s0t1, ur_r1s1t0, ur_r1s1t1,
    /** derivatives of u with respect to s evaluated at the various quadrature points */
    us_r0s0t0, us_r0s0t1, us_r0s1t0, us_r0s1t1, us_r1s0t0, us_r1s0t1, us_r1s1t0, us_r1s1t1,
    /** derivative of u with respect to t evaluated at the tying points (s={1, -1} with r=t=0, r={1, -1} with s=t=0)*/
    ut_r0_sp1_t0, ut_r0_sn1_t0, ut_rp1_s0_t0, ut_rn1_s0_t0;

    /** Cartesian local basis: */
    A2DVec<N, Vec<T, 3>>
        e3_r0s0, e3_r0s1, e3_r1s0, e3_r1s1,
        e1_r0s0t0, e1_r0s0t1, e1_r0s1t0, e1_r0s1t1, e1_r1s0t0, e1_r1s0t1, e1_r1s1t0, e1_r1s1t1,
        e2_r0s0t0, e2_r0s0t1, e2_r0s1t0, e2_r0s1t1, e2_r1s0t0, e2_r1s0t1, e2_r1s1t0, e2_r1s1t1;

    /** Contravariant basis: */
    A2DVec<N, Vec<T, 3>>
        Gr_r0s0t0, Gr_r0s0t1, Gr_r0s1t0, Gr_r0s1t1, Gr_r1s0t0, Gr_r1s0t1, Gr_r1s1t0, Gr_r1s1t1,
        Gs_r0s0t0, Gs_r0s0t1, Gs_r0s1t0, Gs_r0s1t1, Gs_r1s0t0, Gs_r1s0t1, Gs_r1s1t0, Gs_r1s1t1,
        Gt_r0s0t0, Gt_r0s0t1, Gt_r0s1t0, Gt_r0s1t1, Gt_r1s0t0, Gt_r1s0t1, Gt_r1s1t0, Gt_r1s1t1;

    A2DScalar<N, T>
    /** Coefficients for the tying scheme. */
    e_rt_A, e_rt_B, e_st_C, e_st_D,
    /** Covariant transverse shear strains evaluated (using the tying points) at the various quadrature points. NOTE: the
     * values are the same for multiple quadrature points because we're assuming constant covariant transverse shear
     * strain conditions along the edges.*/
    e_rt_rAs0tA, e_rt_rAs1tA, e_st_r0sAtA, e_st_r1sAtA,
    /** Covariant in-plane strain components evaluated at the various quadrature points*/
    e_rr_r0s0t0, e_rr_r0s0t1, e_rr_r0s1t0, e_rr_r0s1t1, e_rr_r1s0t0, e_rr_r1s0t1, e_rr_r1s1t0, e_rr_r1s1t1,
        e_rs_r0s0t0, e_rs_r0s0t1, e_rs_r0s1t0, e_rs_r0s1t1, e_rs_r1s0t0, e_rs_r1s0t1, e_rs_r1s1t0, e_rs_r1s1t1,
        e_ss_r0s0t0, e_ss_r0s0t1, e_ss_r0s1t0, e_ss_r0s1t1, e_ss_r1s0t0, e_ss_r1s0t1, e_ss_r1s1t0, e_ss_r1s1t1;
    /*TODO: move all these to private attributes, and just reset their values whenever this method is called.*/

    // TODO: write code for tying points

    /* Write code for e_rt_r... and e_st_r... values */
    /* Evaluate at s = quad_val_0; */
    auto expression_f_1 = ScalarAxpby((1 + quad_0) * 0.5, e_rt_A, (1 - quad_0) * 0.5, e_rt_B, e_rt_rAs0tA);
    // TODO: make the ScalarAxpby operation (z = a * x + b * y)
    /* Evaluate at r = quad_val_0; */
    auto expression_f_2 = ScalarAxpby((1 + quad_0) * 0.5, e_st_C, (1 - quad_0) * 0.5, e_st_D, e_st_r0sAtA);
    /* Evaluate at s = quad_val_1; */
    auto expression_f_3 = ScalarAxpby((1 + quad_1) * 0.5, e_rt_A, (1 - quad_1) * 0.5, e_rt_B, e_rt_rAs1tA);
    /* Evaluate at r = quad_val_1; */
    auto expression_f_4 = ScalarAxpby((1 + quad_1) * 0.5, e_st_C, (1 - quad_1) * 0.5, e_st_D, e_st_r1sAtA);

    /* write test code for a single quadrature point: r0s0t0 <=> s=r=t=-1/sqrt(3) */
    double r(-1 / sqrt(3)), s(-1 / sqrt(3)), t(-1 / sqrt(3));
    gr_r0s0t0;
    gs_r0s0t0;
    ur_r0s0t0;
    us_r0s0t0;
    A2DVec3DotA2DVecExpr<N, T> expression5 = Vec3Dot(gr_r0s0t0, ur_r0s0t0, e_rr_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> expression6 = Vec3Dot(gs_r0s0t0, us_r0s0t0, e_ss_r0s0t0);
    A2DScalar<N, T> gr_us_r0s0t0, gs_ur_r0s0t0;  // TODO: move declaration
    A2DVec3DotA2DVecExpr<N, T> expression7 = Vec3Dot(gr_r0s0t0, us_r0s0t0, gr_us_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> expression8 = Vec3Dot(gs_r0s0t0, ur_r0s0t0, gs_ur_r0s0t0);
    auto expression9 = ScalarAxpay(0.5, gr_us_r0s0t0, gs_ur_r0s0t0, e_rs_r0s0t0);  // TODO: make this operation
  };

  ShellNodeMITC<N, T>& node1;  /**< The top left node. */
  ShellNodeMITC<N, T>& node2;  /**< The top right node.*/
  ShellNodeMITC<N, T>& node3;  /**< The bottom left node. */
  ShellNodeMITC<N, T>& node4;  /**< The bottom right node.*/

  const LinearIsotropicMaterial<T> material{5.2, 0.5};  /**< I'm using a linear isotropic material assumption here */

 private:
  // TODO

  // replace -1/sqrt(3) and 1/sqrt(3) with quad_0 and quad_1 respectively, do the reverse in comments?
  /** quadrature point values */
  constexpr static const T quad_0{-0.5773502691896257645091488}, quad_1{0.5773502691896257645091488};

  /** shape functions evaluated at various quadrature points */
  /*constexpr static const T  // TODO: remove
      N1_r0s0{(1 - quad_0) * (1 - quad_0) * 0.25},
      N1_r0s1{(1 - quad_0) * (1 - quad_1) * 0.25},
      N1_r1s0{(1 - quad_1) * (1 - quad_0) * 0.25},
      N1_r1s1{(1 - quad_1) * (1 - quad_1) * 0.25},
      N2_r0s0{(1 + quad_0) * (1 - quad_0) * 0.25},
      N2_r0s1{(1 + quad_0) * (1 - quad_1) * 0.25},
      N2_r1s0{(1 + quad_1) * (1 - quad_0) * 0.25},
      N2_r1s1{(1 + quad_1) * (1 - quad_1) * 0.25},
      N3_r0s0{(1 + quad_0) * (1 + quad_0) * 0.25},
      N3_r0s1{(1 + quad_0) * (1 + quad_1) * 0.25},
      N3_r1s0{(1 + quad_1) * (1 + quad_0) * 0.25},
      N3_r1s1{(1 + quad_1) * (1 + quad_1) * 0.25},
      N4_r0s0{(1 - quad_0) * (1 + quad_0) * 0.25},
      N4_r0s1{(1 - quad_0) * (1 + quad_1) * 0.25},
      N4_r1s0{(1 - quad_1) * (1 + quad_0) * 0.25},
      N4_r1s1{(1 - quad_1) * (1 + quad_1) * 0.25};*/
  constexpr static const T Ni_rj_sk_data[16]{
      (1 - quad_0) * (1 - quad_0) * 0.25,  /**< N1 evaluated at r=quad_0, s=quad_0 */
      (1 - quad_0) * (1 - quad_1) * 0.25,  /**< N1 evaluated at r=quad_0, s=quad_1 */
      (1 - quad_1) * (1 - quad_0) * 0.25,  /**< N1 evaluated at r=quad_1, s=quad_0 */
      (1 - quad_1) * (1 - quad_1) * 0.25,  /**< N1 evaluated at r=quad_1, s=quad_1 */
      (1 + quad_0) * (1 - quad_0) * 0.25,  /**< N2 evaluated at r=quad_0, s=quad_0 */
      (1 + quad_0) * (1 - quad_1) * 0.25,  /**< N2 evaluated at r=quad_0, s=quad_1 */
      (1 + quad_1) * (1 - quad_0) * 0.25,  /**< N2 evaluated at r=quad_1, s=quad_0 */
      (1 + quad_1) * (1 - quad_1) * 0.25,  /**< N2 evaluated at r=quad_1, s=quad_1 */
      (1 + quad_0) * (1 + quad_0) * 0.25,  /**< N3 evaluated at r=quad_0, s=quad_0 */
      (1 + quad_0) * (1 + quad_1) * 0.25,  /**< N3 evaluated at r=quad_0, s=quad_1 */
      (1 + quad_1) * (1 + quad_0) * 0.25,  /**< N3 evaluated at r=quad_1, s=quad_0 */
      (1 + quad_1) * (1 + quad_1) * 0.25,  /**< N3 evaluated at r=quad_1, s=quad_1 */
      (1 - quad_0) * (1 + quad_0) * 0.25,  /**< N4 evaluated at r=quad_0, s=quad_0 */
      (1 - quad_0) * (1 + quad_1) * 0.25,  /**< N4 evaluated at r=quad_0, s=quad_1 */
      (1 - quad_1) * (1 + quad_0) * 0.25,  /**< N4 evaluated at r=quad_1, s=quad_0 */
      (1 - quad_1) * (1 + quad_1) * 0.25   /**< N4 evaluated at r=quad_1, s=quad_1 */
  };
  /** getter for shape functions evaluated at the various quadrature points*/  // TODO: doc
  constexpr static T Ni_rj_sk(const int i, const int j, const int k) {
    return Ni_rj_sk_data[i * 4 + j * 2 + k];
  }

  /*/// shape function derivative with respect to r, evaluated at various quadrature points // TODO: remove
  constexpr static const T
      drN1_s0{-0.25 * (1 - quad_0)},
      drN1_s1{-0.25 * (1 - quad_1)},
      drN2_s0{0.25 * (1 - quad_0)},
      drN2_s1{0.25 * (1 - quad_1)},
      drN3_s0{0.25 * (1 + quad_0)},
      drN3_s1{0.25 * (1 + quad_1)},
      drN4_s0{-0.25 * (1 + quad_0)},
      drN4_s1{-0.25 * (1 + quad_1)};

  /// shape function derivative with respect to s, evaluated at various quadrature points
  constexpr static const T
      dsN1_r0{-0.25 * (1 - quad_0)},
      dsN1_r1{-0.25 * (1 - quad_1)},
      dsN2_r0{-0.25 * (1 + quad_0)},
      dsN2_r1{-0.25 * (1 + quad_1)},
      dsN3_r0{0.25 * (1 + quad_0)},
      dsN3_r1{0.25 * (1 + quad_1)},
      dsN4_r0{0.25 * (1 - quad_0)},
      dsN4_r1{0.25 * (1 - quad_1)};*/
  /** shape function derivative with respect to r and s, evaluated at various quadrature points */
  constexpr static const T dNi_d_alpha_alpha_k_data[16]{
      /** shape function derivative with respect to r, evaluated at various quadrature points */
      -0.25 * (1 - quad_0), /**< dN1dr evaluated at s=quad_0 */
      -0.25 * (1 - quad_1), /**< dN1dr evaluated at s=quad_1 */
      0.25 * (1 - quad_0),  /**< dN2dr evaluated at s=quad_0 */
      0.25 * (1 - quad_1),  /**< dN2dr evaluated at s=quad_1 */
      0.25 * (1 + quad_0),  /**< dN3dr evaluated at s=quad_0 */
      0.25 * (1 + quad_1),  /**< dN3dr evaluated at s=quad_1 */
      -0.25 * (1 + quad_0), /**< dN4dr evaluated at s=quad_0 */
      -0.25 * (1 + quad_1), /**< dN4dr evaluated at s=quad_1 */
      /** shape function derivative with respect to s, evaluated at various quadrature points */
      -0.25 * (1 - quad_0), /**< dN1ds evaluated at r=quad_0 */
      -0.25 * (1 - quad_1), /**< dN1ds evaluated at r=quad_1 */
      -0.25 * (1 + quad_0), /**< dN2ds evaluated at r=quad_0 */
      -0.25 * (1 + quad_1), /**< dN2ds evaluated at r=quad_1 */
      0.25 * (1 + quad_0),  /**< dN3ds evaluated at r=quad_0 */
      0.25 * (1 + quad_1),  /**< dN3ds evaluated at r=quad_1 */
      0.25 * (1 - quad_0),  /**< dN4ds evaluated at r=quad_0 */
      0.25 * (1 - quad_1)   /**< dN4ds evaluated at r=quad_1 */
  };
  /** getter for shape function derivatives evaluated at the various quadrature points*/  // TODO: doc
  constexpr static T dNi_d_alpha_alpha_k(const int alpha, const int i, const int k) {
    return dNi_d_alpha_alpha_k_data[alpha * 8 + i * 2 + k];
  }
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

//  x.test;
//  x.Ni_rj_sk(0, 0, 0);

  return 0;
}

}

#endif //A2D_SHELL_DEVELOPMENT_H_
