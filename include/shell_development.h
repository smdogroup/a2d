//
// Created by James on 10/23/2022.
//

#ifndef A2D_SHELL_DEVELOPMENT_H_
#define A2D_SHELL_DEVELOPMENT_H_

#include "a2dtypes.h"
#include "a2dvecops3d.h"
#include "a2dmatops3d.h"

namespace A2D {

// TODO: OPERATOR CLASSES:

template <int N, typename T>  // TODO: make the ScalarMult operation (z = a * b)
class A2DScalarScalarMultExpr {

};


// END TODO: OPERATOR CLASSES

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
     * @param result:           An A2DVec where the resulting g<sub>&alpha</sub> vector should be stored.
     * */
    g_alpha(const int alpha,
            const int n_alpha_var_ind,
            const T& t,
            const ShellElementMITC4<N, T>& element,
            A2DVec<N, Vec<T, 3>>& result) {

      /* Calculations for first node. */
      node1_thickness_scale_expression = ScalarMult(element.node1.thickness, t * 0.5, node1_scaled_thickness);
      node1_addition_expression = Vec3Axpy(node1_scaled_thickness,
                                           element.node1.shell_director,
                                           element.node1.position,
                                           node1_contribution_unscaled);
      node1_scale_expression = Vec3Scale(node1_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 0, n_alpha_var_ind),
                                         node1_contribution);

      /* Calculations for second node. */
      node2_thickness_scale_expression = ScalarMult(element.node2.thickness, t * 0.5, node2_scaled_thickness);
      node2_addition_expression = Vec3Axpy(node2_scaled_thickness,
                                           element.node2.shell_director,
                                           element.node2.position,
                                           node2_contribution_unscaled);
      node2_scale_expression = Vec3Scale(node2_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 1, n_alpha_var_ind),
                                         node2_contribution);

      /* Calculations for third node. */
      node3_thickness_scale_expression = ScalarMult(element.node3.thickness, t * 0.5, node3_scaled_thickness);
      node3_addition_expression = Vec3Axpy(node3_scaled_thickness,
                                           element.node3.shell_director,
                                           element.node3.position,
                                           node3_contribution_unscaled);
      node3_scale_expression = Vec3Scale(node3_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 2, n_alpha_var_ind),
                                         node3_contribution);

      /* Calculations for fourth node. */
      node4_thickness_scale_expression = ScalarMult(element.node4.thickness, t * 0.5, node4_scaled_thickness);
      node4_addition_expression = Vec3Axpy(node4_scaled_thickness,
                                           element.node4.shell_director,
                                           element.node4.position,
                                           node4_contribution_unscaled);
      node4_scale_expression = Vec3Scale(node4_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 3, n_alpha_var_ind),
                                         node4_contribution);

      /* Sum together the components. */
      sum_12_expression = Vec3Axpy(1.0, node1_contribution, node2_contribution, n1c_n2c);
      sum_123_expression = Vec3Axpy(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c);
      sum_1234_expression = Vec3Axpy(1.0, n1c_n2c_n3c, node4_contribution, result);
    }

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
    A2DVec<N, Vec<T, 3>>
        node1_contribution, node1_contribution_unscaled,
        node2_contribution, node2_contribution_unscaled,
        node3_contribution, node3_contribution_unscaled,
        node4_contribution, node4_contribution_unscaled;
    A2DScalar<N, T>
        node1_scaled_thickness,
        node2_scaled_thickness,
        node3_scaled_thickness,
        node4_scaled_thickness;
    A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

    /* Expressions: */

    Vec3VecA2DScalarAxpyExpr<N, T>
        node1_addition_expression,
        node2_addition_expression,
        node3_addition_expression,
        node4_addition_expression;
    A2DVec3ScaleExpr<N, T>
        node1_scale_expression,
        node2_scale_expression,
        node3_scale_expression,
        node4_scale_expression;
    A2DScalarScalarMultExpr<N, T>
        node1_thickness_scale_expression,
        node2_thickness_scale_expression,
        node3_thickness_scale_expression,
        node4_thickness_scale_expression;
    A2DVec3A2DVecScalarAxpyExpr<N, T>
        sum_12_expression,
        sum_123_expression,
        sum_1234_expression;
  };

  class u_alpha {
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
     * @param result:           An A2DVec where the resulting du/d&alpha vector should be stored.
     * */
    u_alpha(const int alpha,
            const int n_alpha_var_ind,
            const T& t,
            const ShellElementMITC4<N, T>& element,
            A2DVec<N, Vec<T, 3>>& result) {

      /* Calculations for first node. */
      node1_phi_expression = Vec3Cross(element.node1.rotation,
                                       element.node1.shell_director,
                                       node1_phi);
      node1_thickness_scale_expression = ScalarMult(element.node1.thickness,
                                                    t * 0.5,
                                                    node1_scaled_thickness);
      node1_addition_expression = Vec3Axpy(node1_scaled_thickness,
                                           node1_phi,
                                           element.node1.displacement,
                                           node1_contribution_unscaled);
      node1_scale_expression = Vec3Scale(node1_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 0, n_alpha_var_ind),
                                         node1_contribution);

      /* Calculations for second node. */
      node2_phi_expression = Vec3Cross(element.node2.rotation,
                                       element.node2.shell_director,
                                       node2_phi);
      node2_thickness_scale_expression = ScalarMult(element.node2.thickness,
                                                    t * 0.5,
                                                    node2_scaled_thickness);
      node2_addition_expression = Vec3Axpy(node2_scaled_thickness,
                                           node2_phi,
                                           element.node2.displacement,
                                           node2_contribution_unscaled);
      node2_scale_expression = Vec3Scale(node2_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 1, n_alpha_var_ind),
                                         node2_contribution);

      /* Calculations for third node. */
      node3_phi_expression = Vec3Cross(element.node3.rotation,
                                       element.node3.shell_director,
                                       node3_phi);
      node3_thickness_scale_expression = ScalarMult(element.node3.thickness,
                                                    t * 0.5,
                                                    node3_scaled_thickness);
      node3_addition_expression = Vec3Axpy(node3_scaled_thickness,
                                           node3_phi,
                                           element.node3.displacement,
                                           node3_contribution_unscaled);
      node3_scale_expression = Vec3Scale(node3_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 2, n_alpha_var_ind),
                                         node3_contribution);

      /* Calculations for fourth node. */
      node4_phi_expression = Vec3Cross(element.node4.rotation,
                                       element.node4.shell_director,
                                       node4_phi);
      node4_thickness_scale_expression = ScalarMult(element.node4.thickness,
                                                    t * 0.5,
                                                    node4_scaled_thickness);
      node4_addition_expression = Vec3Axpy(node4_scaled_thickness,
                                           node4_phi,
                                           element.node4.displacement,
                                           node4_contribution_unscaled);
      node4_scale_expression = Vec3Scale(node4_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 3, n_alpha_var_ind),
                                         node4_contribution);

      /* Sum together the components. */
      sum_12_expression = Vec3Axpy(1.0, node1_contribution, node2_contribution, n1c_n2c);
      sum_123_expression = Vec3Axpy(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c);
      sum_1234_expression = Vec3Axpy(1.0, n1c_n2c_n3c, node4_contribution, result);
    }

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
    A2DVec<N, Vec<T, 3>>
        node1_phi,
        node2_phi,
        node3_phi,
        node4_phi;
    A2DVec<N, Vec<T, 3>>
        node1_contribution, node1_contribution_unscaled,
        node2_contribution, node2_contribution_unscaled,
        node3_contribution, node3_contribution_unscaled,
        node4_contribution, node4_contribution_unscaled;
    A2DScalar<N, T>
        node1_scaled_thickness,
        node2_scaled_thickness,
        node3_scaled_thickness,
        node4_scaled_thickness;
    A2DVec<N, Vec<T, 3>> n1c_n2c, n1c_n2c_n3c;

    /* Expressions: */

    A2DVec3CrossVecExpr<N, T>
        node1_phi_expression,
        node2_phi_expression,
        node3_phi_expression,
        node4_phi_expression;
    Vec3VecA2DScalarAxpyExpr<N, T>
        node1_addition_expression,
        node2_addition_expression,
        node3_addition_expression,
        node4_addition_expression;
    A2DVec3ScaleExpr<N, T>
        node1_scale_expression,
        node2_scale_expression,
        node3_scale_expression,
        node4_scale_expression;
    A2DScalarScalarMultExpr<N, T>
        node1_thickness_scale_expression,
        node2_thickness_scale_expression,
        node3_thickness_scale_expression,
        node4_thickness_scale_expression;
    A2DVec3A2DVecScalarAxpyExpr<N, T>
        sum_12_expression,
        sum_123_expression,
        sum_1234_expression;
  };

  class epsilon_ri_sj_tk {
    // TODO:
    epsilon_ri_sj_tk(const int i, const int j, const int k) {

    }
  };

  void generate_energy() {
    /* Generate the internal energy of the shell element, that will then be used to calculate the global stiffness
     * matrix
     * */

    A2DVec<N, Vec<T, 3>>
    /** gr vector evaluated at the various quadrature points */
    gr_rAs0t0, gr_rAs0t1, gr_rAs1t0, gr_rAs1t1, /*gr_rAs0t0, gr_rAs0t1, gr_rAs1t0, gr_rAs1t1,*/
    /** gs vector evaluated at the various quadrature points */
    gs_r0sAt0, gs_r0sAt1, /*gs_r0sAt0, gs_r0sAt1,*/ gs_r1sAt0, gs_r1sAt1, /*gs_r1sAt0, gs_r1sAt1,*/
    /** gt vector evaluated at the tying points (s={1, -1} with r=t=0, r={1, -1} with s=t=0)*/
    gt_r0_sp1_t0, gt_r0_sn1_t0, gt_rp1_s0_t0, gt_rn1_s0_t0,
    /** derivatives of u with respect to r evaluated at the various quadrature points */
    ur_rAs0t0, ur_rAs0t1, ur_rAs1t0, ur_rAs1t1, /*ur_rAs0t0, ur_rAs0t1, ur_rAs1t0, ur_rAs1t1,*/
    /** derivatives of u with respect to s evaluated at the various quadrature points */
    us_r0sAt0, us_r0sAt1, /*us_r0sAt0, us_r0sAt1,*/ us_r1sAt0, us_r1sAt1, /*us_r1sAt0, us_r1sAt1,*/
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
    /*TODO: move all these to private members, and just reset their values whenever this method is called.*/

    // TODO: write code for tying points
    // TODO: remove non-dependencies *********
    /* e_rt_A calculations: */
    /*auto gr_r0_sp1_t0_expression = ;
    auto gt_r0_sp1_t0_expression = ;
    auto ur_r0_sp1_t0_expression = ;
    auto ut_r0_sp1_t0_expression = ;
    A2DScalar<N, T> gr_ut_A, gt_ur_A;  // TODO: move declaration
    auto gr_ut_A_expression = Vec3Dot(gr_r0_sp1_t0, ut_r0_sp1_t0, gr_ut_A);
    auto gt_ur_A_expression = Vec3Dot(gt_r0_sp1_t0, ur_r0_sp1_t0, gt_ur_A);
    auto e_rt_A_expression = ScalarAxpay(0.5, gr_ut_A, gt_ur_A, e_rt_A);*/
    /* e_rt_B calculations: */
    /*auto gr_r0_sn1_t0_expression = ;
    auto gt_r0_sn1_t0_expression = ;
    auto ur_r0_sn1_t0_expression = ;
    auto ut_r0_sn1_t0_expression = ;
    A2DScalar<N, T> gr_ut_B, gt_ur_B;  // TODO: move declaration
    auto gr_ut_B_expression = Vec3Dot(gr_r0_sn1_t0, ut_r0_sn1_t0, gr_ut_B);
    auto gt_ur_B_expression = Vec3Dot(gt_r0_sn1_t0, ur_r0_sn1_t0, gt_ur_B);
    auto e_rt_B_expression = ScalarAxpay(0.5, gr_ut_B, gt_ur_B, e_rt_B);*/
    /* e_st_C calculations: */
    /*auto gs_rp1_s0_t0_expression = ;
    auto gt_rp1_s0_t0_expression = ;
    auto us_rp1_s0_t0_expression = ;
    auto ut_rp1_s0_t0_expression = ;
    A2DScalar<N, T> gs_ut_C, gt_us_C;  // TODO: move declaration
    auto gs_ut_C_expression = Vec3Dot(gs_rp1_s0_t0, ut_rp1_s0_t0, gs_ut_C);
    auto gt_us_C_expression = Vec3Dot(gt_rp1_s0_t0, us_rp1_s0_t0, gt_us_C);
    auto e_st_C_expression = ScalarAxpay(0.5, gs_ut_C, gt_us_C, e_st_C);*/
    /* e_st_D calculations: */
    /*auto gs_rn1_s0_t0_expression = ;
    auto gt_rn1_s0_t0_expression = ;
    auto us_rn1_s0_t0_expression = ;
    auto ut_rn1_s0_t0_expression = ;
    A2DScalar<N, T> gs_ut_D, gt_us_D;  // TODO: move declaration
    auto gs_ut_D_expression = Vec3Dot(gs_rn1_s0_t0, ut_rn1_s0_t0, gs_ut_D);
    auto gt_us_D_expression = Vec3Dot(gt_rn1_s0_t0, us_rn1_s0_t0, gt_us_D);
    auto e_st_D_expression = ScalarAxpay(0.5, gs_ut_D, gt_us_D, e_st_D);*/


    /* Write code for e_rt_r... and e_st_r... values */
    /* Evaluate at s = quad_0; */
    auto e_rt_rAs0tA_expression = ScalarAxpby((1 + quad_0) * 0.5, e_rt_A, (1 - quad_0) * 0.5, e_rt_B, e_rt_rAs0tA);
    // TODO: make the ScalarAxpby operation (z = a * x + b * y)
    /* Evaluate at r = quad_0; */
    auto e_st_r0sAtA_expression = ScalarAxpby((1 + quad_0) * 0.5, e_st_C, (1 - quad_0) * 0.5, e_st_D, e_st_r0sAtA);
    /* Evaluate at s = quad_1; */
    auto e_rt_rAs1tA_expression = ScalarAxpby((1 + quad_1) * 0.5, e_rt_A, (1 - quad_1) * 0.5, e_rt_B, e_rt_rAs1tA);
    /* Evaluate at r = quad_1; */
    auto e_st_r1sAtA_expression = ScalarAxpby((1 + quad_1) * 0.5, e_st_C, (1 - quad_1) * 0.5, e_st_D, e_st_r1sAtA);

    /* write test code for a single quadrature point: r0s0t0 <=> s=r=t=-1/sqrt(3) <=> s=r=t=quad_0*/
    g_alpha gr_rAs0t0_expression(0, 1, quad_0, this, gr_rAs0t0);
    g_alpha gs_r0sAt0_expression(1, 1, quad_0, this, gs_r0sAt0);
    u_alpha ur_rAs0t0_expression(0, 1, quad_0, this, ur_rAs0t0);
    u_alpha us_r0sAt0_expression(1, 1, quad_0, this, us_r0sAt0);
    A2DVec3DotA2DVecExpr<N, T> e_rr_r0s0t0_expression = Vec3Dot(gr_rAs0t0, ur_rAs0t0, e_rr_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> e_ss_r0s0t0_expression = Vec3Dot(gs_r0sAt0, us_r0sAt0, e_ss_r0s0t0);
    A2DScalar<N, T> gr_us_r0s0t0, gs_ur_r0s0t0;  // TODO: move declaration
    A2DVec3DotA2DVecExpr<N, T> gr_us_r0s0t0_expression = Vec3Dot(gr_rAs0t0, us_r0sAt0, gr_us_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> gs_ur_r0s0t0_expression = Vec3Dot(gs_r0sAt0, ur_rAs0t0, gs_ur_r0s0t0);
    auto e_rs_r0s0t0_expression = ScalarAxpay(0.5, gr_us_r0s0t0, gs_ur_r0s0t0, e_rs_r0s0t0);
    // TODO: make the ScalarAxpay operation (z = a * (x + y))

    /*A2DVec<N, Vec<T, 5>> epsilon_r0s0t0;  // TODO: move declaration
    epsilon_ri_sj_tk epsilon_r0s0t0_expression();*/
  };

  ShellNodeMITC<N, T>& node1;  /**< The top left node. */
  ShellNodeMITC<N, T>& node2;  /**< The top right node.*/
  ShellNodeMITC<N, T>& node3;  /**< The bottom left node. */
  ShellNodeMITC<N, T>& node4;  /**< The bottom right node.*/
  // TODO: make sure this is correct ^

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
  /** getter for shape functions evaluated at the various quadrature points */  // TODO: doc
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
  /** @brief getter for shape function derivatives evaluated at the various quadrature points
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