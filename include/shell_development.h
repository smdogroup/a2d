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

template <int N, typename T>  // TODO: make the Vec3ScaleDiv operation (v = (1/a) * x)
class A2DVec3A2DScaleDivExpr {

};

template <int N, typename T>  // TODO: make the ScalarAxpay operation (z = a * (x + y))
class ScalarA2DScalarA2DScalarAxpayExpr {

};

template <int N, typename T>  // TODO: make the ScalarAxpby operation (z = a * x + b * y)
class ScalarA2DScalarScalarA2DScalarAxpbyExpr {

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
               0, 0, 0, 0, temp3}, D(D_init) {};

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

  // TODO: move documentation inside class (as doc to constructor) after instances are converted to members in this
  //  class.
  // TODO: refactor these helper classes to <class name>Expr

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
  class g_alpha {
    g_alpha(const int alpha,
            const int n_alpha_var_ind,
            const T& t,
            const ShellElementMITC4<N, T>& element,
            A2DVec<N, Vec<T, 3>>& result) {

      /* Calculations for first node. */
      node1_thickness_scale_expression = ScalarMult(element.node1.thickness,
                                                    t * 0.5,
                                                    node1_scaled_thickness);
      node1_addition_expression = Vec3Axpy(node1_scaled_thickness,
                                           element.node1.shell_director,
                                           element.node1.position,
                                           node1_contribution_unscaled);
      node1_scale_expression = Vec3Scale(node1_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 0, n_alpha_var_ind),
                                         node1_contribution);

      /* Calculations for second node. */
      node2_thickness_scale_expression = ScalarMult(element.node2.thickness,
                                                    t * 0.5,
                                                    node2_scaled_thickness);
      node2_addition_expression = Vec3Axpy(node2_scaled_thickness,
                                           element.node2.shell_director,
                                           element.node2.position,
                                           node2_contribution_unscaled);
      node2_scale_expression = Vec3Scale(node2_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 1, n_alpha_var_ind),
                                         node2_contribution);

      /* Calculations for third node. */
      node3_thickness_scale_expression = ScalarMult(element.node3.thickness,
                                                    t * 0.5,
                                                    node3_scaled_thickness);
      node3_addition_expression = Vec3Axpy(node3_scaled_thickness,
                                           element.node3.shell_director,
                                           element.node3.position,
                                           node3_contribution_unscaled);
      node3_scale_expression = Vec3Scale(node3_contribution_unscaled,
                                         dNi_d_alpha_alpha_j(alpha, 2, n_alpha_var_ind),
                                         node3_contribution);

      /* Calculations for fourth node. */
      node4_thickness_scale_expression = ScalarMult(element.node4.thickness,
                                                    t * 0.5,
                                                    node4_scaled_thickness);
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
  class u_alpha {
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
  class g_t {
    g_t(const int r_ind,
        const int s_ind,
        const ShellElementMITC4<N, T>& element,
        A2DVec<N, Vec<T, 3>>& result) {
      /* Calculations for first node contribution */
      node1_thickness_scale_expression = ScalarMult(element.node1.thickness,
                                                    Ni_rj_sk(0, r_ind, s_ind) * 0.5,
                                                    node1_scaled_thickness);
      node1_scale_expression = Vec3Scale(element.node1.shell_director,
                                         node1_scaled_thickness,
                                         node1_contribution);
      /* Calculations for second node contribution */
      node2_thickness_scale_expression = ScalarMult(element.node2.thickness,
                                                    Ni_rj_sk(1, r_ind, s_ind) * 0.5,
                                                    node2_scaled_thickness);
      node2_scale_expression = Vec3Scale(element.node2.shell_director,
                                         node2_scaled_thickness,
                                         node2_contribution);
      /* Calculations for third node contribution */
      node3_thickness_scale_expression = ScalarMult(element.node3.thickness,
                                                    Ni_rj_sk(2, r_ind, s_ind) * 0.5,
                                                    node3_scaled_thickness);
      node3_scale_expression = Vec3Scale(element.node3.shell_director,
                                         node3_scaled_thickness,
                                         node3_contribution);
      /* Calculations for fourth node contribution */
      node4_thickness_scale_expression = ScalarMult(element.node4.thickness,
                                                    Ni_rj_sk(3, r_ind, s_ind) * 0.5,
                                                    node4_scaled_thickness);
      node4_scale_expression = Vec3Scale(element.node4.shell_director,
                                         node4_scaled_thickness,
                                         node4_contribution);

      /* Sum together the components. */
      sum_12_expression = Vec3Axpy(1.0, node1_contribution, node2_contribution, n1c_n2c);
      sum_123_expression = Vec3Axpy(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c);
      sum_1234_expression = Vec3Axpy(1.0, n1c_n2c_n3c, node4_contribution, result);
    };

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

    A2DScalarScalarMultExpr<N, T>
        node1_thickness_scale_expression,
        node2_thickness_scale_expression,
        node3_thickness_scale_expression,
        node4_thickness_scale_expression;
    Vec3A2DScaleExpr<N, T>
        node1_scale_expression,
        node2_scale_expression,
        node3_scale_expression,
        node4_scale_expression;
    A2DVec3A2DVecScalarAxpyExpr<N, T>
        sum_12_expression,
        sum_123_expression,
        sum_1234_expression;
  };

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
  class u_t {
    u_t(const int r_ind,
        const int s_ind,
        const ShellElementMITC4<N, T>& element,
        A2DVec<N, Vec<T, 3>>& result) {
      /* Calculations for first node. */
      node1_phi_expression = Vec3Cross(element.node1.rotation,
                                       element.node1.shell_director,
                                       node1_phi);
      node1_thickness_scale_expression = ScalarMult(element.node1.thickness,
                                                    Ni_rj_sk(0, r_ind, s_ind) * 0.5,
                                                    node1_scaled_thickness);
      node1_scale_expression = Vec3Scale(node1_phi,
                                         node1_scaled_thickness,
                                         node1_contribution);

      /* Calculations for second node. */
      node2_phi_expression = Vec3Cross(element.node2.rotation,
                                       element.node2.shell_director,
                                       node2_phi);
      node2_thickness_scale_expression = ScalarMult(element.node2.thickness,
                                                    Ni_rj_sk(1, r_ind, s_ind) * 0.5,
                                                    node2_scaled_thickness);
      node2_scale_expression = Vec3Scale(node2_phi,
                                         node2_scaled_thickness,
                                         node2_contribution);

      /* Calculations for third node. */
      node3_phi_expression = Vec3Cross(element.node3.rotation,
                                       element.node3.shell_director,
                                       node3_phi);
      node3_thickness_scale_expression = ScalarMult(element.node3.thickness,
                                                    Ni_rj_sk(2, r_ind, s_ind) * 0.5,
                                                    node3_scaled_thickness);
      node3_scale_expression = Vec3Scale(node3_phi,
                                         node3_scaled_thickness,
                                         node3_contribution);

      /* Calculations for fourth node. */
      node4_phi_expression = Vec3Cross(element.node4.rotation,
                                       element.node4.shell_director,
                                       node4_phi);
      node4_thickness_scale_expression = ScalarMult(element.node4.thickness,
                                                    Ni_rj_sk(3, r_ind, s_ind) * 0.5,
                                                    node4_scaled_thickness);
      node4_scale_expression = Vec3Scale(node4_phi,
                                         node4_scaled_thickness,
                                         node4_contribution);

      /* Sum together the components. */
      sum_12_expression = Vec3Axpy(1.0, node1_contribution, node2_contribution, n1c_n2c);
      sum_123_expression = Vec3Axpy(1.0, n1c_n2c, node3_contribution, n1c_n2c_n3c);
      sum_1234_expression = Vec3Axpy(1.0, n1c_n2c_n3c, node4_contribution, result);
    };

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

    A2DVec3CrossVecExpr<N, T>
        node1_phi_expression,
        node2_phi_expression,
        node3_phi_expression,
        node4_phi_expression;
    A2DScalarScalarMultExpr<N, T>
        node1_thickness_scale_expression,
        node2_thickness_scale_expression,
        node3_thickness_scale_expression,
        node4_thickness_scale_expression;
    A2DVec3A2DScaleExpr<N, T>
        node1_scale_expression,
        node2_scale_expression,
        node3_scale_expression,
        node4_scale_expression;
    A2DVec3A2DVecScalarAxpyExpr<N, T>
        sum_12_expression,
        sum_123_expression,
        sum_1234_expression;
  };

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
  class contravariant_basis {
    contravariant_basis(A2DVec<N, Vec<T, 3>>& gr, A2DVec<N, Vec<T, 3>>& gs, A2DVec<N, Vec<T, 3>>& gt,
                        A2DVec<N, Vec<T, 3>>& Gr, A2DVec<N, Vec<T, 3>>& Gs, A2DVec<N, Vec<T, 3>>& Gt) {
      gs_cross_gt_expression = Vec3Cross(gs, gt, gs_cross_gt);
      gt_cross_gr_expression = Vec3Cross(gt, gr, gt_cross_gr);
      gr_cross_gs_expression = Vec3Cross(gr, gs, gr_cross_gs);

      gr_dot_gs_cross_gt_expression = Vec3Dot(gr, gs_cross_gt, gr_dot_gs_cross_gt);

      Gr_expression = Vec3ScaleDiv(gr_dot_gs_cross_gt, gs_cross_gt, Gr);
      Gs_expression = Vec3ScaleDiv(gr_dot_gs_cross_gt, gt_cross_gr, Gs);
      Gt_expression = Vec3ScaleDiv(gr_dot_gs_cross_gt, gr_cross_gs, Gt);
    };

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
  class cartesian_local_basis {
    cartesian_local_basis(A2DVec<N, Vec<T, 3>>& gs, A2DVec<N, Vec<T, 3>>& gt,
                          A2DVec<N, Vec<T, 3>>& e1, A2DVec<N, Vec<T, 3>>& e2, A2DVec<N, Vec<T, 3>>& e3) {
      e3_expression = Vec3Normalize(gt, e3);
      gs_cross_e3_expression = Vec3Cross(gs, e3, gs_cross_e3);
      e1_expression = Vec3Normalize(gs_cross_e3, e1);
      e2_expression = Vec3Cross(e3, e1, e2);
    };

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
    A2DVec<N, Vec<T, 3>> gs_cross_e3;

    /* Expressions: */
    A2DVec3CrossA2DVecExpr<N, T>
        gs_cross_e3_expression,
        e2_expression;
    A2DVec3NormalizeExpr<N, T>
        e3_expression,
        e1_expression;
  };

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
  class in_plane_strain_componentsExpr{
    in_plane_strain_componentsExpr(A2DVec<N, Vec<T, 3>>& gr, A2DVec<N, Vec<T, 3>>& gs,
                                   A2DVec<N, Vec<T, 3>>& ur, A2DVec<N, Vec<T, 3>>& us,
                                   A2DVec<N, Vec<T, 3>>& e_rr, A2DVec<N, Vec<T, 3>>& e_ss, A2DVec<N, Vec<T, 3>>& e_rs){
      e_rr_expression = Vec3Dot(gr, ur, e_rr);
      e_ss_expression = Vec3Dot(gs, us, e_ss);
      gr_us_expression = Vec3Dot(gr, us, gr_us);
      gs_ur_expression = Vec3Dot(gs, ur, gs_ur);
      e_rs_expression = ScalarAxpay(0.5, gr_us, gs_ur, e_rs);
    };

    void reverse(){
      e_rs_expression.reverse();
      gs_ur_expression.reverse();
      gr_us_expression.reverse();
      e_ss_expression.reverse();
      e_rr_expression.reverse();
    };

    void hforward(){
      e_rr_expression.hforward();
      e_ss_expression.hforward();
      gr_us_expression.hforward();
      gs_ur_expression.hforward();
      e_rs_expression.hforward();
    };

    void hreverse(){
      e_rs_expression.hreverse();
      gs_ur_expression.hreverse();
      gr_us_expression.hreverse();
      e_ss_expression.hreverse();
      e_rr_expression.hreverse();
    };

   private:
    A2DScalar<N, T>
        gr_us, gs_ur;

    /* Expressions */
    A2DVec3DotA2DVecExpr<N, T>
        e_rr_expression, e_ss_expression,
        gr_us_expression, gs_ur_expression;
    ScalarA2DScalarA2DScalarAxpayExpr<N, T>
        e_rs_expression;
  };

  class epsilon_ri_sj_tk {
    // TODO:?
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
    /** gr vector evaluated at the tying points (s={1, -1} with t=0)*/
    gr_r0_sp1_t0, gr_r0_sn1_t0,
    /** gs vector evaluated at the tying points (r={1, -1} with t=0)*/
    gs_rp1_s0_t0, gs_rn1_s0_t0,
    /** derivatives of u with respect to r evaluated at the various quadrature points */
    ur_rAs0t0, ur_rAs0t1, ur_rAs1t0, ur_rAs1t1, /*ur_rAs0t0, ur_rAs0t1, ur_rAs1t0, ur_rAs1t1,*/
    /** derivatives of u with respect to s evaluated at the various quadrature points */
    us_r0sAt0, us_r0sAt1, /*us_r0sAt0, us_r0sAt1,*/ us_r1sAt0, us_r1sAt1, /*us_r1sAt0, us_r1sAt1,*/
    /** derivative of u with respect to t evaluated at the tying points (s={1, -1} with r=t=0, r={1, -1} with s=t=0)*/
    ut_r0_sp1_t0, ut_r0_sn1_t0, ut_rp1_s0_t0, ut_rn1_s0_t0,
    /** derivative of u with respect to r evaluated at the tying points (s={1,-1} with t=0)*/
    ur_r0_sp1_t0, ur_r0_sn1_t0,
    /** derivative of u with respect to s evaluated at the tying points (r={1,-1} with t=0)*/
    us_rp1_s0_t0, us_rn1_s0_t0;

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



    // code for tying points
    // TODO: remove non-dependencies *********
    /* e_rt_A calculations: */
    g_alpha gr_r0_sp1_t0_expression(0, 3, 0, this, gr_r0_sp1_t0);
    g_t gt_r0_sp1_t0_expression(2, 4, this, gt_r0_sp1_t0);
    u_alpha ur_r0_sp1_t0_expression(0, 3, 0, this, ur_r0_sp1_t0);
    u_t ut_r0_sp1_t0_expression(2, 4, this, ut_r0_sp1_t0);
    A2DScalar<N, T> gr_ut_A, gt_ur_A;  // TODO: move declaration
    A2DVec3DotA2DVecExpr gr_ut_A_expression = Vec3Dot(gr_r0_sp1_t0, ut_r0_sp1_t0, gr_ut_A);
    A2DVec3DotA2DVecExpr gt_ur_A_expression = Vec3Dot(gt_r0_sp1_t0, ur_r0_sp1_t0, gt_ur_A);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_rt_A_expression = ScalarAxpay(0.5, gr_ut_A, gt_ur_A, e_rt_A);
    /* e_rt_B calculations: */
    g_alpha gr_r0_sn1_t0_expression(0, 0, 0, this, gr_r0_sn1_t0);
    g_t gt_r0_sn1_t0_expression(2, 0, this, gt_r0_sn1_t0);
    u_alpha ur_r0_sn1_t0_expression(0, 0, 0, this, ur_r0_sn1_t0);
    u_t ut_r0_sn1_t0_expression(2, 0, this, ut_r0_sn1_t0);
    A2DScalar<N, T> gr_ut_B, gt_ur_B;  // TODO: move declaration
    A2DVec3DotA2DVecExpr gr_ut_B_expression = Vec3Dot(gr_r0_sn1_t0, ut_r0_sn1_t0, gr_ut_B);
    A2DVec3DotA2DVecExpr gt_ur_B_expression = Vec3Dot(gt_r0_sn1_t0, ur_r0_sn1_t0, gt_ur_B);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_rt_B_expression = ScalarAxpay(0.5, gr_ut_B, gt_ur_B, e_rt_B);
    /* e_st_C calculations: */
    g_alpha gs_rp1_s0_t0_expression(1, 3, 0, this, gs_rp1_s0_t0);
    g_t gt_rp1_s0_t0_expression(4, 2, this, gt_rp1_s0_t0);
    u_alpha us_rp1_s0_t0_expression(1, 3, 0, this, us_rp1_s0_t0);
    u_t ut_rp1_s0_t0_expression(4, 2, this, ut_rp1_s0_t0);
    A2DScalar<N, T> gs_ut_C, gt_us_C;  // TODO: move declaration
    A2DVec3DotA2DVecExpr gs_ut_C_expression = Vec3Dot(gs_rp1_s0_t0, ut_rp1_s0_t0, gs_ut_C);
    A2DVec3DotA2DVecExpr gt_us_C_expression = Vec3Dot(gt_rp1_s0_t0, us_rp1_s0_t0, gt_us_C);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_st_C_expression = ScalarAxpay(0.5, gs_ut_C, gt_us_C, e_st_C);
    /* e_st_D calculations: */
    g_alpha gs_rn1_s0_t0_expression(1, 0, 0, this, gs_rn1_s0_t0);
    g_t gt_rn1_s0_t0_expression(0, 2, this, gt_rn1_s0_t0);
    u_alpha us_rn1_s0_t0_expression(1, 0, 0, this, us_rn1_s0_t0);
    u_t ut_rn1_s0_t0_expression(0, 2, this, ut_rn1_s0_t0);
    A2DScalar<N, T> gs_ut_D, gt_us_D;  // TODO: move declaration
    A2DVec3DotA2DVecExpr gs_ut_D_expression = Vec3Dot(gs_rn1_s0_t0, ut_rn1_s0_t0, gs_ut_D);
    A2DVec3DotA2DVecExpr gt_us_D_expression = Vec3Dot(gt_rn1_s0_t0, us_rn1_s0_t0, gt_us_D);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_st_D_expression = ScalarAxpay(0.5, gs_ut_D, gt_us_D, e_st_D);


    /* Write code for e_rt_r... and e_st_r... values */
    /* Evaluate at s = quad_0; */
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_rt_rAs0tA_expression =
        ScalarAxpby((1 + quad_0) * 0.5, e_rt_A, (1 - quad_0) * 0.5, e_rt_B, e_rt_rAs0tA);
    /* Evaluate at r = quad_0; */
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_st_r0sAtA_expression =
        ScalarAxpby((1 + quad_0) * 0.5, e_st_C, (1 - quad_0) * 0.5, e_st_D, e_st_r0sAtA);
    /* Evaluate at s = quad_1; */
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_rt_rAs1tA_expression =
        ScalarAxpby((1 + quad_1) * 0.5, e_rt_A, (1 - quad_1) * 0.5, e_rt_B, e_rt_rAs1tA);
    /* Evaluate at r = quad_1; */
    ScalarA2DScalarScalarA2DScalarAxpbyExpr<N, T> e_st_r1sAtA_expression =
        ScalarAxpby((1 + quad_1) * 0.5, e_st_C, (1 - quad_1) * 0.5, e_st_D, e_st_r1sAtA);

    /* gr, gs, ur, and us calculations (subset of) */
    g_alpha gr_rAs0t0_expression(0, 1, quad_0, this, gr_rAs0t0);
    g_alpha gs_r0sAt0_expression(1, 1, quad_0, this, gs_r0sAt0);
    u_alpha ur_rAs0t0_expression(0, 1, quad_0, this, ur_rAs0t0);
    u_alpha us_r0sAt0_expression(1, 1, quad_0, this, us_r0sAt0);

    /* code for a single quadrature point: r0s0t0 <=> s=r=t=-1/sqrt(3) <=> s=r=t=quad_0*/
    /*A2DVec3DotA2DVecExpr<N, T> e_rr_r0s0t0_expression = Vec3Dot(gr_rAs0t0, ur_rAs0t0, e_rr_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> e_ss_r0s0t0_expression = Vec3Dot(gs_r0sAt0, us_r0sAt0, e_ss_r0s0t0);
    A2DScalar<N, T> gr_us_r0s0t0, gs_ur_r0s0t0;  // TODO: move declaration
    A2DVec3DotA2DVecExpr<N, T> gr_us_r0s0t0_expression = Vec3Dot(gr_rAs0t0, us_r0sAt0, gr_us_r0s0t0);
    A2DVec3DotA2DVecExpr<N, T> gs_ur_r0s0t0_expression = Vec3Dot(gs_r0sAt0, ur_rAs0t0, gs_ur_r0s0t0);
    ScalarA2DScalarA2DScalarAxpayExpr<N, T> e_rs_r0s0t0_expression =
        ScalarAxpay(0.5, gr_us_r0s0t0, gs_ur_r0s0t0, e_rs_r0s0t0);*/
    in_plane_strain_componentsExpr strain_r0s0t0_expression(gr_rAs0t0, gs_r0sAt0, ur_rAs0t0, us_r0sAt0,
                                                            e_rr_r0s0t0, e_ss_r0s0t0, e_rs_r0s0t0);
    // TODO: add code for other quadrature points


    /* Contravariant basis: */
    contravariant_basis contravariant_basis_r0s0t0_expression(); // TODO

    /* Cartesian local basis: */
    cartesian_local_basis local_basis_r0s0t0_expression(); // TODO


    /* Energy at a quadrature point: */
    A2DScalar<N, T> energy_r0s0t0, energy_r0s0t1 /* etc. */;
    energy_r0s0t0 = ;


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
