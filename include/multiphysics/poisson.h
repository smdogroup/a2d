#ifndef A2D_POISSON_H
#define A2D_POISSON_H

#include "a2dmatops2d.h"
#include "multiphysics/feelement.h"
#include "multiphysics/fespace.h"

/*
  The PDE object
*/
template <typename T, index_t D>
class MixedPoisson2D {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;

  // No data associated with this element
  static const A2D::index_t data_dim = 0;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim> DataSpace;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::L2Space<T, 1, dim>, A2D::HdivSpace<T, D>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1 in 2D
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear elasticity
   *
   * <tau, sigma> - <div(tau), u> + <v, div(sigma)>
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  static T eval_weak_form(T wdetJ, const DataSpace& dobj,
                          const FiniteElementSpace& s,
                          FiniteElementSpace& coef) {
    // Field objects for solution functions
    A2D::L2Space<T, 1, 2>& u = s.template get<0>();
    A2D::Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    A2D::L2Space<T, 1, 2>& v = t.template get<0>();
    A2D::Hdiv2DSpace<T>& tau = t.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // return A2D::VecDot(tau_val, sigma_val) - tau_div * u + v * sigma_div;
    return wdetJ * (tau_val(0) * sigma_val(0) + tau_val(0) * sigma_val(0) -
                    tau_div * u + v * sigma_div);
  }

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  static void eval_weak_coef(T wdetJ, const DataSpace& dobj,
                             const FiniteElementSpace& s,
                             FiniteElementSpace& coef) {
    // Field objects for solution functions
    const A2D::L2Space<T, 1, 2>& u = s.template get<0>();
    const A2D::Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    const A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    const T& sigma_div = sigma.get_div();
    const T& u_val = u.get_value();

    // Test function values
    A2D::L2Space<T, 1, 2>& v = coef.template get<0>();
    A2D::Hdiv2DSpace<T>& tau = coef.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    tau_val(0) = wdetJ * sigma_val(0);
    tau_val(1) = wdetJ * sigma_val(1);

    tau_div = -wdetJ * u_val;

    v_val = wdetJ * sigma_div;
  }

  /**
   * @brief Evaluate the derivative of the weak form coefficients along the
   * solution direction p
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param s The trial solution
   * @param p The direction of the tangent
   * @param coef The coeffcients
   */
  static JVP_t eval_weak_jacobian_vec_product(T wdetJ, const DataSpace& dobj,
                                              const FiniteElementSpace& s,
                                              const FiniteElementSpace& p,
                                              FiniteElementSpace& coef) {
    // Field objects for solution functions
    const A2D::L2Space<T, 1, 2>& u = p.template get<0>();
    const A2D::Hdiv2DSpace<T>& sigma = p.template get<1>();

    // Solution function values
    const A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    const T& sigma_div = sigma.get_div();
    const T& u_val = u.get_value();

    // Test function values
    A2D::L2Space<T, 1, 2>& v = coef.template get<0>();
    A2D::Hdiv2DSpace<T>& tau = coef.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    tau_val(0) = wdetJ * sigma_val(0);
    tau_val(1) = wdetJ * sigma_val(1);

    tau_div = -wdetJ * u_val;

    v_val = wdetJ * sigma_div;
  }
};

#endif  // A2D_POISSON_H
