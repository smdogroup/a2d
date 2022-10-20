#ifndef A2D_POISSON_H
#define A2D_POISSON_H

#include "a2dmatops2d.h"
#include "multiphysics/feelement.h"
#include "multiphysics/fespace.h"

/*
  The PDE object
*/
template <typename T>
class MixedPoisson2D {
 public:
  // Spatial dimension
  static const A2D::index_t dim = 2;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::L1ScalarSpace<T, dim>, A2D::Hdiv2DSpace<T>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1 in 2D
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  /**
   * Evaluate the variational form at a quadrature point
   *
   * <tau, sigma> - <div(tau), u> + <v, div(sigma)>
   */
  static T eval_weak_form(T wdetJ, FiniteElementSpace& s,
                          FiniteElementSpace& t) {
    // Field objects for solution functions
    A2D::L1ScalarSpace<T, 2>& u = s.template get<0>();
    A2D::Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    A2D::L1ScalarSpace<T, 2>& v = t.template get<0>();
    A2D::Hdiv2DSpace<T>& tau = t.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // return A2D::VecDot(tau_val, sigma_val) - tau_div * u + v * sigma_div;
    return wdetJ * (tau_val(0) * sigma_val(0) + tau_val(0) * sigma_val(0) -
                    tau_div * u + v * sigma_div);
  }

  static void eval_weak_coef(T wdetJ, FiniteElementSpace& s,
                             FiniteElementSpace& coef) {
    // Field objects for solution functions
    A2D::L1ScalarSpace<T, 2>& u = s.template get<0>();
    A2D::Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    A2D::L1ScalarSpace<T, 2>& v = coef.template get<0>();
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

  static void eval_weak_jacobian_vec_product(T wdetJ, FiniteElementSpace& s,
                                             FiniteElementSpace& p,
                                             FiniteElementSpace& coef) {
    // Field objects for solution functions
    A2D::L1ScalarSpace<T, 2>& u = p.template get<0>();
    A2D::Hdiv2DSpace<T>& sigma = p.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    A2D::L1ScalarSpace<T, 2>& v = coef.template get<0>();
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
