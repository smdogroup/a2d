#ifndef A2D_POISSON_H
#define A2D_POISSON_H

#include "a2dcore.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fespace.h"

namespace A2D {

/**
 * @brief Regular Poisson problem
 *
 * @tparam T Scalar type for the calculation
 * @tparam D Dimension of the problem
 */
template <typename T, index_t D>
class Poisson {
 public:
  // Spatial dimension
  static const index_t dim = D;

  // No data associated with this element
  static const index_t data_dim = 0;

  // Space for the finite-element data
  using DataSpace = FESpace<T, data_dim>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, 1, dim>>;

  // Space for the element geometry - parametrized by H1 in 2D
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;

  // The type of matrix used to store data at each quadrature point
  using QMatType = SymMat<T, FiniteElementSpace::ncomp>;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& dobj,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    const H1Space<T, 1, dim>& u = s.template get<0>();
    const Vec<T, dim>& u_grad = u.get_grad();

    H1Space<T, 1, dim>& v = coef.template get<0>();
    Vec<T, dim>& v_grad = v.get_grad();

    // Set the terms from the variational statement
    for (index_t k = 0; k < dim; k++) {
      v_grad(k) = wdetJ * u_grad(k);
    }
  }

  /**
   * @brief Compute the Jacobian
   *
   * This functor computes a Jacobian-vector product of the weak form
   *
   * @param integrand The Integrand object for this class
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadrature point
   */
  KOKKOS_FUNCTION void jacobian(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                QMatType& jac) const {}
};

/**
 * @brief Mixed Poisson problem discretization
 *
 * @tparam T Scalar type for the calculation
 * @tparam D Dimension of the problem
 */
template <typename T, index_t D>
class MixedPoisson {
 public:
  // Spatial dimension
  static const index_t dim = D;

  // No data associated with this element
  static const index_t data_dim = 0;

  // Space for the finite-element data
  typedef FESpace<T, data_dim> DataSpace;

  // Finite element space
  typedef FESpace<T, dim, HdivSpace<T, dim>, L2Space<T, 1, dim>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1 in 2D
  typedef FESpace<T, dim, H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef SymMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& dobj,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    // Field objects for solution functions
    const HdivSpace<T, dim>& sigma = s.template get<0>();
    const L2Space<T, 1, dim>& u = s.template get<1>();

    // Solution function values
    const Vec<T, dim>& sigma_val = sigma.get_value();
    const T& sigma_div = sigma.get_div();
    const T& u_val = u.get_value();

    // Test function values
    HdivSpace<T, dim>& tau = coef.template get<0>();
    L2Space<T, 1, dim>& v = coef.template get<1>();

    // Test function values
    Vec<T, dim>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    for (index_t k = 0; k < dim; k++) {
      tau_val(k) = wdetJ * sigma_val(k);
    }

    tau_div = -wdetJ * u_val;
    v_val = wdetJ * sigma_div;
  }

  /**
   * @brief Construct a JacVecProduct functor
   *
   * This functor computes a Jacobian-vector product of the
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class JacVecProduct {
   public:
    KOKKOS_FUNCTION JacVecProduct(const MixedPoisson<T, D>& integrand, T wdetJ,
                                  const DataSpace& data,
                                  const FiniteElementGeometry& geo,
                                  const FiniteElementSpace& s)
        : wdetJ(wdetJ) {}

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const HdivSpace<T, dim>& sigma = p.template get<0>();
      const L2Space<T, 1, dim>& u = p.template get<1>();

      // Solution function values
      const Vec<T, dim>& sigma_val = sigma.get_value();
      const T& sigma_div = sigma.get_div();
      const T& u_val = u.get_value();

      // Test function values
      HdivSpace<T, dim>& tau = Jp.template get<0>();
      L2Space<T, 1, dim>& v = Jp.template get<1>();

      // Test function values
      Vec<T, dim>& tau_val = tau.get_value();
      T& tau_div = tau.get_div();
      T& v_val = v.get_value();

      // Set the terms from the variational statement
      for (index_t k = 0; k < dim; k++) {
        tau_val(k) = wdetJ * sigma_val(k);
      }

      tau_div = -wdetJ * u_val;
      v_val = wdetJ * sigma_div;
    }

   private:
    T wdetJ;
  };
};

}  // namespace A2D

#endif  // A2D_POISSON_H
