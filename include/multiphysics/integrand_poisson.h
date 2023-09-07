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
   * @brief Find the integral of the compliance over the entire domain
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    // Get the value and the gradient of the solution
    const T& u = s.template get<0>().get_value();
    const Vec<T, dim>& grad = s.template get<0>().get_grad();

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    T dot, output;
    VecDot(grad, grad, dot);
    Sum(0.5 * wdetJ, dot, -wdetJ, T(1.0), output);

    return output;
  }

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
    // Get the value and the gradient of the solution
    T u0 = s.template get<0>().get_value();  // Make a copy of the value
    Vec<T, dim> grad0 = s.template get<0>().get_grad();  // Make a copy of grad

    // Grab references to the output values
    T& ub = coef.template get<0>().get_value();
    Vec<T, dim>& gradb = coef.template get<0>().get_grad();

    // Make the objects with references
    ADObj<T&> u(u0, ub);
    ADObj<Vec<T, dim>&> grad(grad0, gradb);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    ADObj<T> dot, output;
    auto stack = MakeStack(VecDot(grad, grad, dot),
                           Sum(0.5 * wdetJ, dot, -wdetJ, T(1.0), output));

    output.bvalue() = 1.0;

    // Reverse the derivatives - values are written directly to the output
    // because we set the references
    stack.reverse();
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
                                QMatType& jac) const {
    FiniteElementSpace in, out;

    // Get the value and the gradient of the solution
    T u0 = s.template get<0>().get_value();     // Make a copy of the value
    T ub;                                       // Value for b
    T& up = in.template get<0>().get_value();   // Set reference for p
    T& uh = out.template get<0>().get_value();  // Set reference for h
    A2DObj<T&> u(u0, ub, up, uh);

    Vec<T, dim> grad0 = s.template get<0>().get_grad();     // Make a copy
    Vec<T, dim> gradb;                                      // Value for b
    Vec<T, dim>& gradp = in.template get<0>().get_grad();   // Set p ref
    Vec<T, dim>& gradh = out.template get<0>().get_grad();  // Set h ref
    A2DObj<Vec<T, dim>&> grad(grad0, gradb, gradp, gradh);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    A2DObj<T> dot, output;
    auto stack = MakeStack(VecDot(grad, grad, dot),
                           Sum(0.5 * wdetJ, dot, -wdetJ, T(1.0), output));

    output.bvalue() = 1.0;

    // Compute the derivatives first..
    stack.reverse();

    // Create data for extracting the Hessian-vector product
    constexpr index_t ncomp = FiniteElementSpace::ncomp;
    auto inters = MakeTieTuple<T, ADseed::h>(dot, output);

    // Extract the matrix
    stack.template hextract<T, ncomp, ncomp>(inters, in, out, jac);
  }
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
