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
 * The mixed poisson problem takes the form
 *
 * inf_{q} /sup_{u} Integral( 1/2 q^{T} q + u * div(q) + f * u )
 *
 * where u is in L2 and q is in H(div).
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
   * @brief Find the integral of the sadle point problem
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    // Extract the variables from the solution
    const Vec<T, dim>& q =
        s.template get<0>().get_value();            // Get the value of q
    const T& divq = s.template get<0>().get_div();  // Get the divergence of q
    const T& u = s.template get<1>().get_value();   // Get the value of u

    T dot, prod, sum, output;
    VecDot(q, q, dot);    // dot = ||q||_{2}^{2}
    Mult(divq, u, prod);  // prod = 0.5 * wdetJ * dot
    Sum(0.5 * wdetJ, dot, wdetJ, prod,
        sum);  // sum = wdetJ * (0.5 * dot + prod)
    Sum(wdetJ, u, T(1.0), sum, output);

    return output;
  }

  /**
   * @brief Evaluate the weak form of the coefficients
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
    // Extract the variables from the solution
    Vec<T, dim> q0 = s.template get<0>().get_value();  // Get the value of q
    T divq0 = s.template get<0>().get_div();  // Get the divergence of q
    T u0 = s.template get<1>().get_value();   // Get the value of u

    // Set references to the coefficients
    Vec<T, dim>& qb = coef.template get<0>().get_value();
    T& divqb = coef.template get<0>().get_div();
    T& ub = coef.template get<1>().get_value();

    ADObj<Vec<T, dim>&> q(q0, qb);
    ADObj<T&> u(u0, ub), divq(divq0, divqb);

    // Intermediate values
    ADObj<T> dot, prod, sum, output;

    auto stack = MakeStack(VecDot(q, q, dot),    // dot = ||q||_{2}^{2}
                           Mult(divq, u, prod),  // prod = 0.5 * wdetJ * dot
                           Sum(0.5 * wdetJ, dot, wdetJ, prod, sum),
                           Sum(wdetJ, u, T(1.0), sum, output));

    output.bvalue() = 1.0;
    stack.reverse();
  }

  /**
   * @brief Evaluate the Jacobian at a quadrature point
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param jac The Jacobian output
   */
  KOKKOS_FUNCTION void jacobian(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                QMatType& jac) const {
    FiniteElementSpace in, out;

    // Extract the variables from the solution
    Vec<T, dim> q0 = s.template get<0>().get_value();  // Get the value of q
    T divq0 = s.template get<0>().get_div();  // Get the divergence of q
    T u0 = s.template get<1>().get_value();   // Get the value of u

    // Set references to the coefficients
    Vec<T, dim> qb;
    T divqb, ub;

    // Extract the variables from the solution
    Vec<T, dim>& qp = in.template get<0>().get_value();  // Get the value of q
    T& divqp = in.template get<0>().get_div();  // Get the divergence of q
    T& up = in.template get<1>().get_value();   // Get the value of u

    // Extract the variables from the solution
    Vec<T, dim>& qh = out.template get<0>().get_value();  // Get the value of q
    T& divqh = out.template get<0>().get_div();  // Get the divergence of q
    T& uh = out.template get<1>().get_value();   // Get the value of u

    A2DObj<Vec<T, dim>&> q(q0, qb, qp, qh);
    A2DObj<T&> u(u0, ub, up, uh), divq(divq0, divqb, divqp, divqh);

    // Intermediate values
    A2DObj<T> dot, prod, sum, output;

    auto stack = MakeStack(VecDot(q, q, dot),    // dot = ||q||_{2}^{2}
                           Mult(divq, u, prod),  // prod = 0.5 * wdetJ * dot
                           Sum(0.5 * wdetJ, dot, wdetJ, prod, sum),
                           Sum(wdetJ, u, T(1.0), sum, output));

    output.bvalue() = 1.0;

    // Compute the derivatives first..
    stack.reverse();

    // Create data for extracting the Hessian-vector product
    constexpr index_t ncomp = FiniteElementSpace::ncomp;
    auto inters = MakeTieTuple<T, ADseed::h>(dot, prod, sum);

    // Extract the matrix
    stack.template hextract<T, ncomp, ncomp>(inters, in, out, jac);
  }
};

}  // namespace A2D

#endif  // A2D_POISSON_H
