#ifndef A2D_HEAT_CONDUCTION_H
#define A2D_HEAT_CONDUCTION_H

#include "a2dcore.h"
#include "multiphysics/femapping.h"

namespace A2D {

/*
  Heat conduction problem

  The governing equation is

  -k div(u) = g in Omega
  u = u0 on S

  The weak form is

  Integral( k dot(grad(u), grad(v)) - g * v ) dOmega = 0

  for any test function v in the FE space
*/
template <typename T, index_t D>
class HeatConduction {
 public:
  HeatConduction(T kappa, T q, T heat_source)
      : kappa(kappa), q(q), heat_source(heat_source) {}

  // Spatial dimension
  static const index_t dim = D;

  // No data associated with this element
  static const index_t data_dim = 1;

  // Space for the finite-element data
  typedef FESpace<T, data_dim, L2Space<T, data_dim, dim>> DataSpace;

  // Finite element space
  typedef FESpace<T, dim, H1Space<T, 1, dim>> FiniteElementSpace;

  // Space for the element geometry - parametrized by H1
  typedef FESpace<T, dim, H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef SymMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  // Data for the element
  T kappa;        // Thermal conductivity
  T q;            // The RAMP penalty parameter
  T heat_source;  // the constant heat source term

  /**
   * @brief Find the integral of the thermal compliance over the entire domain
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

    T rho(data[0]);
    T k, dot, output;

    // Compute the RAMP penalty
    k = kappa / (1.0 + q * (1.0 - rho));
    VecDot(grad, grad, dot);  // dot = (grad, grad)
    output = wdetJ * (0.5 * dot - heat_source * u);

    return output;
  }

  /**
   * @brief Evaluate the derivatives of the weak form w.r.t. test function
   * components
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& data,
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

    // Make temporaries directly
    ADObj<T> rho(data[0]);
    ADObj<T> k, dot, output;

    auto stack = MakeStack(Eval(kappa / (1.0 + q * (1.0 - rho)), k),
                           VecDot(grad, grad, dot),  // dot = (grad, grad)
                           Eval(wdetJ * (0.5 * dot - heat_source * u), output));

    // Set the seed value
    output.bvalue() = 1.0;

    // Reverse the derivatives - values are written directly to the output
    // because we set the references
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

    // Make temporaries directly
    A2DObj<T> rho(data[0]);
    A2DObj<T> k, dot, output;

    // Make the stack
    auto stack = MakeStack(Eval(kappa / (1.0 + q * (1.0 - rho)), k),
                           VecDot(grad, grad, dot),  // dot = (grad, grad)
                           Eval(wdetJ * (0.5 * dot - heat_source * u), output));

    // Set the seed value
    output.bvalue() = 1.0;

    // Compute the derivatives first..
    stack.reverse();

    // Create data for extracting the Hessian-vector product
    constexpr index_t ncomp = FiniteElementSpace::ncomp;
    auto inters = MakeTieTuple<T, ADseed::h>(k);

    // Extract the matrix
    stack.template hextract<T, ncomp, ncomp>(inters, in, out, jac);
  }

  /**
   * @brief Compute the derivative of the adjoint-residual product with respect
   * to the data
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadrature point
   * @param adj The adjoint solution at the quadrature point
   */
  KOKKOS_FUNCTION void data_adjoint_product(T wdetJ, const DataSpace& data,
                                            const FiniteElementGeometry& geo,
                                            const FiniteElementSpace& s,
                                            const FiniteElementSpace& adj,
                                            DataSpace& dfdx) const {
    // Get the value and the gradient of the solution
    T u0 = s.template get<0>().get_value();    // Make a copy of the value
    T up = adj.template get<0>().get_value();  // Set reference for p
    T ub, uh;
    A2DObj<T> u(u0, ub, up, uh);

    Vec<T, dim> grad0 = s.template get<0>().get_grad();    // Make a copy
    Vec<T, dim> gradp = adj.template get<0>().get_grad();  // Set p ref
    Vec<T, dim> gradb, gradh;
    A2DObj<Vec<T, dim>&> grad(grad0, gradb, gradp, gradh);

    // Make temporaries directly
    A2DObj<T> rho(data[0]);
    A2DObj<T> k, dot, output;

    // Make the stack
    auto stack = MakeStack(Eval(kappa / (1.0 + q * (1.0 - rho)), k),
                           VecDot(grad, grad, dot),  // dot = (grad, grad)
                           Eval(wdetJ * (0.5 * dot - heat_source * u), output));

    // Set the seed value
    output.bvalue() = 1.0;

    // Compute the derivatives first..
    stack.reverse();

    // Perform the second derivatives to get the adjoint-vector product
    stack.hforward();
    stack.hreverse();

    // Get the derivative
    dfdx[0] = rho.hvalue();
  }
};

/*
  Mixed heat conduction problem

  The governing equations for the problem are

  q / kappa + grad(t) = 0        (heat flux)
  div(q) - f = 0                 (conservation of energy)

  Multiplying the first and second equations by the test functions qb and tb,
  respectively, and integrating over the domain gives

  (q / kappa + grad(t), qb) + (div(kappa * q) - f, tb) = 0

  Expanding this gives

  (q / kappa, qb) + (grad(t), qb) + (div(q), tb) = (f, tb)

  Since we have

  div(t * qb) = grad(t) * qb + t div(qb)  =>
  grad(t) * qb = div(t * qb) - t div(qb)

  Then (grad(t), qb) = - (t, div(qb)) so

  (q / kappa, qb) - (t, div(qb)) + (div(q), tb) = (f, tb)

  The weak coefficients must be provided by the function "weak".
*/
template <typename T, index_t D>
class MixedHeatConduction {
 public:
  // Spatial dimension
  static const index_t dim = D;

  // No data associated with this element
  static const index_t data_dim = 1;

  // Space for the finite-element data
  typedef FESpace<T, data_dim, L2Space<T, data_dim, dim>> DataSpace;

  // Finite element space
  typedef FESpace<T, dim, HdivSpace<T, dim>, L2Space<T, 1, dim>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1
  typedef FESpace<T, dim, H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef SymMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical
  // element
  using SolutionMapping = InteriorMapping<T, dim>;

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear
   elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    // Field objects for solution functions
    const Vec<T, dim>& q = s.template get<0>().get_value();
    const T& div = s.template get<0>().get_div();
    const T& temp = s.template get<1>().get_value();

    // Test function values
    Vec<T, dim>& qb = coef.template get<0>().get_value();
    T& divb = coef.template get<0>().get_div();
    T& tempb = coef.template get<1>().get_value();

    T kappa = data[0];

    // <q, qb> - (div(q), tempb) + (temp, div(qb)) = 0
    for (index_t k = 0; k < dim; k++) {
      qb(k) = wdetJ * q(k) / kappa;
    }

    divb = -wdetJ * kappa * temp;
    tempb = wdetJ * kappa * div;
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
    KOKKOS_FUNCTION JacVecProduct(const MixedHeatConduction<T, D>& integrand,
                                  T wdetJ, const DataSpace& data,
                                  const FiniteElementGeometry& geo,
                                  const FiniteElementSpace& s)
        : wdetJ(wdetJ), kappa(data[0]) {}

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const Vec<T, dim>& q = p.template get<0>().get_value();
      const T& div = p.template get<0>().get_div();
      const T& temp = p.template get<1>().get_value();

      // Test function values
      Vec<T, dim>& qb = Jp.template get<0>().get_value();
      T& divb = Jp.template get<0>().get_div();
      T& tempb = Jp.template get<1>().get_value();

      // <q, qb> - (div(q), tempb) + (temp, div(qb)) = 0
      for (index_t k = 0; k < dim; k++) {
        qb(k) = wdetJ * q(k) / kappa;
      }

      divb = -wdetJ * kappa * temp;
      tempb = wdetJ * kappa * div;
    }

   private:
    T wdetJ;
    T kappa;
  };
};

}  // namespace A2D

#endif  //  A2D_HEAT_CONDUCTION_H