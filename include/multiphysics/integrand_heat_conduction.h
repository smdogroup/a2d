#ifndef A2D_HEAT_CONDUCTION_H
#define A2D_HEAT_CONDUCTION_H

#include "multiphysics/femapping.h"

namespace A2D {

/*
  Heat conduction problem

  The governing equation is

  -k∆u = g in Ω
  u = u0 on ∂Ω

  The weak form is

  ∫ k ∇u ∇v dΩ = ∫gv dΩ for test function v

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
    return wdetJ * (s.template get<0>()).get_value() * heat_source;
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
  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    // Field objects for solution functions
    const Vec<T, dim>& tx = s.template get<0>().get_grad();
    Vec<T, dim>& cx = coef.template get<0>().get_grad();
    T& c = coef.template get<0>().get_value();
    c = -wdetJ * heat_source;  // Add a constant heat source here

    // Set the thermal conductivity coefficient
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    for (index_t k = 0; k < dim; k++) {
      cx(k) = wdetJ * kappa * penalty * rho * tx(k);
    }
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
    KOKKOS_FUNCTION JacVecProduct(const HeatConduction<T, D>& integrand,
                                  T wdetJ, const DataSpace& data,
                                  const FiniteElementGeometry& geo,
                                  const FiniteElementSpace& s)
        : wdetJ(wdetJ),
          kappa(integrand.kappa),
          rho(data[0]),
          penalty(1.0 / (1.0 + integrand.q * (1.0 - rho))) {}

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const Vec<T, dim>& tx = p.template get<0>().get_grad();
      Vec<T, dim>& cx = Jp.template get<0>().get_grad();

      for (index_t k = 0; k < dim; k++) {
        cx(k) = wdetJ * kappa * penalty * rho * tx(k);
      }
    }

   private:
    T wdetJ;
    T kappa;
    T rho;
    T penalty;
  };

  /**
   * @brief Construct a JacVecProduct functor
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class AdjVecProduct {
   public:
    KOKKOS_FUNCTION AdjVecProduct(const HeatConduction<T, D>& integrand,
                                  T wdetJ, const DataSpace& data,
                                  const FiniteElementGeometry& geo,
                                  const FiniteElementSpace& s)
        : wdetJ(wdetJ),
          kappa(integrand.kappa),
          rho(data[0]),
          q(integrand.q),
          penalty(1.0 / (1.0 + integrand.q * (1.0 - rho))),
          tx(s.template get<0>().get_grad()) {}

    KOKKOS_FUNCTION void operator()(const FiniteElementSpace& p,
                                    DataSpace& dfdx) {
      T denom = (1.0 + q * (1.0 - rho));
      T dprhodrho = (q + 1.0) / (denom * denom);

      const Vec<T, dim>& psi = p.template get<0>().get_grad();

      for (index_t k = 0; k < dim; k++) {
        dfdx[0] += wdetJ * kappa * dprhodrho * tx(k) * psi(k);
      }
    }

   private:
    T wdetJ;
    T kappa;
    T rho;
    T q;
    T penalty;
    const Vec<T, dim>& tx;
  };
};

/**
 * @brief Compute the right-hand-side of the adjoint equation for thermal
 * compliance function.
 */
template <typename T, index_t D>
class AdjRHS {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename HeatConduction<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename HeatConduction<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace = typename HeatConduction<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  T heat_source;

  AdjRHS(T heat_source) : heat_source(heat_source) {}

  KOKKOS_FUNCTION void weak(T wdetJ, const DataSpace& data,
                            const FiniteElementGeometry& geo,
                            const FiniteElementSpace& s,
                            FiniteElementSpace& coef) const {
    T& c = coef.template get<0>().get_value();
    c = -wdetJ * heat_source;
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