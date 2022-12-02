#ifndef A2D_HEAT_CONDUCTION_H
#define A2D_HEAT_CONDUCTION_H

namespace A2D {

/*
  Heat conduction problem
*/
template <typename T, index_t D>
class HeatConduction {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;

  // No data associated with this element
  static const A2D::index_t data_dim = 1;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim, A2D::H1Space<T, data_dim, dim>> DataSpace;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, 1, dim>> FiniteElementSpace;

  // Space for the element geometry - parametrized by H1
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef A2D::SymmMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = A2D::VolumeMapping<T, dim>;

  /**
   * @brief Evaluate the weak form of the coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    // Field objects for solution functions
    const A2D::Vec<T, dim>& tx = s.template get<0>().get_grad();
    A2D::Vec<T, dim>& cx = coef.template get<0>().get_grad();

    // Set the thermal conductivity coefficient
    T kappa = data[0];

    for (index_t k = 0; k < dim; k++) {
      cx(k) = wdetJ * kappa * tx(k);
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
    A2D_INLINE_FUNCTION JacVecProduct(const HeatConduction<T, D>& pde, T wdetJ,
                                      const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        : wdetJ(wdetJ), kappa(data[0]) {}

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const A2D::Vec<T, dim>& tx = p.template get<0>().get_grad();
      A2D::Vec<T, dim>& cx = Jp.template get<0>().get_grad();

      for (index_t k = 0; k < dim; k++) {
        cx(k) = wdetJ * kappa * tx(k);
      }
    }

   private:
    T wdetJ;
    T kappa;
  };
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
template <typename T, A2D::index_t D>
class MixedHeatConduction {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;

  // No data associated with this element
  static const A2D::index_t data_dim = 1;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim, A2D::H1Space<T, data_dim, dim>> DataSpace;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::HdivSpace<T, dim>, A2D::L2Space<T, 1, dim>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef A2D::SymmMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = A2D::VolumeMapping<T, dim>;

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
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    // Field objects for solution functions
    const A2D::Vec<T, dim>& q = s.template get<0>().get_value();
    const T& div = s.template get<0>().get_div();
    const T& temp = s.template get<1>().get_value();

    // Test function values
    A2D::Vec<T, dim>& qb = coef.template get<0>().get_value();
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
    A2D_INLINE_FUNCTION JacVecProduct(const MixedHeatConduction<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        : wdetJ(wdetJ), kappa(data[0]) {}

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const A2D::Vec<T, dim>& q = p.template get<0>().get_value();
      const T& div = p.template get<0>().get_div();
      const T& temp = p.template get<1>().get_value();

      // Test function values
      A2D::Vec<T, dim>& qb = Jp.template get<0>().get_value();
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