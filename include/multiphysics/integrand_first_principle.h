#ifndef A2D_POISSON_H
#define A2D_POISSON_H

#include "a2dcore.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fespace.h"

namespace A2D {

/**
 * @brief Calculate First Principle Strain
 *
 * @tparam T Scalar type for the calculation
 * @tparam D Dimension of the problem
 */
template <typename T, index_t D>
class FirstPrinciple {
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

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      FEVarSelect<wrt, DataSpace
      , FiniteElementGeometry, FiniteElementSpace>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian =
      FESymMatSelect<of, wrt, T, DataSpace::ncomp, FiniteElementGeometry::ncomp,
                     FiniteElementSpace::ncomp>;

  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    T detJ, dot, output;
    FiniteElementSpace s;

    // Get the value and the gradient of the solution
    const T& u = get_value<0>(s);
    const Vec<T, dim>& grad = get_grad<0>(s);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    RefElementTransform(geo, sref, detJ, s);
    VecDot(grad, grad, dot);
    output = weight * detJ * (0.5 * dot - u);

    return output;
  }

  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementVar<wrt>& res) const {
    ADObj<FiniteElementSpace> sref(sref0);
    ADObj<FiniteElementGeometry> geo(geo0);

    // Intermediate values
    ADObj<T> detJ, dot, output;
    ADObj<FiniteElementSpace> s;

    // Grab references to the input values
    ADObj<T&> u = get_value<0>(s);
    ADObj<Vec<T, dim>&> grad = get_grad<0>(s);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    auto stack = MakeStack(RefElementTransform(geo, sref, detJ, s),
                           VecDot(grad, grad, dot),
                           Eval(weight * detJ * (0.5 * dot - u), output));

    output.bvalue() = 1.0;
    stack.reverse();

    if constexpr (wrt == FEVarType::DATA) {
    } else if constexpr (wrt == FEVarType::GEOMETRY) {
      res.copy(geo.bvalue());
    } else if constexpr (wrt == FEVarType::STATE) {
      res.copy(sref.bvalue());
    }
  }

  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian_product(T weight, const DataSpace& data0,
                                        const FiniteElementGeometry& geo0,
                                        const FiniteElementSpace& sref0,
                                        const FiniteElementVar<wrt>& p,
                                        FiniteElementVar<of>& res) const {
    A2DObj<DataSpace> data(data0);
    A2DObj<FiniteElementSpace> sref(sref0);
    A2DObj<FiniteElementGeometry> geo(geo0);

    // Intermediate values
    A2DObj<T> detJ, dot, output;
    A2DObj<FiniteElementSpace> s;

    // Grab references to the input values
    A2DObj<T&> u = get_value<0>(s);
    A2DObj<Vec<T, dim>&> grad = get_grad<0>(s);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    auto stack = MakeStack(RefElementTransform(geo, sref, detJ, s),
                           VecDot(grad, grad, dot),
                           Eval(weight * detJ * (0.5 * dot - u), output));

    output.bvalue() = 1.0;

    // Compute the Jacobian-vector product
    JacobianProduct<of, wrt>(stack, data, geo, sref, p, res);
  }

  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementJacobian<of, wrt>& jac) const {
    A2DObj<DataSpace> data(data0);
    A2DObj<FiniteElementSpace> sref(sref0);
    A2DObj<FiniteElementGeometry> geo(geo0);

    // Intermediate values
    A2DObj<T> detJ, dot, output;
    A2DObj<FiniteElementSpace> s;

    // Grab references to the input values
    A2DObj<T&> u = get_value<0>(s);
    A2DObj<Vec<T, dim>&> grad = get_grad<0>(s);

    // Compute wdetJ * (0.5 * || grad ||_{2}^{2} - u)
    auto stack = MakeStack(RefElementTransform(geo, sref, detJ, s),
                           VecDot(grad, grad, dot),
                           Eval(weight * detJ * (0.5 * dot - u), output));

    output.bvalue() = 1.0;

    // Extract the Jacobian
    ExtractJacobian<of, wrt>(stack, data, geo, sref, jac);
  }
};

template <class Impl, index_t degree>
class HexPoissonElement
    : public ElementIntegrand<
          Impl, Poisson<typename Impl::type, 3>, HexGaussQuadrature<degree>,
          FEBasis<typename Impl::type>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 3, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 1, degree>>> {
 public:
  using T = typename Impl::type;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;

  HexPoissonElement(Poisson<T, 3> integrand,
                    std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                    std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                    std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const Poisson<T, 3>& get_integrand() { return integrand; }

 private:
  Poisson<T, 3> integrand;
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

    T dot, output;
    VecDot(q, q, dot);  // dot = ||q||_{2}^{2}
    output = wdetJ * (0.5 * dot + divq * u + u);

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
    ADObj<T> dot, output;

    auto stack = MakeStack(VecDot(q, q, dot),  // dot = ||q||_{2}^{2}
                           Eval(wdetJ * (0.5 * dot + divq * u + u), output));

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
    A2DObj<T> dot, output;

    auto stack = MakeStack(VecDot(q, q, dot),  // dot = ||q||_{2}^{2}
                           Eval(wdetJ * (0.5 * dot + divq * u + u), output));

    output.bvalue() = 1.0;

    // Compute the derivatives first..
    stack.reverse();

    // Create data for extracting the Hessian-vector product
    constexpr index_t ncomp = FiniteElementSpace::ncomp;
    auto inters = MakeTieTuple<T, ADseed::h>(dot);

    // Extract the matrix
    stack.template hextract<T, ncomp, ncomp>(inters, in, out, jac);
  }
};

}  // namespace A2D

#endif  // A2D_POISSON_H
