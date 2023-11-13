#ifndef A2D_COMPLIANCE_H
#define A2D_COMPLIANCE_H

#include "a2dcore.h"
#include "multiphysics/febase.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/fespace.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/lagrange_hypercube_basis.h"

namespace A2D {

template <typename T, index_t D,
          GreenStrainType etype = GreenStrainType::LINEAR>
class Compliance {
 public:
  Compliance(T E, T nu, T q) : q(q) {
    // Needs to be different for 2D and 3D
    mu0 = 0.5 * E * (1.0 - nu) / (1.0 - nu * nu);
    lambda0 = E * nu / (1.0 - nu * nu);
  }

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename TopoElasticityIntegrand<T, D, etype>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename TopoElasticityIntegrand<T, D, etype>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename TopoElasticityIntegrand<T, D, etype>::FiniteElementSpace;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      typename TopoElasticityIntegrand<T, D,
                                       etype>::template FiniteElementVar<wrt>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian = typename TopoElasticityIntegrand<
      T, D, etype>::template FiniteElementJacobian<of, wrt>;

  // Material parameters
  T mu0;
  T lambda0;

  // The RAMP penalty parameter
  T q;

  /**
   * @brief Find the integral of the compliance over the entire domain
   *
   * @param weight The quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param sref The solution at the quadurature point
   * @return T The integrand contribution
   */
  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    // Input values
    T rho = data[0];

    // Intermediate variables
    FiniteElementSpace s;
    SymMat<T, dim> E, S;
    T detJ, energy;

    // Extract the physical coordinate derivatives
    const Mat<T, dim, dim>& Ux = get_grad<0>(s);

    RefElementTransform(geo, sref, detJ, s);
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;
    MatGreenStrain<etype>(Ux, E);
    SymIsotropic(mu, lambda, E, S);
    SymMatMultTrace(E, S, energy);
    T output = 0.5 * weight * detJ * energy;

    return output;
  }

  /**
   * @brief Compute the contribution to the residual
   *
   * @tparam wrt Variable type (DATA, GEOMETRY, STATE)
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param res Residual contribution
   */
  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementVar<wrt>& res) const {
    ADObj<DataSpace> data(data0);
    ADObj<FiniteElementSpace> sref(sref0);
    ADObj<FiniteElementGeometry> geo(geo0);

    // Intermediate variables
    ADObj<T> detJ, penalty, mu, lambda, energy, output;
    ADObj<FiniteElementSpace> s;
    ADObj<SymMat<T, dim>> E, S;

    // Set the derivative of the solution
    ADObj<T&> rho = get_value<0>(data);
    ADObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

    // Make a stack of the operations
    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),       // transform
        Eval(1.0 / (1.0 + q * (1.0 - rho)), penalty),  // penalty parameter
        Eval(penalty * mu0, mu), Eval(penalty * lambda0, lambda),
        MatGreenStrain<etype>(Ux, E),
        SymIsotropic(mu, lambda, E, S),               // Evaluate the stress
        SymMatMultTrace(E, S, energy),                // Compute the energy
        Eval(0.5 * weight * detJ * energy, output));  // Compute the output

    output.bvalue() = 1.0;
    stack.reverse();

    if constexpr (wrt == FEVarType::DATA) {
      res[0] = rho.bvalue();
    } else if constexpr (wrt == FEVarType::GEOMETRY) {
      res.copy(geo.bvalue());
    } else if constexpr (wrt == FEVarType::STATE) {
      res.copy(sref.bvalue());
    }
  }
};

template <class Impl, GreenStrainType etype, index_t degree>
class HexTopoCompliance
    : public IntegralFunctional<
          Impl, Compliance<typename Impl::type, 3, etype>,
          HexGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 3, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 3, degree>>> {
 public:
  using T = typename Impl::type;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;

  HexTopoCompliance(Compliance<T, 3, etype> integrand,
                    std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                    std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                    std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const Compliance<T, 3, etype>& get_integrand() { return integrand; }

 private:
  Compliance<T, 3, etype> integrand;
};

template <class Impl, GreenStrainType etype, index_t degree>
class QuadTopoCompliance
    : public IntegralFunctional<
          Impl, Compliance<typename Impl::type, 2, etype>,
          QuadGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>> {
 public:
  using T = typename Impl::type;
  using DataBasis = FEBasis<T, LagrangeH1QuadBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;

  QuadTopoCompliance(Compliance<T, 2, etype> integrand,
                     std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                     std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                     std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const Compliance<T, 2, etype>& get_integrand() { return integrand; }

 private:
  Compliance<T, 2, etype> integrand;
};

}  // namespace A2D
#endif  // A2D_COMPLIANCE_H
