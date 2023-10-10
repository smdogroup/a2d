#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "a2dcore.h"
#include "multiphysics/febase.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/fespace.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hypercube_basis.h"

namespace A2D {

template <typename T, index_t D,
          GreenStrainType etype = GreenStrainType::LINEAR>
class TopoElasticityIntegrand {
 public:
  TopoElasticityIntegrand(T E, T nu, T q) : q(q) {
    mu0 = 0.5 * E / (1.0 + nu);
    lambda0 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  // Data for the element
  T mu0;      // Second Lame parameter
  T lambda0;  // First Lame parameter
  T q;        // The RAMP penalty parameter

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, dim, H1Space<T, data_dim, dim>>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      FEVarSelect<wrt, DataSpace, FiniteElementGeometry, FiniteElementSpace>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian =
      FESymMatSelect<of, wrt, T, DataSpace::ncomp, FiniteElementGeometry::ncomp,
                     FiniteElementSpace::ncomp>;

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

  /**
   * @brief Compute the Jacobian-vector product
   *
   * @tparam of The residual that we're taking a derivative of
   * @tparam wrt The derivative that we're taking
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param p Direction for Jacobian-vector product
   * @param res Output product
   */
  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian_product(T weight, const DataSpace& data0,
                                        const FiniteElementGeometry& geo0,
                                        const FiniteElementSpace& sref0,
                                        const FiniteElementVar<wrt>& p,
                                        FiniteElementVar<of>& res) const {
    A2DObj<DataSpace> data(data0);
    A2DObj<FiniteElementSpace> sref(sref0);
    A2DObj<FiniteElementGeometry> geo(geo0);

    // Intermediate variables
    A2DObj<T> detJ, penalty, mu, lambda, energy, output;
    A2DObj<FiniteElementSpace> s;
    A2DObj<SymMat<T, dim>> E, S;

    // Set the derivative of the solution
    A2DObj<T&> rho = get_value<0>(data);
    A2DObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

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

    // Compute the Jacobian-vector product
    JacobianProduct<of, wrt>(stack, data, geo, sref, p, res);
  }

  /**
   * @brief Compute the Jacobian at a quadrature point
   *
   * @tparam of The residual that we're taking a derivative of
   * @tparam wrt The derivative that we're taking
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param jac The Jacobian output
   */
  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementJacobian<of, wrt>& jac) const {
    A2DObj<DataSpace> data(data0);
    A2DObj<FiniteElementSpace> sref(sref0);
    A2DObj<FiniteElementGeometry> geo(geo0);

    // Intermediate variables
    A2DObj<T> detJ, penalty, mu, lambda, energy, output;
    A2DObj<FiniteElementSpace> s;
    A2DObj<SymMat<T, dim>> E, S;

    // Set the derivative of the solution
    A2DObj<T&> rho = get_value<0>(data);
    A2DObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

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

    // Extract the Jacobian
    ExtractJacobian<of, wrt>(stack, data, geo, sref, jac);
  }
};

template <class Impl, GreenStrainType etype, index_t degree>
class HexTopoElement
    : public ElementIntegrand<
          Impl, TopoElasticityIntegrand<typename Impl::type, 3, etype>,
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

  HexTopoElement(TopoElasticityIntegrand<T, 3, etype> integrand,
                 std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                 std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                 std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoElasticityIntegrand<T, 3, etype>& get_integrand() {
    return integrand;
  }

 private:
  TopoElasticityIntegrand<T, 3, etype> integrand;
};

template <class Impl, GreenStrainType etype, index_t degree>
class QuadTopoElement
    : public ElementIntegrand<
          Impl, TopoElasticityIntegrand<typename Impl::type, 2, etype>,
          QuadGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>> {
 public:
  using T = typename Impl::type;
  using Integrand = TopoElasticityIntegrand<T, 2, etype>;
  using DataBasis = FEBasis<T, LagrangeH1QuadBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Vec_t = typename Impl::Vec_t;
  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  QuadTopoElement(TopoElasticityIntegrand<T, 2, etype> integrand,
                  std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                  std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoElasticityIntegrand<T, 2, etype>& get_integrand() {
    return integrand;
  }

  void to_vtk(Vec_t& data, Vec_t& geo, Vec_t& sol, const std::string filename) {
    const int num_out = 2;  // Number of outputs to the file
    ElementVector<DataBasis> elem_data(*this->data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*this->geo_mesh, geo);
    ElementVector<Basis> elem_sol(*this->sol_mesh, sol);

    write_quad_to_vtk<num_out, degree, T, DataBasis, GeoBasis, Basis,
                      Integrand>(
        elem_data, elem_geo, elem_sol, filename,
        [](index_t k, typename Integrand::DataSpace& d,
           typename Integrand::FiniteElementGeometry& g,
           typename Integrand::FiniteElementSpace& s) {
          return get_value<0>(s)[k];
        });
  }

 private:
  TopoElasticityIntegrand<T, 2, etype> integrand;
};

/*
  Evalute the KS functional of the stress, given the constitutive class
*/
template <typename T, index_t D,
          GreenStrainType etype = GreenStrainType::LINEAR>
class TopoVonMisesKS {
 public:
  TopoVonMisesKS(T E, T nu, T q, T design_stress, T ks_penalty)
      : q(q), design_stress(design_stress), ks_penalty(ks_penalty) {
    mu = 0.5 * E / (1.0 + nu);
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    max_failure_index = 1.0;
    failure_index_integral = 1.0;
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
  T mu;
  T lambda;

  // The RAMP penalty parameter
  T q;

  // Design stress value used in the constraint
  T design_stress;

  // The KS penalty parameter
  T ks_penalty;

  // Offset value - should be the maximum failure index anywhere in the domain
  T max_failure_index;

  // Integral of e^{ks_penalty * (failure_index - offset)}
  T failure_index_integral;

  /**
   * @brief Set the maximum failure index value (or approximate value)
   *
   * @param max_failure_index_ Maximum value of the failure index anywhere in
   * the domain
   */
  KOKKOS_FUNCTION void set_max_failure_index(T max_failure_index_) {
    max_failure_index = max_failure_index_;
  }

  /**
   * @brief Evaluate the functional value based on the max failure index in the
   * domain and the failure index integral
   *
   * @param failure_index_integral_ Integral of the failure index
   * @return T The failure value
   */
  KOKKOS_FUNCTION T evaluate_functional(T failure_index_integral_) {
    failure_index_integral = failure_index_integral_;
    return max_failure_index + log(failure_index_integral) / ks_penalty;
  }

  /**
   * @brief Compute the failure index at a quadrature point
   *
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  KOKKOS_FUNCTION T max(const DataSpace& data, const FiniteElementGeometry& geo,
                        const FiniteElementSpace& sref) const {
    T rho = data[0];
    FiniteElementSpace s;
    SymMat<T, dim> E, S;
    T detJ, trS, trSS, trS2;

    // Extract the physical coordinate derivatives
    const Mat<T, dim, dim>& Ux = get_grad<0>(s);

    RefElementTransform(geo, sref, detJ, s);  // Transform to physical
    MatGreenStrain<etype>(Ux, E);             // E = E(Ux)
    SymIsotropic(mu, lambda, E, S);           // S = S(mu, lambda, E)
    MatTrace(S, trS);                         // trS = tr(S)
    SymMatMultTrace(S, S, trSS);              // trSS = tr(S * S)
    T vm = sqrt(1.5 * trSS - 0.5 * trS * trS);
    T relaxed_stress = (vm * ((q + 1.0) / (q * rho + 1.0)));
    T failure_index = relaxed_stress / design_stress;

    return failure_index;
  }

  /**
   * @brief Compute the integrand for this functional
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    T rho = data[0];
    FiniteElementSpace s;
    SymMat<T, dim> E, S;
    T detJ, trS, trSS, trS2;

    // Extract the physical coordinate derivatives
    const Mat<T, dim, dim>& Ux = get_grad<0>(s);

    RefElementTransform(geo, sref, detJ, s);  // Transform to physical
    MatGreenStrain<etype>(Ux, E);             // E = E(Ux)
    SymIsotropic(mu, lambda, E, S);           // S = S(mu, lambda, E)
    MatTrace(S, trS);                         // trS = tr(S)
    SymMatMultTrace(S, S, trSS);              // trSS = tr(S * S)
    T vm = sqrt(1.5 * trSS - 0.5 * trS * trS);
    T relaxed_stress = (vm * ((q + 1.0) / (q * rho + 1.0)));
    T failure_index = relaxed_stress / design_stress;

    // Compute the integrand for the KS function
    T output =
        weight * detJ * exp(ks_penalty * (failure_index - max_failure_index));

    return output;
  }

  /**
   * @brief Evaluate the weak form coefficients for linear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data,
                                const FiniteElementGeometry& geo_,
                                const FiniteElementSpace& sref_,
                                FiniteElementVar<wrt>& res) const {
    ADObj<T> rho(data[0]);
    ADObj<FiniteElementSpace> sref(sref_);
    ADObj<FiniteElementGeometry> geo(geo_);

    ADObj<FiniteElementSpace> s;
    ADObj<SymMat<T, dim>> E, S;
    ADObj<T> trS, trSS;
    ADObj<T> output, detJ, vm, relaxed_stress, failure_index;

    // Set the derivative of the solution in the physical coordinates
    ADObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

    auto stack = MakeStack(
        // Compute the strain and stress
        RefElementTransform(geo, sref, detJ, s),  // Transform to physical
        MatGreenStrain<etype>(Ux, E),             // Compute strain
        SymIsotropic(mu, lambda, E, S),           // Compute stress
        MatTrace(S, trS), SymMatMultTrace(S, S, trSS),
        // Evaluate the von Mises stress output
        Eval(sqrt(1.5 * trSS - 0.5 * trS * trS), vm),
        Eval((vm * ((q + 1.0) / (q * rho + 1.0))), relaxed_stress),
        Eval(relaxed_stress / design_stress, failure_index),
        Eval(weight * detJ *
                 exp(ks_penalty * (failure_index - max_failure_index)),
             output));

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
class HexTopoVonMises
    : public IntegralFunctional<
          Impl, TopoVonMisesKS<typename Impl::type, 3, etype>,
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

  HexTopoVonMises(TopoVonMisesKS<T, 3, etype> integrand,
                  std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                  std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoVonMisesKS<T, 3, etype>& get_integrand() { return integrand; }

 private:
  TopoVonMisesKS<T, 3, etype> integrand;
};

template <class Impl, GreenStrainType etype, index_t degree>
class QuadTopoVonMises
    : public IntegralFunctional<
          Impl, TopoVonMisesKS<typename Impl::type, 2, etype>,
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

  QuadTopoVonMises(TopoVonMisesKS<T, 2, etype> integrand,
                   std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                   std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                   std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoVonMisesKS<T, 2, etype>& get_integrand() { return integrand; }

 private:
  TopoVonMisesKS<T, 2, etype> integrand;
};

/*
  Evaluate the volume of the structure, given the constitutive class
*/
template <typename T, index_t D,
          class Integrand = TopoElasticityIntegrand<T, D>>
class TopoVolume {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename Integrand::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry = typename Integrand::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace = typename Integrand::FiniteElementSpace;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar = typename Integrand::template FiniteElementVar<wrt>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian =
      typename Integrand::template FiniteElementJacobian<of, wrt>;

  TopoVolume() = default;

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
    const Mat<T, dim, dim>& J = get_grad<0>(geo);
    T detJ;
    MatDet(J, detJ);
    return weight * detJ * data[0];
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
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data,
                                const FiniteElementGeometry& geo_,
                                const FiniteElementSpace& sref_,
                                FiniteElementVar<wrt>& res) const {
    if constexpr (wrt == FEVarType::DATA) {
      const Mat<T, dim, dim>& J = get_grad<0>(geo_);
      T detJ;
      MatDet(J, detJ);
      res[0] = weight * detJ;
    } else if constexpr (wrt == FEVarType::GEOMETRY) {
      ADObj<FiniteElementGeometry> geo(geo_);
      ADObj<Mat<T, dim, dim>&> J = get_grad<0>(geo);
      ADObj<T> detJ;
      auto stack = MakeStack(MatDet(J, detJ));
      detJ.bvalue() = weight * data[0];
      stack.reverse();
      res.copy(geo.bvalue());
    } else {
      res.zero();
    }
  }
};

template <class Impl, index_t degree>
class QuadTopoVolume
    : public IntegralFunctional<
          Impl, TopoVolume<typename Impl::type, 2>,
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

  QuadTopoVolume(TopoVolume<T, 2> integrand,
                 std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                 std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                 std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoVolume<T, 2>& get_integrand() { return integrand; }

 private:
  TopoVolume<T, 2> integrand;
};

template <typename T, index_t D>
class TopoBodyForceIntegrand {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename TopoElasticityIntegrand<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename TopoElasticityIntegrand<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename TopoElasticityIntegrand<T, D>::FiniteElementSpace;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      typename TopoElasticityIntegrand<T, D>::template FiniteElementVar<wrt>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian = typename TopoElasticityIntegrand<
      T, D>::template FiniteElementJacobian<of, wrt>;

  KOKKOS_FUNCTION TopoBodyForceIntegrand(T q, const T tx[]) : q(q), tx(tx) {}

  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    T rho(data[0]);
    const Vec<T, dim>& U = get_value<0>(sref);
    const Mat<T, dim, dim>& J = get_grad<0>(geo);
    T detJ, dot;
    MatDet(J, detJ);
    VecDot(U, tx, dot);
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);
    T energy = -penalty * weight * detJ * dot;
    return energy;
  }

  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementVar<wrt>& res) const {
    ADObj<DataSpace> data(data0);
    ADObj<FiniteElementSpace> sref(sref0);
    ADObj<FiniteElementGeometry> geo(geo0);

    ADObj<T&> rho = get_value<0>(data);
    ADObj<Vec<T, dim>&> U = get_value<0>(sref);
    ADObj<Mat<T, dim, dim>&> J = get_grad<0>(geo);
    ADObj<T> detJ, dot, penalty, energy;

    auto stack = MakeStack(MatDet(J, detJ), VecDot(U, tx, dot),
                           Eval((q + 1.0) * rho / (q * rho + 1.0), penalty),
                           Eval(-penalty * weight * detJ * dot, energy));

    energy.bvalue() = 1.0;
    stack.reverse();

    if constexpr (wrt == FEVarType::DATA) {
      res[0] = rho.bvalue();
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

    A2DObj<T&> rho = get_value<0>(data);
    A2DObj<Vec<T, dim>&> U = get_value<0>(sref);
    A2DObj<Mat<T, dim, dim>&> J = get_grad<0>(geo);
    A2DObj<T> detJ, dot, penalty, energy;

    auto stack = MakeStack(MatDet(J, detJ), VecDot(U, tx, dot),
                           Eval((q + 1.0) * rho / (q * rho + 1.0), penalty),
                           Eval(-penalty * weight * detJ * dot, energy));
    energy.bvalue() = 1.0;

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

    A2DObj<T&> rho = get_value<0>(data);
    A2DObj<Vec<T, dim>&> U = get_value<0>(sref);
    A2DObj<Mat<T, dim, dim>&> J = get_grad<0>(geo);
    A2DObj<T> detJ, dot, penalty, energy;

    auto stack = MakeStack(MatDet(J, detJ), VecDot(U, tx, dot),
                           Eval((q + 1.0) * rho / (q * rho + 1.0), penalty),
                           Eval(-penalty * weight * detJ * dot, energy));
    energy.bvalue() = 1.0;

    // Extract the Jacobian
    ExtractJacobian<of, wrt>(stack, data, geo, sref, jac);
  }

 private:
  T q;             // RAMP parameter
  Vec<T, dim> tx;  // body force values
};

template <class Impl, index_t degree>
class QuadBodyForceTopoElement
    : public ElementIntegrand<
          Impl, TopoBodyForceIntegrand<typename Impl::type, 2>,
          QuadGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>> {
 public:
  using T = typename Impl::type;
  using Integrand = TopoBodyForceIntegrand<T, 2>;
  using DataBasis = FEBasis<T, LagrangeH1QuadBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Vec_t = typename Impl::Vec_t;
  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  QuadBodyForceTopoElement(TopoBodyForceIntegrand<T, 2> integrand,
                           std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                           std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                           std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const TopoBodyForceIntegrand<T, 2>& get_integrand() { return integrand; }

 private:
  TopoBodyForceIntegrand<T, 2> integrand;
};

/**
 * @brief Apply surface traction and/or surface torque.
 */
template <typename T, index_t D>
class SurfaceTractionIntegrand {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, dim>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim - 1>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, dim, dim - 1>>;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      FEVarSelect<wrt, DataSpace, FiniteElementGeometry, FiniteElementSpace>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian =
      FESymMatSelect<of, wrt, T, DataSpace::ncomp, FiniteElementGeometry::ncomp,
                     FiniteElementSpace::ncomp>;

  // Data used for information
  Vec<T, dim> x0;  // Position vector
  Vec<T, dim> tx;  // Surface traction
  typename std::conditional<dim == 2, T, Vec<T, dim>>::type torx;

  KOKKOS_FUNCTION SurfaceTractionIntegrand(const T tx_[] = nullptr,
                                           const T torx_[] = nullptr,
                                           const T x0_[] = nullptr) {
    if (tx_) {
      for (index_t i = 0; i < dim; i++) {
        tx[i] = tx_[i];
      }
    }

    if (x0_) {
      for (index_t i = 0; i < dim; i++) {
        x0[i] = x0_[i];
      }
    }

    if constexpr (dim == 2) {
      if (torx_) {
        torx = torx_[0];
      } else {
        torx = 0.0;  // must initialize to zero
      }
    } else if constexpr (dim == 3) {
      if (torx_) {
        for (index_t i = 0; i < dim; i++) {
          torx[i] = torx_[i];
        }
      }
    }
  }

  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    FiniteElementSpace s;
    Vec<T, dim> d, tqx, fx;
    T dot, detJ, output;

    Vec<T, dim>& u = get_value<0>(s);
    const Vec<T, dim>& x = get_value<0>(geo);

    BoundaryElementTransform(geo, sref, detJ, s);
    VecSum(T(1.0), x, T(-1.0), x0, d);
    VecCross(d, torx, tqx);
    VecSum(tx, tqx, fx);
    VecDot(u, fx, dot);
    output = -weight * detJ * dot;

    return output;
  }

  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementVar<wrt>& res) const {
    if constexpr (wrt == FEVarType::STATE) {
      ADObj<FiniteElementSpace> sref(sref0), s;
      ADObj<T> dot, output;
      ADObj<Vec<T, dim>&> u = get_value<0>(s);

      T detJ;
      const Vec<T, dim>& x = get_value<0>(geo0);

      Vec<T, dim> d, tqx, fx;
      VecSum(T(1.0), x, T(-1.0), x0, d);
      VecCross(d, torx, tqx);
      VecSum(tx, tqx, fx);

      auto stack =
          MakeStack(BoundaryElementTransform(geo0, sref, detJ, s),
                    VecDot(u, fx, dot), Eval(-weight * detJ * dot, output));

      output.bvalue() = 1.0;
      stack.reverse();

      res.copy(sref.bvalue());
    }
  }

  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian_product(T weight, const DataSpace& data0,
                                        const FiniteElementGeometry& geo0,
                                        const FiniteElementSpace& sref0,
                                        const FiniteElementVar<wrt>& p,
                                        FiniteElementVar<of>& res) const {}

  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian(T weight, const DataSpace& data0,
                                const FiniteElementGeometry& geo0,
                                const FiniteElementSpace& sref0,
                                FiniteElementJacobian<of, wrt>& jac) const {}
};

template <class Impl, index_t degree>
class QuadSurfTraction
    : public ElementIntegrand<
          Impl, SurfaceTractionIntegrand<typename Impl::type, 2>,
          QuadGaussQuadrature<degree + 1>, FEBasis<typename Impl::type>,
          FEBasis<typename Impl::type,
                  LagrangeH1LineBasis<typename Impl::type, 2, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1LineBasis<typename Impl::type, 2, degree>>> {
 public:
  using T = typename Impl::type;
  using Integrand = SurfaceTractionIntegrand<T, 2>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1LineBasis<T, 2, degree>>;
  using Basis = FEBasis<T, LagrangeH1LineBasis<T, 2, degree>>;
  using Vec_t = typename Impl::Vec_t;
  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  QuadSurfTraction(SurfaceTractionIntegrand<T, 2> integrand,
                   std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                   std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                   std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const SurfaceTractionIntegrand<T, 2>& get_integrand() { return integrand; }

 private:
  SurfaceTractionIntegrand<T, 2> integrand;
};

}  // namespace A2D

#endif  // A2D_ELASTICITY_H
