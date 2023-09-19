#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "a2dcore.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fespace.h"

namespace A2D {

template <typename T, index_t D>
class TopoElasticity {
 public:
  TopoElasticity(T E, T nu, T q) : q(q) {
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
  using DataSpace = FESpace<T, data_dim, L2Space<T, data_dim, dim>>;

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
  T integrand(T weight, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& sref) const {
    // Input values
    T rho = data[0];

    // Intermediate variables
    FiniteElementSpace s;
    SymMat<T, dim> E, S;
    T detJ, energy;

    // Extract the pointer to the physical coordinate derivatives
    const Mat<T, dim, dim>& Ux = get_grad<0>(s);

    // Get the constitutive data at the points
    RefElementTransform(geo, sref, detJ, s);
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;
    MatGreenStrain<GreenStrain::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);
    SymMatMultTrace(E, S, energy);
    T output = 0.5 * weight * detJ * energy;

    return output;
  }

  template <FEVarType wrt>
  void residual(T weight, const DataSpace& data,
                const FiniteElementGeometry& geo_,
                const FiniteElementSpace& sref_,
                FiniteElementVar<wrt>& res) const {
    ADObj<T> rho(data[0]);
    ADObj<FiniteElementSpace> sref(sref_);
    ADObj<FiniteElementGeometry> geo(geo_);

    // Intermediate variables
    ADObj<T> detJ, penalty, mu, lambda, energy, output;
    ADObj<FiniteElementSpace> s;
    ADObj<SymMat<T, dim>> E, S;

    // Set the derivative of the solution
    ADObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

    // Make a stack of the operations
    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),       // transform
        Eval(1.0 / (1.0 + q * (1.0 - rho)), penalty),  // penalty parameter
        Eval(penalty * mu0, mu), Eval(penalty * lambda0, lambda),
        MatGreenStrain<GreenStrain::LINEAR>(Ux, E),
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

  template <FEVarType of, FEVarType wrt>
  void jacobian_product(T weight, const DataSpace& data,
                        const FiniteElementGeometry& geo_,
                        const FiniteElementSpace& sref_,
                        const FiniteElementVar<wrt>& p,
                        FiniteElementVar<of>& res) const {
    A2DObj<T> rho(data[0]);
    A2DObj<FiniteElementSpace> sref(sref_);
    A2DObj<FiniteElementGeometry> geo(geo_);

    // Intermediate variables
    A2DObj<T> detJ, penalty, mu, lambda, energy, output;
    A2DObj<FiniteElementSpace> s;
    A2DObj<SymMat<T, dim>> E, S;

    // Set the derivative of the solution
    A2DObj<Mat<T, dim, dim>&> Ux = get_grad<0>(s);

    if constexpr (wrt == FEVarType::DATA) {
      rho.pvalue() = p[0];
    } else if constexpr (wrt == FEVarType::GEOMETRY) {
      geo.pvalue().copy(p);
    } else if constexpr (wrt == FEVarType::STATE) {
      sref.pvalue().copy(p);
    }

    // Make a stack of the operations
    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),       // transform
        Eval(1.0 / (1.0 + q * (1.0 - rho)), penalty),  // penalty parameter
        Eval(penalty * mu0, mu), Eval(penalty * lambda0, lambda),
        MatGreenStrain<GreenStrain::LINEAR>(Ux, E),
        SymIsotropic(mu, lambda, E, S),               // Evaluate the stress
        SymMatMultTrace(E, S, energy),                // Compute the energy
        Eval(0.5 * weight * detJ * energy, output));  // Compute the output

    output.bvalue() = 1.0;
    stack.reverse();
    stack.hforward();
    stack.hreverse();

    if constexpr (of == FEVarType::DATA) {
      res[0] = rho.hvalue();
    } else if constexpr (of == FEVarType::GEOMETRY) {
      res.copy(geo.hvalue());
    } else if constexpr (of == FEVarType::STATE) {
      res.copy(sref.hvalue());
    }
  }

  template <FEVarType of, FEVarType wrt>
  void jacobian(T weight, const DataSpace& data,
                const FiniteElementGeometry& geo,
                const FiniteElementSpace& sref,
                FiniteElementJacobian<of, wrt>& jac) const {}
};

template <typename T, index_t D>
class IntegrandTopoLinearElasticity {
 public:
  IntegrandTopoLinearElasticity(T E, T nu, T q) : q(q) {
    mu0 = 0.5 * E / (1.0 + nu);
    lambda0 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, data_dim, L2Space<T, data_dim, dim>>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  // The type of matrix used to store data at each quadrature point
  static const index_t ncomp = FiniteElementSpace::ncomp;

  // Finite-element matrix at a quadrature point
  using QMatType = SymMat<T, ncomp>;

  // Data for the element
  T mu0;      // Second Lame parameter
  T lambda0;  // First Lame parameter
  T q;        // The RAMP penalty parameter

  /**
   * @brief Find the integral of the compliance over the entire domain
   *
   * @param weight The integration weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    // Input values
    T rho = data[0];
    Mat<T, dim, dim> Ux = (s.template get<0>()).get_grad();

    // The output scalar
    T output;

    // The Green-Langrange strain terms
    SymMat<T, dim> E, S;

    // Get the constitutive data at the points
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;
    MatGreenStrain<GreenStrain::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);
    SymMatMultTrace(E, S, output);

    return 0.5 * wdetJ * output;
  }

  /**
   * @brief Evaluate the weak form coefficients for the residual
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Create the values needed in the computation
    ADObj<Mat<T, dim, dim>> Ux((s.template get<0>()).get_grad());
    ADObj<SymMat<T, dim>> E, S;

    // The strain energy output
    ADObj<T> output;

    auto stack =
        MakeStack(MatGreenStrain<GreenStrain::LINEAR>(Ux, E),  // E = E(Ux)
                  SymIsotropic(mu, lambda, E, S),              // S = S(E)
                  SymMatMultTrace(E, S, output));  // output = tr(E * S)

    // Seed the output value with the wdetJ
    output.bvalue() = 0.5 * wdetJ;

    // Reverse the derivatives through the code
    stack.reverse();

    // Set the derivative output value
    Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    Uxb.copy(Ux.bvalue());
  }

  // Evaluate the second order derivatives of the integral
  KOKKOS_FUNCTION void jacobian(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                QMatType& jac) const {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    FiniteElementSpace in, out;

    // Set up the pointers into the input/output matrices
    Mat<T, dim, dim> Ux0 = s.template get<0>().get_grad();
    Mat<T, dim, dim> Uxb;
    Mat<T, dim, dim>& Uxp = in.template get<0>().get_grad();
    Mat<T, dim, dim>& Uxh = out.template get<0>().get_grad();

    // Extract displacement gradient
    A2DObj<Mat<T, dim, dim>&> Ux(Ux0, Uxb, Uxp, Uxh);

    // The Green-Lagrange strain terms
    A2DObj<SymMat<T, dim>> E, S;

    // The strain energy output
    A2DObj<T> output;

    auto stack =
        MakeStack(MatGreenStrain<GreenStrain::LINEAR>(Ux, E),  // E = E(Ux)
                  SymIsotropic(mu, lambda, E, S),              // S = S(E)
                  SymMatMultTrace(E, S, output));  // output = tr(E * S)

    // Seed the output value with the wdetJ
    output.bvalue() = 0.5 * wdetJ;

    // Reverse the derivatives through the code
    stack.reverse();

    // Create data for extracting the Hessian-vector product
    constexpr index_t ncomp = FiniteElementSpace::ncomp;
    auto inters = MakeTieTuple<T, ADseed::h>(S, E);

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
    A2DObj<T> rho(data[0]);

    // Extract displacement gradient
    A2DObj<Mat<T, dim, dim>> Ux(s.template get<0>().get_grad());
    Ux.pvalue().copy(adj.template get<0>().get_grad());

    // Set the displacement
    A2DObj<T> penalty, mu, lambda;

    // The Green-Lagrange strain terms
    A2DObj<SymMat<T, dim>> E, S;

    // The strain energy output
    A2DObj<T> output;

    auto stack =
        MakeStack(Eval(1.0 / (1.0 + q * (1.0 - rho)), penalty),
                  Eval(mu0 * penalty, mu), Eval(lambda0 * penalty, lambda),
                  MatGreenStrain<GreenStrain::LINEAR>(Ux, E),  // E = E(Ux)
                  SymIsotropic(mu, lambda, E, S),              // S = S(E)
                  SymMatMultTrace(E, S, output));  // output = tr(E * S)

    // Seed the output value with the wdetJ
    output.bvalue() = 0.5 * wdetJ;

    // Reverse the derivatives through the code
    stack.reverse();

    // Forward/reverse through the code
    stack.hforward();
    stack.hreverse();

    // Set the derivative
    dfdx[0] = rho.hvalue();
  }
};

/*
  Evaluate the volume of the structure, given the constitutive class
*/
template <typename T, index_t C, index_t D, class Integrand>
class IntegrandTopoVolume {
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

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = typename Integrand::SolutionMapping;

  IntegrandTopoVolume() = default;

  /**
   * @brief Compute the integrand for this functional
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) const {
    return wdetJ * data[0];
  }

  /**
   * @brief Derivative of the integrand with respect to the data
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @param dfdx The output derivative value
   */
  void data_derivative(T wdetJ, const DataSpace& data,
                       const FiniteElementGeometry& geo,
                       const FiniteElementSpace& s, DataSpace& dfdx) const {
    dfdx.zero();
    dfdx[0] = wdetJ;
  }
};

template <typename T, index_t D>
class IntegrandTopoBodyForce {
 public:
  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename IntegrandTopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = InteriorMapping<T, dim>;

  KOKKOS_FUNCTION IntegrandTopoBodyForce(T q, const T tx_[]) : q(q) {
    for (index_t i = 0; i < dim; i++) {
      tx[i] = tx_[i];
    }
  }

  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    T rho = data[0];
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // Add body force components
    Vec<T, dim>& Ub = (coef.template get<0>()).get_value();
    for (index_t i = 0; i < dim; i++) {
      Ub(i) = wdetJ * penalty * tx[i];
    }
  }

  KOKKOS_FUNCTION void data_adjoint_product(T wdetJ, const DataSpace& data,
                                            const FiniteElementGeometry& geo,
                                            const FiniteElementSpace& s,
                                            const FiniteElementSpace& adj,
                                            DataSpace& dfdx) const {
    const Vec<T, dim>& Uadj = (adj.template get<0>()).get_value();

    T rho = data[0];
    T dpdrho = (q + 1.0) / ((q * rho + 1.0) * (q * rho + 1.0));

    for (index_t i = 0; i < dim; i++) {
      dfdx[0] += wdetJ * dpdrho * Uadj(i) * tx[i];
    }
  }

 private:
  T q;        // RAMP parameter
  T tx[dim];  // body force values
};

/*
  Evalute the KS functional of the stress, given the constitutive class
*/
template <typename T, index_t D>
class IntegrandTopoVonMisesKS {
 public:
  IntegrandTopoVonMisesKS(T E, T nu, T q, T design_stress, T ks_penalty)
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
  using DataSpace = typename IntegrandTopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename IntegrandTopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping =
      typename IntegrandTopoLinearElasticity<T, D>::SolutionMapping;

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
                        const FiniteElementSpace& s) const {
    // Inputs
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    T rho = data[0];

    // Intermediaries
    SymMat<T, dim> E, S;
    T trS, trSS, trS2;

    // Compute the strain and stress
    MatGreenStrain<GreenStrain::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);

    // Compute the von Mises stress = sqrt(1.5 * tr(S * S) - 0.5 * tr(S)**2)
    MatTrace(S, trS);
    SymMatMultTrace(S, S, trSS);
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
  KOKKOS_FUNCTION T integrand(T wdetJ, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& s) const {
    // Inputs
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    T rho = data[0];

    // Intermediaries
    SymMat<T, dim> E, S;
    T trS, trSS, trS2;

    // Compute the strain and stress
    MatGreenStrain<GreenStrain::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);

    // Compute the von Mises stress = sqrt(1.5 * tr(S * S) - 0.5 * tr(S)**2)
    MatTrace(S, trS);
    SymMatMultTrace(S, S, trSS);
    T vm = sqrt(1.5 * trSS - 0.5 * trS * trS);
    T relaxed_stress = (vm * ((q + 1.0) / (q * rho + 1.0)));
    T failure_index = relaxed_stress / design_stress;

    // Compute the integrand for the KS function
    T output = exp(ks_penalty * (failure_index - max_failure_index));

    return wdetJ * output;
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
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();

    // Inputs/outputs
    ADObj<Mat<T, dim, dim>&> Ux(Ux0, Uxb);
    T rho(data[0]);

    // Intermediaries
    ADObj<SymMat<T, dim>> E, S;
    ADObj<T> trS, trSS;
    ADObj<T> vm, relaxed_stress, failure_index;

    auto stack = MakeStack(
        // Compute the strain and stress
        MatGreenStrain<GreenStrain::LINEAR>(Ux, E),
        SymIsotropic(mu, lambda, E, S), MatTrace(S, trS),
        SymMatMultTrace(S, S, trSS),

        // Evaluate the von Mises stress output
        Eval(sqrt(1.5 * trSS - 0.5 * trS * trS), vm),
        Eval((vm * ((q + 1.0) / (q * rho + 1.0))), relaxed_stress),
        Eval(relaxed_stress / design_stress, failure_index));

    // Compute the integrand for the KS function
    T ks_factor =
        exp(ks_penalty * (failure_index.value() - max_failure_index)) /
        failure_index_integral;

    failure_index.bvalue() = wdetJ * ks_factor;

    stack.reverse();
  }

  // A2D::operator*<double, A2D::ADObj<std::complex<double>>>

  /**
   * @brief Derivative of the integrand with respect to the data
   *
   * @param wdetJ The determinant of the Jacobian times the quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @param dfdx The output derivative value
   */
  KOKKOS_FUNCTION void data_derivative(T wdetJ, const DataSpace& data,
                                       const FiniteElementGeometry& geo,
                                       const FiniteElementSpace& s,
                                       DataSpace& dfdx) const {
    const Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();

    SymMat<T, dim> E, S;
    T trS, trSS;

    // Compute the strain and stress
    MatGreenStrain<GreenStrain::LINEAR>(Ux, E);
    SymIsotropic(mu, lambda, E, S);
    MatTrace(S, trS);
    SymMatMultTrace(S, S, trSS);

    // Evaluate the von Mises stress output
    T vm = sqrt(1.5 * trSS - 0.5 * trS * trS);

    // Inputs
    ADObj<T> rho(data[0]);

    // Intermediaries
    ADObj<T> relaxed_stress, failure_index;

    auto stack = MakeStack(
        Eval((vm * ((T(q + 1.0)) / (q * rho + T(1.0)))), relaxed_stress),
        Eval(relaxed_stress / design_stress, failure_index));

    // Compute the integrand for the KS function
    T ks_factor =
        exp(ks_penalty * (failure_index.value() - max_failure_index)) /
        failure_index_integral;

    failure_index.bvalue() = wdetJ * ks_factor;

    stack.reverse();

    dfdx[0] = rho.bvalue();
  }
};

/**
 * @brief Apply surface traction and/or surface torque.
 */
template <typename T, index_t D>
class IntegrandTopoSurfaceTraction {
 public:
  KOKKOS_FUNCTION IntegrandTopoSurfaceTraction(const T tx_[] = nullptr,
                                               const T torx_[] = nullptr,
                                               const T x0_[] = nullptr) {
    has_traction = false;
    has_torque = false;

    if (tx_) {
      for (index_t i = 0; i < dim; i++) {
        tx[i] = tx_[i];
      }
      has_traction = true;
    }

    if (torx_ && x0_) {
      for (index_t i = 0; i < dim; i++) {
        x0[i] = x0_[i];
        if constexpr (dim == 3) {
          torx[i] = torx_[i];
        }
      }
      if constexpr (dim == 2) {
        torx[0] = torx_[0];
      }
      has_torque = true;
    }
  }

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

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = SurfaceMapping<T, dim>;

  T tx[dim];  // surface traction vector
  T torx[conditional_value<index_t, dim == 3, 3, 1>::value];  // surface torque
                                                              // vector
  T x0[dim];                                                  // torque origin
  bool has_traction;
  bool has_torque;

  /**
   * @brief Evaluate the weak form coefficients for linear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  KOKKOS_FUNCTION void residual(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) const {
    // Extract the solution
    Vec<T, dim>& U = (coef.template get<0>()).get_value();
    for (index_t i = 0; i < dim; i++) {
      U(i) = 0.0;
    }

    if (has_traction) {
      for (index_t i = 0; i < dim; i++) {
        U(i) -= wdetJ * tx[i];
      }
    }

    if (has_torque) {
      // Extract location
      const Vec<T, dim>& x = (geo.template get<0>()).get_value();

      if constexpr (dim == 2) {
        // Force at this point is (x - x0) cross torque
        U(0) += wdetJ * torx[0] * (x(1) - x0[1]);
        U(1) += -wdetJ * torx[0] * (x(0) - x0[0]);
      } else {  // dim == 3
        // Force at this point is (x - x0) cross torque
        U(0) += wdetJ * ((x(1) - x0[1]) * torx[2] - (x(2) - x0[2]) * torx[1]);
        U(1) += wdetJ * ((x(2) - x0[2]) * torx[0] - (x(0) - x0[0]) * torx[2]);
        U(2) += wdetJ * ((x(0) - x0[0]) * torx[1] - (x(1) - x0[1]) * torx[0]);
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_ELASTICITY_H
