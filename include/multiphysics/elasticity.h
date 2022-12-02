#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fespace.h"

namespace A2D {

template <typename T, A2D::index_t D>
class TopoLinearElasticity {
 public:
  TopoLinearElasticity(T E, T nu, T q) : q(q) {
    mu0 = 0.5 * E / (1.0 + nu);
    lambda0 = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  // Number of dimensions
  static const A2D::index_t dim = D;

  // Number of data dimensions
  static const A2D::index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = A2D::FESpace<T, data_dim, A2D::L2Space<T, data_dim, dim>>;

  // Space for the element geometry
  using FiniteElementGeometry = A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>>;

  // Finite element space
  using FiniteElementSpace = A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>>;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = A2D::VolumeMapping<T, dim>;

  // The type of matrix used to store data at each quadrature point
  static const A2D::index_t ncomp = FiniteElementSpace::ncomp;
  using QMatType = A2D::SymmMat<T, ncomp>;

  // Data for the element
  T mu0;      // First Lame parameter
  T lambda0;  // Second Lame parameter
  T q;        // The RAMP penalty parameter

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
              const FiniteElementSpace& s) {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Extract the solution
    A2D::Mat<T, dim, dim> Ux = (s.template get<0>()).get_grad();

    // The Green-Langrange strain terms
    A2D::SymmMat<T, dim> E;

    T output;
    A2D::MatLinearGreenStrain(Ux, E);
    A2D::SymmIsotropicEnergy(mu, lambda, E, output);

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
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    T rho = data[0];
    T penalty = 1.0 / (1.0 + q * (1.0 - rho));

    // Get the constitutive data at the points
    T mu = penalty * mu0;
    T lambda = penalty * lambda0;

    // Extract the trial solution gradient and the coefficient terms. Here
    // Uxb is the output computed as the derivative of the strain energy
    // w.r.t. Ux
    A2D::Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    A2D::Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    A2D::ADMat<A2D::Mat<T, dim, dim>> Ux(Ux0, Uxb);

    // The Green-Langrange strain terms
    A2D::SymmMat<T, dim> E0, Eb;
    A2D::ADMat<A2D::SymmMat<T, dim>> E(E0, Eb);

    // The strain energy output
    A2D::ADScalar<T> output;

    auto strain = A2D::MatLinearGreenStrain(Ux, E);
    auto energy = A2D::SymmIsotropicEnergy(mu, lambda, E, output);

    // Seed the output value with the wdetJ
    output.bvalue = wdetJ;

    // Reverse the derivatives through the code
    energy.reverse();
    strain.reverse();
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
    A2D_INLINE_FUNCTION JacVecProduct(const TopoLinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          rho(data[0]),
          penalty(1.0 / (1.0 + pde.q * (1.0 - rho))),
          mu(penalty * pde.mu0),
          lambda(penalty * pde.lambda0),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      Ux.get_hvalue((Jp.template get<0>()).get_grad());
    }

   private:
    T rho, penalty;
    T mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatLinearGreenStrain(Ux, E)) strain;
    decltype(A2D::SymmIsotropicEnergy(mu, lambda, E, output)) energy;
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
    A2D_INLINE_FUNCTION AdjVecProduct(const TopoLinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          rho(data[0]),
          q(pde.q),
          penalty(1.0 / (1.0 + q * (1.0 - rho))),
          mu0(pde.mu0),
          lambda0(pde.lambda0),
          mu(penalty * mu0),
          lambda(penalty * lambda0),

          // Initialize the displacement gradient
          Ux(s.template get<0>().get_grad()),

          // Compute the strain from the displacement gradient
          strain(Ux, E),

          // Compute the strain energy from the strain
          energy(mu, lambda, E, output) {
      // Set the seed on the derivative
      output.bvalue = wdetJ;

      // Reverse mode for the first derivative
      energy.reverse();
      strain.reverse();
    }

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        DataSpace& dfdx) {
      Ux.set_pvalue((p.template get<0>()).get_grad());

      strain.hforward();
      energy.hreverse();
      strain.hreverse();

      T denom = (1.0 + q * (1.0 - rho));
      dfdx[0] +=
          (mu0 * mu.hvalue + lambda0 * lambda.hvalue) * q / (denom * denom);
    }

   private:
    T rho;
    T q;
    T penalty;
    T mu0, lambda0;
    A2D::A2DScalar<T> mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatLinearGreenStrain(Ux, E)) strain;
    decltype(A2D::SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };
};

/*
  Evaluate the volume of the structure, given the constitutive class
*/
template <typename T, A2D::index_t D>
class TopoVolume {
 public:
  // Number of dimensions
  static const A2D::index_t dim = D;

  // Number of data dimensions
  static const A2D::index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename TopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename TopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename TopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = typename TopoLinearElasticity<T, D>::SolutionMapping;

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
              const FiniteElementSpace& s) {
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
                       const FiniteElementSpace& s, DataSpace& dfdx) {
    dfdx.zero();
    dfdx[0] = wdetJ;
  }
};

/*
  Evalute the KS functional of the stress, given the constitutive class
*/
template <typename T, A2D::index_t D>
class TopoVonMisesAggregation {
 public:
  TopoVonMisesAggregation(T E, T nu, T q, T design_stress, T ks_penalty)
      : q(q), design_stress(design_stress), ks_penalty(ks_penalty) {
    mu = 0.5 * E / (1.0 + nu);
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    max_failure_index = 1.0;
    failure_index_integral = 1.0;
  }

  // Number of dimensions
  static const A2D::index_t dim = D;

  // Number of data dimensions
  static const A2D::index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = typename TopoLinearElasticity<T, D>::DataSpace;

  // Space for the element geometry
  using FiniteElementGeometry =
      typename TopoLinearElasticity<T, D>::FiniteElementGeometry;

  // Finite element space
  using FiniteElementSpace =
      typename TopoLinearElasticity<T, D>::FiniteElementSpace;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = typename TopoLinearElasticity<T, D>::SolutionMapping;

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
  void set_max_failure_index(T max_failure_index_) {
    max_failure_index = max_failure_index_;
  }

  /**
   * @brief Evaluate the functional value based on the max failure index in the
   * domain and the failure index integral
   *
   * @param failure_index_integral_ Integral of the failure index
   * @return T The failure value
   */
  T evaluate_functional(T failure_index_integral_) {
    failure_index_integral = failure_index_integral_;
    return max_failure_index + A2D::log(failure_index_integral) / ks_penalty;
  }

  /**
   * @brief Compute the failure index at a quadrature point
   *
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadurature point
   * @return T The integrand contribution
   */
  T max(const DataSpace& data, const FiniteElementGeometry& geo,
        const FiniteElementSpace& s) {
    const A2D::Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    A2D::SymmMat<T, dim> E, S;
    T trS, trSS;

    A2D::MatLinearGreenStrain(Ux, E);
    A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
    A2D::SymmTrace(S, trS);
    A2D::SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * rho/(q * rho + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises^2 = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

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
  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) {
    const A2D::Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    A2D::SymmMat<T, dim> E, S;
    T trS, trSS;

    A2D::MatLinearGreenStrain(Ux, E);
    A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
    A2D::SymmTrace(S, trS);
    A2D::SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * rho/(q * rho + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises^2 = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    return wdetJ * A2D::exp(ks_penalty * (failure_index - max_failure_index));
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
  void weak(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
            const FiniteElementSpace& s, FiniteElementSpace& coef) {
    A2D::Mat<T, dim, dim> Ux0 = (s.template get<0>()).get_grad();
    A2D::Mat<T, dim, dim>& Uxb = (coef.template get<0>()).get_grad();
    A2D::SymmMat<T, dim> E0, Eb;
    A2D::SymmMat<T, dim> S0, Sb;

    A2D::ADMat<A2D::Mat<T, dim, dim>> Ux(Ux0, Uxb);
    A2D::ADMat<A2D::SymmMat<T, dim>> E(E0, Eb);
    A2D::ADMat<A2D::SymmMat<T, dim>> S(S0, Sb);
    A2D::ADScalar<T> trS, trSS;

    auto strain = A2D::MatLinearGreenStrain(Ux, E);
    auto cons = A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
    auto trace1 = A2D::SymmTrace(S, trS);
    auto trace2 = A2D::SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS.value - 0.5 * trS.value * trS.value) /
           design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    // Compute the exponential contribution
    T ks_exp = A2D::exp(ks_penalty * (failure_index - max_failure_index));

    T scale = 0.5 * wdetJ * penalty * ks_exp /
              (vm * design_stress * design_stress * failure_index_integral);

    trSS.bvalue = 1.5 * scale;
    trS.bvalue = -trS.value * scale;

    trace2.reverse();
    trace1.reverse();
    cons.reverse();
    strain.reverse();
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
                       const FiniteElementSpace& s, DataSpace& dfdx) {
    const A2D::Mat<T, dim, dim>& Ux = (s.template get<0>()).get_grad();
    A2D::SymmMat<T, dim> E, S;
    T trS, trSS;

    A2D::MatLinearGreenStrain(Ux, E);
    A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
    A2D::SymmTrace(S, trS);
    A2D::SymmSymmMultTrace(S, S, trSS);

    // Extract the design density value
    T rho = data[0];

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T penalty = (q + 1.0) * rho / (q * rho + 1.0);

    // Compute the penalty = (q + 1) * x/(q * x + 1)
    T denom = (q * rho + 1.0) * (q * rho + 1.0);
    T dpenalty = (q + 1.0) / denom;

    // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
    T vm = std::sqrt(1.5 * trSS - 0.5 * trS * trS) / design_stress;

    // Compute the failure index
    T failure_index = penalty * vm;

    // Compute the exponential contribution
    T ks_exp = A2D::exp(ks_penalty * (failure_index - max_failure_index));

    T scale = wdetJ * vm * ks_exp / (failure_index_integral);

    dfdx[0] = scale * dpenalty;
  }
};

}  // namespace A2D

#endif  // A2D_ELASTICITY_H
