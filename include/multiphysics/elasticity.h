#ifndef A2D_ELASTICITY_H
#define A2D_ELASTICITY_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "multiphysics/fespace.h"

template <typename T, A2D::index_t D>
class NonlinearElasticity {
 public:
  // Number of dimensions
  static const A2D::index_t dim = D;

  // Number of data dimensions
  static const A2D::index_t data_dim = 2;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim, A2D::H1Space<T, data_dim, dim>> DataSpace;

  // Space for the element geometry
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementSpace;

  // The type of matrix used to store data at each quadrature point
  typedef A2D::SymmMat<T, FiniteElementSpace::ncomp> QMatType;

  /**
   * @brief Evaluate the weak form coefficients for nonlinear elasticity
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The trial solution
   * @param coef Output weak form coefficients of the test space
   */
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& data,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    // Get the constitutive data at the points
    T mu = data[0];
    T lambda = data[1];

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

    auto strain = A2D::MatGreenStrain(Ux, E);
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
    A2D_INLINE_FUNCTION JacVecProduct(const NonlinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          mu(data[0]),
          lambda(data[1]),

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
    T mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatGreenStrain(Ux, E)) strain;
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
    A2D_INLINE_FUNCTION AdjVecProduct(const NonlinearElasticity<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        :  // Initialize constitutive data
          mu(data[0]),
          lambda(data[1]),

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

      dfdx.set_value(0, mu.hvalue);
      dfdx.set_value(1, lambda.hvalue);
    }

   private:
    A2D::A2DScalar<T> mu, lambda;
    A2D::A2DMat<A2D::Mat<T, dim, dim>> Ux;
    A2D::A2DMat<A2D::SymmMat<T, dim>> E;
    A2D::A2DScalar<T> output;

    // Declare types of the operators
    decltype(A2D::MatGreenStrain(Ux, E)) strain;
    decltype(A2D::SymmIsotropicEnergy(mu, lambda, E, output)) energy;
  };
};

// template <typename T>
// class TopoIsoConstitutive : public ConstitutiveBase<T> {
//  public:
//   TopoIsoConstitutive(
//       std::shared_ptr<ElementBasis<
//           I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
//           element,
//       T q, T E, T nu, T density, T design_stress, T beta = 20.0,
//       T xoffset = 0.5)
//       : q(q),
//         xoffset(xoffset),
//         beta(beta),
//         E(E),
//         nu(nu),
//         density(density),
//         design_stress(design_stress),
//         element(element) {
//     Timer t("TopoIsoConstitutive::TopoIsoConstitutive");
//     xe = ElemDesignArray("xe", element->nelems);
//     xq = QuadDesignArray("xq", element->nelems);
//     mu = 0.5 * E / (1.0 + nu);
//     lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
//   }

//   // Penalization value
//   const T q;

//   // Heaviside filter approximation
//   const T xoffset;
//   const T beta;

//   // Constitutive data
//   const T E;
//   const T nu;
//   const T density;
//   const T design_stress;

//   // Get the Lame parameters
//   void get_lame_parameters(T& mu_, T& lambda_) {
//     mu_ = mu;
//     lambda_ = lambda;
//   }

//   /*
//     Set the design variables values into the element object
//   */
//   void set_design_vars(
//       typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::DesignArray&
//       x) {
//     // Set the design variable values
//     // x -> xe -> xq
//     auto conn = element->get_conn();
//     VecElementScatter(conn, x, xe);
//     BasisOps::template interp<dvs_per_point>(xe, xq);

//     auto data = element->get_quad_data();
//     for (I i = 0; i < data.extent(0); i++) {
//       for (I j = 0; j < data.extent(1); j++) {
//         T rho_exp = A2D::exp(-beta * (xq(i, j, 0) - xoffset));
//         T rho = 1.0 / (1.0 + rho_exp);
//         T penalty = rho / (1.0 + q * (1.0 - rho));

//         data(i, j, 0) = mu * penalty;
//         data(i, j, 1) = lambda * penalty;
//       }
//     }
//   }

//   /*
//     Compute the derivative of the adjoint-residual product data w.r.t. x
//   */
//   void add_adjoint_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM,
//   I,
//                                                    T>::SolutionArray& psi,
//                         typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM,
//                         I,
//                                                    T>::DesignArray& dfdx) {
//     typename ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//     T>,
//                           BasisOps>::QuadDataArray dfddata("dfddata",
//                                                            element->nelems);

//     // Compute the product of the adjoint with the derivatives of
//     // the residuals w.r.t. the element data
//     element->add_adjoint_dfddata(psi, dfddata);

//     // Set the result into the dfdxq array
//     QuadDesignArray dfdxq("dfdxq", element->nelems);
//     for (I i = 0; i < dfddata.extent(0); i++) {
//       for (I j = 0; j < dfddata.extent(1); j++) {
//         T rho_exp = A2D::exp(-beta * (xq(i, j, 0) - xoffset));
//         T rho = 1.0 / (1.0 + rho_exp);
//         T denom = (1.0 + q * (1.0 - rho));
//         T dpenalty = (q + 1.0) / (denom * denom);
//         dpenalty *= beta * rho_exp * rho * rho;

//         dfdxq(i, j, 0) =
//             dpenalty * (mu * dfddata(i, j, 0) + lambda * dfddata(i, j, 1));
//       }
//     }

//     ElemDesignArray dfdxe("dfdxe", element->nelems);
//     BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

//     auto conn = element->get_conn();
//     VecElementGatherAdd(conn, dfdxe, dfdx);
//   }

//   std::shared_ptr<ElementBasis<
//       I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
//   get_element() {
//     return element;
//   }

//   ElemDesignArray& get_elem_design() { return xe; }
//   QuadDesignArray& get_quad_design() { return xq; }

//  private:
//   // Reference to the element class
//   std::shared_ptr<ElementBasis<
//       I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
//       element;

//   // Design variable views
//   ElemDesignArray xe;
//   QuadDesignArray xq;

//   // Parameter value
//   T mu, lambda;
// };

// /*
//   Evaluate the volume of the structure, given the constitutive class
// */
// template <typename I, typename T, class BasisOps>
// class TopoVolume
//     : public ElementFunctional<I, T,
//                                ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//                                T>> {
//  public:
//   static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
//   static const index_t vars_per_node =
//       ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::vars_per_node;
//   static const index_t dvs_per_point =
//       ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;

//   TopoVolume(std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con)
//       : con(con) {}

//   T eval_functional() {
//     auto element = con->get_element();
//     auto detJ = element->get_detJ();
//     auto xq = con->get_quad_design();

//     T integral = BasisOps::template integrate<T>(
//         detJ, A2D_LAMBDA(index_t i, index_t j, T wdetJ)->T {
//           return xq(i, j, 0) * wdetJ;
//         });

//     return integral;
//   }
//   void add_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//                                            T>::DesignArray& dfdx) {
//     auto element = con->get_element();
//     auto detJ = element->get_detJ();
//     auto xq = con->get_quad_design();

//     typename TopoIsoConstitutive<I, T, BasisOps>::QuadDesignArray dfdxq(
//         "dfdxq", element->nelems);

//     T integral = BasisOps::template integrate<T>(
//         detJ, A2D_LAMBDA(index_t i, index_t j, T wdetJ)->T {
//           dfdxq(i, j, 0) = wdetJ;
//           return xq(i, j, 0) * wdetJ;
//         });

//     typename TopoIsoConstitutive<I, T, BasisOps>::ElemDesignArray dfdxe(
//         "dfdxe", element->nelems);
//     BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

//     auto conn = element->get_conn();
//     VecElementGatherAdd(conn, dfdxe, dfdx);
//   }

//  private:
//   std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con;
// };

// /*
//   Evalute the KS functional of the stress, given the constitutive class
// */
// template <typename I, typename T, class BasisOps>
// class TopoVonMisesAggregation
//     : public ElementFunctional<I, T,
//                                ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//                                T>> {
//  public:
//   static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
//   static const int vars_per_node =
//       ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::vars_per_node;
//   static const index_t dvs_per_point =
//       ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;

//   TopoVonMisesAggregation(
//       std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con,
//       T weight = 100.0)
//       : weight(weight), con(con) {
//     offset = 0.0;
//     integral = 1.0;
//   }

//   // The KS aggregation weight
//   const T weight;

//   /*
//     Reset the maximum value
//   */
//   void compute_offset() {
//     auto element = con->get_element();
//     auto data = element->get_quad_data();
//     auto detJ = element->get_detJ();
//     auto Jinv = element->get_Jinv();
//     auto Uxi = element->get_quad_gradient();
//     auto xq = con->get_quad_design();

//     T ys = con->design_stress;
//     T mu, lambda;
//     con->get_lame_parameters(mu, lambda);
//     T qval = con->q;

//     // Compute the maximum value over all quadrature points
//     offset = BasisOps::template maximum<T, vars_per_node>(
//         detJ, Jinv, Uxi,
//         A2D_LAMBDA(index_t i, index_t j, T wdetJ,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
//             ->T {
//               A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
//               A2D::SymmMat<T, SPATIAL_DIM> E, S;
//               T vm, trS, trSS;

//               A2D::MatMatMult(Uxi0, Jinv0, Ux);
//               A2D::MatGreenStrain(Ux, E);
//               A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
//               A2D::SymmTrace(S, trS);
//               A2D::SymmSymmMultTrace(S, S, trSS);

//               // Compute the penalty = (q + 1) * x/(q * x + 1)
//               T penalty =
//                   (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

//               // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
//               vm = penalty * (1.5 * trSS - 0.5 * trS * trS) / ys;

//               return vm;
//             });
//   }

//   T eval_functional() {
//     auto element = con->get_element();
//     auto data = element->get_quad_data();
//     auto detJ = element->get_detJ();
//     auto Jinv = element->get_Jinv();
//     auto Uxi = element->get_quad_gradient();
//     auto xq = con->get_quad_design();

//     T ys = con->design_stress;
//     T mu, lambda;
//     con->get_lame_parameters(mu, lambda);
//     T qval = con->q;
//     T off = offset;
//     T wgt = weight;

//     // Compute the maximum value over all quadrature points
//     integral = BasisOps::template maximum<T, vars_per_node>(
//         detJ, Jinv, Uxi,
//         A2D_LAMBDA(index_t i, index_t j, T wdetJ,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
//             ->T {
//               A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
//               A2D::SymmMat<T, SPATIAL_DIM> E, S;
//               T vm, trS, trSS;

//               A2D::MatMatMult(Uxi0, Jinv0, Ux);
//               A2D::MatGreenStrain(Ux, E);
//               A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
//               A2D::SymmTrace(S, trS);
//               A2D::SymmSymmMultTrace(S, S, trSS);

//               // Compute the penalty = (q + 1) * x/(q * x + 1)
//               T penalty =
//                   (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

//               // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
//               vm = penalty * (1.5 * trSS - 0.5 * trS * trS) / ys;

//               return wdetJ * A2D::exp(wgt * (vm - off));
//             });

//     return offset + A2D::log(integral) / weight;
//   }

//   void add_dfdu(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//                                            T>::SolutionArray& dfdu) {
//     auto element = con->get_element();
//     typename ElementBasis<I, T, ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//     T>,
//                           BasisOps>::ElemResArray elem_dfdu("elem_dfdu",
//                                                             element->nelems);

//     auto data = element->get_quad_data();
//     auto detJ = element->get_detJ();
//     auto Jinv = element->get_Jinv();
//     auto Uxi = element->get_quad_gradient();
//     auto xq = con->get_quad_design();

//     T ys = con->design_stress;
//     T mu, lambda;
//     con->get_lame_parameters(mu, lambda);
//     T qval = con->q;
//     T off = offset;
//     T wgt = weight;
//     T intgrl = integral;

//     BasisOps::template residuals<T, vars_per_node>(
//         detJ, Jinv, Uxi,
//         A2D_LAMBDA(index_t i, index_t j, T wdetJ,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxib)
//             ->void {
//               A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux0, Uxb;
//               A2D::SymmMat<T, SPATIAL_DIM> E0, Eb;
//               A2D::SymmMat<T, SPATIAL_DIM> S0, Sb;

//               A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>> Uxi(Uxi0,
//               Uxib); A2D::ADMat<A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>>
//               Ux(Ux0, Uxb); A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> E(E0,
//               Eb); A2D::ADMat<A2D::SymmMat<T, SPATIAL_DIM>> S(S0, Sb);
//               A2D::ADScalar<T> trS, trSS;

//               auto mult = A2D::MatMatMult(Uxi, Jinv0, Ux);
//               auto strain = A2D::MatGreenStrain(Ux, E);
//               auto cons = A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
//               auto trace1 = A2D::SymmTrace(S, trS);
//               auto trace2 = A2D::SymmSymmMultTrace(S, S, trSS);

//               // Compute the penalty = (q + 1) * x/(q * x + 1)
//               T penalty =
//                   (qval + 1.0) * xq(i, j, 0) / (qval * xq(i, j, 0) + 1.0);

//               // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
//               T vm = penalty *
//                      (1.5 * trSS.value - 0.5 * trS.value * trS.value) / ys;

//               T scale =
//                   wdetJ * penalty * A2D::exp(wgt * (vm - off)) / (ys *
//                   intgrl);

//               trSS.bvalue = 1.5 * scale;
//               trS.bvalue = -trS.value * scale;

//               trace2.reverse();
//               trace1.reverse();
//               cons.reverse();
//               strain.reverse();
//               mult.reverse();
//             },
//         elem_dfdu);

//     VecElementGatherAdd(element->get_conn(), elem_dfdu, dfdu);
//   }

//   void add_dfdx(typename ElasticityPDEInfo<BasisOps::SPATIAL_DIM, I,
//                                            T>::DesignArray& dfdx) {
//     auto element = con->get_element();
//     auto data = element->get_quad_data();
//     auto detJ = element->get_detJ();
//     auto Jinv = element->get_Jinv();
//     auto Uxi = element->get_quad_gradient();
//     auto xq = con->get_quad_design();

//     T ys = con->design_stress;
//     T mu, lambda;
//     con->get_lame_parameters(mu, lambda);
//     T qval = con->q;
//     T off = offset;
//     T wgt = weight;
//     T intgrl = integral;

//     typename TopoIsoConstitutive<I, T, BasisOps>::QuadDesignArray dfdxq(
//         "dfdxq", element->nelems);

//     BasisOps::template maximum<T, vars_per_node>(
//         detJ, Jinv, Uxi,
//         A2D_LAMBDA(index_t i, index_t j, T wdetJ,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Jinv0,
//                    A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> & Uxi0)
//             ->T {
//               A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Ux;
//               A2D::SymmMat<T, SPATIAL_DIM> E, S;
//               T trS, trSS;

//               A2D::MatMatMult(Uxi0, Jinv0, Ux);
//               A2D::MatGreenStrain(Ux, E);
//               A2D::SymmIsotropicConstitutive(mu, lambda, E, S);
//               A2D::SymmTrace(S, trS);
//               A2D::SymmSymmMultTrace(S, S, trSS);

//               // Compute the penalty = (q + 1) * x/(q * x + 1)
//               T denom = (qval * xq(i, j, 0) + 1.0) * (qval * xq(i, j, 0)
//               + 1.0); T dpenalty = (qval + 1.0) / denom;

//               // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
//               T vm = (1.5 * trSS - 0.5 * trS * trS) / ys;

//               T scale = wdetJ * vm * A2D::exp(wgt * (vm - off)) / intgrl;

//               dfdxq(i, j, 0) += scale * dpenalty;

//               return wdetJ * A2D::exp(wgt * (vm - off));
//             });

//     typename TopoIsoConstitutive<I, T, BasisOps>::ElemDesignArray dfdxe(
//         "dfdxe", element->nelems);
//     BasisOps::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

//     auto conn = element->get_conn();
//     VecElementGatherAdd(conn, dfdxe, dfdx);
//   }

//  private:
//   std::shared_ptr<TopoIsoConstitutive<I, T, BasisOps>> con;
//   T offset;    // Offset value for computing the KS function value
//   T integral;  // Integral: int_{Omega} e^{weight*(vm - offset)} dOmega
// };

// }  // namespace A2D

#endif  // A2D_ELASTICITY_H
