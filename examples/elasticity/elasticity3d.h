#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "constitutive.h"
#include "element.h"
#include "functional.h"
#include "multiarray.h"

namespace A2D {

/*
  Store all the type information about the PDE in one place
*/
template <typename I, typename T>
class ElasticityPDE {
 public:
  static const index_t spatial_dim = 3;
  static const index_t vars_per_node = 3;
  static const index_t null_space_dim = 6;
  static const index_t data_per_point = 2;
  static const index_t dvs_per_point = 1;

  // Layout for the boundary conditions
  typedef A2D::CLayout<2> BCsLayout;
  typedef A2D::MultiArray<I, BCsLayout> BCsArray;

  // Layout for the nodes
  typedef A2D::CLayout<spatial_dim> NodeLayout;
  typedef A2D::MultiArray<T, NodeLayout> NodeArray;

  // Layout for the solution
  typedef A2D::CLayout<vars_per_node> SolutionLayout;
  typedef A2D::MultiArray<T, SolutionLayout> SolutionArray;

  // Layout for the design variables
  typedef A2D::CLayout<dvs_per_point> DesignLayout;
  typedef A2D::MultiArray<T, DesignLayout> DesignArray;

  // Near null space layout - for the AMG preconditioner
  typedef A2D::CLayout<vars_per_node, null_space_dim> NullSpaceLayout;
  typedef A2D::MultiArray<T, NullSpaceLayout> NullSpaceArray;

  // Jacobian matrix
  typedef A2D::BSRMat<I, T, vars_per_node, vars_per_node> SparseMat;

  // Sparse matrix multigrid type
  typedef A2D::BSRMatAmg<I, T, vars_per_node, null_space_dim> SparseAmg;

  static void compute_null_space(NodeArray& X, NullSpaceArray& B) {
    B.zero();
    for (I i = 0; i < B.extent(0); i++) {
      B(i, 0, 0) = 1.0;
      B(i, 1, 1) = 1.0;
      B(i, 2, 2) = 1.0;

      // Rotation about the x-axis
      B(i, 1, 3) = X(i, 2);
      B(i, 2, 3) = -X(i, 1);

      // Rotation about the y-axis
      B(i, 0, 4) = X(i, 2);
      B(i, 2, 4) = -X(i, 0);

      // Rotation about the z-axis
      B(i, 0, 5) = X(i, 1);
      B(i, 1, 5) = -X(i, 0);
    }
  }
};

template <typename I, typename T, class Basis>
class NonlinElasticityElement3D
    : public ElementBasis<I, T, ElasticityPDE<I, T>, Basis> {
 public:
  // Finite-element basis class
  static const index_t NUM_VARS = 3;  // Number of variables per node
  static const index_t NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef ElementBasis<I, T, ElasticityPDE<I, T>, Basis> base;

  NonlinElasticityElement3D(const index_t nelems)
      : ElementBasis<I, T, ElasticityPDE<I, T>, Basis>(nelems) {}

  template <typename IdxType>
  NonlinElasticityElement3D(const index_t nelems, const IdxType conn_[])
      : ElementBasis<I, T, ElasticityPDE<I, T>, Basis>(nelems, conn_) {}

  T energy() {
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    T engry = Basis::template integrate<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0) -> T {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux;
          A2D::SymmMat<T, 3> E;
          T output;

          A2D::Mat3x3MatMult(Uxi0, Jinv0, Ux);
          A2D::Mat3x3GreenStrain(Ux, E);
          A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          return wdetJ * output;
        });

    return engry;
  }

  void add_residual(typename ElasticityPDE<I, T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());
    elem_res.zero();

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    Basis::template residuals<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib) -> void {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          A2D::ADMat<A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::ADMat<A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::ADMat<A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::ADScalar<T> output;

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3GreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();
        },
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  void add_jacobian(typename ElasticityPDE<I, T>::SparseMat& J) {
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());
    elem_jac.zero();

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    Basis::template jacobians<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib,
                A2D::SymmTensor<T, 3, 3>& jac) -> void {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          const int N = 9;
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::A2DMat<N, A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::A2DScalar<N, T> output;

          // Set up the seed values
          for (int k = 0; k < N; k++) {
            A2D::Mat<T, 3, 3>& Up = Uxi.pvalue(k);
            Up(k / 3, k % 3) = 1.0;
          }

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3GreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();

          mult.hforward();
          strain.hforward();
          energy.hreverse();
          strain.hreverse();
          mult.hreverse();

          for (int k = 0; k < N; k++) {
            A2D::Mat<T, 3, 3>& Uxih = Uxi.hvalue(k);
            for (int i = 0; i < 3; i++) {
              for (int j = 0; j < 3; j++) {
                jac(i, j, k / 3, k % 3) = Uxih(i, j);
              }
            }
          }
        },
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  void add_adjoint_dfddata(typename ElasticityPDE<I, T>::SolutionArray& psi,
                           typename base::QuadDataArray& dfdx) {
    // Retrieve the element data
    auto conn = this->get_conn();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();
    auto data = this->get_quad_data();

    // Compute the element adjoint data
    typename base::ElemSolnArray psie(this->get_elem_solution_layout());
    typename base::QuadGradLayout psixi(this->get_quad_gradient_layout());

    VecElementScatter(conn, psi, psie);
    Basis::template gradient<T, NUM_VARS>(psie, psixi);

    // Compute the product
    Basis::template adjoint_product<T, NUM_VARS>(
        detJ, Jinv, Uxi, psixi,
        [&data, &dfdx](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                       A2D::Mat<T, 3, 3> Uxi0,
                       A2D::Mat<T, 3, 3> Psixi0) -> void {
          A2D::Mat<T, 3, 3> Uxib, Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          const int N = 1;
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::A2DMat<N, A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::A2DScalar<N, T> output;
          A2D::A2DScalar<N, T> mu(data(i, j, 0)), lambda(data(i, j, 1));

          // Set the seed values
          A2D::Mat<T, 3, 3>& Psi = Uxi.pvalue(0);
          for (int k2 = 0; k2 < 3; k2++) {
            for (int k1 = 0; k1 < 3; k1++) {
              Psi(k1, k2) = Psixi0(k1, k2);
            }
          }

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3GreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();

          mult.hforward();
          strain.hforward();
          energy.hreverse();

          dfdx(i, j, 0) = mu.hvalue[0];
          dfdx(i, j, 1) = lambda.hvalue[0];
        });
  }
};

template <typename I, typename T, class Basis>
class LinElasticityElement3D
    : public ElementBasis<I, T, ElasticityPDE<I, T>, Basis> {
 public:
  // Finite-element basis class
  static const index_t NUM_VARS = 3;  // Number of variables per node
  static const index_t NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef ElementBasis<I, T, ElasticityPDE<I, T>, Basis> base;

  LinElasticityElement3D(const index_t nelems)
      : ElementBasis<I, T, ElasticityPDE<I, T>, Basis>(nelems) {}

  template <typename IdxType>
  LinElasticityElement3D(const index_t nelems, const IdxType conn_[])
      : ElementBasis<I, T, ElasticityPDE<I, T>, Basis>(nelems, conn_) {}

  T energy() {
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    T elem_energy = Basis::template integrate<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3> Uxi0) -> T {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux;
          A2D::SymmMat<T, 3> E;
          T output;

          A2D::Mat3x3MatMult(Uxi0, Jinv0, Ux);
          A2D::Mat3x3LinearGreenStrain(Ux, E);
          A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          return wdetJ * output;
        });

    return elem_energy;
  }

  void add_residual(typename ElasticityPDE<I, T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());
    elem_res.zero();

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    Basis::template residuals<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib) -> void {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          A2D::ADMat<A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::ADMat<A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::ADMat<A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::ADScalar<T> output;

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();
        },
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  void add_jacobian(typename ElasticityPDE<I, T>::SparseMat& J) {
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());
    elem_jac.zero();

    // Retrieve the element data
    auto data = this->get_quad_data();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();

    Basis::template jacobians<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib,
                A2D::SymmTensor<T, 3, 3>& jac) -> void {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          const int N = 9;
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::A2DMat<N, A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::A2DScalar<N, T> output;

          // Set up the seed values
          for (int k = 0; k < N; k++) {
            A2D::Mat<T, 3, 3>& Up = Uxi.pvalue(k);
            Up(k / 3, k % 3) = 1.0;
          }

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();

          mult.hforward();
          strain.hforward();
          energy.hreverse();
          strain.hreverse();
          mult.hreverse();

          for (int k = 0; k < N; k++) {
            A2D::Mat<T, 3, 3>& Uxih = Uxi.hvalue(k);
            for (int i = 0; i < 3; i++) {
              for (int j = 0; j < 3; j++) {
                jac(i, j, k / 3, k % 3) = Uxih(i, j);
              }
            }
          }
        },
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  void add_adjoint_dfddata(typename ElasticityPDE<I, T>::SolutionArray& psi,
                           typename base::QuadDataArray& dfdx) {
    // Retrieve the element data
    auto conn = this->get_conn();
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uxi = this->get_quad_gradient();
    auto data = this->get_quad_data();

    // Compute the element adjoint data
    typename base::ElemSolnArray psie(this->get_elem_solution_layout());
    typename base::QuadGradArray psixi(this->get_quad_gradient_layout());

    VecElementScatter(conn, psi, psie);
    Basis::template gradient<T, NUM_VARS>(psie, psixi);

    // Compute the product
    Basis::template adjoint_product<T, NUM_VARS>(
        detJ, Jinv, Uxi, psixi,
        [&data, &dfdx](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                       A2D::Mat<T, 3, 3> Uxi0,
                       A2D::Mat<T, 3, 3> Psixi0) -> void {
          A2D::Mat<T, 3, 3> Uxib, Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;

          const int N = 1;
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::A2DMat<N, A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::A2DMat<N, A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::A2DScalar<N, T> output;
          A2D::A2DScalar<N, T> mu(data(i, j, 0)), lambda(data(i, j, 1));

          // Set the seed values
          A2D::Mat<T, 3, 3>& Psi = Uxi.pvalue(0);
          for (int k2 = 0; k2 < 3; k2++) {
            for (int k1 = 0; k1 < 3; k1++) {
              Psi(k1, k2) = Psixi0(k1, k2);
            }
          }

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
          auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

          output.bvalue = wdetJ;

          energy.reverse();
          strain.reverse();
          mult.reverse();

          mult.hforward();
          strain.hforward();
          energy.hreverse();

          dfdx(i, j, 0) = mu.hvalue[0];
          dfdx(i, j, 1) = lambda.hvalue[0];
        });
  }
};

template <typename I, typename T, class Basis>
class RAMPIsoConstitutive : public Constitutive<I, T, ElasticityPDE<I, T>> {
 public:
  static const index_t dvs_per_point = ElasticityPDE<I, T>::dvs_per_point;
  static const index_t nodes_per_elem = Basis::NUM_NODES;
  static const index_t quad_pts_per_elem = Basis::quadrature::NUM_QUAD_PTS;

  typedef A2D::CLayout<nodes_per_elem, dvs_per_point> ElemDesignLayout;
  typedef A2D::CLayout<quad_pts_per_elem, dvs_per_point> QuadDesignLayout;

  typedef A2D::MultiArray<T, ElemDesignLayout> ElemDesignArray;
  typedef A2D::MultiArray<T, QuadDesignLayout> QuadDesignArray;

  RAMPIsoConstitutive(ElementBasis<I, T, ElasticityPDE<I, T>, Basis>& element,
                      T q, T E, T nu)
      : element(element),
        q(q),
        E(E),
        nu(nu),
        elem_design_layout(element.nelems),
        quad_design_layout(element.nelems),
        xe(elem_design_layout),
        xq(quad_design_layout) {
    mu = 0.5 * E / (1.0 + nu);
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  }

  /*
    Set the design variables values into the element object
  */
  void set_design_vars(typename ElasticityPDE<I, T>::DesignArray& x) {
    // Set the design variable values
    // x -> xe -> xq
    auto conn = element.get_conn();
    VecElementScatter(conn, x, xe);
    Basis::template interp<dvs_per_point>(xe, xq);

    auto data = element.get_quad_data();
    for (I i = 0; i < data.extent(0); i++) {
      for (I j = 0; j < data.extent(1); j++) {
        T penalty = xq(i, j, 0) / (1.0 + q * (1.0 - xq(i, j, 0)));
        data(i, j, 0) = mu * penalty;
        data(i, j, 1) = lambda * penalty;
      }
    }
  }

  /*
    Compute the derivative of the adjoint-residual product data w.r.t. x
  */
  void add_adjoint_dfdx(typename ElasticityPDE<I, T>::SolutionArray& psi,
                        typename ElasticityPDE<I, T>::DesignArray& dfdx) {
    typename ElementBasis<I, T, ElasticityPDE<I, T>, Basis>::QuadDataArray
        dfddata(element.get_quad_data_layout());

    dfddata.zero();
    element.add_adjoint_dfddata(psi, dfddata);

    // Set the result into the dfdxq array
    QuadDesignArray dfdxq(quad_design_layout);
    for (I i = 0; i < dfddata.extent(0); i++) {
      for (I j = 0; j < dfddata.extent(1); j++) {
        T denom = (1.0 + q * (1.0 - xq(i, j, 0)));
        T dpenalty = 1.0 / (denom * denom);

        dfdxq(i, j) =
            dpenalty * (mu * dfddata(i, j, 0) + lambda * dfddata(i, j, 1));
      }
    }

    ElemDesignArray dfdxe(elem_design_layout);
    Basis::template interpReverseAdd<dvs_per_point>(dfdxq, dfdxe);

    auto conn = element.get_conn();
    VecElementGatherAdd(conn, dfdxe, dfdx);
  }

 private:
  // Reference to the element class
  ElementBasis<I, T, ElasticityPDE<I, T>, Basis>& element;

  // Design variable views
  ElemDesignLayout elem_design_layout;
  QuadDesignLayout quad_design_layout;
  ElemDesignArray xe;
  QuadDesignArray xq;

  // RAMP penalization
  T q;

  // Constitutive data
  T E, nu;

  // Parameter value
  T mu, lambda;
};

template <typename I, typename T, class Basis>
class StressIntegral3D : public ElementFunctional<I, T, ElasticityPDE<I, T>> {
 public:
  typedef ElementBasis<I, T, ElasticityPDE<I, T>, Basis> base;
  typedef ElasticityPDE<I, T> PDE;

  StressIntegral3D(ElementBasis<I, T, ElasticityPDE<I, T>, Basis>& element,
                   T yield_stress)
      : element(element), yield_stress(yield_stress) {}

  static const int NUM_VARS = 3;

  T eval_functional() {
    auto data = element.get_quad_data();
    auto detJ = element.get_detJ();
    auto Jinv = element.get_Jinv();
    auto Uxi = element.get_quad_gradient();
    T ys = yield_stress;

    T integral = Basis::template integrate<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data, ys](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                    A2D::Mat<T, 3, 3>& Uxi0) -> T {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux;
          A2D::SymmMat<T, 3> E, S;
          T output, trS, trSS;

          A2D::Mat3x3MatMult(Uxi0, Jinv0, Ux);
          A2D::Mat3x3GreenStrain(Ux, E);
          A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
          A2D::Symm3x3Trace(S, trS);
          A2D::Symm3x3SymmMultTrace(S, S, trSS);

          // von Mises = 1.5 * tr(S * S) - 0.5 * tr(S)**2;
          output = 1.5 * trSS - 0.5 * trS * trS;

          return wdetJ * output;
        });

    return integral;
  }

  void add_dfdu(typename PDE::SolutionArray& dfdu) {
    typename base::ElemResArray elem_dfdu(element.get_elem_res_layout());
    elem_dfdu.zero();

    auto data = element.get_quad_data();
    auto detJ = element.get_detJ();
    auto Jinv = element.get_Jinv();
    auto Uxi = element.get_quad_gradient();
    T ys = yield_stress;

    Basis::template residuals<T, NUM_VARS>(
        detJ, Jinv, Uxi,
        [&data](index_t i, index_t j, T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib) -> void {
          T mu(data(i, j, 0)), lambda(data(i, j, 1));
          A2D::Mat<T, 3, 3> Ux0, Uxb;
          A2D::SymmMat<T, 3> E0, Eb;
          A2D::SymmMat<T, 3> S0, Sb;

          A2D::ADMat<A2D::Mat<T, 3, 3>> Uxi(Uxi0, Uxib);
          A2D::ADMat<A2D::Mat<T, 3, 3>> Ux(Ux0, Uxb);
          A2D::ADMat<A2D::SymmMat<T, 3>> E(E0, Eb);
          A2D::ADMat<A2D::SymmMat<T, 3>> S(S0, Sb);
          A2D::ADScalar<T> output, trS, trSS;

          auto mult = A2D::Mat3x3MatMult(Uxi, Jinv0, Ux);
          auto strain = A2D::Mat3x3GreenStrain(Ux, E);
          auto cons = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
          auto trace1 = A2D::Symm3x3Trace(S, trS);
          auto trace2 = A2D::Symm3x3SymmMultTrace(S, S, trSS);

          trSS.bvalue = 1.5 * wdetJ;
          trS.bvalue = -trS.value * wdetJ;

          trace2.reverse();
          trace1.reverse();
          cons.reverse();
          strain.reverse();
          mult.reverse();
        },
        elem_dfdu);

    VecElementGatherAdd(element.get_conn(), elem_dfdu, dfdu);
  }

 private:
  ElementBasis<I, T, ElasticityPDE<I, T>, Basis>& element;
  T yield_stress;
};

}  // namespace A2D

#endif  // ELASTICITY_3D_H
