#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "model.h"
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

  // Layout for the boundary conditions
  typedef A2D::CLayout<2> BCsLayout;
  typedef A2D::MultiArray<I, BCsLayout> BCsArray;

  // Layout for the nodes
  typedef A2D::CLayout<spatial_dim> NodeLayout;
  typedef A2D::MultiArray<T, NodeLayout> NodeArray;

  // Layout for the solution
  typedef A2D::CLayout<vars_per_node> SolutionLayout;
  typedef A2D::MultiArray<T, SolutionLayout> SolutionArray;

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

  T energy() {
    T engry = 0.0;
    Basis::template energy<T, NonlinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), engry);
    return engry;
  }

  void add_residual(typename ElasticityPDE<I, T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());
    elem_res.zero();

    Basis::template residuals<T, NonlinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), elem_res);

    element_gather_add(this->get_conn(), elem_res, res);
  }

  void add_jacobian(typename ElasticityPDE<I, T>::SparseMat& J) {
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());
    elem_jac.zero();

    Basis::template jacobians<T, NonlinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  class Impl {
   public:
    static const index_t NUM_VARS = 3;

    template <class QuadPointData>
    static T compute_energy(I i, I j, QuadPointData& data, T wdetJ,
                            A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux;
      SymmMat3x3 E;
      T output;

      A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      A2D::Mat3x3GreenStrain(Ux, E);
      A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

      return output * wdetJ;
    }

    template <class QuadPointData>
    static void compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 3, 3>& Uxi0,
                                 A2D::Mat<T, 3, 3>& Uxib) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;

      A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
      A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
      A2D::ADMat<SymmMat3x3> E(E0, Eb);
      A2D::ADScalar<T> output;

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3GreenStrain(Ux, E);
      auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

      output.bvalue = wdetJ;

      energy.reverse();
      strain.reverse();
      mult.reverse();
    }

    template <class QuadPointData>
    static void compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 3, 3>& Uxi0,
                                 A2D::Mat<T, 3, 3>& Uxib,
                                 A2D::SymmTensor<T, 3, 3>& jac) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;

      const int N = 9;
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;

      // Set up the seed values
      for (int k = 0; k < N; k++) {
        Mat3x3& Up = Uxi.pvalue(k);
        Up(k / 3, k % 3) = 1.0;
      }

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
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
        Mat3x3& Uxih = Uxi.hvalue(k);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            jac(i, j, k / 3, k % 3) = Uxih(i, j);
          }
        }
      }
    }
  };
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

  T energy() {
    T engry = 0.0;
    Basis::template energy<T, LinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), engry);
    return engry;
  }

  void add_residual(typename ElasticityPDE<I, T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());
    elem_res.zero();

    Basis::template residuals<T, LinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), elem_res);

    element_gather_add(this->get_conn(), elem_res, res);
  }

  void add_jacobian(typename ElasticityPDE<I, T>::SparseMat& J) {
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());
    elem_jac.zero();

    Basis::template jacobians<T, LinElasticityElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  class Impl {
   public:
    static const index_t NUM_VARS = 3;

    template <class QuadPointData>
    static T compute_energy(I i, I j, QuadPointData& data, T wdetJ,
                            A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux;
      SymmMat3x3 E;
      T output;

      A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      A2D::Mat3x3LinearGreenStrain(Ux, E);
      A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

      return output * wdetJ;
    }

    template <class QuadPointData>
    static void compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 3, 3>& Uxi0,
                                 A2D::Mat<T, 3, 3>& Uxib) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;

      A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
      A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
      A2D::ADMat<SymmMat3x3> E(E0, Eb);
      A2D::ADScalar<T> output;

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
      auto energy = A2D::Symm3x3IsotropicEnergy(mu, lambda, E, output);

      output.bvalue = wdetJ;

      energy.reverse();
      strain.reverse();
      mult.reverse();
    }

    template <class QuadPointData>
    static void compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 3, 3>& Uxi0,
                                 A2D::Mat<T, 3, 3>& Uxib,
                                 A2D::SymmTensor<T, 3, 3>& jac) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;

      const int N = 9;
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;

      // Set up the seed values
      for (int k = 0; k < N; k++) {
        Mat3x3& Up = Uxi.pvalue(k);
        Up(k / 3, k % 3) = 1.0;
      }

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
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
        Mat3x3& Uxih = Uxi.hvalue(k);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            jac(i, j, k / 3, k % 3) = Uxih(i, j);
          }
        }
      }
    }
  };
};

}  // namespace A2D

#endif  // ELASTICITY_3D_H
