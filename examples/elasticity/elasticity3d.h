#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "model.h"
#include "multiarray.h"

template <class IdxType, class ScalarType, class Basis>
class NonlinearElasticity3D
    : public PDEModel<IdxType, ScalarType, Basis, 3, 2> {
 public:
  // Finite-element basis class
  static const int NUM_VARS = 3;  // Number of variables per node
  static const int NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef PDEModel<IdxType, ScalarType, Basis, 3, 2> base;

  NonlinearElasticity3D(const int nelems, const int nnodes)
      : PDEModel<IdxType, ScalarType, Basis, 3, 2>(nelems, nnodes) {}

  ScalarType energy() {
    ScalarType engry;
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();

    Basis::template energy<
        ScalarType, NonlinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, engry);

    return engry;
  }

  void add_residuals(typename base::SolutionArray& res) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemResArray& elem_res = this->get_elem_res();

    Basis::template residuals<
        ScalarType, NonlinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_res);

    typename base::ConnArray& conn = this->get_conn();
    element_gather_add(conn, elem_res, res);
  }

  void add_jacobians() {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemJacArray& elem_jac = this->get_elem_jac();

    Basis::template jacobians<
        ScalarType, NonlinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 3;

    template <typename T, class I, class QuadPointData>
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

    template <typename T, class I, class QuadPointData>
    static T compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi0,
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

      return output.value;
    }

    template <typename T, class I, class QuadPointData>
    static T compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi0,
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

      return output.value;
    }
  };
};

template <class IdxType, class ScalarType, class Basis>
class LinearElasticity3D : public PDEModel<IdxType, ScalarType, Basis, 3, 2> {
 public:
  // Finite-element basis class
  static const int NUM_VARS = 3;  // Number of variables per node
  static const int NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef PDEModel<IdxType, ScalarType, Basis, 3, 2> base;

  LinearElasticity3D(const int nelems, const int nnodes)
      : PDEModel<IdxType, ScalarType, Basis, 3, 2>(nelems, nnodes) {}

  ScalarType energy() {
    ScalarType engry;
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();

    Basis::template energy<
        ScalarType, LinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, engry);

    return engry;
  }

  void add_residuals(typename base::SolutionArray& res) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemResArray& elem_res = this->get_elem_res();

    Basis::template residuals<
        ScalarType, LinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_res);

    typename base::ConnArray& conn = this->get_conn();
    element_gather_add(conn, elem_res, res);
  }

  void add_jacobians() {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemJacArray& elem_jac = this->get_elem_jac();

    Basis::template jacobians<
        ScalarType, LinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 3;

    template <typename T, class I, class QuadPointData>
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

    template <typename T, class I, class QuadPointData>
    static T compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi0,
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

      return output.value;
    }

    template <typename T, class I, class QuadPointData>
    static T compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi0,
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

      return output.value;
    }
  };
};

#endif  // ELASTICITY_3D_H
