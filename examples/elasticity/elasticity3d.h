#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "model.h"
#include "multiarray.h"

namespace A2D {

template <class IdxType, class ScalarType, class Basis>
class NonlinearElasticity3D
    : public PDEModel<IdxType, ScalarType, Basis, 3, 2, 6> {
 public:
  // Finite-element basis class
  static const int NUM_VARS = 3;  // Number of variables per node
  static const int NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef PDEModel<IdxType, ScalarType, Basis, 3, 2, 6> base;

  NonlinearElasticity3D(const int nelems, const int nnodes, const int nbcs)
      : PDEModel<IdxType, ScalarType, Basis, 3, 2, 6>(nelems, nnodes, nbcs) {}

  void reset_nodes() {
    typename base::ConnArray& conn = this->get_conn();
    typename base::NodeArray& X = this->get_nodes();
    typename base::ElemNodeArray& Xe = this->get_elem_nodes();
    typename base::QuadNodeArray& Xq = this->get_quad_nodes();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::NullSpaceArray& B = this->get_null_space();

    element_scatter(conn, X, Xe);
    Basis::template interp<base::spatial_dim>(Xe, Xq);
    Basis::template compute_jtrans<ScalarType>(Xe, detJ, Jinv);

    B.zero();
    for (IdxType i = 0; i < this->nnodes; i++) {
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

    typename base::BCsArray& bcs = this->get_bcs();
    A2D::VecZeroBCRows(bcs, B);
  }

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

  void add_residual(typename base::SolutionArray& res) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemResArray& elem_res = this->get_elem_res();

    Basis::template residuals<
        ScalarType, NonlinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_res);

    typename base::ConnArray& conn = this->get_conn();
    typename base::BCsArray& bcs = this->get_bcs();
    element_gather_add(conn, elem_res, res);
    A2D::VecZeroBCRows(bcs, res);
  }

  void add_jacobian(typename base::SparseMat& J) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemJacArray& elem_jac = this->get_elem_jac();

    Basis::template jacobians<
        ScalarType, NonlinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_jac);

    typename base::ConnArray& conn = this->get_conn();
    typename base::BCsArray& bcs = this->get_bcs();
    A2D::BSRMatAddElementMatrices(conn, elem_jac, J);
    A2D::BSRMatZeroBCRows(bcs, J);
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

    template <typename T, class I, class QuadPointData>
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

    template <typename T, class I, class QuadPointData>
    static void compute_res_adjoint_data(I i, I j, QuadPointData& data, T wdetJ,
                                         A2D::Mat<T, 3, 3>& Jinv,
                                         A2D::Mat<T, 3, 3>& Uxi0,
                                         A2D::Mat<T, 3, 3>& Psi,
                                         QuadPointData& dfdx) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      Mat3x3 Ux0, Uxb, Uxib;
      SymmMat3x3 E0, Eb;

      const int N = 1;
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;
      A2D::A2DScalar<N, T> mu(data(i, j, 0)), lambda(data(i, j, 1));

      // Copy over the Psi values
      Uxi.Ap[0] = Psi;

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

      dfdx(i, j, 0) += mu.hvalue[0];
      dfdx(i, j, 1) += lambda.hvalue[0];
    }

    template <typename T, class I, class QuadPointData>
    static void compute_res_adjoint_nodes(I i, I j, QuadPointData& data,
                                          T wdetJ, A2D::Mat<T, 3, 3>& Jinv0,
                                          A2D::Mat<T, 3, 3>& Uxi0,
                                          A2D::Mat<T, 3, 3>& Psi, T dfdwdetJ,
                                          A2D::Mat<T, 3, 3>& dfdJinv) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      Mat3x3 Ux0, Uxb, Uxib, Jinvb;
      SymmMat3x3 E0, Eb;

      const int N = 1;
      A2D::A2DMat<N, Mat3x3> Jinv(Jinv0, Jinvb);
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;
      T mu(data(i, j, 0)), lambda(data(i, j, 1));

      // Copy over the Psi values
      Uxi.Ap[0] = Psi;

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

      dfdwdetJ = output.value;
      dfdJinv = Jinv.Ah;
    }
  };
};

template <class IdxType, class ScalarType, class Basis>
class LinearElasticity3D
    : public PDEModel<IdxType, ScalarType, Basis, 3, 2, 6> {
 public:
  // Finite-element basis class
  static const int NUM_VARS = 3;  // Number of variables per node
  static const int NUM_DATA = 2;  // Data points per quadrature point

  // Short cut for base class name
  typedef PDEModel<IdxType, ScalarType, Basis, 3, 2, 6> base;

  LinearElasticity3D(const int nelems, const int nnodes, const int nbcs)
      : PDEModel<IdxType, ScalarType, Basis, 3, 2, 6>(nelems, nnodes, nbcs) {}

  void reset_nodes() {
    typename base::ConnArray& conn = this->get_conn();
    typename base::NodeArray& X = this->get_nodes();
    typename base::ElemNodeArray& Xe = this->get_elem_nodes();
    typename base::QuadNodeArray& Xq = this->get_quad_nodes();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::NullSpaceArray& B = this->get_null_space();

    element_scatter(conn, X, Xe);
    Basis::template interp<base::spatial_dim>(Xe, Xq);
    Basis::template compute_jtrans<ScalarType>(Xe, detJ, Jinv);

    B.zero();
    for (IdxType i = 0; i < this->nnodes; i++) {
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

    typename base::BCsArray& bcs = this->get_bcs();
    A2D::VecZeroBCRows(bcs, B);
  }

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

  void add_residual(typename base::SolutionArray& res) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemResArray& elem_res = this->get_elem_res();

    Basis::template residuals<
        ScalarType, LinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_res);

    typename base::ConnArray& conn = this->get_conn();
    typename base::BCsArray& bcs = this->get_bcs();
    element_gather_add(conn, elem_res, res);
    A2D::VecZeroBCRows(bcs, res);
  }

  void add_jacobian(typename base::SparseMat& J) {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemJacArray& elem_jac = this->get_elem_jac();

    Basis::template jacobians<
        ScalarType, LinearElasticity3D<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_jac);

    typename base::ConnArray& conn = this->get_conn();
    typename base::BCsArray& bcs = this->get_bcs();
    A2D::BSRMatAddElementMatrices(conn, elem_jac, J);
    A2D::BSRMatZeroBCRows(bcs, J);
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

    template <typename T, class I, class QuadPointData>
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
