#ifndef ELASTICITY_3D_H
#define ELASTICITY_3D_H

#include <cstddef>

#include "a2dtmp.h"
#include "basis3d.h"

template <class Basis>
class NonlinearElasticity3D {
 public:
  static const int SPATIAL_DIM = 3;
  static const int NUM_DATA = 2;
  static const int NUM_VARS = 3;

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray>
  static void energy(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                     QuadPointJacobianArray& Jinv, QuadPointGradientArray& Uxi,
                     T& energy) {
    Basis::template energy<T, NonlinearElasticity3D<Basis>::Impl>(
        data, detJ, Jinv, Uxi, energy);
  }

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& res) {
    Basis::template residuals<T, NonlinearElasticity3D<Basis>::Impl>(
        data, detJ, Jinv, Uxi, res);
  }

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void jacobians(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& jac) {
    Basis::template jacobians<T, NonlinearElasticity3D<Basis>::Impl>(
        data, detJ, Jinv, Uxi, jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 3;

    template <typename T, class IdxType, class QuadPointData>
    static T compute_energy(IdxType i, IdxType j, QuadPointData& data, T wdetJ,
                            A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux;
      SymmMat3x3 E, S;
      T output;

      A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      A2D::Mat3x3GreenStrain(Ux, E);
      A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      A2D::Symm3x3SymmMultTrace(S, E, output);

      return 0.5 * output * wdetJ;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_residual(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 3, 3>& Uxi0,
                              A2D::Mat<T, 3, 3>& Uxib) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;
      SymmMat3x3 S0, Sb;

      A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
      A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
      A2D::ADMat<SymmMat3x3> S(S0, Sb);
      A2D::ADMat<SymmMat3x3> E(E0, Eb);
      A2D::ADScalar<T> output;

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3GreenStrain(Ux, E);
      auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

      output.bvalue = 0.5 * wdetJ;

      trace.reverse();
      constitutive.reverse();
      strain.reverse();
      mult.reverse();

      return output.value;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_jacobian(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib,
                              A2D::SymmTensor<T, 3, 3>& jac) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb, Uxp, Uxh;
      Mat3x3 Uxip, Uxih;
      SymmMat3x3 E0, Eb, Ep, Eh;
      SymmMat3x3 S0, Sb, Sp, Sh;

      const int N = 9;
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> S(S0, Sb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;

      // Set up the seed values
      for (int k = 0; k < N; k++) {
        Mat3x3& Up = Uxi.pvalue(k);
        Up(k / 3, k % 3) = 1.0;
      }

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3GreenStrain(Ux, E);
      auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

      output.bvalue = 0.5 * wdetJ;

      trace.reverse();
      constitutive.reverse();
      strain.reverse();
      mult.reverse();

      mult.hforward();
      strain.hforward();
      constitutive.hforward();
      trace.hreverse();
      constitutive.hreverse();
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

template <class Basis>
class LinearElasticity3D {
 public:
  static const int SPATIAL_DIM = 3;
  static const int NUM_DATA = 2;
  static const int NUM_VARS = 3;

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray>
  static void energy(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                     QuadPointJacobianArray& Jinv, QuadPointGradientArray& Uxi,
                     T& energy) {
    Basis::template energy<T, LinearElasticity3D<Basis>::Impl>(data, detJ, Jinv,
                                                               Uxi, energy);
  }

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& res) {
    Basis::template residuals<T, LinearElasticity3D<Basis>::Impl>(
        data, detJ, Jinv, Uxi, res);
  }

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void jacobians(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& jac) {
    Basis::template jacobians<T, LinearElasticity3D<Basis>::Impl>(
        data, detJ, Jinv, Uxi, jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 3;

    template <typename T, class IdxType, class QuadPointData>
    static T compute_energy(IdxType i, IdxType j, QuadPointData& data, T wdetJ,
                            A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 3, 3>& Uxi) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux;
      SymmMat3x3 E, S;
      T output;

      A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      A2D::Mat3x3LinearGreenStrain(Ux, E);
      A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      A2D::Symm3x3SymmMultTrace(S, E, output);

      return 0.5 * output * wdetJ;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_residual(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 3, 3>& Uxi0,
                              A2D::Mat<T, 3, 3>& Uxib) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb;
      SymmMat3x3 E0, Eb;
      SymmMat3x3 S0, Sb;

      A2D::ADMat<Mat3x3> Uxi(Uxi0, Uxib);
      A2D::ADMat<Mat3x3> Ux(Ux0, Uxb);
      A2D::ADMat<SymmMat3x3> S(S0, Sb);
      A2D::ADMat<SymmMat3x3> E(E0, Eb);
      A2D::ADScalar<T> output;

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
      auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

      output.bvalue = 0.5 * wdetJ;

      trace.reverse();
      constitutive.reverse();
      strain.reverse();
      mult.reverse();

      return output.value;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_jacobian(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 3, 3>& Uxi0, A2D::Mat<T, 3, 3>& Uxib,
                              A2D::SymmTensor<T, 3, 3>& jac) {
      typedef A2D::SymmMat<T, 3> SymmMat3x3;
      typedef A2D::Mat<T, 3, 3> Mat3x3;

      T mu(data(i, j, 0)), lambda(data(i, j, 1));
      Mat3x3 Ux0, Uxb, Uxp, Uxh;
      Mat3x3 Uxip, Uxih;
      SymmMat3x3 E0, Eb, Ep, Eh;
      SymmMat3x3 S0, Sb, Sp, Sh;

      const int N = 9;
      A2D::A2DMat<N, Mat3x3> Uxi(Uxi0, Uxib);
      A2D::A2DMat<N, Mat3x3> Ux(Ux0, Uxb);
      A2D::A2DMat<N, SymmMat3x3> S(S0, Sb);
      A2D::A2DMat<N, SymmMat3x3> E(E0, Eb);
      A2D::A2DScalar<N, T> output;

      // Set up the seed values
      for (int k = 0; k < N; k++) {
        Mat3x3& Up = Uxi.pvalue(k);
        Up(k / 3, k % 3) = 1.0;
      }

      auto mult = A2D::Mat3x3MatMult(Uxi, Jinv, Ux);
      auto strain = A2D::Mat3x3LinearGreenStrain(Ux, E);
      auto constitutive = A2D::Symm3x3IsotropicConstitutive(mu, lambda, E, S);
      auto trace = A2D::Symm3x3SymmMultTrace(S, E, output);

      output.bvalue = 0.5 * wdetJ;

      trace.reverse();
      constitutive.reverse();
      strain.reverse();
      mult.reverse();

      mult.hforward();
      strain.hforward();
      constitutive.hforward();
      trace.hreverse();
      constitutive.hreverse();
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
