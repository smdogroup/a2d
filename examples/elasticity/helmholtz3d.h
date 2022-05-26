#ifndef HELMHOLTZ_3D_H
#define HELMHOLTZ_3D_H

#include "a2dtmp.h"
#include "basis3d.h"

template <class Basis>
class HelmholtzPDE {
 public:
  static const int SPATIAL_DIM = 3;
  static const int NUM_DATA = 1;
  static const int NUM_VARS = 1;

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointSolutionArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv, QuadPointSolutionArray& U,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& res) {
    Basis::template residuals<T, HelmholtzPDE<Basis>::Impl>(data, detJ, Jinv,
                                                            Uxi, res);
    Basis::template residuals<T, HelmholtzPDE<Basis>::Impl>(data, detJ, U, res);
  }

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointSolutionArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void jacobians(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv, QuadPointSolutionArray& U,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& jac) {
    Basis::template jacobians<T, HelmholtzPDE<Basis>::Impl>(data, detJ, Jinv,
                                                            Uxi, jac);
    Basis::template jacobians<T, HelmholtzPDE<Basis>::Impl>(data, detJ, U, jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 1;

    template <typename T, class IdxType, class QuadPointData>
    static T compute_residual(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub) {
      Ub(0) = wdetJ * U0(0);

      return 0.0;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_residual(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 1, 3>& Uxi, A2D::Mat<T, 1, 3>& Uxib) {
      T r0 = data(i, j, 0);

      A2D::Vec<T, 3> Ux;
      // Ux = Uxi * Jinv
      Ux(0) = Uxi(0, 0) * Jinv(0, 0) + Uxi(0, 1) * Jinv(1, 0) +
              Uxi(0, 2) * Jinv(2, 0);
      Ux(1) = Uxi(0, 0) * Jinv(0, 1) + Uxi(0, 1) * Jinv(1, 1) +
              Uxi(0, 2) * Jinv(2, 1);
      Ux(2) = Uxi(0, 0) * Jinv(0, 2) + Uxi(0, 1) * Jinv(1, 2) +
              Uxi(0, 2) * Jinv(2, 2);

      A2D::Vec<T, 3> Uxb;
      Uxb(0) = wdetJ * r0 * r0 * Ux(0);
      Uxb(1) = wdetJ * r0 * r0 * Ux(1);
      Uxb(2) = wdetJ * r0 * r0 * Ux(2);

      // Ux = Uxi * Jinv
      // Uxb^{T} dot{Ux} = Uxb^{T} * dot{Uxi} * Jinv = Jinv * Uxb^{T} * dot{Uxi}
      // => Uxib^{T} = Jinv * Uxb^{T}
      // => Uxib = Uxb * Jinv^{T}

      Uxib(0, 0) =
          Uxb(0) * Jinv(0, 0) + Uxb(1) * Jinv(0, 1) + Uxb(2) * Jinv(0, 2);
      Uxib(0, 1) =
          Uxb(0) * Jinv(1, 0) + Uxb(1) * Jinv(1, 1) + Uxb(2) * Jinv(1, 2);
      Uxib(0, 2) =
          Uxb(0) * Jinv(2, 0) + Uxb(1) * Jinv(2, 1) + Uxb(2) * Jinv(2, 2);

      return 0.0;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_jacobian(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub,
                              A2D::Mat<T, 1, 1>& jac) {
      jac(0, 0) = wdetJ;

      return 0.0;
    }

    template <typename T, class IdxType, class QuadPointData>
    static T compute_jacobian(IdxType i, IdxType j, QuadPointData& data,
                              T wdetJ, A2D::Mat<T, 3, 3>& Jinv,
                              A2D::Mat<T, 1, 3>& Uxi0, A2D::Mat<T, 1, 3>& Uxib,
                              A2D::SymmTensor<T, 1, 3>& jac) {
      T r0 = data(i, j, 0);
      T r2 = r0 * r0 * wdetJ;

      // Uxib = Uxb * Jinv^{T}
      // Uxb = r0 * r0 * Uxb = r0 * r0 * Uxi * Jinv
      // Uxib = r0 * r0 * Jinv * Jinv^{T}

      jac(0, 0, 0, 0) =
          r2 * (Jinv(0, 0) * Jinv(0, 0) + Jinv(0, 1) * Jinv(0, 1) +
                Jinv(0, 2) * Jinv(0, 2));
      jac(0, 1, 0, 1) =
          r2 * (Jinv(1, 0) * Jinv(1, 0) + Jinv(1, 1) * Jinv(1, 1) +
                Jinv(1, 2) * Jinv(1, 2));
      jac(0, 2, 0, 2) =
          r2 * (Jinv(2, 0) * Jinv(2, 0) + Jinv(2, 1) * Jinv(2, 1) +
                Jinv(2, 2) * Jinv(2, 2));

      jac(0, 0, 0, 1) =
          r2 * (Jinv(0, 0) * Jinv(1, 0) + Jinv(0, 1) * Jinv(1, 1) +
                Jinv(0, 2) * Jinv(1, 2));
      jac(0, 0, 0, 2) =
          r2 * (Jinv(0, 0) * Jinv(2, 0) + Jinv(0, 1) * Jinv(2, 1) +
                Jinv(0, 2) * Jinv(2, 2));
      jac(0, 1, 0, 2) =
          r2 * (Jinv(1, 0) * Jinv(2, 0) + Jinv(1, 1) * Jinv(2, 1) +
                Jinv(1, 2) * Jinv(2, 2));

      return 0.0;
    }
  };
};

#endif  // HELMHOLTZ_3D_H