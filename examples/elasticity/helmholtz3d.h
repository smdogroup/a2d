#ifndef HELMHOLTZ_3D_H
#define HELMHOLTZ_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "model.h"
#include "multiarray.h"

template <class IdxType, class ScalarType, class Basis>
class HelmholtzPDE : public PDEModel<IdxType, ScalarType, Basis, 1, 1> {
 public:
  // Finite-element basis class
  static const int NUM_VARS = 1;  // Number of variables per node
  static const int NUM_DATA = 1;  // Data points per quadrature point

  // Short cut for base class name
  typedef PDEModel<IdxType, ScalarType, Basis, 1, 1> base;

  HelmholtzPDE(const int nelems, const int nnodes)
      : PDEModel<IdxType, ScalarType, Basis, 1, 1>(nelems, nnodes) {}

  template <class ResArray>
  void add_residuals(ResArray& res) {
    // Get data for computing the residuals
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadSolnArray& Uq = this->get_quad_solution();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemResArray& elem_res = this->get_elem_res();

    Basis::template residuals<ScalarType,
                              HelmholtzPDE<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_res);
    Basis::template residuals<ScalarType,
                              HelmholtzPDE<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Uq, elem_res);

    typename base::ConnArray& conn = this->get_conn();
    element_gather_add(conn, elem_res, res);
  }

  void add_jacobians() {
    typename base::QuadDataArray& data = this->get_quad_data();
    typename base::QuadDetArray& detJ = this->get_detJ();
    typename base::QuadJtransArray& Jinv = this->get_Jinv();
    typename base::QuadSolnArray& Uq = this->get_quad_solution();
    typename base::QuadGradArray& Uxi = this->get_quad_gradient();
    typename base::ElemJacArray& elem_jac = this->get_elem_jac();

    Basis::template jacobians<ScalarType,
                              HelmholtzPDE<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Jinv, Uxi, elem_jac);
    Basis::template jacobians<ScalarType,
                              HelmholtzPDE<IdxType, ScalarType, Basis>::Impl>(
        data, detJ, Uq, elem_jac);
  }

  class Impl {
   public:
    static const int NUM_VARS = 1;

    template <typename T, class I, class QuadPointData>
    static T compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub) {
      Ub(0) = wdetJ * U0(0);

      return 0.0;
    }

    template <typename T, class I, class QuadPointData>
    static T compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 1, 3>& Uxi,
                              A2D::Mat<T, 1, 3>& Uxib) {
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

    template <typename T, class I, class QuadPointData>
    static T compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub,
                              A2D::Mat<T, 1, 1>& jac) {
      jac(0, 0) = wdetJ;

      return 0.0;
    }

    template <typename T, class I, class QuadPointData>
    static T compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                              A2D::Mat<T, 3, 3>& Jinv, A2D::Mat<T, 1, 3>& Uxi0,
                              A2D::Mat<T, 1, 3>& Uxib,
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