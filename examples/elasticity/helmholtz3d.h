#ifndef HELMHOLTZ_3D_H
#define HELMHOLTZ_3D_H

#include "a2dtmp.h"
#include "basis3d.h"
#include "model.h"
#include "multiarray.h"

namespace A2D {

/*
  Store all the type information about the PDE in one place
*/
template <typename I, typename T>
class HelmholtzPDE {
 public:
  static const index_t spatial_dim = 3;
  static const index_t vars_per_node = 1;
  static const index_t null_space_dim = 1;
  static const index_t data_per_point = 1;

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
    B.fill(1.0);
  }
};

template <typename I, typename T, class Basis>
class HelmholtzElement3D
    : public ElementBasis<I, T, HelmholtzPDE<I, T>, Basis> {
 public:
  // Finite-element basis class
  static const index_t NUM_VARS = 1;  // Number of variables per node
  static const index_t NUM_DATA = 1;  // Data points per quadrature point

  // Short cut for base class name
  typedef ElementBasis<I, T, HelmholtzPDE<I, T>, Basis> base;

  HelmholtzElement3D(const int nelems)
      : ElementBasis<I, T, HelmholtzPDE<I, T>, Basis>(nelems) {}

  // Add the residual contribution
  void add_residual(typename HelmholtzPDE<I, T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());
    elem_res.zero();

    Basis::template residuals<T, HelmholtzElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_Jinv(),
        this->get_quad_gradient(), elem_res);
    Basis::template residuals<T, HelmholtzElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_quad_solution(),
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  // Add the element Jacobian contribution
  void add_jacobian(typename HelmholtzPDE<I, T>::SparseMat& J) {
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());
    elem_jac.zero();

    Basis::template jacobians<T, HelmholtzElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_quad_gradient(),
        elem_jac);

    Basis::template jacobians<T, HelmholtzElement3D<I, T, Basis>::Impl>(
        this->get_quad_data(), this->get_detJ(), this->get_quad_solution(),
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  class Impl {
   public:
    static const index_t NUM_VARS = 1;

    template <class QuadPointData>
    static void compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub) {
      Ub(0) = wdetJ * U0(0);
    }

    template <class QuadPointData>
    static void compute_residual(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 1, 3>& Uxi,
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
      // Uxb^{T} dot{Ux} = Uxb^{T} * dot{Uxi} * Jinv
      //                 = Jinv * Uxb^{T} * dot{Uxi}
      // => Uxib^{T} = Jinv * Uxb^{T}
      // => Uxib = Uxb * Jinv^{T}

      Uxib(0, 0) =
          Uxb(0) * Jinv(0, 0) + Uxb(1) * Jinv(0, 1) + Uxb(2) * Jinv(0, 2);
      Uxib(0, 1) =
          Uxb(0) * Jinv(1, 0) + Uxb(1) * Jinv(1, 1) + Uxb(2) * Jinv(1, 2);
      Uxib(0, 2) =
          Uxb(0) * Jinv(2, 0) + Uxb(1) * Jinv(2, 1) + Uxb(2) * Jinv(2, 2);
    }

    template <class QuadPointData>
    static void compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Vec<T, 1>& U0, A2D::Vec<T, 1>& Ub,
                                 A2D::Mat<T, 1, 1>& jac) {
      jac(0, 0) = wdetJ;
    }

    template <class QuadPointData>
    static void compute_jacobian(I i, I j, QuadPointData& data, T wdetJ,
                                 A2D::Mat<T, 3, 3>& Jinv,
                                 A2D::Mat<T, 1, 3>& Uxi0,
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
    }
  };
};

}  // namespace A2D

#endif  // HELMHOLTZ_3D_H