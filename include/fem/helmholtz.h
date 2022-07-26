#ifndef A2D_HELMHOLTZ_H
#define A2D_HELMHOLTZ_H

#include "a2dtmp3d.h"
#include "basis.h"
#include "constitutive.h"
#include "element.h"
#include "multiarray.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/*
  Store all the type information about the PDE in one place
*/
template <index_t spatial_dim, typename I, typename T>
class HelmholtzPDEInfo {
 public:
  static_assert(spatial_dim == 3 or spatial_dim == 2,
                "spatial_dim must be 2 or 3");
  static const index_t SPATIAL_DIM = spatial_dim;
  static const index_t vars_per_node = 1;
  static const index_t null_space_dim = 1;
  static const index_t data_per_point = 1;  // Right-hand-side data
  static const index_t dvs_per_point = 1;   // Same as the solution space

  // Layout for the boundary conditions
  typedef A2D::CLayout<2> BCsLayout;
  typedef A2D::MultiArray<I, BCsLayout> BCsArray;

  // Layout for the nodes
  typedef A2D::CLayout<SPATIAL_DIM> NodeLayout;
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
    B.fill(1.0);
  }
};

template <typename I, typename T, class BasisOps>
class HelmholtzElement
    : public ElementBasis<I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps> {
 public:
  // Finite-element basis class
  static const index_t SPATIAL_DIM = BasisOps::SPATIAL_DIM;
  static const index_t vars_per_node =
      HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I,
                       T>::vars_per_node;  // Number of variables per node

  // Short cut for base class name
  typedef ElementBasis<I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                       BasisOps>
      base;

  HelmholtzElement(const int nelems, double r0)
      : ElementBasis<I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems),
        r0(r0) {}

  template <typename IdxType>
  HelmholtzElement(const index_t nelems, const IdxType conn_[], double r0)
      : ElementBasis<I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                     BasisOps>(nelems, conn_),
        r0(r0) {}

  // The radius (global value)
  const double r0;

  // Add the residual contribution
  void add_residual(typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I,
                                              T>::SolutionArray& res) {
    // Allocate the element residual
    typename base::ElemResArray elem_res(this->get_elem_res_layout());

    // Retrieve the element data
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uq = this->get_quad_solution();
    auto Uxi = this->get_quad_gradient();
    auto data = this->get_quad_data();

    double r2 = r0 * r0;

    BasisOps::template residuals<T, vars_per_node>(
        detJ, Uq,
        [&data](index_t i, index_t j, T wdetJ, A2D::Vec<T, 1>& U0,
                A2D::Vec<T, 1>& Ub) -> void {
          Ub(0) = wdetJ * (U0(0) - data(i, j, 0));
        },
        elem_res);

    BasisOps::template residuals<T, vars_per_node>(
        detJ, Jinv, Uxi,
        [r2](index_t i, index_t j, T wdetJ,
             A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Jinv,
             A2D::Mat<T, 1, SPATIAL_DIM>& Uxi,
             A2D::Mat<T, 1, SPATIAL_DIM>& Uxib) -> void {
          A2D::Vec<T, SPATIAL_DIM> Ux;
          // Ux = Uxi * Jinv
          for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
            Ux(idim) = 0.0;
            for (index_t jdim = 0; jdim < SPATIAL_DIM; jdim++) {
              Ux(idim) += Uxi(0u, jdim) * Jinv(jdim, idim);
            }
          }

          A2D::Vec<T, SPATIAL_DIM> Uxb;
          for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
            Uxb(idim) = wdetJ * r2 * Ux(idim);
          }

          // Ux = Uxi * Jinv
          // Uxb^{T} dot{Ux} = Uxb^{T} * dot{Uxi} * Jinv
          //                 = Jinv * Uxb^{T} * dot{Uxi}
          // => Uxib^{T} = Jinv * Uxb^{T}
          // => Uxib = Uxb * Jinv^{T}

          for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
            Uxib(0u, idim) = 0.0;
            for (index_t jdim = 0; jdim < SPATIAL_DIM; jdim++) {
              Uxib(0u, idim) += Uxb(jdim) * Jinv(idim, jdim);
            }
          }
        },
        elem_res);

    VecElementGatherAdd(this->get_conn(), elem_res, res);
  }

  // Add the element Jacobian contribution
  void add_jacobian(
      typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>::SparseMat& J) {
    Timer timer("HelmholtzElement::add_jacobian()");
    typename base::ElemJacArray elem_jac(this->get_elem_jac_layout());

    // Retrieve the element data
    auto detJ = this->get_detJ();
    auto Jinv = this->get_Jinv();
    auto Uq = this->get_quad_solution();
    auto Uxi = this->get_quad_gradient();

    double r2 = r0 * r0;

    BasisOps::template jacobians<T, vars_per_node>(
        detJ, Uq,
        [](index_t i, index_t j, T wdetJ, A2D::Vec<T, 1>& U0,
           A2D::SymmMat<T, 1>& jac) -> void { jac(0, 0) = wdetJ; },
        elem_jac);

    BasisOps::template jacobians<T, vars_per_node>(
        detJ, Jinv, Uxi,
        [r2](index_t i, index_t j, T wdetJ,
             A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM>& Jinv,
             A2D::Mat<T, 1, SPATIAL_DIM>& Uxi,
             A2D::Mat<T, 1, SPATIAL_DIM>& Uxib,
             A2D::SymmTensor<T, 1, SPATIAL_DIM>& jac) -> void {
          T wr2 = r2 * wdetJ;

          // Uxib = Uxb * Jinv^{T}
          // Uxb = r0 * r0 * Uxb = r0 * r0 * Uxi * Jinv
          // Uxib = r0 * r0 * Jinv * Jinv^{T}

          for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
            for (index_t jdim = idim; jdim < SPATIAL_DIM; jdim++) {
              jac(0u, idim, 0u, jdim) = 0.0;
              for (index_t kdim = 0; kdim < SPATIAL_DIM; kdim++) {
                jac(0u, idim, 0u, jdim) += Jinv(idim, kdim) * Jinv(jdim, kdim);
              }
              jac(0u, idim, 0u, jdim) *= wr2;
            }
          }
        },
        elem_jac);

    A2D::BSRMatAddElementMatrices(this->get_conn(), elem_jac, J);
  }

  void add_adjoint_dfddata(typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                     T>::SolutionArray& psi,
                           typename base::QuadDataArray& dfdx) {
    // Retrieve the element data
    auto conn = this->get_conn();
    auto detJ = this->get_detJ();
    auto Uq = this->get_quad_solution();

    // Compute the element adjoint data
    typename base::ElemSolnArray pe(this->get_elem_solution_layout());
    typename base::QuadSolnArray psiq(this->get_quad_solution_layout());

    VecElementScatter(conn, psi, pe);
    BasisOps::template interp<vars_per_node>(pe, psiq);

    // Compute the product
    BasisOps::template adjoint_product<T, vars_per_node>(
        detJ, Uq, psiq,
        [&dfdx](index_t i, index_t j, T wdetJ, A2D::Vec<T, 1>& U0,
                A2D::Vec<T, 1>& Psi) { dfdx(i, j, 0) -= wdetJ * Psi(0); });
  }
};

template <typename I, typename T, class BasisOps>
class HelmholtzConstitutive
    : public ConstitutiveBase<I, T,
                              HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>> {
 public:
  static const index_t dvs_per_point =
      HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>::dvs_per_point;
  static const index_t nodes_per_elem = BasisOps::NUM_NODES;

  typedef A2D::CLayout<nodes_per_elem, dvs_per_point> ElemDesignLayout;
  typedef A2D::MultiArray<T, ElemDesignLayout> ElemDesignArray;

  HelmholtzConstitutive(
      std::shared_ptr<ElementBasis<
          I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
          element)
      : element(element),
        elem_design_layout(element->nelems),
        xe(elem_design_layout) {}

  /*
    Set the design variables values into the element object
  */
  void set_design_vars(
      typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>::DesignArray& x) {
    // Set the design variable values
    // x -> xe -> xq
    auto conn = element->get_conn();
    auto data = element->get_quad_data();
    VecElementScatter(conn, x, xe);
    BasisOps::template interp<dvs_per_point>(xe, data);
  }

  /*
    Compute the derivative of the adjoint-residual product data w.r.t. x
  */
  void add_adjoint_dfdx(typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                  T>::SolutionArray& psi,
                        typename HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I,
                                                  T>::DesignArray& dfdx) {
    typename ElementBasis<I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>,
                          BasisOps>::QuadDataArray
        dfddata(element->get_quad_data_layout());

    // Compute the product of the adjoint with the derivatives of
    // the residuals w.r.t. the element data
    element->add_adjoint_dfddata(psi, dfddata);

    ElemDesignArray dfdxe(elem_design_layout);
    BasisOps::template interpReverseAdd<dvs_per_point>(dfddata, dfdxe);

    auto conn = element->get_conn();
    VecElementGatherAdd(conn, dfdxe, dfdx);
  }

 private:
  // Reference to the element class
  std::shared_ptr<ElementBasis<
      I, T, HelmholtzPDEInfo<BasisOps::SPATIAL_DIM, I, T>, BasisOps>>
      element;

  // Design variable views
  ElemDesignLayout elem_design_layout;
  ElemDesignArray xe;
};

}  // namespace A2D

#endif  // A2D_HELMHOLTZ_H