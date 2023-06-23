#ifndef A2D_ELEMENT_H
#define A2D_ELEMENT_H

#include "a2dobjs.h"
#include "sparse/sparse_amg.h"
#include "sparse/sparse_matrix.h"

namespace A2D {

/**
 * @brief Element base class.
 *
 * This defines an element that is compatible with the given PDEInfo. You can
 * set the node locations into the element, add the non-zero-pattern of the
 * Jacobian matrix via the add_node_set as well as adding the residual and
 * Jacobian contributions
 *
 * @tparam I index type
 * @tparam T data type
 * @tparam PDEInfo a class that stores type information for problem's PDE
 */
template <typename I, typename T, class PDEInfo>
class ElementBase {
 public:
  virtual ~ElementBase() {}
  virtual void set_nodes(typename PDEInfo::NodeArray& X) = 0;
  virtual void add_node_set(Kokkos::UnorderedMap<COO<I>, void>& node_set) = 0;
  virtual void set_solution(typename PDEInfo::SolutionArray& U) = 0;
  virtual T energy() { return T(0.0); }
  virtual void add_residual(typename PDEInfo::SolutionArray& res) = 0;
  virtual void add_jacobian(typename PDEInfo::SparseMat& jac) = 0;
  virtual void add_adjoint_dfdnodes(typename PDEInfo::SolutionArray& psi,
                                    typename PDEInfo::NodeArray& dfdx) {}
};

/**
 * @brief A partial implementation of the element class with a prescribed basis
 *        function
 *
 * @tparam I index type
 * @tparam T data type
 * @tparam PDEInfo a PDEInfo class that stores type information
 * @tparam BasisOps the collection of basis operations
 */
template <typename I, typename T, class PDEInfo, class BasisOps>
class ElementBasis : public ElementBase<I, T, PDEInfo> {
 public:
  static_assert(PDEInfo::SPATIAL_DIM == BasisOps::SPATIAL_DIM,
                "Check consistency");
  static const index_t spatial_dim = PDEInfo::SPATIAL_DIM;
  static const index_t vars_per_node = PDEInfo::vars_per_node;
  static const index_t data_per_point = PDEInfo::data_per_point;
  static const index_t nodes_per_elem = BasisOps::NUM_NODES;
  static const index_t quad_pts_per_elem = BasisOps::quadrature::NUM_QUAD_PTS;

  // Connectivity array
  using ConnArray = A2D::MultiArrayNew<I* [nodes_per_elem]>;

  // Element-node arrays
  using ElemNodeArray = A2D::MultiArrayNew<T* [nodes_per_elem][spatial_dim]>;
  using ElemSolnArray = A2D::MultiArrayNew<T* [nodes_per_elem][vars_per_node]>;

  // Element-quadrature arrays
  using QuadNodeArray = A2D::MultiArrayNew<T* [quad_pts_per_elem][spatial_dim]>;
  using QuadSolnArray =
      A2D::MultiArrayNew<T* [quad_pts_per_elem][vars_per_node]>;
  using QuadGradArray =
      A2D::MultiArrayNew<T* [quad_pts_per_elem][vars_per_node][spatial_dim]>;
  using QuadDetArray = A2D::MultiArrayNew<T* [quad_pts_per_elem]>;
  using QuadJtransArray =
      A2D::MultiArrayNew<T* [quad_pts_per_elem][spatial_dim][spatial_dim]>;
  using QuadDataArray =
      A2D::MultiArrayNew<T* [quad_pts_per_elem][data_per_point]>;

  // Residual data
  using ElemResArray = A2D::MultiArrayNew<T* [nodes_per_elem][vars_per_node]>;
  using ElemJacArray = A2D::MultiArrayNew<
      T* [nodes_per_elem][nodes_per_elem][vars_per_node][vars_per_node]>;

  ElementBasis(const index_t nelems) : nelems(nelems) {
    Timer t("ElementBasis::ElementBasis(1)");
    conn = ConnArray("conn", nelems);
    Xe = ElemNodeArray("Xe", nelems);
    Ue = ElemSolnArray("Ue", nelems);
    Xq = QuadNodeArray("Xq", nelems);
    Uq = QuadSolnArray("Uq", nelems);
    Uxi = QuadGradArray("Uxi", nelems);
    detJ = QuadDetArray("detJ", nelems);
    Jinv = QuadJtransArray("Jinv", nelems);
    data = QuadDataArray("data", nelems);
  }

  template <typename IdxType>
  ElementBasis(const index_t nelems, const IdxType conn_[]) : nelems(nelems) {
    Timer t("ElementBasis::ElementBasis(2)");
    conn = ConnArray("conn", nelems);
    Xe = ElemNodeArray("Xe", nelems);
    Ue = ElemSolnArray("Ue", nelems);
    Xq = QuadNodeArray("Xq", nelems);
    Uq = QuadSolnArray("Uq", nelems);
    Uxi = QuadGradArray("Uxi", nelems);
    detJ = QuadDetArray("detJ", nelems);
    Jinv = QuadJtransArray("Jinv", nelems);
    data = QuadDataArray("data", nelems);
    // Set the connectivity
    for (I i = 0; i < nelems; i++) {
      for (I j = 0; j < nodes_per_elem; j++) {
        conn(i, j) = conn_[nodes_per_elem * i + j];
      }
    }
  }
  virtual ~ElementBasis() {}

  const index_t nelems;  // Number of elements in the model

  // Get the connectivity array
  ConnArray& get_conn() { return conn; }

  // Get the data associated with the quadrature points
  QuadDataArray& get_quad_data() { return data; }

  // Add the node set to the connectivity, use Kokkos unordered set
  void add_node_set(Kokkos::UnorderedMap<COO<I>, void>& node_set) {
    BSRMatAddConnectivity(conn, node_set);
  }

  // Set the nodes from the array
  void set_nodes(typename PDEInfo::NodeArray& X) {
    VecElementScatter(conn, X, Xe);
    BasisOps::template interp<spatial_dim>(Xe, Xq);
    BasisOps::template compute_jtrans<T>(Xe, detJ, Jinv);
  }

  // Set the solution
  void set_solution(typename PDEInfo::SolutionArray& U) {
    VecElementScatter(conn, U, Ue);
    BasisOps::template interp<vars_per_node>(Ue, Uq);
    BasisOps::template gradient<T, vars_per_node>(Ue, Uxi);
  }

  // Pure virtual member functions that need an override
  // virtual void add_residual(typename PDEInfo::SolutionArray& res) = 0;
  // virtual void add_jacobian(typename PDEInfo::SparseMat& jac) = 0;

  virtual void add_adjoint_dfddata(typename PDEInfo::SolutionArray& psi,
                                   QuadDataArray& dfdx) {}

  // Expose the underlying element data
  ElemNodeArray& get_elem_nodes() { return Xe; }
  QuadNodeArray& get_quad_nodes() { return Xq; }
  QuadDetArray& get_detJ() { return detJ; }
  QuadJtransArray& get_Jinv() { return Jinv; }
  ElemSolnArray& get_elem_solution() { return Ue; }
  QuadSolnArray& get_quad_solution() { return Uq; }
  QuadGradArray& get_quad_gradient() { return Uxi; }

 private:
  // Instances of the multi-dimensional arrays
  ConnArray conn;

  // Element-node arrays
  ElemNodeArray Xe;
  ElemSolnArray Ue;

  // Quadrature point arrays
  QuadNodeArray Xq;
  QuadSolnArray Uq;
  QuadGradArray Uxi;
  QuadDetArray detJ;
  QuadJtransArray Jinv;
  QuadDataArray data;
};

}  // namespace A2D

#endif  // A2D_ELEMENT_H