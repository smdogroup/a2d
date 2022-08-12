#ifndef A2D_ELEMENT_H
#define A2D_ELEMENT_H

#include "a2dlayout.h"
#include "multiarray.h"
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
  virtual void add_node_set(
      std::unordered_set<std::pair<I, I>, pair_hash_fun<I>>& node_set) = 0;
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

  // Set layout

  // Connectivity layout
  typedef A2D_Layout<nodes_per_elem> ConnLayout;

  // Element-node layouts
  typedef A2D_Layout<nodes_per_elem, spatial_dim> ElemNodeLayout;
  typedef A2D_Layout<nodes_per_elem, vars_per_node> ElemSolnLayout;

  // Element-quadrature layouts
  typedef A2D_Layout<quad_pts_per_elem, spatial_dim> QuadNodeLayout;
  typedef A2D_Layout<quad_pts_per_elem, vars_per_node> QuadSolnLayout;
  typedef A2D_Layout<quad_pts_per_elem, vars_per_node, spatial_dim>
      QuadGradLayout;
  typedef A2D_Layout<quad_pts_per_elem> QuadDetLayout;
  typedef A2D_Layout<quad_pts_per_elem, spatial_dim, spatial_dim>
      QuadJtransLayout;
  typedef A2D_Layout<quad_pts_per_elem, data_per_point> QuadDataLayout;

  // Residual/Jacobian layouts
  typedef A2D_Layout<nodes_per_elem, vars_per_node> ElemResLayout;
  typedef A2D_Layout<nodes_per_elem, nodes_per_elem, vars_per_node,
                     vars_per_node>
      ElemJacLayout;

  // Connectivity array
  typedef A2D::MultiArray<I, ConnLayout> ConnArray;

  // Element-node arrays
  typedef A2D::MultiArray<T, ElemNodeLayout> ElemNodeArray;
  typedef A2D::MultiArray<T, ElemSolnLayout> ElemSolnArray;

  // Element-quadrature arrays
  typedef A2D::MultiArray<T, QuadNodeLayout> QuadNodeArray;
  typedef A2D::MultiArray<T, QuadSolnLayout> QuadSolnArray;
  typedef A2D::MultiArray<T, QuadGradLayout> QuadGradArray;
  typedef A2D::MultiArray<T, QuadDetLayout> QuadDetArray;
  typedef A2D::MultiArray<T, QuadJtransLayout> QuadJtransArray;
  typedef A2D::MultiArray<T, QuadDataLayout> QuadDataArray;

  // Residual data
  typedef A2D::MultiArray<T, ElemResLayout> ElemResArray;
  typedef A2D::MultiArray<T, ElemJacLayout> ElemJacArray;

  ElementBasis(const index_t nelems)
      : nelems(nelems),
        conn_layout(nelems),
        elem_node_layout(nelems),
        elem_soln_layout(nelems),
        quad_node_layout(nelems),
        quad_soln_layout(nelems),
        quad_grad_layout(nelems),
        quad_detJ_layout(nelems),
        quad_jtrans_layout(nelems),
        quad_data_layout(nelems),
        elem_res_layout(nelems),
        elem_jac_layout(nelems) {
    Timer t("ElementBasis::ElementBasis(1)");
    conn = ConnArray(conn_layout);
    Xe = ElemNodeArray(elem_node_layout);
    Ue = ElemSolnArray(elem_soln_layout);
    Xq = QuadNodeArray(quad_node_layout);
    Uq = QuadSolnArray(quad_soln_layout);
    Uxi = QuadGradArray(quad_grad_layout);
    detJ = QuadDetArray(quad_detJ_layout);
    Jinv = QuadJtransArray(quad_jtrans_layout);
    data = QuadDataArray(quad_data_layout);
  }

  template <typename IdxType>
  ElementBasis(const index_t nelems, const IdxType conn_[])
      : nelems(nelems),
        conn_layout(nelems),
        elem_node_layout(nelems),
        elem_soln_layout(nelems),
        quad_node_layout(nelems),
        quad_soln_layout(nelems),
        quad_grad_layout(nelems),
        quad_detJ_layout(nelems),
        quad_jtrans_layout(nelems),
        quad_data_layout(nelems),
        elem_res_layout(nelems),
        elem_jac_layout(nelems) {
    Timer t("ElementBasis::ElementBasis(2)");
    conn = ConnArray(conn_layout);
    Xe = ElemNodeArray(elem_node_layout);
    Ue = ElemSolnArray(elem_soln_layout);
    Xq = QuadNodeArray(quad_node_layout);
    Uq = QuadSolnArray(quad_soln_layout);
    Uxi = QuadGradArray(quad_grad_layout);
    detJ = QuadDetArray(quad_detJ_layout);
    Jinv = QuadJtransArray(quad_jtrans_layout);
    data = QuadDataArray(quad_data_layout);
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

  // Get the layout associated with the quadrature point data
  QuadDataLayout& get_quad_data_layout() { return quad_data_layout; }

  // Get the layout associated with the element residuals
  ElemResLayout& get_elem_res_layout() { return elem_res_layout; }

  // Get the layout associated with the element jacobians
  ElemJacLayout& get_elem_jac_layout() { return elem_jac_layout; }

  // Add the node set to the connectivity
  void add_node_set(
      std::unordered_set<std::pair<I, I>, pair_hash_fun<I>>& node_set) {
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

  // Get the layouts
  ElemNodeLayout& get_elem_node_layout() { return elem_node_layout; }
  QuadNodeLayout& get_quad_node_layout() { return quad_node_layout; }
  QuadDetLayout& get_detJ_layout() { return quad_detJ_layout; }
  QuadJtransLayout& get_jinv_layout() { return quad_jtrans_layout; }
  ElemSolnLayout& get_elem_solution_layout() { return elem_soln_layout; }
  QuadSolnLayout& get_quad_solution_layout() { return quad_soln_layout; }
  QuadGradLayout& get_quad_gradient_layout() { return quad_grad_layout; }

 private:
  // Connectivity layout
  ConnLayout conn_layout;

  // Element-node layouts
  ElemNodeLayout elem_node_layout;
  ElemSolnLayout elem_soln_layout;

  // Element-quadrature layouts
  QuadNodeLayout quad_node_layout;
  QuadSolnLayout quad_soln_layout;
  QuadGradLayout quad_grad_layout;
  QuadDetLayout quad_detJ_layout;
  QuadJtransLayout quad_jtrans_layout;
  QuadDataLayout quad_data_layout;

  // Layout instances for the residuals/jacobians
  ElemResLayout elem_res_layout;
  ElemJacLayout elem_jac_layout;

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