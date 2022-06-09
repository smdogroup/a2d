#ifndef A2D_MODEL_H
#define A2D_MODEL_H

#include "a2dtmp.h"
#include "multiarray.h"
#include "sparse_amg.h"
#include "sparse_matrix.h"
#include "sparse_numeric.h"
#include "sparse_symbolic.h"

namespace A2D {

/*
  Scatter the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray, class NodeArray, class ElementArray>
void element_scatter(ConnArray& conn, NodeArray& X, ElementArray& Xe) {
  // Loop over the elements
  for (A2D::index_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (A2D::index_t j = 0; j < conn.extent(1); j++) {
      const A2D::index_t index = conn(i, j);

      // Loop over the variables
      for (A2D::index_t k = 0; k < X.extent(1); k++) {
        Xe(i, j, k) = X(index, k);
      }
    }
  }
}

/*
  Gather the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray, class ElementArray, class NodeArray>
void element_gather_add(ConnArray& conn, ElementArray& Xe, NodeArray& X) {
  // Loop over the elements
  for (A2D::index_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (A2D::index_t j = 0; j < conn.extent(1); j++) {
      const A2D::index_t index = conn(i, j);

      // Loop over the variables
      for (A2D::index_t k = 0; k < X.extent(1); k++) {
        X(index, k) += Xe(i, j, k);
      }
    }
  }
}

template <typename I, typename T, class Basis, int vars_per_node,
          int data_per_point, int null_space_size>
class PDEModel {
 public:
  PDEModel(const A2D::index_t nelems, const A2D::index_t nnodes,
           const A2D::index_t nbcs)
      : nelems(nelems),
        nnodes(nnodes),
        nbcs(nbcs),
        conn_layout(nelems),
        node_layout(nnodes),
        bcs_layout(nbcs),
        null_space_layout(nnodes),
        solution_layout(nnodes),
        elem_node_layout(nelems),
        elem_soln_layout(nelems),
        quad_data_layout(nelems),
        quad_node_layout(nelems),
        quad_soln_layout(nelems),
        quad_grad_layout(nelems),
        quad_detJ_layout(nelems),
        quad_jtrans_layout(nelems),
        res_layout(nelems),
        jac_layout(nelems) {
    // Allocate the connectivity
    conn_ = new ConnArray(conn_layout);
    bcs_ = new BCsArray(bcs_layout);
    B_ = new NullSpaceArray(null_space_layout);

    // Node location information
    X_ = new NodeArray(node_layout);
    Xe_ = new ElemNodeArray(elem_node_layout);
    Xq_ = new QuadNodeArray(quad_node_layout);
    detJ_ = new QuadDetArray(quad_detJ_layout);
    Jinv_ = new QuadJtransArray(quad_jtrans_layout);

    // Solution information
    U_ = new SolutionArray(solution_layout);
    Ue_ = new ElemSolnArray(elem_soln_layout);
    Uq_ = new QuadSolnArray(quad_soln_layout);
    Uxi_ = new QuadGradArray(quad_grad_layout);

    // Material data
    data_ = new QuadDataArray(quad_data_layout);

    // Residual data
    res_ = new ElemResArray(res_layout);
    jac_ = new ElemJacArray(jac_layout);
  }
  ~PDEModel() {
    delete conn_;
    delete bcs_;
    delete B_;
    delete X_;
    delete Xe_;
    delete Xq_;
    delete detJ_;
    delete Jinv_;
    delete U_;
    delete Ue_;
    delete Uq_;
    delete Uxi_;
    delete data_;
    delete res_;
    delete jac_;
  }

  const index_t nelems;  // Number of elements in the model
  const index_t nnodes;  // Number of nodes in the model
  const index_t nbcs;    // Number of nodes with Dirichlet bcs

  static const index_t vpn = vars_per_node;
  static const index_t dpp = data_per_point;
  static const index_t nss = null_space_size;

  static const index_t spatial_dim = Basis::SPATIAL_DIM;
  static const index_t nodes_per_elem = Basis::NUM_NODES;
  static const index_t quad_pts_per_elem = Basis::quadrature::NUM_QUAD_PTS;

  // Basic layouts
  typedef A2D::CLayout<nodes_per_elem> ConnLayout;
  typedef A2D::CLayout<2> BCsLayout;
  typedef A2D::CLayout<spatial_dim> NodeLayout;
  typedef A2D::CLayout<vars_per_node> SolutionLayout;

  // Near null space layout - for the AMG preconditioner
  typedef A2D::CLayout<vars_per_node, null_space_size> NullSpaceLayout;

  // Layouts with element view
  typedef A2D::CLayout<nodes_per_elem, spatial_dim> ElemNodeLayout;
  typedef A2D::CLayout<nodes_per_elem, vars_per_node> ElemSolnLayout;

  // Layouts with quadrature-point view
  typedef A2D::CLayout<quad_pts_per_elem, data_per_point> QuadDataLayout;
  typedef A2D::CLayout<quad_pts_per_elem, spatial_dim> QuadNodeLayout;
  typedef A2D::CLayout<quad_pts_per_elem, vars_per_node> QuadSolnLayout;
  typedef A2D::CLayout<quad_pts_per_elem> QuadDetLayout;
  typedef A2D::CLayout<quad_pts_per_elem, spatial_dim, spatial_dim>
      QuadJtransLayout;
  typedef A2D::CLayout<quad_pts_per_elem, vars_per_node, spatial_dim>
      QuadGradLayout;

  // Residual/Jacobian layouts
  typedef A2D::CLayout<nodes_per_elem, vars_per_node> ResLayout;
  typedef A2D::CLayout<nodes_per_elem, nodes_per_elem, vars_per_node,
                       vars_per_node>
      JacLayout;

  // Jacobian matrix
  typedef A2D::BSRMat<I, T, vars_per_node, vars_per_node> SparseMat;

  // Sparse matrix multigrid type
  typedef A2D::BSRMatAmg<I, T, vars_per_node, null_space_size> SparseAmg;

  // Typedef the array types
  typedef A2D::MultiArray<I, ConnLayout> ConnArray;
  typedef A2D::MultiArray<I, BCsLayout> BCsArray;
  typedef A2D::MultiArray<T, NodeLayout> NodeArray;
  typedef A2D::MultiArray<T, ElemNodeLayout> ElemNodeArray;
  typedef A2D::MultiArray<T, NullSpaceLayout> NullSpaceArray;
  typedef A2D::MultiArray<T, QuadNodeLayout> QuadNodeArray;
  typedef A2D::MultiArray<T, QuadDetLayout> QuadDetArray;
  typedef A2D::MultiArray<T, QuadJtransLayout> QuadJtransArray;

  // Solution information
  typedef A2D::MultiArray<T, SolutionLayout> SolutionArray;
  typedef A2D::MultiArray<T, ElemSolnLayout> ElemSolnArray;
  typedef A2D::MultiArray<T, QuadSolnLayout> QuadSolnArray;
  typedef A2D::MultiArray<T, QuadGradLayout> QuadGradArray;

  // Material data
  typedef A2D::MultiArray<T, QuadDataLayout> QuadDataArray;

  // Residual data
  typedef A2D::MultiArray<T, ResLayout> ElemResArray;
  typedef A2D::MultiArray<T, JacLayout> ElemJacArray;

  // Create a new solution/residual vector
  SolutionArray* new_solution() { return new SolutionArray(solution_layout); }
  ElemSolnArray* new_elem_solution() {
    return new ElemSolnArray(elem_soln_layout);
  }

  ConnArray& get_conn() { return *conn_; }
  NodeArray& get_nodes() { return *X_; }
  BCsArray& get_bcs() { return *bcs_; }
  NullSpaceArray& get_null_space() { return *B_; }
  virtual void reset_nodes() {
    element_scatter(*conn_, *X_, *Xe_);
    Basis::template interp<spatial_dim>(*Xe_, *Xq_);
    Basis::template compute_jtrans<T>(*Xe_, *detJ_, *Jinv_);
  }

  SolutionArray& get_solution() { return *U_; }
  void reset_solution() {
    element_scatter(*conn_, *U_, *Ue_);
    Basis::template interp<vars_per_node>(*Ue_, *Uq_);
    Basis::template gradient<T, vars_per_node>(*Ue_, *Uxi_);
  }

  QuadDataArray& get_quad_data() { return *data_; }

  // Get a new matrix
  SparseMat* new_matrix() {
    return A2D::BSRMatFromConnectivity<I, T, vars_per_node>(*conn_);
  }

  // With a matrix, create a preconditioner. Note that the entries
  // in the matrix must be filled at this point, e.g. after a call to
  // add_jacobian
  SparseAmg* new_amg(int num_levels, double omega, SparseMat* mat,
                     bool print_info = false) {
    NullSpaceArray& B = get_null_space();
    return new SparseAmg(num_levels, omega, mat, &B, print_info);
  }

 protected:
  ElemResArray& get_elem_res() { return *res_; }
  ElemJacArray& get_elem_jac() { return *jac_; }
  ElemNodeArray& get_elem_nodes() { return *Xe_; }
  QuadNodeArray& get_quad_nodes() { return *Xq_; }
  QuadDetArray& get_detJ() { return *detJ_; }
  QuadJtransArray& get_Jinv() { return *Jinv_; }
  ElemSolnArray& get_elem_solution() { return *Ue_; }
  QuadSolnArray& get_quad_solution() { return *Uq_; }
  QuadGradArray& get_quad_gradient() { return *Uxi_; }

 private:
  // Basic layout instances
  ConnLayout conn_layout;
  NodeLayout node_layout;
  SolutionLayout solution_layout;
  BCsLayout bcs_layout;
  NullSpaceLayout null_space_layout;

  // Layout instances with element view
  ElemNodeLayout elem_node_layout;
  ElemSolnLayout elem_soln_layout;

  // Layout instances with quadrature-point view
  QuadDataLayout quad_data_layout;
  QuadNodeLayout quad_node_layout;
  QuadSolnLayout quad_soln_layout;
  QuadGradLayout quad_grad_layout;

  // Data derived from the nodes
  QuadDetLayout quad_detJ_layout;
  QuadJtransLayout quad_jtrans_layout;

  // Layout instances for the residuals/jacobians
  ResLayout res_layout;
  JacLayout jac_layout;

  // Instances of the multi-dimensional arrays
  ConnArray* conn_;
  BCsArray* bcs_;

  // Near null space array
  NullSpaceArray* B_;

  // Node location information
  NodeArray* X_;
  ElemNodeArray* Xe_;
  QuadNodeArray* Xq_;
  QuadDetArray* detJ_;
  QuadJtransArray* Jinv_;

  // Solution information
  SolutionArray* U_;
  ElemSolnArray* Ue_;
  QuadSolnArray* Uq_;
  QuadGradArray* Uxi_;

  // Material data
  QuadDataArray* data_;

  // Residual data
  ElemResArray* res_;
  ElemJacArray* jac_;
};

}  // namespace A2D

#endif  // A2D_MODEL_H
