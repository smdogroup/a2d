#ifndef A2D_MODEL_H
#define A2D_MODEL_H

#include "a2dtmp.h"
#include "multiarray.h"

/*
  Scatter the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray, class NodeArray, class ElementArray>
void element_scatter(ConnArray& conn, NodeArray& X, ElementArray& Xe) {
  // Loop over the elements
  for (std::size_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (std::size_t j = 0; j < conn.extent(1); j++) {
      const std::size_t index = conn(i, j);

      // Loop over the variables
      for (std::size_t k = 0; k < X.extent(1); k++) {
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
  for (std::size_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (std::size_t j = 0; j < conn.extent(1); j++) {
      const std::size_t index = conn(i, j);

      // Loop over the variables
      for (std::size_t k = 0; k < X.extent(1); k++) {
        X(index, k) += Xe(i, j, k);
      }
    }
  }
}

template <class IdxType, class ScalarType, class Basis, int vars_per_node,
          int data_per_point>
class PDEModel {
 public:
  PDEModel(const std::size_t nelems, const std::size_t nnodes)
      : nelems(nelems),
        nnodes(nnodes),
        conn_layout(nelems),
        node_layout(nnodes),
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

  const int nelems;  // Number of elements in the model
  const int nnodes;  // Number of nodes in the model

  // static const int vars_per_node = Model::NUM_VARS;
  // static const int data_per_point = Model::NUM_DATA;

  static const int spatial_dim = Basis::SPATIAL_DIM;
  static const int nodes_per_elem = Basis::NUM_NODES;
  static const int quad_pts_per_elem = Basis::quadrature::NUM_QUAD_PTS;

  // Basic layouts
  typedef A2D::CLayout<nodes_per_elem> ConnLayout;
  typedef A2D::CLayout<spatial_dim> NodeLayout;
  typedef A2D::CLayout<vars_per_node> SolutionLayout;

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

  // Typedef the array types
  typedef A2D::MultiArray<IdxType, ConnLayout> ConnArray;

  typedef A2D::MultiArray<ScalarType, NodeLayout> NodeArray;
  typedef A2D::MultiArray<ScalarType, ElemNodeLayout> ElemNodeArray;
  typedef A2D::MultiArray<ScalarType, QuadNodeLayout> QuadNodeArray;
  typedef A2D::MultiArray<ScalarType, QuadDetLayout> QuadDetArray;
  typedef A2D::MultiArray<ScalarType, QuadJtransLayout> QuadJtransArray;

  // Solution information
  typedef A2D::MultiArray<ScalarType, SolutionLayout> SolutionArray;
  typedef A2D::MultiArray<ScalarType, ElemSolnLayout> ElemSolnArray;
  typedef A2D::MultiArray<ScalarType, QuadSolnLayout> QuadSolnArray;
  typedef A2D::MultiArray<ScalarType, QuadGradLayout> QuadGradArray;

  // Material data
  typedef A2D::MultiArray<ScalarType, QuadDataLayout> QuadDataArray;

  // Residual data
  typedef A2D::MultiArray<ScalarType, ResLayout> ElemResArray;
  typedef A2D::MultiArray<ScalarType, JacLayout> ElemJacArray;

  SolutionArray* new_solution() { return new SolutionArray(solution_layout); }
  ElemSolnArray* new_elem_solution() {
    return new ElemSolnArray(elem_soln_layout);
  }

  ConnArray& get_conn() { return *conn_; }
  NodeArray& get_nodes() { return *X_; }
  void reset_nodes() {
    element_scatter(*conn_, *X_, *Xe_);
    Basis::template interp<spatial_dim>(*Xe_, *Xq_);
    Basis::template compute_jtrans<ScalarType>(*Xe_, *detJ_, *Jinv_);
  }

  SolutionArray& get_solution() { return *U_; }
  void reset_solution() {
    element_scatter(*conn_, *U_, *Ue_);
    Basis::template interp<vars_per_node>(*Ue_, *Uq_);
    Basis::template gradient<ScalarType, vars_per_node>(*Ue_, *Uxi_);
  }

  QuadDataArray& get_quad_data() { return *data_; }
  ElemResArray& get_elem_res() { return *res_; }
  ElemJacArray& get_elem_jac() { return *jac_; }

 protected:
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

#endif  // A2D_MODEL_H