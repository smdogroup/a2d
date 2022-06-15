#ifndef A2D_MODEL_H
#define A2D_MODEL_H

#include <list>

#include "a2dtmp.h"
#include "multiarray.h"
#include "sparse_amg.h"
#include "sparse_matrix.h"
#include "sparse_numeric.h"
#include "sparse_symbolic.h"

namespace A2D {

/*
  Element base class.

  This defines an element that is compatible with the given PDE. You can
  set the node locations into the element, add the non-zero-pattern of the
  Jacobian matrix via the add_node_set as well as adding the residual and
  Jacobian contributions
*/
template <typename I, typename T, class PDE>
class Element {
 public:
  virtual ~Element() {}
  virtual void set_nodes(typename PDE::NodeArray& X) = 0;
  virtual void add_node_set(std::set<std::pair<I, I>>& node_set) = 0;
  virtual void set_solution(typename PDE::SolutionArray& U) = 0;
  virtual T energy() { return T(0.0); }
  virtual void add_residual(typename PDE::SolutionArray& res) = 0;
  virtual void add_jacobian(typename PDE::SparseMat& jac) = 0;
  virtual void add_adjoint_dfdnodes(typename PDE::SolutionArray& psi,
                                    typename PDE::NodeArray& dfdx) {}
};

/*
  Constitutive base class.

  This is an optional layer of classes that define a relationship between
  the design variables and the material data stored in each class.

  The constitutive object implementation is responsible for storing the
  associated element and defining that relationship for a given type of PDE.
*/
template <typename I, typename T, class PDE>
class Constitutive {
 public:
  virtual ~Constitutive() {}
  virtual void set_design_vars(typename PDE::DesignArray& x) = 0;
  virtual void add_dfdx(typename PDE::DesignArray& dfdx) = 0;
  virtual void add_adjoint_dfdx(typename PDE::SolutionArray& psi,
                                typename PDE::NodeArray& dfdx) = 0;
};

/*
  ElementFunction base class.

  This class enables the computation of a function that sums up over all
  the elements within the element.
*/
template <typename I, typename T, class PDE>
class ElementFunctional {
 public:
  virtual ~ElementFunctional() {}
  virtual T eval_functional() = 0;
  virtual void add_dfdu(typename PDE::SolutionArray& dfdu) = 0;
  virtual void add_dfdx(typename PDE::NodeArray& dfdx) = 0;
};

/*
  Container class for all functionals.

  The functional must be the result of a sum over the ElementFunctions in the
  container.
*/
template <typename I, typename T, class PDE>
class Functional {
 public:
  Functional() { num_functionals = 0; }
  virtual ~Functional() {}

  void add_functional(ElementFunctional<I, T, PDE>* functional) {
    functionals[num_functionals] = functional;
    num_functionals++;
  }

  index_t get_num_functionals() { return num_functionals; }
  ElementFunctional<I, T, PDE>* get_functional(index_t index) {
    return functionals[index];
  }

  /*
    Evaluate the functional by summing the values over all components
  */
  T eval_functional(Functional& functional) {
    T value = 0.0;
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      value += functional.eval_functional();
    }
  }

  /*
    Compute the derivative of the functional w.r.t. state variables
  */
  void eval_dfdu(Functional& functional, typename PDE::SolutionArray& dfdu) {
    dfdu.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdu(dfdu);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. design variables
  */
  void eval_dfdx(Functional& functional, typename PDE::DesignArray& dfdx) {
    dfdx.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdx(dfdx);
    }
  }

  /*
    Compute the derivative of the functional w.r.t. nodes
  */
  void eval_dfdnodes(Functional& functional, typename PDE::NodeArray& dfdx) {
    dfdx.zero();
    for (index_t index = 0; index < functional.get_num_functionals(); index++) {
      ElementFunctional<I, T, PDE>* elem = functional.get_functional(index);
      functional.add_dfdx(dfdx);
    }
  }

 private:
  index_t num_functionals;
  ElementFunctional<I, T, PDE>* functionals[10];
};

/*
  The FE Model base class.

  This class holds all the elements and constitutive objects in the
  model. It is used to compute the residual, Jacobian and derivatives
  needed for adjoint-based gradient evaluation
*/
template <typename I, typename T, class PDE>
class FEModel {
 public:
  FEModel(const index_t nnodes, const index_t nbcs)
      : nnodes(nnodes),
        nbcs(nbcs),
        bcs_layout(nbcs),
        node_layout(nnodes),
        solution_layout(nnodes),
        null_space_layout(nnodes),
        bcs(bcs_layout),
        X(node_layout),
        U(solution_layout),
        B(null_space_layout) {
    num_elements = 0;
    num_constitutive = 0;
  }
  ~FEModel() {}

  /*
    Set the design variables into the elements -

    The constitutive objects serve as a mapping from the design variables
    x to the material parameter data stored in each element array.
  */
  void set_design_vars(typename PDE::DesignArray& x) {
    for (I i = 0; i < num_constitutive; i++) {
      constitutive[i]->set_design_vars(x);
    }
  }

  // void add_element(Element<I, T, PDE>* element) {
  // elements.push_back(element); }

  void add_constitutive(Constitutive<I, T, PDE>* con) {
    constitutive[num_constitutive] = con;
    num_constitutive++;
  }

  void add_element(Element<I, T, PDE>* element) {
    elements[num_elements] = element;
    num_elements++;
  }

  // Create a new solution vector
  typename PDE::SolutionArray* new_solution() {
    return new typename PDE::SolutionArray(solution_layout);
  }

  typename PDE::NodeArray& get_nodes() { return X; }
  typename PDE::BCsArray& get_bcs() { return bcs; }

  void set_nodes(typename PDE::NodeArray& X) {
    // for (auto it = elements.begin(); it != elements.end(); it++) {
    //   it->set_nodes(X);
    // }
    for (I i = 0; i < num_elements; i++) {
      elements[i]->set_nodes(X);
    }
  }

  void set_solution(typename PDE::SolutionArray& U) {
    // for (auto it = elements.begin(); it != elements.end(); it++) {
    //   it->set_solution(U);
    // }
    for (I i = 0; i < num_elements; i++) {
      elements[i]->set_solution(U);
    }
  }

  // Zero the dirichlet boundary conditions
  void zero_bcs(typename PDE::SolutionArray& U) { A2D::VecZeroBCRows(bcs, U); }

  T energy() {
    T value = 0.0;
    for (I i = 0; i < num_elements; i++) {
      value += elements[i]->energy();
    }
    return value;
  }

  void residual(typename PDE::SolutionArray& res) {
    res.zero();
    // for (auto it = elements.begin(); it != elements.end(); it++) {
    //   it->add_residual(res);
    // }
    for (I i = 0; i < num_elements; i++) {
      elements[i]->add_residual(res);
    }
    A2D::VecZeroBCRows(bcs, res);
  }

  void jacobian(typename PDE::SparseMat& jac) {
    jac.zero();
    // for (auto it = elements.begin(); it != elements.end(); it++) {
    //   it->add_jacobian(jac);
    // }
    for (I i = 0; i < num_elements; i++) {
      elements[i]->add_jacobian(jac);
    }
    A2D::BSRMatZeroBCRows(bcs, jac);
  }

  // Get a new matrix
  typename PDE::SparseMat* new_matrix() {
    std::set<std::pair<I, I>> node_set;
    // for (auto it = elements.begin(); it != elements.end(); it++) {
    //   it->add_node_set(node_set);
    // }
    for (I i = 0; i < num_elements; i++) {
      elements[i]->add_node_set(node_set);
    }
    return A2D::BSRMatFromNodeSet<I, T, PDE::vars_per_node>(nnodes, node_set);
  }

  // With a matrix, create a preconditioner. Note that the entries
  // in the matrix must be filled at this point, e.g. after a call to
  // add_jacobian
  typename PDE::SparseAmg* new_amg(int num_levels, double omega,
                                   typename PDE::SparseMat* mat,
                                   bool print_info = false) {
    PDE::compute_null_space(X, B);
    A2D::VecZeroBCRows(bcs, B);
    return new typename PDE::SparseAmg(num_levels, omega, mat, &B, print_info);
  }

 private:
  // std::list<Element<I, T, PDE>*> elements;
  index_t num_elements;
  Element<I, T, PDE>* elements[10];

  index_t num_constitutive;
  Constitutive<I, T, PDE>* constitutive;

  const index_t nnodes;  // Number of nodes in the model
  const index_t nbcs;    // Number of nodes with Dirichlet bcs

  typename PDE::BCsLayout bcs_layout;
  typename PDE::NodeLayout node_layout;
  typename PDE::SolutionLayout solution_layout;
  typename PDE::NullSpaceLayout null_space_layout;

  typename PDE::BCsArray bcs;
  typename PDE::NodeArray X;
  typename PDE::SolutionArray U;
  typename PDE::NullSpaceArray B;
};

/*
  A partial implementation of the element class with a prescribed basis function
*/
template <typename I, typename T, class PDE, class Basis>
class ElementBasis : public Element<I, T, PDE> {
 public:
  static const index_t spatial_dim = PDE::spatial_dim;
  static const index_t vars_per_node = PDE::vars_per_node;
  static const index_t data_per_point = PDE::data_per_point;

  static const index_t nodes_per_elem = Basis::NUM_NODES;
  static const index_t quad_pts_per_elem = Basis::quadrature::NUM_QUAD_PTS;

  // Connectivity layout
  typedef A2D::CLayout<nodes_per_elem> ConnLayout;

  // Element-node layouts
  typedef A2D::CLayout<nodes_per_elem, spatial_dim> ElemNodeLayout;
  typedef A2D::CLayout<nodes_per_elem, vars_per_node> ElemSolnLayout;

  // Element-quadrature layouts
  typedef A2D::CLayout<quad_pts_per_elem, spatial_dim> QuadNodeLayout;
  typedef A2D::CLayout<quad_pts_per_elem, vars_per_node> QuadSolnLayout;
  typedef A2D::CLayout<quad_pts_per_elem, vars_per_node, spatial_dim>
      QuadGradLayout;
  typedef A2D::CLayout<quad_pts_per_elem> QuadDetLayout;
  typedef A2D::CLayout<quad_pts_per_elem, spatial_dim, spatial_dim>
      QuadJtransLayout;
  typedef A2D::CLayout<quad_pts_per_elem, data_per_point> QuadDataLayout;

  // Residual/Jacobian layouts
  typedef A2D::CLayout<nodes_per_elem, vars_per_node> ElemResLayout;
  typedef A2D::CLayout<nodes_per_elem, nodes_per_elem, vars_per_node,
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
        elem_jac_layout(nelems),
        conn(conn_layout),
        Xe(elem_node_layout),
        Ue(elem_soln_layout),
        Xq(quad_node_layout),
        Uq(quad_soln_layout),
        Uxi(quad_grad_layout),
        detJ(quad_detJ_layout),
        Jinv(quad_jtrans_layout),
        data(quad_data_layout) {}
  virtual ~ElementBasis() {}

  const index_t nelems;  // Number of elements in the model

  // Get the connectivity array
  ConnArray& get_conn() { return conn; }

  // Get the data associated with the quadrature points
  QuadDataArray& get_quad_data() { return data; }

  // Get the layout associated with the element residuals
  ElemResLayout& get_elem_res_layout() { return elem_res_layout; }

  // Get the layout associated with the element jacobians
  ElemJacLayout& get_elem_jac_layout() { return elem_jac_layout; }

  // Add the node set to the connectivity
  void add_node_set(std::set<std::pair<I, I>>& node_set) {
    BSRMatAddConnectivity(conn, node_set);
  }

  // Set the nodes from the array
  void set_nodes(typename PDE::NodeArray& X) {
    VecElementScatter(conn, X, Xe);
    Basis::template interp<spatial_dim>(Xe, Xq);
    Basis::template compute_jtrans<T>(Xe, detJ, Jinv);
  }

  // Set the solution
  void set_solution(typename PDE::SolutionArray& U) {
    VecElementScatter(conn, U, Ue);
    Basis::template interp<vars_per_node>(Ue, Uq);
    Basis::template gradient<T, vars_per_node>(Ue, Uxi);
  }

  // Pure virtual member functions that need an override
  // virtual void add_residual(typename PDE::SolutionArray& res) = 0;
  // virtual void add_jacobian(typename PDE::SparseMat& jac) = 0;

  // Expose the underlying element data
  ElemNodeArray& get_elem_nodes() { return Xe; }
  QuadNodeArray& get_quad_nodes() { return Xq; }
  QuadDetArray& get_detJ() { return detJ; }
  QuadJtransArray& get_Jinv() { return Jinv; }
  ElemSolnArray& get_elem_solution() { return Ue; }
  QuadSolnArray& get_quad_solution() { return Uq; }
  QuadGradArray& get_quad_gradient() { return Uxi; }

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

/*
  A minimal set of data needed for the ElementFunctional class
*/
template <typename I, typename T, class PDE, class Basis>
class FunctionalBasis : public ElementFunctional<I, T, PDE> {
 public:
  FunctionalBasis(ElementBasis<I, T, PDE, Basis>& element) : element(element) {}

  ElementBasis<I, T, PDE, Basis>& element;
};

}  // namespace A2D

#endif  // A2D_MODEL_H
