#ifndef A2D_MODEL_H
#define A2D_MODEL_H

#include <list>
#include <memory>

#include "a2dtmp3d.h"
#include "constitutive.h"
#include "element.h"
#include "multiarray.h"
#include "sparse/sparse_amg.h"
#include "sparse/sparse_matrix.h"
#include "sparse/sparse_numeric.h"
#include "sparse/sparse_symbolic.h"

namespace A2D {

/*
  The FE Model base class.

  This class holds all the elements and constitutive objects in the
  model. It is used to compute the residual, Jacobian and derivatives
  needed for adjoint-based gradient evaluation
*/
template <typename I, typename T, class PDEInfo>
class FEModel {
 public:
  static const index_t SPATIAL_DIM = PDEInfo::SPATIAL_DIM;
  FEModel(const index_t nnodes, const index_t nbcs)
      : nnodes(nnodes),
        nbcs(nbcs),
        bcs_layout(nbcs),
        node_layout(nnodes),
        solution_layout(nnodes),
        null_space_layout(nnodes),
        bcs(bcs_layout),
        X(node_layout),
        U(solution_layout) {
    B = std::make_shared<typename PDEInfo::NullSpaceArray>(null_space_layout);
  }
  template <typename Ttype, typename IdxType>
  FEModel(const index_t nnodes, const Ttype X_[], const index_t nbcs,
          const IdxType bcs_[])
      : nnodes(nnodes),
        nbcs(nbcs),
        bcs_layout(nbcs),
        node_layout(nnodes),
        solution_layout(nnodes),
        null_space_layout(nnodes),
        bcs(bcs_layout),
        X(node_layout),
        U(solution_layout) {
    B = std::make_shared<typename PDEInfo::NullSpaceArray>(null_space_layout);
    // Copy the x values
    for (I i = 0; i < nnodes; i++) {
      for (I j = 0; j < SPATIAL_DIM; j++) {
        X(i, j) = X_[SPATIAL_DIM * i + j];
      }
    }

    // Copy the bcs values
    for (I i = 0; i < nbcs; i++) {
      for (I j = 0; j < 2; j++) {
        bcs(i, j) = bcs_[2 * i + j];
      }
    }
  }
  ~FEModel() {}

  const index_t nnodes;  // Number of nodes in the model
  const index_t nbcs;    // Number of nodes with Dirichlet bcs

  /*
    Add an element object to the model
  */
  void add_element(std::shared_ptr<ElementBase<I, T, PDEInfo>> element) {
    elements.push_back(element);
  }

  /*
    Add a constitutive object to the model
  */
  void add_constitutive(std::shared_ptr<ConstitutiveBase<I, T, PDEInfo>> con) {
    constitutive.push_back(con);
  }

  /*
    Perform initialization tasks after nodes, connectivities and elements have
    been set into the model class
  */
  void init() {
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->set_nodes(X);
    }
  }

  /*
    Create a new solution vector
  */
  std::shared_ptr<typename PDEInfo::SolutionArray> new_solution() {
    return std::make_shared<typename PDEInfo::SolutionArray>(solution_layout);
  }

  /*
    Create a new node vector
  */
  std::shared_ptr<typename PDEInfo::NodeArray> new_nodes() {
    return std::make_shared<typename PDEInfo::NodeArray>(solution_layout);
  }

  /*
    Get the node locations
  */
  typename PDEInfo::NodeArray& get_nodes() { return X; }

  /*
    Get the boundary conditions
  */
  typename PDEInfo::BCsArray& get_bcs() { return bcs; }

  /*
    Get the solution
  */
  typename PDEInfo::SolutionArray& get_solution() { return U; }

  /*
    Set new node locations for each of the elements
  */
  void set_nodes(std::shared_ptr<typename PDEInfo::NodeArray> Xnew) {
    X.copy(*Xnew);
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->set_nodes(X);
    }
  }

  /*
    Set the solution into the vector
  */
  void set_solution(std::shared_ptr<typename PDEInfo::SolutionArray> Unew) {
    U.copy(*Unew);
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->set_solution(U);
    }
  }

  /*
    Zero the dirichlet boundary conditions in the vector
  */
  void zero_bcs(std::shared_ptr<typename PDEInfo::SolutionArray> U0) {
    A2D::VecZeroBCRows(bcs, *U0);
  }

  /*
    Compute the energy from all the elements, if they define an energy
    functional
  */
  T energy() {
    T value = 0.0;
    for (auto it = elements.begin(); it != elements.end(); it++) {
      value += (*it)->energy();
    }
    return value;
  }

  /*
    Compute the residual
  */
  void residual(std::shared_ptr<typename PDEInfo::SolutionArray> res) {
    res->zero();
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->add_residual(*res);
    }
    A2D::VecZeroBCRows(bcs, *res);
  }

  /*
    Compute the Jacobian matrix
  */
  void jacobian(std::shared_ptr<typename PDEInfo::SparseMat> jac) {
    jac->zero();
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->add_jacobian(*jac);
    }
    A2D::BSRMatZeroBCRows(bcs, *jac);
  }

  /*
    Set the design variables
  */
  void set_design_vars(std::shared_ptr<typename PDEInfo::DesignArray> x) {
    for (auto it = constitutive.begin(); it != constitutive.end(); it++) {
      (*it)->set_design_vars(*x);
    }
  }

  /*
    Add the derivative of the adjoint-residual product
  */
  void add_adjoint_dfdx(std::shared_ptr<typename PDEInfo::SolutionArray> psi,
                        std::shared_ptr<typename PDEInfo::DesignArray> dfdx) {
    for (auto it = constitutive.begin(); it != constitutive.end(); it++) {
      (*it)->add_adjoint_dfdx(*psi, *dfdx);
    }
  }

  /*
    Create a new matrix
  */
  std::shared_ptr<typename PDEInfo::SparseMat> new_matrix() {
    std::set<std::pair<I, I>> node_set;
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->add_node_set(node_set);
    }
    return std::shared_ptr<typename PDEInfo::SparseMat>(
        A2D::BSRMatFromNodeSet<I, T, PDEInfo::vars_per_node>(nnodes, node_set));
  }

  // With a matrix, create a preconditioner. Note that the entries
  // in the matrix must be filled at this point, e.g. after a call to
  // add_jacobian
  std::shared_ptr<typename PDEInfo::SparseAmg> new_amg(
      int num_levels, double omega, double epsilon,
      std::shared_ptr<typename PDEInfo::SparseMat> mat,
      bool print_info = false) {
    PDEInfo::compute_null_space(X, *B);
    A2D::VecZeroBCRows(bcs, *B);
    return std::make_shared<typename PDEInfo::SparseAmg>(
        num_levels, omega, epsilon, mat, B, print_info);
  }

 private:
  std::list<std::shared_ptr<ElementBase<I, T, PDEInfo>>> elements;
  std::list<std::shared_ptr<ConstitutiveBase<I, T, PDEInfo>>> constitutive;

  typename PDEInfo::BCsLayout bcs_layout;
  typename PDEInfo::NodeLayout node_layout;
  typename PDEInfo::SolutionLayout solution_layout;
  typename PDEInfo::NullSpaceLayout null_space_layout;

  typename PDEInfo::BCsArray bcs;
  typename PDEInfo::NodeArray X;
  typename PDEInfo::SolutionArray U;
  std::shared_ptr<typename PDEInfo::NullSpaceArray> B;
};

}  // namespace A2D

#endif  // A2D_MODEL_H
