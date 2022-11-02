#ifndef A2D_MODEL_H
#define A2D_MODEL_H

#include <list>
#include <memory>

#include "febase.h"

namespace A2D {

/**
 * @brief The Finite-element model base class. This class holds all the elements
 * and constitutive objects in the model. It is used to compute the residual,
 * Jacobian and derivatives needed for adjoint-based gradient evaluation
 */
template <typename T>
class FEModel {
 public:
  // static const index_t SPATIAL_DIM = PDEInfo::SPATIAL_DIM;

  // FEModel(const index_t nnodes, const index_t nbcs)
  //     : nnodes(nnodes), nbcs(nbcs) {
  //   Timer t("FEModel::FEModel(1)");
  //   bcs = typename PDEInfo::BCsArray("bcs", nbcs);
  //   X = typename PDEInfo::NodeArray("X", nnodes);
  //   U = typename PDEInfo::SolutionArray("U", nnodes);
  //   B = typename PDEInfo::NullSpaceArray("B", nnodes);
  // }
  // template <typename Ttype, typename IdxType>
  // FEModel(const index_t nnodes, const Ttype X_[], const index_t nbcs,
  //         const IdxType bcs_[])
  //     : nnodes(nnodes), nbcs(nbcs) {
  //   Timer t("FEModel::FEModel(2)");
  //   bcs = typename PDEInfo::BCsArray("bcs", nbcs);
  //   X = typename PDEInfo::NodeArray("X", nnodes);
  //   U = typename PDEInfo::SolutionArray("U", nnodes);
  //   B = typename PDEInfo::NullSpaceArray("B", nnodes);
  //   // Copy the x values
  //   for (I i = 0; i < nnodes; i++) {
  //     for (I j = 0; j < SPATIAL_DIM; j++) {
  //       X(i, j) = X_[SPATIAL_DIM * i + j];
  //     }
  //   }

  //   // Copy the bcs values
  //   for (I i = 0; i < nbcs; i++) {
  //     for (I j = 0; j < 2; j++) {
  //       bcs(i, j) = bcs_[2 * i + j];
  //     }
  //   }
  // }
  ~FEModel() {}

  /**
   * @brief Add an element object to the model
   *
   * @param element The element object
   */
  void add_element(std::shared_ptr<ElementBase<T>> element) {
    elements.push_back(element);
  }

  /**
   * @brief  Add a constitutive object to the model
   *
   * @param con The constitutive object
   */
  void add_constitutive(std::shared_ptr<ConstitutiveBase<T>> con) {
    constitutive.push_back(con);
  }

  /**
   * @brief Perform initialization tasks after nodes, connectivities and
   * elements have been set into the model class
   */
  void init() {
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->set_geo(X);
    }
  }

  /**
   * @brief Create a new solution vector
   */
  std::shared_ptr<typename PDEInfo::SolutionArray> new_solution() {
    Timer timer("FEModel::new_solution()");
    return std::make_shared<typename PDEInfo::SolutionArray>("new_solution",
                                                             nnodes);
  }

  /**
   * @brief Create a new node vector
   */
  std::shared_ptr<typename PDEInfo::NodeArray> new_nodes() {
    return std::make_shared<typename PDEInfo::NodeArray>("new_nodes", nnodes);
  }

  // /**
  //  * @brief Get the node locations
  //  */
  // typename PDEInfo::NodeArray& get_nodes() { return X; }

  // /**
  //  * @brief Get the bcs object
  //  */
  // typename PDEInfo::BCsArray& get_bcs() { return bcs; }

  // /**
  //  * @brief Get the solution
  //  */
  // typename PDEInfo::SolutionArray& get_solution() { return U; }

  /**
   * @brief Set node locations for each of the elements
   *
   * @param Xnew
   */
  void set_nodes(std::shared_ptr<SolutionVector<T>> X) {
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->set_nodes(X);
    }
  }

  /*
    Set the solution into the vector
  */
  void set_solution(std::shared_ptr<typename PDEInfo::SolutionArray> Unew) {
    A2D::BLAS::copy(U, *Unew);
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
    Compute the residual
  */
  void residual(std::shared_ptr<typename PDEInfo::SolutionArray> res) {
    A2D::BLAS::zero(*res);
    for (auto it = elements.begin(); it != elements.end(); it++) {
      (*it)->add_residual(*res);
    }
    A2D::VecZeroBCRows(bcs, *res);
  }

  // /*
  //   Compute the Jacobian matrix
  // */
  // void jacobian(std::shared_ptr<typename PDEInfo::SparseMat> jac) {
  //   Timer timer("FEModel::jacobian()");
  //   jac->zero();
  //   for (auto it = elements.begin(); it != elements.end(); it++) {
  //     (*it)->add_jacobian(*jac);
  //   }
  //   A2D::BSRMatZeroBCRows(bcs, *jac);
  // }

  // /*
  //   Set the design variables
  // */
  // void set_design_vars(std::shared_ptr<typename PDEInfo::DesignArray> x) {
  //   for (auto it = constitutive.begin(); it != constitutive.end(); it++) {
  //     (*it)->set_design_vars(*x);
  //   }
  // }

  // /*
  //   Add the derivative of the adjoint-residual product
  // */
  // void add_adjoint_dfdx(std::shared_ptr<typename PDEInfo::SolutionArray> psi,
  //                       std::shared_ptr<typename PDEInfo::DesignArray> dfdx)
  //                       {
  //   for (auto it = constitutive.begin(); it != constitutive.end(); it++) {
  //     (*it)->add_adjoint_dfdx(*psi, *dfdx);
  //   }
  // }

  // /*
  //   Create a new matrix
  // */
  // std::shared_ptr<typename PDEInfo::SparseMat> new_matrix() {
  //   Timer t1("FEModel::new_matrix()");
  //   Kokkos::UnorderedMap<COO<I>, void> node_set;

  //   for (auto it = elements.begin(); it != elements.end(); it++) {
  //     (*it)->add_node_set(node_set);
  //   }
  //   return std::shared_ptr<typename PDEInfo::SparseMat>(
  //       A2D::BSRMatFromNodeSet<I, T, PDEInfo::vars_per_node>(nnodes,
  //       node_set));
  // }

  // // With a matrix, create a preconditioner. Note that the entries
  // // in the matrix must be filled at this point, e.g. after a call to
  // // add_jacobian
  // std::shared_ptr<typename PDEInfo::SparseAmg> new_amg(
  //     int num_levels, double omega, double epsilon,
  //     std::shared_ptr<typename PDEInfo::SparseMat> mat,
  //     bool print_info = false) {
  //   Timer timer("FEModel::new_amg()");
  //   PDEInfo::compute_null_space(X, B);
  //   A2D::VecZeroBCRows(bcs, B);
  //   return std::make_shared<typename PDEInfo::SparseAmg>(
  //       num_levels, omega, epsilon, mat, B, print_info);
  // }

 private:
  std::list<std::shared_ptr<ElementBase<T>>> elements;
  std::list<std::shared_ptr<ConstitutiveBase<T>>> constitutive;

  // typename PDEInfo::BCsArray bcs;
  // X;
  // typename PDEInfo::SolutionArray U;
  // typename PDEInfo::NullSpaceArray B;
};

}  // namespace A2D

#endif  // A2D_MODEL_H