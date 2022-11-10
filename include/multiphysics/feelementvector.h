
#ifndef A2D_FE_ELEMENT_VECTOR_H
#define A2D_FE_ELEMENT_VECTOR_H

#include "array.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"

namespace A2D {

/*
  This file contains objects that enable a vector-centric view of the
  finite-element degrees of freedom.

  The element vector class must implement the following:

  1. A lightweight object ElementVector::FEDof that is designed to access
  the degrees of freedom associated with an element, without exposing the
  details of the underlying element vector implementation. This class must
  contain the following members:

  1.1 The FEDof constructor must take the element index and a reference to
  the element vector object itself as arguments

  1.2 The FEDof must be indexable via operator[](const A2D::index_t)

  2. get_element_values(elem, dof)

  This function gets the element degrees of freedom associated with the index
  elem, and sets them into the FEDof object.  The dof object may store a local
  copy of them or just contain a pointer to the array of element values.

  2. add_element_values(elem, dof)

  This function adds the values to the element vector. This function may be
  empty if the implementation uses references to the values.

  3. init_values()

  Initialize the values from the solution vector into any temporary or local
  storage in the vector.

  4. init_zero_values()

  Initialize and zero values from the temporary or local storage only. Note
  that this does not zero values in the source vector.

  5. add_values()

  Add values into the source vector from any local storage. This may be an
  empty function if the values are stored directly.
*/

/*
  In-place element vector implementation
*/
template <typename T, class Basis>
class ElementVector_Serial {
 public:
  ElementVector_Serial(A2D::ElementMesh<Basis>& mesh, A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {}

  // Required DOF container object (different for each element vector
  // implementation)
  class FEDof {
   public:
    FEDof(A2D::index_t elem, ElementVector_Serial& elem_vec) {
      std::fill(dof, dof + Basis::ndof, T(0.0));
    }

    /**
     * @brief Get a reference to the underlying element data
     *
     * @return A reference to the degree of freedom
     */
    T& operator[](const A2D::index_t index) { return dof[index]; }
    const T& operator[](const A2D::index_t index) const { return dof[index]; }

   private:
    // Variables for all the basis functions
    T dof[Basis::ndof];
  };

  /**
   * @brief Get the number of elements
   */
  A2D::index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Initialize the element vector values
   *
   * This function may be called once before element values are accessed.
   */
  void init_values() {}
  /**
   * @brief Initialize any local element vector values to zero
   *
   * This function may be called before element values are added.
   */
  void init_zero_values() {}

  /**
   * @brief Finish adding values to the element vector
   *
   * Add any values from the element vector into the source vector.
   */
  void add_values() {}

  /**
   * @brief Get the element values from the object and store them in the FEDof
   *
   * @param elem the element index
   * @param dof the object that stores a reference to the degrees of freedom
   */
  // Get values for this element from the vector
  void get_element_values(A2D::index_t elem, FEDof& dof) {
    get_element_values_<Basis::nbasis - 1>(elem, dof);
  }

  /**
   * @brief Add the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   *
   * If FEDof contains a pointer to data, this function may do nothing
   */
  // Add values for this element to the vector
  void add_element_values(A2D::index_t elem, const FEDof& dof) {
    add_element_values_<Basis::nbasis - 1>(elem, dof);
  }

 private:
  template <A2D::index_t basis>
  void get_element_values_(A2D::index_t elem, FEDof& dof) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      dof[i + Basis::template get_dof_offset<basis>()] = sign * vec[dof_index];
    }
    if constexpr (basis > 0) {
      get_element_values_<basis - 1>(elem, dof);
    }
  }

  template <A2D::index_t basis>
  void add_element_values_(A2D::index_t elem, const FEDof& dof) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      vec[dof_index] += sign * dof[i + Basis::template get_dof_offset<basis>()];
    }
    if constexpr (basis > 0) {
      add_element_values_<basis - 1>(elem, dof);
    }
  }

  A2D::ElementMesh<Basis>& mesh;
  A2D::SolutionVector<T>& vec;
};

template <typename T, class Basis>
class ElementMat_Serial {
 public:
  ElementMat_Serial(A2D::ElementMesh<Basis>& mesh, A2D::BSRMat<T>& mat) : mesh(mesh), mat(mat) {}

  // Required DOF container object (different for each element vector
  // implementation)
  class FEMat {
    FEMat : FEDof(A2D::index_t elem, ElementMat_Serial& elem_mat) {
      std::fill(A, A + Basis::ndof * Basis::ndof, T(0.0));
    }

    /**
     * @brief Get a reference to the underlying element data
     *
     * @return A reference to the degree of freedom
     */
    T& operator()(const A2D::index_t i, const A2D::index_t j) { return A[i * Basis::ndof + j]; }
    const T& operator()(const A2D::index_t i, const A2D::index_t j) const {
      return A[i * Basis::ndof + j];
    }

   private:
    // Variables for all the basis functions
    T A[Basis::ndof * Basis::ndof];
  };

  /**
   * @brief Get the number of elements
   */
  A2D::index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Initialize the element vector values
   *
   * This function may be called once before element values are accessed.
   */
  void init_values() {}
  /**
   * @brief Initialize any local element vector values to zero
   *
   * This function may be called before element values are added.
   */
  void init_zero_values() {}

  /**
   * @brief Finish adding values to the element vector
   *
   * Add any values from the element vector into the source vector.
   */
  void add_values() {}

  /**
   * @brief Get the element values from the object and store them in the FEDof
   *
   * @param elem the element index
   * @param dof the object that stores a reference to the degrees of freedom
   */
  // Get values for this element from the vector
  void get_element_values(A2D::index_t elem, FEDof& dof) {
    get_element_values_<Basis::nbasis - 1>(elem, dof);
  }

  /**
   * @brief Add the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   *
   * If FEDof contains a pointer to data, this function may do nothing
   */
  // Add values for this element to the vector
  void add_element_values(A2D::index_t elem, const FEDof& dof) {
    add_element_values_<Basis::nbasis - 1>(elem, dof);
  }

 private:
  template <A2D::index_t basis>
  void get_element_values_(A2D::index_t elem, FEDof& dof) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      dof[i + Basis::template get_dof_offset<basis>()] = sign * vec[dof_index];
    }
    if constexpr (basis > 0) {
      get_element_values_<basis - 1>(elem, dof);
    }
  }

  template <A2D::index_t basis>
  void add_element_values_(A2D::index_t elem, const FEDof& dof) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      vec[dof_index] += sign * dof[i + Basis::template get_dof_offset<basis>()];
    }
    if constexpr (basis > 0) {
      add_element_values_<basis - 1>(elem, dof);
    }
  }

  A2D::ElementMesh<Basis>& mesh;
  A2D::SolutionVector<T>& vec;
};

/*
   This class allocates a heavy-weight 2-dimensional array to store
   (potentially duplicated) local degrees of freedom to realize a
   between-element parallelization.

   global to local dof population is done in parallel, local to global dof add
   is done by atomic operation to resolve write conflicts.
 */
template <typename T, class Basis>
class ElementVector_Parallel {
 public:
  using ElemVecArray_t = A2D::MultiArrayNew<T * [Basis::ndof]>;

  ElementVector_Parallel(A2D::ElementMesh<Basis>& mesh, A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec), elem_vec_array("elem_vec_array", mesh.get_num_elements()) {}

  // Required DOF container object (different for each element vector
  // implementation)
  class FEDof {
   public:
    FEDof(A2D::index_t elem, ElementVector_Parallel& elem_vec)
        : elem(elem), elem_vec_array(elem_vec.elem_vec_array) {}

    T& operator[](const int index) { return elem_vec_array(elem, index); }

    const T& operator[](const int index) const { return elem_vec_array(elem, index); }

   private:
    const A2D::index_t elem;
    ElemVecArray_t& elem_vec_array;
  };

  /**
   * @brief Get number of elements
   */
  A2D::index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Populate local dofs from global dof
   */
  void init_values() {
    for (A2D::index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
      _get_element_values<Basis::nbasis>(elem);
    }
  }

  /**
   * @brief Initialize local dof values to zero
   */
  void init_zero_values() { A2D::BLAS::zero(elem_vec_array); }

  /**
   * @brief Add local dof to global dof
   */
  void add_values() {
    for (A2D::index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
      _add_element_values<Basis::nbasis>(elem);
    }
  }

  /**
   * @brief Does nothing for this parallel implementation
   */
  void get_element_values(A2D::index_t elem, FEDof& dof) {}

  /**
   * @brief Does nothing for this parallel implementation
   */
  void add_element_values(A2D::index_t elem, const FEDof& dof) {}

 private:
  /**
   * @brief Populate local dof for a single element
   *
   * @tparam nbasis number of function spaces for the element, at least 1
   * @param elem_idx element index
   */
  template <A2D::index_t nbasis>
  void _get_element_values(const A2D::index_t& elem_idx) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<nbasis - 1>(); i++) {
      const int& sign = mesh.get_global_dof_sign(elem_idx, nbasis - 1, i);
      const A2D::index_t& dof_index = mesh.get_global_dof(elem_idx, nbasis - 1, i);
      const A2D::index_t dof_idx = i + Basis::template get_dof_offset<nbasis - 1>();
      elem_vec_array(elem_idx, dof_idx) = sign * vec[dof_index];
    }
    if constexpr (nbasis > 1) {
      _get_element_values<nbasis - 2>(elem_idx);
    }
    return;
  }

  /**
   * @brief Add local dof to global dof for a single element
   *
   * @tparam nbasis number of function spaces for the element, at least 1
   */
  template <A2D::index_t nbasis>
  void _add_element_values(const A2D::index_t& elem_idx) {
    for (A2D::index_t i = 0; i < Basis::template get_ndof<nbasis - 1>(); i++) {
      const int& sign = mesh.get_global_dof_sign(elem_idx, nbasis - 1, i);
      const A2D::index_t& dof_index = mesh.get_global_dof(elem_idx, nbasis - 1, i);
      const A2D::index_t dof_idx = i + Basis::template get_dof_offset<nbasis - 1>();
      Kokkos::atomic_add(&vec[dof_index], sign * elem_vec_array(elem_idx, dof_idx));
    }
    if constexpr (nbasis > 1) {
      _add_element_values<nbasis - 2>(elem_idx);
    }
    return;
  }

  A2D::ElementMesh<Basis>& mesh;
  A2D::SolutionVector<T>& vec;
  ElemVecArray_t elem_vec_array;  // The heavy-weight storage
};

}  // namespace A2D

#endif  //  A2D_FE_ELEMENT_VECTOR_H
