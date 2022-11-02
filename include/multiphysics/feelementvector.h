
#ifndef A2D_FE_ELEMENT_VECTOR_H
#define A2D_FE_ELEMENT_VECTOR_H

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

  1.1 Constructor with no arguments
  1.2 T* get() and const T* get const () that retrieves a pointer to the
  element data

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
  ElementVector_Serial(A2D::ElementMesh<Basis>& mesh,
                       A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {}

  // Required DOF container object (different for each element vector
  // implementation)
  class FEDof {
   public:
    FEDof(A2D::index_t elem, ElementVector_Serial& elem_vec) {
      std::fill(dof, dof + Basis::ndof, T(0.0));
    }

    template <A2D::index_t index>
    A2D_INLINE_FUNCTION T& operator[](const int i) {
      return dof[Basis::template get_dof_offset<index>() + i];
    }

    template <A2D::index_t index>
    A2D_INLINE_FUNCTION const T& operator[](const int i) const {
      return dof[Basis::template get_dof_offset<index>() + i];
    }

    /**
     * @brief Get the values associated with the given basis
     *
     * @return A pointer to the degrees of freedom
     */
    template <A2D::index_t index>
    T* get() {
      return &dof[Basis::template get_dof_offset<index>()];
    }

    /**
     * @brief Get the values associated with the given basis
     *
     * @return A pointer to the degrees of freedom
     */
    template <A2D::index_t index>
    const T* get() const {
      return &dof[Basis::template get_dof_offset<index>()];
    }

   private:
    // Variables for all the basis functions
    T dof[Basis::ndof];
  };

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
   * freeom
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
    T* values = dof.template get<basis>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      values[i] = sign * vec[dof_index];
    }
    if constexpr (basis > 0) {
      get_element_values_<basis - 1>(elem, dof);
    }
  }

  template <A2D::index_t basis>
  void add_element_values_(A2D::index_t elem, const FEDof& dof) {
    const T* values = dof.template get<basis>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      vec[dof_index] += sign * values[i];
    }
    if constexpr (basis > 0) {
      add_element_values_<basis - 1>(elem, dof);
    }
  }

  A2D::ElementMesh<Basis>& mesh;
  A2D::SolutionVector<T>& vec;
};

/*
   This class allocates a heavy-weight 2-dimensional array to store (potentially
   duplicated) local degrees of freedom to realize a between-element
   parallelization.

   global to local dof population is done in parallel, local to global dof add
   is done by atomic operation to resolve write conflicts.
 */
template <typename T, class Basis>
class ElementVector_Parallel {
 public:
  ElementVector_Parallel(A2D::ElementMesh<Basis>& mesh,
                         A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec), elem_vec_array("elem_vec_array", Basis::ndof) {
    // Populate the array
    init_values();
  }

  // Required DOF container object (different for each element vector
  // implementation)
  class FEDof {
   public:
    FEDof(A2D::index_t elem, ElementVector_Parallel& elem_vec_array)
        : elem(elem), elem_vec_array(elem_vec_array) {}

    template <A2D::index_t index>
    A2D_INLINE_FUNCTION T& operator[](const int i) {
      return elem_vec_array(elem, Basis::template get_dof_offset<index>() + i);
    }

    template <A2D::index_t index>
    A2D_INLINE_FUNCTION const T& operator[](const int i) const {
      return elem_vec_array(elem, Basis::template get_dof_offset<index>() + i);
    }

    /**
     * @brief Get the values associated with the given basis
     *
     * @return A pointer to the degrees of freedom
     */
    template <A2D::index_t index>
    T& operator[](int i) T* get() {
      return &dof[Basis::template get_dof_offset<index>()];
    }

    /**
     * @brief Get the values associated with the given basis
     *
     * @return A pointer to the degrees of freedom
     */
    template <A2D::index_t index>
    const T* get() const {
      return &dof[Basis::template get_dof_offset<index>()];
    }

   private:
    ElemVecArray_t& elem_vec_array;
    const A2D::index_t elem;
  };

  /**
   * @brief Add local dof to global dof
   */
  void add_values() {
    for (A2D::index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
      add_element_values<Basis::nbasis>(elem);
    }
  }

 private:
  /**
   * @brief Populate local dofs from global dof
   */
  void init_values() {
    for (A2D::index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
      populate_element_values<Basis::nbasis>(elem);
    }
  }

  /**
   * @brief Populate local dof for a single element
   *
   * @tparam nbasis number of function spaces for the element, at least 1
   * @param elem_idx element index
   */
  template <A2D::index_t nbasis>
  void populate_element_values(const A2D::index_t& elem_idx) {
    for (A2D::index_t i = 0; i < Basis::get_ndof<nbasis - 1>; i++) {
      const int& sign = mesh.get_global_dof_sign(elem_idx, nbasis - 1, i);
      const A2D::index_t& dof_index =
          mesh.get_global_dof(elem_idx, nbasis - 1, i);
      element_vec_array(elem_idx, dof_index) = sign * vec[dof_index];
    }
    if constexpr (nbasis > 1) {
      populate_element_values<nbasis - 2>(elem_idx);
    }
    return;
  }

  /**
   * @brief Add local dof to global dof for a single element
   *
   * @tparam nbasis number of function spaces for the element, at least 1
   */
  template <A2D::index_t nbasis>
  void add_element_values(const A2D::index_t& elem_idx) {
    for (A2D::index_t i = 0; i < Basis::get_ndof<nbasis - 1>; i++) {
      const int& sign = mesh.get_global_dof_sign(elem_idx, nbasis - 1, i);
      const A2D::index_t& dof_index =
          mesh.get_global_dof(elem_idx, nbasis - 1, i);
      Kokkos::atomic_add(&vec[dof_index],
                         sign * element_vec_array(elem_idx, dof_index));
    }
    if constexpr (nbasis > 1) {
      populate_element_values<nbasis - 2>(elem_idx);
    }
    return;
  }

  using ElemVecArray_t = A2D::MultiArrayNew<T*>;
  A2D::ElementMesh<Basis>& mesh;
  A2D::SolutionVector<T>& vec;
  ElemVecArray_t elem_vec_array;  // The heavy-weight storage
};

}  // namespace A2D

#endif  //  A2D_FE_ELEMENT_VECTOR_H
