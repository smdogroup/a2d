
#ifndef A2D_FE_ELEMENT_VECTOR_H
#define A2D_FE_ELEMENT_VECTOR_H

#include "array.h"
#include "multiphysics/femesh.h"
#include "utils/complex_math.h"

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

  1.2 The FEDof must be indexable via operator[](const index_t)

  2. get_element_values(elem, dof)

  This function gets the element degrees of freedom associated with the index
  elem, and sets them into the FEDof object.  The dof object may store a local
  copy of them or just contain a pointer to the array of element values.

  2. add_element_values(elem, dof)

  This function adds the values to the element vector. This function may be
  empty if the implementation uses references to the values.

  3. get_zero_values()

  Initialize and zero values from the temporary or local storage only. Note
  that this does not zero values in the source vector.

  4. get_values()

  Initialize the values from the solution vector into any temporary or local
  storage in the vector.

  5. add_values()

  Add values into the source vector from any local storage. This may be an
  empty function if the values are stored directly.
*/

enum class ElemVecType { Serial, Parallel, Empty };

/**
 * @brief Base class for the vector-centric view of the degrees of freedom.
 *
 * This class uses CRTP to achieve compile-time polymorphism to eliminate the
 * overhead of vtable lookup. To ``derive'' from this base class, implementation
 * class should have the following inheritance syntax:
 *
 * template <class... Args>
 * class ElementVectorImpl
 *     : public ElementVectorBase<ElemVecType::[Serial,Parallel],
 *                                ElementVectorImpl<Args..>> {...};
 *
 * @tparam _evtype, type of element vector implementation, serial or parallel
 * @tparam ElementVectorImpl the actual implementation class
 */
template <ElemVecType _evtype, class ElementVectorImpl>
class ElementVectorBase {
 public:
  static constexpr ElemVecType evtype = _evtype;

  /**
   * @brief Placeholder for a helper lightweight object to access DOFs
   * associated with the element(s)
   *
   * @tparam T numeric type
   * @tparam FEDofImpl the actual implementation class
   */
  template <typename T, class FEDofImpl>
  class FEDofBase {
   public:
    /**
     * @brief The following two operators should be implemented in derived FEDof
     *
     * @param index dof index local to an element
     */
    T& operator[](const index_t index) { return static_cast<FEDofImpl*>(this)->operator[](index); }
    const T& operator[](const index_t index) const {
      return static_cast<const FEDofImpl*>(this)->operator[](index);
    }
  };

  /// @brief Return number of elements
  index_t get_num_elements() const {
    return static_cast<const ElementVectorImpl*>(this)->get_num_elements();
  }

  /**
   * @brief ElementVectorImpl should implement the following three methods if
   *        ElementVectorImpl::evtype == ElemVecType::Serial
   *
   * @param elem_index the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   *            freedom
   */
  /// @brief Get the element values from the object and store them in the FEDof
  template <typename T, class FEDofImpl>
  void get_element_values(index_t elem_index, FEDofBase<T, FEDofImpl>& dof) {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Serial) {
      return static_cast<ElementVectorImpl*>(this)->get_element_values(elem_index, dof);
    }
  }
  /// @brief Add the degree of freedom values to the element vector
  template <typename T, class FEDofImpl>
  void add_element_values(index_t elem_index, const FEDofBase<T, FEDofImpl>& dof) {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Serial) {
      return static_cast<ElementVectorImpl*>(this)->add_element_values(elem_index, dof);
    }
  }
  /// @brief Set the degree of freedom values to the element vector
  template <typename T, class FEDofImpl>
  void set_element_values(index_t elem_index, const FEDofBase<T, FEDofImpl>& dof) {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Serial) {
      return static_cast<ElementVectorImpl*>(this)->set_element_values(elem_index, dof);
    }
  }

  /**
   * @brief ElementVectorImpl should implement the following three methods if
   *        ElementVectorImpl::evtype == ElemVecType::Parallel
   */
  /// @brief Initialize local dof values to zero
  void get_zero_values() {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Parallel) {
      return static_cast<ElementVectorImpl*>(this)->get_zero_values();
    }
  }
  /// @brief Populate element-view data from global data for all elements
  void get_values() {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Parallel) {
      return static_cast<ElementVectorImpl*>(this)->get_values();
    }
  }
  /// @brief Set global data from element-view data for all elements.
  /// Note: Data consistency is assumed
  void set_values() {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Parallel) {
      return static_cast<ElementVectorImpl*>(this)->set_values();
    }
  }
  /// @brief Add global data from element-view data for all elements
  void add_values() {
    if constexpr (ElementVectorImpl::evtype == ElemVecType::Parallel) {
      return static_cast<ElementVectorImpl*>(this)->add_values();
    }
  }
};

/*
  Element vector implementation for empty values - does nothing
*/
class ElementVector_Empty {
 public:
  static constexpr ElemVecType evtype = ElemVecType::Empty;
  ElementVector_Empty() = default;

  // All implementations should have FEDof and get_num_elements
  class FEDof {
   public:
    FEDof(index_t elem, ElementVector_Empty& elem_vec) {}
  };
  index_t get_num_elements() const { return 0; }

  // Parallel implementation should have these four methods
  void get_zero_values() {}
  void get_values() {}
  void add_values() {}
  void set_values() {}

  // Serial implementation should have these three methods
  template <class Dof>
  void get_element_values(index_t elem, Dof& dof) const {}
  template <class Dof>
  void add_element_values(index_t elem, const Dof& dof) const {}
  template <class Dof>
  void set_element_values(index_t elem, const Dof& dof) const {}
};

/**
 * @brief In-place element vector implementation
 *
 * @tparam Basis type of the basis, e.g. FEBasis<...>
 * @tparam VecType type of the solution vector, e.g. SolutionVector<...>
 */
template <typename T, class Basis, class VecType>
class ElementVector_Serial : public ElementVector_Empty {
 public:
  static constexpr ElemVecType evtype = ElemVecType::Serial;
  ElementVector_Serial(ElementMesh<Basis>& mesh, VecType& vec) : mesh(mesh), vec(vec) {}

  // Required DOF container object
  class FEDof {
   public:
    FEDof(index_t elem, const ElementVector_Serial& elem_vec) {
      std::fill(dof, dof + Basis::ndof, T(0.0));
    }

    /**
     * @brief Get a reference to the underlying element data
     *
     * @return A reference to the degree of freedom
     */
    T& operator[](const index_t index) { return dof[index]; }
    const T& operator[](const index_t index) const { return dof[index]; }

   private:
    // Variables for all the basis functions
    T dof[Basis::ndof];
  };

 public:
  /**
   * @brief Get the number of elements
   */
  index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Get the element values from the object and store them in the FEDof
   *
   * @param elem the element index
   * @param dof the object that stores a reference to the degrees of freedom
   */
  void get_element_values(index_t elem, FEDof& dof) {
    if constexpr (Basis::nbasis > 0) {
      operate_element_values<ELEM_VALS_OP::GET, 0>(elem, dof);
    }
  }

  /**
   * @brief Add the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   */
  void add_element_values(index_t elem, const FEDof& dof) {
    if constexpr (Basis::nbasis > 0) {
      operate_element_values<ELEM_VALS_OP::ADD, 0>(elem, dof);
    }
  }

  /**
   * @brief Set the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   */
  void set_element_values(index_t elem, const FEDof& dof) {
    if constexpr (Basis::nbasis > 0) {
      operate_element_values<ELEM_VALS_OP::SET, 0>(elem, dof);
    }
  }

 private:
  enum class ELEM_VALS_OP { GET, ADD, SET };

  template <ELEM_VALS_OP op, index_t basis>
  void operate_element_values(
      index_t elem,
      typename std::conditional<op == ELEM_VALS_OP::GET, FEDof, const FEDof>::type& dof) {
    for (index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.template get_global_dof_sign<basis>(elem, i);
      const index_t dof_index = mesh.template get_global_dof<basis>(elem, i);
      if constexpr (op == ELEM_VALS_OP::GET) {
        dof[i + Basis::template get_dof_offset<basis>()] = sign * vec[dof_index];
      } else if constexpr (op == ELEM_VALS_OP::ADD) {
        vec[dof_index] += sign * dof[i + Basis::template get_dof_offset<basis>()];
      } else if constexpr (op == ELEM_VALS_OP::SET) {
        vec[dof_index] = sign * dof[i + Basis::template get_dof_offset<basis>()];
      }
    }
    if constexpr (basis + 1 < Basis::nbasis) {
      operate_element_values<op, basis + 1>(elem, dof);
    }
  }

  ElementMesh<Basis>& mesh;
  VecType& vec;
};

/**
 * @brief Parallel element vector implementation
 *
 * This class allocates a heavy-weight 2-dimensional array to store (potentially
 * duplicated) local degrees of freedom to achieve parallelization among
 * elements.
 *
 * global to local dof population is done in parallel, local to global dof add
 * is done by atomic operation to resolve write conflicts.
 *
 * @tparam Basis type of the basis, e.g. FEBasis<...>
 * @tparam VecType type of the solution vector, e.g. SolutionVector<...>
 */
template <typename T, class Basis, class VecType>
class ElementVector_Parallel : public ElementVector_Empty {
  using ElemVecArray_t = MultiArrayNew<T * [Basis::ndof]>;

 public:
  static constexpr ElemVecType evtype = ElemVecType::Parallel;
  ElementVector_Parallel(ElementMesh<Basis>& mesh, VecType& vec)
      : mesh(mesh), vec(vec), elem_vec_array("elem_vec_array", mesh.get_num_elements()) {}

  class FEDof {
   public:
    FEDof(index_t elem, const ElementVector_Parallel& elem_vec)
        : elem(elem), elem_vec_array(elem_vec.elem_vec_array) {}

    T& operator[](const int index) { return elem_vec_array(elem, index); }
    const T& operator[](const int index) const { return elem_vec_array(elem, index); }

   private:
    const index_t elem;
    ElemVecArray_t elem_vec_array;
  };

 public:
  /**
   * @brief Get number of elements
   */
  index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Initialize local dof values to zero
   */
  void get_zero_values() { BLAS::zero(elem_vec_array); }

  /**
   * @brief Populate element-view data from global data for all elements
   */
  void get_values() {
    auto loop_body = A2D_LAMBDA(const index_t elem) {
      operate_element_values<ELEM_VALS_OP::GET, Basis::nbasis>(elem);
    };

    // for (index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
    //   loop_body(elem);
    // }

    Kokkos::parallel_for("get_values", mesh.get_num_elements(), loop_body);
  }

  /**
   * @brief Set global data from element-view data for all elements.
   * Note: Data consistency is assumed
   */
  void set_values() {
    auto loop_body = A2D_LAMBDA(const index_t elem) {
      operate_element_values<ELEM_VALS_OP::SET, Basis::nbasis>(elem);
    };

    // for (index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
    //   loop_body(elem);
    // }

    Kokkos::parallel_for("set_values", mesh.get_num_elements(), loop_body);
  }

  /**
   * @brief Add global data from element-view data for all elements
   */
  void add_values() {
    auto loop_body = A2D_LAMBDA(const index_t elem) {
      operate_element_values<ELEM_VALS_OP::ADD, Basis::nbasis>(elem);
    };

    // for (index_t elem = 0; elem < mesh.get_num_elements(); elem++) {
    //   loop_body(elem);
    // }
    Kokkos::parallel_for("add_values", mesh.get_num_elements(), loop_body);
  }

 private:
  enum class ELEM_VALS_OP { GET, ADD, SET };

  /**
   * @brief Populate element data for a single element
   *
   * @tparam nbasis number of function spaces for the element, at least 1
   * @param elem_idx element index
   */
  template <ELEM_VALS_OP op, index_t nbasis>
  void operate_element_values(const index_t& elem_idx) {
    for (index_t i = 0; i < Basis::template get_ndof<nbasis - 1>(); i++) {
      const int& sign = mesh.template get_global_dof_sign<nbasis - 1>(elem_idx, i);
      const index_t& dof_index = mesh.template get_global_dof<nbasis - 1>(elem_idx, i);
      const index_t dof_idx = i + Basis::template get_dof_offset<nbasis - 1>();
      if constexpr (op == ELEM_VALS_OP::GET) {
        elem_vec_array(elem_idx, dof_idx) = sign * vec[dof_index];
      } else if constexpr (op == ELEM_VALS_OP::ADD) {
        Kokkos::atomic_add(&vec[dof_index], sign * elem_vec_array(elem_idx, dof_idx));

      } else if constexpr (op == ELEM_VALS_OP::SET) {
        // Note: we encounter race condition here but it's (assumed) fine for
        // now as we assume data consistency across elements
        vec[dof_index] = sign * elem_vec_array(elem_idx, dof_idx);
      }
    }
    if constexpr (nbasis > 1) {
      operate_element_values<op, nbasis - 2>(elem_idx);
    }
    return;
  }

  ElementMesh<Basis>& mesh;
  VecType& vec;
  ElemVecArray_t elem_vec_array;  // The heavy-weight storage
};

}  // namespace A2D

#endif  //  A2D_FE_ELEMENT_VECTOR_H
