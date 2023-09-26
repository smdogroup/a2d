#ifndef A2D_FE_ELEMENT_MAT_H
#define A2D_FE_ELEMENT_MAT_H

#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "utils/complex_math.h"

namespace A2D {

template <typename T, class Basis, class MatType>
class ElementMat_Serial {
 public:
  ElementMat_Serial(ElementMesh<Basis>& mesh, MatType& mat)
      : mesh(mesh), mat(mat) {}

  // Required DOF container object (different for each element vector
  // implementation)
  class FEMat {
   public:
    static const index_t size = Basis::ndof * Basis::ndof;

    FEMat(index_t elem, ElementMat_Serial<T, Basis, MatType>& elem_mat)
        : A(size, T(0.0)) {}

    /**
     * @brief Get a reference to the underlying element data
     *
     * @return A reference to the degree of freedom
     */
    T& operator()(const index_t i, const index_t j) {
      return A[i * Basis::ndof + j];
    }
    const T& operator()(const index_t i, const index_t j) const {
      return A[i * Basis::ndof + j];
    }

   private:
    // Variables for all the basis functions
    std::vector<T> A;
  };

  /**
   * @brief Get the number of elements
   */
  index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Add the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   *
   * If FEDof contains a pointer to data, this function may do nothing
   */
  void add_element_values(index_t elem, FEMat& elem_mat) {
    index_t dof[Basis::ndof];
    int sign[Basis::ndof];
    if constexpr (Basis::nbasis > 0) {
      get_dof<0>(elem, dof, sign);
    }

    for (index_t i = 0; i < Basis::ndof; i++) {
      for (index_t j = 0; j < Basis::ndof; j++) {
        elem_mat(i, j) = sign[i] * sign[j] * elem_mat(i, j);
      }
    }

    mat.add_values(Basis::ndof, dof, Basis::ndof, dof, elem_mat);
  }

 private:
  template <index_t basis>
  void get_dof(index_t elem, index_t dof[], int sign[]) {
    for (index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      sign[i + Basis::template get_dof_offset<basis>()] =
          mesh.template get_global_dof_sign<basis>(elem, i);
      dof[i + Basis::template get_dof_offset<basis>()] =
          mesh.template get_global_dof<basis>(elem, i);
    }
    if constexpr (basis + 1 < Basis::nbasis) {
      get_dof<basis + 1>(elem, dof, sign);
    }
  }

  ElementMesh<Basis>& mesh;
  MatType& mat;
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_MAT_H