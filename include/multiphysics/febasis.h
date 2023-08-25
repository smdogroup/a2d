#ifndef A2D_FE_BASIS_H
#define A2D_FE_BASIS_H

#include "a2dobjs.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/fespace.h"

namespace A2D {

/**
 * @brief The type of basis function that a given basis class implements
 *
 */
enum BasisType { H1, HCURL, HDIV, L2 };

/*
  The Basis class type is a class with all const static member data and
  static member functions.

  The basis classes provide:

  1. static const ndof: The number of degrees of freedom for the basis function.

  2. static const ncomp: The number of components in the function space output.

  3. static function interp: This function evaluates a function space object at
  all the quadrature points. The interpolation computes both values and any
  spatial derivatives that may be used in the weak form of the governing
  equations. This can be written as

  u = N(pt) * dof

  where u may be the solution as well as the derivatives of the solution, the
  divergence or curl. Here pt is a vector of the quadrature points in the
  reference element space.

  4. static function add: This function adds the residual contribution from the
  corresponding function space object to the degrees of freedom. This can be
  written as

  dof += N(pt)^T * u

  where the dimension of u and N are the same as above.

  5. stride: Indicates that the basis functions are arranged in a special
  structure.

  6. basis: This function provides the full interpolation matrix that takes the
  degrees of freedom and computes the solution. To save space, the structure of
  the matrix may be leveraged for computational efficiency. For instance when
  the static member stride = 2, then the full interpolation matrix N is

  dof = [ u1    v1    u2     v2    u3    v3    u4    v4  ]

  u   = [ N1    0     N2     0     N3    0     N4    0   ]
  u,x = [ N1,x  0     N2,x   0     N3,x  0     N4,x  0   ]
  u,y = [ N1,y  0     N2,y   0     N3,y  0     N4,y  0   ]
  u,z = [ N1,z  0     N2,z   0     N3,y  0     N4,z  0   ]
  u   = [ 0     N1    0      N2   0      N3    0     N4  ]
  u,x = [ 0     N1,x  0      N2,x 0      N3,x  0     N4,x]
  u,y = [ 0     N1,y  0      N2,y 0      N3,y  0     N4,y]
  u,z = [ 0     N1,z  0      N2,z 0      N3,y  0     N4,z]

  Note this may be stored as a 4 x 4 matrix, rather than an 8 x 8 matrix.

  In this case, a call to the function basis returns:

  N =
  [ N1,  N1,x  N1,y  N1,z
    N2,  N2,x  N2,y  N2,z
    N3,  N3,x  N3,y  N3,z
    N4,  N4,x  N4,y  N4,z ]

  Note that the matrix is in column-major order, not row major order!
  This enables more efficient memory access.

  For an H(div) space in 2D, this would look something like this:

  dof = [ dof1   dof2   dof3   dof4   dof5]

  u   = [ N1    N4    N7    N10    N13 ]
  v   = [ N2    N5    N8    N11    N14 ]
  div = [ N3    N6    N9    N12    N15 ]

  where N1 through N15 are each separate functions.

  In this case, stride = 1.

  Note: The entries in the matrix must be consistent in ordering with the
  set_value()/get_value() calls from the function space objects.
*/

/*
  Template pattern for computing the number of dof in a basis
*/
template <class... Basis>
struct __count_basis_dof;

template <>
struct __count_basis_dof<> {
  static const index_t ndof = 0;
};

template <class First, class... Remain>
struct __count_basis_dof<First, Remain...> {
  static const index_t ndof = First::ndof + __count_basis_dof<Remain...>::ndof;
};

/*
  Template pattern for computing the total number of components
*/
template <class... Basis>
struct __count_basis_ncomp;

template <>
struct __count_basis_ncomp<> {
  static const index_t ncomp = 0;
};

template <class First, class... Remain>
struct __count_basis_ncomp<First, Remain...> {
  static const index_t ncomp =
      First::ncomp + __count_basis_ncomp<Remain...>::ncomp;
};

/*
  Template pattern for computing the size of the interpolation matrices
*/
template <class... Basis>
struct __count_basis_size;

template <>
struct __count_basis_size<> {
  static const index_t basis_size = 0;
};

template <class First, class... Remain>
struct __count_basis_size<First, Remain...> {
  static const index_t basis_size =
      First::basis_size + __count_basis_size<Remain...>::basis_size;
};

/*
  Check if all dimensions are the same
*/
template <class... Bases>
struct have_same_dim;

template <>
struct have_same_dim<> {
  static constexpr bool value = true;
  static constexpr index_t dim = 0;
};

template <class B1>
struct have_same_dim<B1> {
  static constexpr bool value = true;
  static constexpr index_t dim = B1::dim;
};

template <class B1, class B2>
struct have_same_dim<B1, B2> {
  static constexpr bool value = B1::dim == B2::dim;
  static constexpr index_t dim = B1::dim;
};

template <class B1, class B2, class... rest>
struct have_same_dim<B1, B2, rest...> {
  static constexpr bool value =
      have_same_dim<B1, B2>::value and have_same_dim<B2, rest...>::value;
  static constexpr index_t dim = B1::dim;
};

/*
  The finite element basis class.

  This class stores a collection of basis function objects
*/
template <typename T, class... Basis>
class FEBasis {
 public:
  using BasisSpace = std::tuple<Basis...>;
  static_assert(have_same_dim<Basis...>::value,
                "Bases should have same dimensions");

  static constexpr index_t dim = have_same_dim<Basis...>::dim;

  // Use the definitions from the element types
  using ET = ElementTypes;

  /**
   * @brief Number of basis function objects
   */
  static constexpr index_t nbasis = std::tuple_size<std::tuple<Basis...>>();

  /**
   * @brief Count the total number of components from all the basis functions
   */
  static constexpr index_t ncomp = __count_basis_ncomp<Basis...>::ncomp;

  /**
   * @brief Count up the total number of degrees of freedom for this set of
   * basis functions
   */
  static constexpr index_t ndof = __count_basis_dof<Basis...>::ndof;

  /**
   * @brief Count up the total basis size required to store all of the
   * interpolation matrices for this set of basis functions
   */
  static constexpr index_t basis_size =
      __count_basis_size<Basis...>::basis_size;

  /**
   * @brief Get the number of degrees of freedom associated with the given basis
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_ndof() {
    if constexpr (nbasis == 0) {
      return 0;
    } else {
      return std::tuple_element<index, BasisSpace>::type::ndof;
    }
  }

  /**
   * @brief Get the cumulative number of degrees of freedom before the given
   * basis index
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_dof_offset() {
    return get_dof_offset_<0, index, Basis...>();
  }

  /**
   * @brief Get the basis_size offset - the cumulative size of the basis
   * function matrix up to the given basis index
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_basis_size_offset() {
    return get_basis_size_offset_<0, index, Basis...>();
  }

  /**
   * @brief Get the component offset for the given basis
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_comp_offset() {
    return get_comp_offset_<0, index, Basis...>();
  }

  /**
   * @brief Get the basis function type
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr BasisType get_basis_type() {
    if constexpr (nbasis == 0) {
      return H1;
    } else {
      return get_basis_type_<0, index, Basis...>();
    }
  }

  /**
   * @brief Get the stride for a given basis
   */
  template <index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_stride() {
    if constexpr (nbasis == 0) {
      return 0;
    } else {
      return get_stride_<0, index, Basis...>();
    }
  }

  /**
   * @brief Get the stride of the given basis index
   *
   * @param basis
   * @return index_t
   */
  KOKKOS_FUNCTION static index_t get_stride(index_t basis) {
    if constexpr (nbasis == 0) {
      return 0;
    } else {
      return get_stride_<0, Basis...>(basis);
    }
  }

  /**
   * @brief Interpolate using the basis functions at the quadrature point
   *
   * @param dof The degree of freedom object
   * @param s The output finite element space object
   */
  template <class Quadrature, class FEDof, class FiniteElementSpace>
  KOKKOS_FUNCTION static void interp(
      const FEDof& dof, QptSpace<Quadrature, FiniteElementSpace>& s) {
    interp_<Quadrature, FEDof, FiniteElementSpace, 0, Basis...>(dof, s);
  }

  /**
   * @brief Add values to the degree of freedom using the interpolation
   * functions
   *
   * @param s The finite element space objec
   * @param dof The degree of freedom object that values are added to
   */
  template <class Quadrature, class FiniteElementSpace, class FEDof>
  KOKKOS_FUNCTION static void add(
      const QptSpace<Quadrature, FiniteElementSpace>& s, FEDof& dof) {
    add_<Quadrature, FiniteElementSpace, FEDof, 0, Basis...>(s, dof);
  }

  /**
   * @brief Interpolate with the basis functions. This can be used to check for
   * consistency
   *
   * @param pt The interpolation point index
   * @param dof The degree of freedom object
   * @param s The output finite element space object
   */
  template <class Quadrature, class FEDof, class FiniteElementSpace>
  KOKKOS_FUNCTION static void interp_basis(index_t pt, const FEDof& dof,
                                           FiniteElementSpace& s) {
    // Evaluate the basis functions
    double N[basis_size];
    eval_basis<Quadrature, Basis...>(pt, N);

    interp_basis_<FEDof, FiniteElementSpace, 0, Basis...>(N, dof, s);
  }

  /**
   * @brief Add values to the degree of freedom using the interpolation
   * functions
   *
   * @param pt The interpolation point index
   * @param s The finite element space objec
   * @param dof The degree of freedom object that values are added to
   */
  template <class Quadrature, class FiniteElementSpace, class FEDof>
  KOKKOS_FUNCTION static void add_basis(index_t pt, const FiniteElementSpace& s,
                                        FEDof& dof) {
    // Evaluate the basis functions
    double N[basis_size];
    eval_basis<Quadrature, Basis...>(pt, N);

    add_basis_<FiniteElementSpace, FEDof, 0, Basis...>(N, s, dof);
  }

  /**
   * @brief Add the outer product of the interpolation matrix with the Jacobian
   * of the components at the quadrature point
   *
   * @tparam Quadrature The quadrature point object
   * @tparam QMat The Jacobian matrix at the quadrature point
   * @tparam Mat The matrix type
   * @param pt The quadrature point index
   * @param jac The Jacobian at the quadrature point
   * @param mat The element Jacobian matrix
   */
  template <class Quadrature, class QMat, class Mat>
  KOKKOS_FUNCTION static void add_outer(index_t pt, const QMat& jac, Mat& mat) {
    // Evaluate the basis functions
    double N[basis_size];
    eval_basis<Quadrature, Basis...>(pt, N);

    // Add the outer-product
    add_outer_<QMat, Mat, 0, Basis...>(N, jac, mat);
  }

  /**
   * @brief Get the number of degrees of freedom for the entity
   *
   * @param basis The basis index
   * @param entity The entity type defined in the ElementTypes object
   * @param index The index of the entity e.g. vertex index
   * @return The number of degrees of freedom
   */
  KOKKOS_FUNCTION static index_t get_entity_ndof(index_t basis,
                                                 ET::ElementEntity entity,
                                                 index_t index) {
    if constexpr (sizeof...(Basis) == 0) {
      return 0;
    } else {
      return get_entity_ndof<0, Basis...>(basis, entity, index);
    }
  }

  /**
   * @brief Get the number of degrees of freedom for the entity
   *
   * @param basis The basis index
   * @param entity The entity type defined in the ElementTypes object
   * @param index The index of the entity e.g. vertex index
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */
  template <class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_entity_dof(index_t basis,
                                             ET::ElementEntity entity,
                                             index_t index,
                                             const ElemDof& element_dof,
                                             EntityDof& entity_dof) {
    if constexpr (sizeof...(Basis) == 0) {
      return;
    } else {
      get_entity_dof<0, ElemDof, EntityDof, Basis...>(basis, entity, index,
                                                      element_dof, entity_dof);
    }
  }

  /**
   * @brief Get the number of degrees of freedom for the entity
   *
   * @param basis The basis index
   * @param entity The entity type defined in the ElementTypes object
   * @param index The index of the entity e.g. vertex index
   * @param orient Orientation flag indicating the relative orientation
   * @param entity_dof Entity DOF in the global orientation
   * @param element_dof Degrees of freedom for this element
   */
  template <class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_entity_dof(index_t basis,
                                             ET::ElementEntity entity,
                                             index_t index, index_t orient,
                                             const EntityDof& entity_dof,
                                             ElemDof& element_dof) {
    if constexpr (sizeof...(Basis) == 0) {
      return;
    } else {
      set_entity_dof<0, EntityDof, ElemDof, Basis...>(
          basis, entity, index, orient, entity_dof, element_dof);
    }
  }

  /**
   * @brief Get the sign of the local dof relative to the reference
   *
   * @param basis The basis index
   * @param entity The entity type defined in the ElementTypes object
   * @param index The index of the entity e.g. vertex index
   * @param orient Orientation flag indicating the relative orientation
   * @param elem_signs Degrees of freedom for this element
   */
  KOKKOS_FUNCTION static void set_entity_signs(index_t basis,
                                               ET::ElementEntity entity,
                                               index_t index, index_t orient,
                                               int elem_signs[]) {
    if constexpr (sizeof...(Basis) == 0) {
      return;
    } else {
      set_entity_signs<0, Basis...>(basis, entity, index, orient, elem_signs);
    }
  }

  /**
   * @brief Get the number of low-order elements
   */
  KOKKOS_FUNCTION constexpr static index_t get_num_lorder_elements() {
    if constexpr (sizeof...(Basis) == 0) {
      return 0;
    } else {
      return get_num_lorder_elements_<Basis...>();
    }
  }

  /**
   * @brief Get the degrees of freedom from an array associated with the
   * high-order dof
   *
   * @tparam HOrderDof High-order degree of freedom array type
   * @tparam LOrderDof Low-order degree of freedom array type
   * @param n Index of the low-order element
   * @param hdof High-order array input
   * @param ldof Low-order array output
   */
  template <class HOrderDof, class LOrderDof>
  KOKKOS_FUNCTION static void get_lorder_dof(const index_t n,
                                             const HOrderDof& hdof,
                                             LOrderDof& ldof) {
    if constexpr (sizeof...(Basis) > 0) {
      get_lorder_dof_<0, 0, HOrderDof, LOrderDof, Basis...>(n, hdof, ldof);
    }
  }

  /**
   * @brief Get the low order degree of freedom signs relative to the high-order
   * ones
   *
   * @param n Index of the low-order element
   * @param horder_signs The signs from the high-order element
   * @param lorder_signs The output signs for the low-order element
   */
  KOKKOS_FUNCTION static void get_lorder_signs(const index_t n,
                                               const int horder_signs[],
                                               int lorder_signs[]) {
    if constexpr (sizeof...(Basis) > 0) {
      get_lorder_signs_<0, 0, Basis...>(n, horder_signs, lorder_signs);
    }
  }

  /**
   * @brief Get the parametric point location associated with the given
   * degree of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  KOKKOS_FUNCTION static void get_dof_point(index_t index, double pt[]) {
    if constexpr (sizeof...(Basis) > 0) {
      get_dof_point_<0, Basis...>(index, pt);
    }
  }

 private:
  template <class Quadrature, class FEDof, class FiniteElementSpace,
            index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static void interp_(
      const FEDof& dof, QptSpace<Quadrature, FiniteElementSpace>& s) {
    // Interpolate
    First::template interp<index, Quadrature, FiniteElementSpace,
                           get_dof_offset<index>()>(dof, s);

    // Do the next solution space, if any...
    interp_<Quadrature, FEDof, FiniteElementSpace, index + 1, Remain...>(dof,
                                                                         s);
  }

  template <class Quadrature, class FEDof, class FiniteElementSpace,
            index_t index>
  KOKKOS_FUNCTION static void interp_(
      const FEDof& dof, QptSpace<Quadrature, FiniteElementSpace>& s) {}

  template <class Quadrature, class FiniteElementSpace, class FEDof,
            index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static void add_(
      const QptSpace<Quadrature, FiniteElementSpace>& s, FEDof& dof) {
    // Add the interpolation
    First::template add<index, Quadrature, FiniteElementSpace,
                        get_dof_offset<index>()>(s, dof);

    // Do the next solution space, if any...
    add_<Quadrature, FiniteElementSpace, FEDof, index + 1, Remain...>(s, dof);
  }

  template <class Quadrature, class FiniteElementSpace, class FEDof,
            index_t index>
  KOKKOS_FUNCTION static void add_(
      const QptSpace<Quadrature, FiniteElementSpace>& s, FEDof& dof) {}

  template <class Quadrature, class First, class... Remain>
  KOKKOS_FUNCTION static void eval_basis(const index_t pt, double N[]) {
    First::template basis<Quadrature>(pt, N);

    eval_basis<Quadrature, Remain...>(pt, &N[First::basis_size]);
  }

  template <class Quadrature>
  KOKKOS_FUNCTION static void eval_basis(const index_t pt, double N[]) {}

  // Get the offset recursively
  template <index_t r, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static constexpr index_t get_dof_offset_() {
    if (r == index) {
      return 0;
    }
    return First::ndof + get_dof_offset_<r + 1, index, Remain...>();
  }

  template <index_t r, index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_dof_offset_() {
    return 0;
  }

  // Get the component offset recursively
  template <index_t r, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static constexpr index_t get_comp_offset_() {
    if (r == index) {
      return 0;
    }
    return First::ncomp + get_comp_offset_<r + 1, index, Remain...>();
  }

  template <index_t r, index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_comp_offset_() {
    return 0;
  }

  // Get the basis size offset recursively
  template <index_t r, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static constexpr index_t get_basis_size_offset_() {
    if (r == index) {
      return 0;
    }
    return First::basis_size +
           get_basis_size_offset_<r + 1, index, Remain...>();
  }

  template <index_t r, index_t index>
  KOKKOS_FUNCTION static constexpr index_t get_basis_size_offset_() {
    return 0;
  }

  // Interpolate with the basis functions evaluated
  template <class FEDof, class FiniteElementSpace, index_t index, class First,
            class... Remain>
  KOKKOS_FUNCTION static void interp_basis_(const double N[], const FEDof& dof,
                                            FiniteElementSpace& s) {
    T value[First::ncomp];
    for (index_t icomp = 0; icomp < First::ncomp; icomp++) {
      value[icomp] = 0.0;
    }

    index_t idof = get_dof_offset<index>();
    for (index_t idx = 0; idx < First::ndof_per_stride; idx++) {
      T* v = value;
      for (index_t istride = 0; istride < First::stride; istride++, idof++) {
        for (index_t icomp = 0; icomp < First::ncomp_per_stride; icomp++) {
          v[0] += N[icomp] * dof[idof];
          v++;
        }
      }
      N += First::ncomp_per_stride;
    }

    for (index_t icomp = 0; icomp < First::ncomp; icomp++) {
      s[get_comp_offset<index>() + icomp] = value[icomp];
    }

    interp_basis_<FEDof, FiniteElementSpace, index + 1, Remain...>(N, dof, s);
  }

  template <class FEDof, class FiniteElementSpace, index_t index>
  KOKKOS_FUNCTION static void interp_basis_(const double N[], const FEDof& dof,
                                            FiniteElementSpace& s) {}

  // Interpolate with the basis functions evaluated
  template <class FiniteElementSpace, class FEDof, index_t index, class First,
            class... Remain>
  KOKKOS_FUNCTION static void add_basis_(const double N[],
                                         const FiniteElementSpace& s,
                                         FEDof& dof) {
    T values[First::ncomp];
    for (index_t icomp = 0; icomp < First::ncomp; icomp++) {
      values[icomp] = s[get_comp_offset<index>() + icomp];
    }

    index_t idof = get_dof_offset<index>();
    for (index_t idx = 0; idx < First::ndof_per_stride; idx++) {
      T* v = values;
      for (index_t istride = 0; istride < First::stride; istride++, idof++) {
        for (index_t icomp = 0; icomp < First::ncomp_per_stride; icomp++) {
          dof[idof] += N[icomp] * v[0];
          v++;
        }
      }
      N += First::ncomp_per_stride;
    }

    add_basis_<FiniteElementSpace, FEDof, index + 1, Remain...>(N, s, dof);
  }

  template <class FiniteElementSpace, class FEDof, index_t index>
  KOKKOS_FUNCTION static void add_basis_(const double N[],
                                         const FiniteElementSpace& s,
                                         FEDof& dof) {}

  template <class QMat, class Mat, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static void add_outer_(const double N0[], const QMat& jac,
                                         Mat& mat) {
    const double* N = &N0[get_basis_size_offset<index>()];

    index_t idof = get_dof_offset<index>();
    for (index_t idx = 0; idx < First::ndof_per_stride; idx++) {
      for (index_t istride = 0; istride < First::stride; istride++, idof++) {
        // Compute the values in a row vector values = N_{i}^{T} * jac
        const index_t offset =
            istride * First::ncomp_per_stride + get_comp_offset<index>();

        T values[ncomp];
        for (index_t jcomp = 0; jcomp < ncomp; jcomp++) {
          values[jcomp] = 0.0;
          for (index_t icomp = 0; icomp < First::ncomp_per_stride; icomp++) {
            values[jcomp] += N[icomp] * jac(icomp + offset, jcomp);
          }
        }

        // Make the call to complete the product with the row
        add_outer_row_<Mat, 0, Basis...>(idof, N0, values, mat);
      }

      N += First::ncomp_per_stride;
    }

    add_outer_<QMat, Mat, index + 1, Remain...>(N0, jac, mat);
  }

  template <class QMat, class Mat, index_t index>
  KOKKOS_FUNCTION static void add_outer_(const double N0[], const QMat& jac,
                                         Mat& mat) {}

  template <class Mat, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static void add_outer_row_(const index_t idof,
                                             const double N[], const T values[],
                                             Mat& mat) {
    index_t jdof = get_dof_offset<index>();
    for (index_t idx = 0; idx < First::ndof_per_stride; idx++) {
      const T* v = &values[get_comp_offset<index>()];

      for (index_t jstride = 0; jstride < First::stride; jstride++, jdof++) {
        T val = 0.0;
        for (index_t jcomp = 0; jcomp < First::ncomp_per_stride; jcomp++) {
          val += N[jcomp] * v[0];
          v++;
        }

        mat(idof, jdof) += val;
      }

      N += First::ncomp_per_stride;
    }

    add_outer_row_<Mat, index + 1, Remain...>(idof, N, values, mat);
  }

  template <class Mat, index_t index>
  KOKKOS_FUNCTION static void add_outer_row_(const index_t idof,
                                             const double N[], const T values[],
                                             Mat& mat) {}

  template <index_t r, class First, class... Remain>
  KOKKOS_FUNCTION static index_t get_entity_ndof(index_t basis,
                                                 ET::ElementEntity entity,
                                                 index_t index) {
    if (basis == r) {
      return First::get_entity_ndof(entity, index);
    }
    if constexpr (sizeof...(Remain) == 0) {
      return First::get_entity_ndof(entity, index);
    } else {
      return get_entity_ndof<r + 1, Remain...>(basis, entity, index);
    }
  }

  template <index_t r, class ElemDof, class EntityDof, class First,
            class... Remain>
  KOKKOS_FUNCTION static void get_entity_dof(index_t basis,
                                             ET::ElementEntity entity,
                                             index_t index,
                                             const ElemDof& element_dof,
                                             EntityDof& entity_dof) {
    if (basis == r) {
      First::template get_entity_dof<get_dof_offset<r>(), ElemDof, EntityDof>(
          entity, index, element_dof, entity_dof);
    }
    if constexpr (sizeof...(Remain) == 0) {
      First::template get_entity_dof<get_dof_offset<r>(), ElemDof, EntityDof>(
          entity, index, element_dof, entity_dof);
    } else {
      get_entity_dof<r + 1, ElemDof, EntityDof, Remain...>(
          basis, entity, index, element_dof, entity_dof);
    }
  }

  template <index_t r, class EntityDof, class ElemDof, class First,
            class... Remain>
  KOKKOS_FUNCTION static void set_entity_dof(index_t basis,
                                             ET::ElementEntity entity,
                                             index_t index, index_t orient,
                                             const EntityDof& entity_dof,
                                             ElemDof& element_dof) {
    if (basis == r) {
      First::template set_entity_dof<get_dof_offset<r>(), EntityDof, ElemDof>(
          entity, index, orient, entity_dof, element_dof);
    }
    if constexpr (sizeof...(Remain) == 0) {
      First::template set_entity_dof<get_dof_offset<r>(), EntityDof, ElemDof>(
          entity, index, orient, entity_dof, element_dof);
    } else {
      set_entity_dof<r + 1, EntityDof, ElemDof, Remain...>(
          basis, entity, index, orient, entity_dof, element_dof);
    }
  }

  template <index_t r, class First, class... Remain>
  KOKKOS_FUNCTION static void set_entity_signs(index_t basis,
                                               ET::ElementEntity entity,
                                               index_t index, index_t orient,
                                               int elem_signs[]) {
    if (basis == r) {
      First::template set_entity_signs<get_dof_offset<r>()>(entity, index,
                                                            orient, elem_signs);
    }
    if constexpr (sizeof...(Remain) == 0) {
      First::template set_entity_signs<get_dof_offset<r>()>(entity, index,
                                                            orient, elem_signs);
    } else {
      set_entity_signs<r + 1, Remain...>(basis, entity, index, orient,
                                         elem_signs);
    }
  }

  template <index_t r, class First, class... Remain>
  KOKKOS_FUNCTION static void get_dof_point_(index_t index, double pt[]) {
    if (index < get_ndof<r>()) {
      return First::get_dof_point(index, pt);
    } else if constexpr (sizeof...(Remain) > 0) {
      get_dof_point_<r + 1, Remain...>(index - get_ndof<r>(), pt);
    }
  }

  template <index_t r, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static constexpr BasisType get_basis_type_() {
    if (r == index) {
      return First::get_basis_type();
    }
    if constexpr (sizeof...(Remain) == 0) {
      return First::get_basis_type();
    } else {
      return get_basis_type_<r + 1, index, Remain...>();
    }
  }

  template <index_t r, index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static constexpr index_t get_stride_() {
    if (r == index) {
      return First::stride;
    }
    if constexpr (sizeof...(Remain) == 0) {
      return First::stride;
    } else {
      return get_stride_<r + 1, index, Remain...>();
    }
  }

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION static index_t get_stride_(const index_t basis) {
    if (index == basis) {
      return First::stride;
    }
    if constexpr (sizeof...(Remain) == 0) {
      return 0;
    } else {
      return get_stride_<index + 1, Remain...>(basis);
    }
  }

  /**
   * @brief Get the number of low-order elements
   */
  template <class First, class... Remain>
  KOKKOS_FUNCTION constexpr static index_t get_num_lorder_elements_() {
    if constexpr (sizeof...(Remain) == 0) {
      return First::get_num_lorder_elements();
    } else {
      constexpr index_t size = get_num_lorder_elements_<Remain...>();
      static_assert(size == First::get_num_lorder_elements(),
                    "Number of low-order elements is not compatible");
      return size;
    }
  }

  template <index_t hoffset, index_t loffset, class HOrderDof, class LOrderDof,
            class First, class... Remain>
  KOKKOS_FUNCTION static void get_lorder_dof_(const index_t n,
                                              const HOrderDof& hdof,
                                              LOrderDof& ldof) {
    First::template get_lorder_dof<hoffset, loffset, HOrderDof, LOrderDof>(
        n, hdof, ldof);
    if constexpr (sizeof...(Remain) > 0) {
      get_lorder_dof_<hoffset + First::ndof, loffset + First::LOrderBasis::ndof,
                      HOrderDof, LOrderDof, Remain...>(n, hdof, ldof);
    }
  }

  template <index_t hoffset, index_t loffset, class First, class... Remain>
  KOKKOS_FUNCTION static void get_lorder_signs_(const index_t n,
                                                const int horder_signs[],
                                                int lorder_signs[]) {
    First::template get_lorder_signs<hoffset, loffset>(n, horder_signs,
                                                       lorder_signs);
    if constexpr (sizeof...(Remain) > 0) {
      get_lorder_signs_<hoffset + First::ndof,
                        loffset + First::LOrderBasis::ndof, Remain...>(
          n, horder_signs, lorder_signs);
    }
  }
};

}  // namespace A2D

#endif  // A2D_FE_BASIS_H