#ifndef A2D_LAGRANGE_HEX_BASIS_H
#define A2D_LAGRANGE_HEX_BASIS_H

#include "multiphysics/febasis.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

template <typename T, index_t C, index_t degree>
class LagrangeH1HexBasis {
 public:
  using ET = ElementTypes;

  static const index_t dim = 3;             // Spatial dimension
  static const index_t order = degree + 1;  // Number of nodes along each edge

  static const index_t ndof =
      C * order * order * order;  // Total number of degrees of freedom

  // Number of components
  static const index_t ncomp = H1Space<T, C, dim>::ncomp;

  /**
   * @brief Degree of freedom handling on the vertices, edges, faces and volume
   *
   * @param entity The type of topological entity (vertex, edge, face or volume)
   * @param index The index of the topological entity (e.g. edge index)
   * @return The number of degrees of freedom
   */
  static index_t get_entity_ndof(ET::ElementEntity entity, index_t index) {
    if (entity == ET::VERTEX) {
      return C;
    } else if (entity == ET::EDGE) {
      return C * (order - 2);
    } else if (entity == ET::FACE) {
      return C * (order - 2) * (order - 2);
    } else if (entity == ET::VOLUME) {
      return C * (order - 2) * (order - 2) * (order - 2);
    }
    return 0;
  }

  /**
   * @brief Get the entity DOF from the element DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of entity (vertex, edge, face or volume)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */
  template <index_t offset, class ElemDof, class EntityDof>
  static void get_entity_dof(ET::ElementEntity entity, index_t index,
                             const ElemDof& element_dof,
                             EntityDof& entity_dof) {
    if (entity == ET::VERTEX) {
      ET::get_hex_vert_dof<offset, C, order, order, order, ElemDof, EntityDof>(
          index, element_dof, entity_dof);
    } else if (entity == ET::EDGE) {
      const bool endp = false;
      ET::get_hex_edge_dof<offset, endp, C, order, order, order, ElemDof,
                           EntityDof>(index, element_dof, entity_dof);
    } else if (entity == ET::FACE) {
      const bool endp = false;
      ET::get_hex_face_dof<offset, endp, C, order, order, order, ElemDof,
                           EntityDof>(index, element_dof, entity_dof);
    } else {
      const bool endp = false;
      ET::get_hex_volume_dof<offset, endp, C, order, order, order, ElemDof,
                             EntityDof>(element_dof, entity_dof);
    }
  }

  /**
   * @brief Set the element DOF and signs from the entity DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of entity (vertex, edge, face or volume)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param entity_dof Entity DOF in the global orientation
   * @param element_dof Degrees of freedom for this element
   */
  template <index_t offset, class EntityDof, class ElemDof>
  static void set_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const EntityDof& entity_dof,
                             ElemDof& element_dof) {
    if (entity == ET::VERTEX) {
      ET::set_hex_vert_dof<offset, C, order, order, order, EntityDof, ElemDof>(
          index, entity_dof, element_dof);
    } else if (entity == ET::EDGE) {
      const bool endp = false;
      ET::set_hex_edge_dof<offset, endp, C, order, order, order, EntityDof,
                           ElemDof>(index, orient, entity_dof, element_dof);
    } else if (entity == ET::FACE) {
      const bool endp = false;
      ET::set_hex_face_dof<offset, endp, C, order, order, order, EntityDof,
                           ElemDof>(index, orient, entity_dof, element_dof);
    } else {
      const bool endp = false;
      ET::set_hex_volume_dof<offset, endp, C, order, order, order, EntityDof,
                             ElemDof>(entity_dof, element_dof);
    }
  }

  /**
   * @brief Set the signs for the entity
   *
   * @tparam offset Offset into the element array
   * @param entity Geometric entity type
   * @param index Index of the entity
   * @param orient Orientation of the entity relative to the reference
   * @param signs Array of sign values
   */
  template <index_t offset>
  static void set_entity_signs(ET::ElementEntity entity, index_t index,
                               index_t orient, int signs[]) {
    int sgns[ndof];
    const index_t entity_ndof = get_entity_ndof(entity, index);
    for (index_t i = 0; i < entity_ndof; i++) {
      sgns[i] = 1;
    }
    set_entity_dof<offset>(entity, index, orient, sgns, signs);
  }

  /**
   * @brief Get the parametric point location associated with the given degree
   * of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  static void get_dof_point(index_t index, double pt[]) {
    constexpr const double* pts = get_gauss_lobatto_pts<order>();

    index_t n = index / C;
    pt[0] = pts[n % order];
    pt[1] = pts[(n % order * order) / order];
    pt[2] = pts[n / (order * order)];
  }

  /**
   * @brief Interpolate the degrees of freedom to obtain the values in the space
   * object
   *
   * @tparam Quadrature The quadrature object
   * @tparam offset Degree of freedom offset into the array
   * @tparam SolnType Solution array type
   * @param n The quadrature point index
   * @param sol The solution array
   * @param out The finite element space output object
   */
  template <class Quadrature, index_t offset, class SolnType>
  static void interp(index_t n, const SolnType sol, H1Space<T, C, dim>& out) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    Vec<T, C>& u = out.get_value();
    Mat<T, C, dim>& grad = out.get_grad();

    u.zero();
    grad.zero();

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    double d1[order], d2[order], d3[order];
    lagrange_basis<order>(pt[0], n1, d1);
    lagrange_basis<order>(pt[1], n2, d2);
    lagrange_basis<order>(pt[2], n3, d3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t i = 0; i < C; i++) {
            u(i) +=
                n1[j1] * n2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];

            grad(i, 0) +=
                d1[j1] * n2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];

            grad(i, 1) +=
                n1[j1] * d2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];

            grad(i, 2) +=
                n1[j1] * n2[j2] * d3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];
          }
        }
      }
    }
  }

  /**
   * @brief Add the derivative contained in the solution space to the output
   * residual object
   *
   * @tparam Quadrature The quadrature object
   * @tparam offset Degree of freedom offset into the array
   * @tparam SolnType Solution array type
   * @param n The quadrature point index
   * @param in The finite element space output object
   * @param res The residual array - same shape as the solution array
   */
  template <class Quadrature, index_t offset, class SolnType>
  static void add(index_t n, const H1Space<T, C, dim>& in, SolnType res) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    const Vec<T, C>& u = in.get_value();
    const Mat<T, C, dim>& grad = in.get_grad();

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    double d1[order], d2[order], d3[order];
    lagrange_basis<order>(pt[0], n1, d1);
    lagrange_basis<order>(pt[1], n2, d2);
    lagrange_basis<order>(pt[2], n3, d3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t i = 0; i < C; i++) {
            res[offset + C * (j1 + j2 * order + j3 * order * order) + i] +=
                n1[j1] * n2[j2] * n3[j3] * u(i) +
                d1[j1] * n2[j2] * n3[j3] * grad(i, 0) +
                n1[j1] * d2[j2] * n3[j3] * grad(i, 1) +
                n1[j1] * n2[j2] * d3[j3] * grad(i, 2);
          }
        }
      }
    }
  }

  // Set the matrix stride
  static const index_t stride = C;

  // Set the basis size
  static const index_t basis_size = (dim + 1) * order * order * order;

  // Set the derived quantities - number of dof for each stride
  static const index_t ndof_per_stride = ndof / stride;

  // Number of components per stride
  static const index_t ncomp_per_stride = ncomp / stride;

  /**
   * @brief Evaluate the full set of basis functions for this object
   *
   * @tparam Quadrature The quadrature object
   * @tparam BasisType The type of the basis function array
   * @param n The quadrature point index
   * @param N The basis functions
   */
  template <class Quadrature, class BasisType>
  static void basis(index_t n, BasisType N) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    double d1[order], d2[order], d3[order];
    lagrange_basis<order>(pt[0], n1, d1);
    lagrange_basis<order>(pt[1], n2, d2);
    lagrange_basis<order>(pt[2], n3, d3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          const index_t node = j1 + j2 * order + j3 * order * order;
          N[(dim + 1) * node] = n1[j1] * n2[j2] * n3[j3];
          N[(dim + 1) * node + 1] = d1[j1] * n2[j2] * n3[j3];
          N[(dim + 1) * node + 2] = n1[j1] * d2[j2] * n3[j3];
          N[(dim + 1) * node + 3] = n1[j1] * n2[j2] * d3[j3];
        }
      }
    }
  }
};

template <typename T, index_t C, index_t degree>
class LagrangeL2HexBasis {
 public:
  using ET = ElementTypes;

  static const index_t dim = 3;             // Spatial dimension
  static const index_t order = degree + 1;  // Number of nodes along each edge

  static const index_t ndof =
      C * order * order * order;  // Total number of degrees of freedom

  // Number of components
  static const index_t ncomp = L2Space<T, C, dim>::ncomp;

  /**
   * @brief Degree of freedom handling on the vertices, edges, faces and volume
   *
   * @param entity The type of topological entity (vertex, edge, face or volume)
   * @param index The index of the topological entity (e.g. edge index)
   * @return The number of degrees of freedom
   */
  static index_t get_entity_ndof(ET::ElementEntity entity, index_t index) {
    if (entity == ET::VOLUME) {
      return ndof;
    }
    return 0;
  }

  /**
   * @brief Get the entity DOF from the element DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of entity (vertex, edge, face or volume)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */
  template <index_t offset, class ElemDof, class EntityDof>
  static void get_entity_dof(ET::ElementEntity entity, index_t index,
                             const ElemDof& element_dof,
                             EntityDof& entity_dof) {
    if (entity == ET::VOLUME) {
      for (index_t i = 0; i < ndof; i++) {
        entity_dof[i] = element_dof[offset + i];
      }
    }
  }

  /**
   * @brief Set the element DOF and signs from the entity DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of entity (vertex, edge, face or volume)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param entity_dof Entity DOF in the global orientation
   * @param element_dof Degrees of freedom for this element
   * @param element_sign Sign indices for each degree of freedom
   */
  template <index_t offset, class EntityDof, class ElemDof>
  static void set_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const EntityDof& entity_dof,
                             ElemDof& element_dof) {
    if (entity == ET::VOLUME) {
      for (index_t i = 0; i < ndof; i++) {
        element_dof[offset + i] = entity_dof[i];
      }
    }
  }

  /**
   * @brief Set the signs for the entity
   *
   * @tparam offset Offset into the element array
   * @param entity Geometric entity type
   * @param index Index of the entity
   * @param orient Orientation of the entity relative to the reference
   * @param signs Array of sign values
   */
  template <index_t offset>
  static void set_entity_signs(ET::ElementEntity entity, index_t index,
                               index_t orient, int signs[]) {
    int sgns[ndof];
    const index_t entity_ndof = get_entity_ndof(entity, index);
    for (index_t i = 0; i < entity_ndof; i++) {
      sgns[i] = 1;
    }
    set_entity_dof<offset>(entity, index, orient, sgns, signs);
  }

  /**
   * @brief Get the parametric point location associated with the given degree
   * of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  static void get_dof_point(index_t index, double pt[]) {
    constexpr const double* pts = get_gauss_lobatto_pts<order>();

    index_t n = index / C;
    pt[0] = pts[n % order];
    pt[1] = pts[(n % order * order) / order];
    pt[2] = pts[n / (order * order)];
  }

  template <class Quadrature, index_t offset, class SolnType>
  static void interp(index_t n, const SolnType sol, L2Space<T, C, dim>& out) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    typename L2Space<T, C, dim>::VarType u = out.get_value();
    if constexpr (C == 1) {
      u = 0.0;
    } else {
      u.zero();
    }

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    lagrange_basis<order>(pt[0], n1);
    lagrange_basis<order>(pt[1], n2);
    lagrange_basis<order>(pt[2], n3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          if constexpr (C == 1) {
            u += n1[j1] * n2[j2] * n3[j3] *
                 sol[offset + C * (j1 + j2 * order + j3 * order * order)];
          } else {
            for (index_t i = 0; i < C; i++) {
              u(i) +=
                  n1[j1] * n2[j2] * n3[j3] *
                  sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];
            }
          }
        }
      }
    }
  }

  template <class Quadrature, index_t offset, class SolnType>
  static void add(index_t n, const L2Space<T, C, dim>& in, SolnType res) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    const typename L2Space<T, C, dim>::VarType u = in.get_value();

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    lagrange_basis<order>(pt[0], n1);
    lagrange_basis<order>(pt[1], n2);
    lagrange_basis<order>(pt[2], n3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          if constexpr (C == 1) {
            res[offset + C * (j1 + j2 * order + j3 * order * order)] +=
                n1[j1] * n2[j2] * n3[j3] * u;
          } else {
            for (index_t i = 0; i < C; i++) {
              res[offset + C * (j1 + j2 * order + j3 * order * order) + i] +=
                  n1[j1] * n2[j2] * n3[j3] * u(i);
            }
          }
        }
      }
    }
  }

  // Set the matrix stride
  static const index_t stride = C;

  // Set the basis size
  static const index_t basis_size = order * order * order;

  // Set the derived quantities - number of dof for each stride
  static const index_t ndof_per_stride = ndof / stride;

  // Number of components per stride
  static const index_t ncomp_per_stride = ncomp / stride;

  // Compute the full matrix of basis functions
  template <class Quadrature, class BasisType>
  static void basis(index_t n, BasisType N) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    lagrange_basis<order>(pt[0], n1);
    lagrange_basis<order>(pt[1], n2);
    lagrange_basis<order>(pt[2], n3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          const index_t node = j1 + j2 * order + j3 * order * order;
          N[node] = n1[j1] * n2[j2] * n3[j3];
        }
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_LAGRANGE_HEX_BASIS_H