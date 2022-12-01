#ifndef A2D_QHDIV_HEX_BASIS_H
#define A2D_QHDIV_HEX_BASIS_H

#include "multiphysics/febasis.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

template <typename T, index_t degree>
class QHdivHexBasis {
 public:
  using ET = ElementTypes;

  static const index_t dim = 3;             // Spatial dimension
  static const index_t order = degree + 1;  // Number of nodes along each edge

  // Total number of degrees of freedom
  static const index_t ndof = 3 * order * (order - 1) * (order - 1);

  // Number of components
  static const index_t ncomp = HdivSpace<T, dim>::ncomp;

  // Get the type of basis class implemented
  static constexpr BasisType get_basis_type() { return HDIV; }

  // Define the equivalent low-order basis class if any
  using LOrderBasis = QHdivHexBasis<T, 1>;

  /**
   * @brief Degree of freedom handling on the vertices, edges, faces and volume
   *
   * @param entity The type of topological entity (vertex, edge, face or volume)
   * @param index The index of the topological entity (e.g. edge index)
   * @return The number of degrees of freedom
   */
  static index_t get_entity_ndof(ET::ElementEntity entity, index_t index) {
    if (entity == ET::FACE) {
      return (order - 1) * (order - 1);
    } else if (entity == ET::VOLUME) {
      return 3 * (order - 2) * (order - 1) * (order - 1);
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
    if (entity == ET::FACE) {
      if (index < 2) {
        const bool endp = true;
        const index_t ndof = 1;
        ET::get_hex_face_dof<offset, endp, ndof, order, order - 1, order - 1,
                             ElemDof, EntityDof>(index, element_dof,
                                                 entity_dof);
      } else if (index < 4) {
        const index_t off = offset + order * (order - 1) * (order - 1);
        const bool endp = true;
        const index_t ndof = 1;
        ET::get_hex_face_dof<off, endp, ndof, order - 1, order, order - 1,
                             ElemDof, EntityDof>(index, element_dof,
                                                 entity_dof);
      } else {
        const index_t off = offset + 2 * order * (order - 1) * (order - 1);
        const bool endp = true;
        const index_t ndof = 1;
        ET::get_hex_face_dof<off, endp, ndof, order - 1, order - 1, order,
                             ElemDof, EntityDof>(index, element_dof,
                                                 entity_dof);
      }
    } else if (entity == ET::VOLUME) {
      index_t enode = 0;
      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 1; j2 < order - 1; j1++, enode++) {
            index_t node1 = offset + j1 + order * j2 + (order - 1) * order * j3;
            entity_dof[enode] = element_dof[node1];
          }
        }
      }

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 1; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j2 < order - 1; j1++, enode++) {
            index_t node2 = (offset + order * (order - 1) * (order - 1)) +
                            (j1 + (order - 1) * j2 + order * (order - 1) * j3);
            entity_dof[enode] = element_dof[node2];
          }
        }
      }

      for (index_t j3 = 0; j3 < order; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j2 < order - 1; j1++, enode++) {
            index_t node3 =
                (offset + 2 * order * (order - 1) * (order - 1)) +
                (j1 + (order - 1) * j2 + (order - 1) * (order - 1) * j3);
            entity_dof[enode] = element_dof[node3];
          }
        }
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
   */
  template <index_t offset, class EntityDof, class ElemDof>
  static void set_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const EntityDof& entity_dof,
                             ElemDof& element_dof) {
    if (entity == ET::FACE) {
      if (index < 2) {
        const bool endp = true;
        const index_t ndof = 1;
        ET::set_hex_face_dof<offset, endp, ndof, order, order - 1, order - 1,
                             EntityDof, ElemDof>(index, orient, entity_dof,
                                                 element_dof);
      } else if (index < 4) {
        const index_t off = offset + order * (order - 1) * (order - 1);
        const bool endp = true;
        const index_t ndof = 1;
        ET::set_hex_face_dof<off, endp, ndof, order - 1, order, order - 1,
                             EntityDof, ElemDof>(index, orient, entity_dof,
                                                 element_dof);
      } else {
        const index_t off = offset + 2 * order * (order - 1) * (order - 1);
        const bool endp = true;
        const index_t ndof = 1;
        ET::set_hex_face_dof<off, endp, ndof, order - 1, order - 1, order,
                             EntityDof, ElemDof>(index, orient, entity_dof,
                                                 element_dof);
      }
    } else if (entity == ET::VOLUME) {
      index_t enode = 0;

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 1; j1 < order - 1; j1++, enode++) {
            index_t node1 = offset + j1 + order * j2 + (order - 1) * order * j3;
            element_dof[node1] = entity_dof[enode];
          }
        }
      }

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 1; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, enode++) {
            index_t node2 = (offset + order * (order - 1) * (order - 1)) +
                            (j1 + (order - 1) * j2 + order * (order - 1) * j3);
            element_dof[node2] = entity_dof[enode];
          }
        }
      }

      for (index_t j3 = 1; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, enode++) {
            index_t node3 =
                (offset + 2 * order * (order - 1) * (order - 1)) +
                (j1 + (order - 1) * j2 + (order - 1) * (order - 1) * j3);
            element_dof[node3] = entity_dof[enode];
          }
        }
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
    for (index_t k = 0; k < entity_ndof; k++) {
      sgns[k] = 1;
    }

    if (entity == ET::FACE) {
      if (orient >= 4) {
        for (index_t k = 0; k < entity_ndof; k++) {
          sgns[k] = -1;
        }
      }
    }

    set_entity_dof<offset>(entity, index, orient, sgns, signs);
  }

  /**
   * @brief Get the number of low order elements defined by this basis
   */
  constexpr static index_t get_num_lorder_elements() {
    return degree * degree * degree;
  }

  /**
   * @brief Get the low order degrees of freedom associated with this element
   *
   * @tparam horder_offset Offset into the high-order dof
   * @tparam lorder_offset Offset into the low-order dof
   * @tparam HOrderDof High-order dof array type
   * @tparam LOrderDof Low-order dof array type
   * @param n Index of the low order element to be projected
   * @param hdof High-order dof input
   * @param ldof Low-order dof output
   */
  template <index_t horder_offset, index_t lorder_offset, class HOrderDof,
            class LOrderDof>
  static void get_lorder_dof(const index_t n, const HOrderDof& hdof,
                             LOrderDof& ldof) {
    const index_t i = n % degree;
    const index_t j = (n % (degree * degree)) / degree;
    const index_t k = n / (degree * degree);

    // 0-direction
    ldof[lorder_offset] =
        hdof[horder_offset + i + order * (j + (order - 1) * k)];
    ldof[lorder_offset + 1] =
        hdof[horder_offset + i + 1 + order * (j + (order - 1) * k)];

    // 1-direction
    ldof[lorder_offset + 2] =
        hdof[horder_offset + i + (order - 1) * (j + (order - 1) * k) +
             order * (order - 1) * (order - 1)];
    ldof[lorder_offset + 3] =
        hdof[horder_offset + i + (order - 1) * (j + 1 + (order - 1) * k) +
             order * (order - 1) * (order - 1)];

    // 2-direction
    ldof[lorder_offset + 4] =
        hdof[horder_offset + i + (order - 1) * (j + (order - 1) * k) +
             2 * order * (order - 1) * (order - 1)];
    ldof[lorder_offset + 5] =
        hdof[horder_offset + i + (order - 1) * (j + (order - 1) * (k + 1)) +
             2 * order * (order - 1) * (order - 1)];
  }

  /**
   * @brief Get the low-order signs relative to the high-order degrees of
   * freedom
   *
   * @tparam horder_offset Offset into the high-order dof
   * @tparam lorder_offset Offset into the low-order dof
   * @param n Index of the low order element to be projected
   * @param horder_signs The signs for the high-order degrees of freedom
   * @param signs The signs relative to the high-order element
   */
  template <index_t horder_offset, index_t lorder_offset>
  static void get_lorder_signs(const index_t n, const int horder_signs[],
                               int signs[]) {
    const index_t i = n % degree;
    const index_t j = (n % (degree * degree)) / degree;
    const index_t k = n / (degree * degree);

    // 0-direction
    signs[lorder_offset] =
        horder_signs[horder_offset + i + order * (j + (order - 1) * k)];
    if (i > 0) {
      signs[lorder_offset] *= -1;
    }
    signs[lorder_offset + 1] =
        horder_signs[horder_offset + i + 1 + order * (j + (order - 1) * k)];

    // 1-direction
    signs[lorder_offset + 2] =
        horder_signs[horder_offset + i + (order - 1) * (j + (order - 1) * k) +
                     order * (order - 1) * (order - 1)];
    if (j > 0) {
      signs[lorder_offset + 2] *= -1;
    }
    signs[lorder_offset + 3] =
        horder_signs[horder_offset + i +
                     (order - 1) * (j + 1 + (order - 1) * k) +
                     order * (order - 1) * (order - 1)];

    // 2-direction
    signs[lorder_offset + 4] =
        horder_signs[horder_offset + i + (order - 1) * (j + (order - 1) * k) +
                     2 * order * (order - 1) * (order - 1)];
    if (k > 0) {
      signs[lorder_offset + 4] *= -1;
    }
    signs[lorder_offset + 5] =
        horder_signs[horder_offset + i +
                     (order - 1) * (j + (order - 1) * (k + 1)) +
                     2 * order * (order - 1) * (order - 1)];
  }

  /**
   * @brief Get the parametric point location associated with the given degree
   * of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  static void get_dof_point(index_t index, double pt[]) {
    const double* knots = get_gauss_quadrature_pts<degree>();
    constexpr const double* pts = get_gauss_lobatto_pts<order>();

    if (index < order * (order - 1) * (order - 1)) {
      pt[0] = pts[index % order];
      pt[1] = knots[(index % (order * (order - 1))) / order];
      pt[2] = knots[index / (order * (order - 1))];
    } else if (index < 2 * order * (order - 1) * (order - 1)) {
      index = index - order * (order - 1) * (order - 1);
      pt[0] = knots[index % (order - 1)];
      pt[1] = pts[(index % ((order - 1) * order)) / (order - 1)];
      pt[2] = knots[index / ((order - 1) * order)];
    } else {
      index = index - order * (order - 1) * (order - 1);
      pt[0] = knots[index % (order - 1)];
      pt[1] = knots[(index % ((order - 1) * (order - 1))) / (order - 1)];
      pt[2] = pts[index / ((order - 1) * (order - 1))];
    }
  }

  /**
   * @brief Interpolate over all quadrature points on the element
   *
   * @tparam space The finite element space index
   * @tparam Quadrature The quadrature scheme
   * @tparam FiniteElementSpace The finite element space object
   * @tparam offset Offset index into the solution degrees of freedom
   * @tparam SolnType
   * @param sol The solution array
   * @param out The finite element space object at all quadrature points
   */
  template <index_t space, class Quadrature, class FiniteElementSpace,
            index_t offset, class SolnType>
  static void interp(const SolnType& sol,
                     QptSpace<Quadrature, FiniteElementSpace>& out) {
    for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
      // Get the quadrature point
      double pt[dim];
      Quadrature::get_point(q, pt);

      FiniteElementSpace& s = out.get(q);
      HdivSpace<T, dim>& hdiv = s.template get<space>();
      Vec<T, dim>& u = hdiv.get_value();
      T& div = hdiv.get_div();

      u.zero();
      div = 0.0;

      const double* knots = get_gauss_quadrature_pts<degree>();

      // Evaluate the basis functions
      double dx[order];
      double n1[order], n2[order], n3[order];
      lagrange_basis<order>(pt[0], n1, dx);
      lagrange_basis<order - 1>(knots, pt[1], n2);
      lagrange_basis<order - 1>(knots, pt[2], n3);

      // Flip the first basis function on the negative face
      n1[0] *= -1.0;
      dx[0] *= -1.0;

      index_t node = offset;
      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order; j1++, node++) {
            u(0) += n1[j1] * n2[j2] * n3[j3] * sol[node];
            div += dx[j1] * n2[j2] * n3[j3] * sol[node];
          }
        }
      }

      lagrange_basis<order - 1>(knots, pt[0], n1);
      lagrange_basis<order>(pt[1], n2, dx);
      lagrange_basis<order - 1>(knots, pt[2], n3);

      // Flip the first basis function on the negative face
      n2[0] *= -1.0;
      dx[0] *= -1.0;

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
            u(1) += n1[j1] * n2[j2] * n3[j3] * sol[node];
            div += n1[j1] * dx[j2] * n3[j3] * sol[node];
          }
        }
      }

      lagrange_basis<order - 1>(knots, pt[0], n1);
      lagrange_basis<order - 1>(knots, pt[1], n2);
      lagrange_basis<order>(pt[2], n3, dx);

      // Flip the first basis function on the negative face
      n3[0] *= -1.0;
      dx[0] *= -1.0;

      for (index_t j3 = 0; j3 < order; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
            u(2) += n1[j1] * n2[j2] * n3[j3] * sol[node];
            div += n1[j1] * n2[j2] * dx[j3] * sol[node];
          }
        }
      }
    }
  }

  /**
   * @brief Add the derivative contained in the solution space to the output
   * residual object
   *
   * @tparam space The finite element space index
   * @tparam Quadrature The quadrature object
   * @tparam FiniteElementSpace
   * @tparam offset Degree of freedom offset into the array
   * @tparam SolnType Solution array type
   * @param in The finite element space output object
   * @param res The residual array - same shape as the solution array
   */
  template <index_t space, class Quadrature, class FiniteElementSpace,
            index_t offset, class SolnType>
  static void add(const QptSpace<Quadrature, FiniteElementSpace>& in,
                  SolnType& res) {
    for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
      // Get the quadrature point
      double pt[dim];
      Quadrature::get_point(q, pt);

      const FiniteElementSpace& s = in.get(q);
      const HdivSpace<T, dim>& hdiv = s.template get<space>();
      const Vec<T, dim>& u = hdiv.get_value();
      const T& div = hdiv.get_div();

      const double* knots = get_gauss_quadrature_pts<degree>();

      // Evaluate the basis functions
      double dx[order];
      double n1[order], n2[order], n3[order];
      lagrange_basis<order>(pt[0], n1, dx);
      lagrange_basis<order - 1>(knots, pt[1], n2);
      lagrange_basis<order - 1>(knots, pt[2], n3);

      // Flip the first basis function on the negative face
      n1[0] *= -1.0;
      dx[0] *= -1.0;

      index_t node = offset;
      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order; j1++, node++) {
            res[node] += n1[j1] * n2[j2] * n3[j3] * u(0);
            res[node] += dx[j1] * n2[j2] * n3[j3] * div;
          }
        }
      }

      lagrange_basis<order - 1>(knots, pt[0], n1);
      lagrange_basis<order>(pt[1], n2, dx);
      lagrange_basis<order - 1>(knots, pt[2], n3);

      // Flip the first basis function on the negative face
      n2[0] *= -1.0;
      dx[0] *= -1.0;

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
            res[node] += n1[j1] * n2[j2] * n3[j3] * u(1);
            res[node] += n1[j1] * dx[j2] * n3[j3] * div;
          }
        }
      }

      lagrange_basis<order - 1>(knots, pt[0], n1);
      lagrange_basis<order - 1>(knots, pt[1], n2);
      lagrange_basis<order>(pt[2], n3, dx);

      // Flip the first basis function on the negative face
      n3[0] *= -1.0;
      dx[0] *= -1.0;

      for (index_t j3 = 0; j3 < order; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
            res[node] += n1[j1] * n2[j2] * n3[j3] * u(2);
            res[node] += n1[j1] * n2[j2] * dx[j3] * div;
          }
        }
      }
    }
  }

  // Set the matrix stride
  static const index_t stride = 1;

  // Set the basis size
  static const index_t basis_size = (dim + 1) * ndof;

  // Set the derived quantities - number of dof for each stride
  static const index_t ndof_per_stride = ndof / stride;

  // Number of components per stride
  static const index_t ncomp_per_stride = ncomp / stride;

  // Compute the full matrix of basis functions
  template <class Quadrature, class BasisType>
  static void basis(index_t n, BasisType N) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    const double* knots = get_gauss_quadrature_pts<degree>();

    // Evaluate the basis functions
    double dx[order];
    double n1[order], n2[order], n3[order];
    lagrange_basis<order>(pt[0], n1, dx);
    lagrange_basis<order - 1>(knots, pt[1], n2);
    lagrange_basis<order - 1>(knots, pt[2], n3);

    // Flip the first basis function on the negative face
    n1[0] *= -1.0;
    dx[0] *= -1.0;

    index_t node = 0;
    for (index_t j3 = 0; j3 < order - 1; j3++) {
      for (index_t j2 = 0; j2 < order - 1; j2++) {
        for (index_t j1 = 0; j1 < order; j1++, node++) {
          N[4 * node] = n1[j1] * n2[j2] * n3[j3];
          N[4 * node + 1] = 0.0;
          N[4 * node + 2] = 0.0;
          N[4 * node + 3] = dx[j1] * n2[j2] * n3[j3];
        }
      }
    }

    lagrange_basis<order - 1>(knots, pt[0], n1);
    lagrange_basis<order>(pt[1], n2, dx);
    lagrange_basis<order - 1>(knots, pt[2], n3);

    // Flip the first basis function on the negative face
    n2[0] *= -1.0;
    dx[0] *= -1.0;

    for (index_t j3 = 0; j3 < order - 1; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
          N[4 * node] = 0.0;
          N[4 * node + 1] = n1[j1] * n2[j2] * n3[j3];
          N[4 * node + 2] = 0.0;
          N[4 * node + 3] = n1[j1] * dx[j2] * n3[j3];
        }
      }
    }

    lagrange_basis<order - 1>(knots, pt[0], n1);
    lagrange_basis<order - 1>(knots, pt[1], n2);
    lagrange_basis<order>(pt[2], n3, dx);

    // Flip the first basis function on the negative face
    n3[0] *= -1.0;
    dx[0] *= -1.0;

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order - 1; j2++) {
        for (index_t j1 = 0; j1 < order - 1; j1++, node++) {
          N[4 * node] = 0.0;
          N[4 * node + 1] = 0.0;
          N[4 * node + 2] = n1[j1] * n2[j2] * n3[j3];
          N[4 * node + 3] = n1[j1] * dx[j2] * n3[j3];
        }
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_QHDIV_HEX_BASIS_H