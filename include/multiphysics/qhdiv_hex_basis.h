#ifndef A2D_QHDIV_HEX_BASIS_H
#define A2D_QHDIV_HEX_BASIS_H

#include "multiphysics/febasis.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

static const double QHdivKnots1[] = {0.0};

static const double QHdivKnots2[] = {-0.5, 0.5};

static const double QHdivKnots3[] = {-2.0 / 3.0, 0.0, 2.0 / 3.0};

static const double QHdivKnots4[] = {-0.75, -0.25, 0.25, 0.75};

template <index_t degree>
inline const double* get_qhdiv_knots() {
  if (degree == 1) {
    return QHdivKnots1;
  } else if (degree == 2) {
    return QHdivKnots2;
  } else if (degree == 3) {
    return QHdivKnots3;
  } else if (degree == 4) {
    return QHdivKnots4;
  }
  return NULL;
}

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
  template <index_t offset>
  static void get_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const index_t element_dof[],
                             index_t entity_dof[]) {
    if (entity == ET::FACE) {
      // Loop over the reference face and transform to the local face
      if (index < 2) {
        for (index_t j3 = 0; j3 < order - 1; j3++) {
          for (index_t j2 = 0; j2 < order - 1; j2++) {
            index_t j1 = (index % 2) * (order - 1);

            index_t i2, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j2, j3, &i2,
                                               &i3);

            // Get the node on the face
            index_t node1 = offset + j1 + order * j2 + (order - 1) * order * j3;
            index_t enode = i2 + (order - 1) * i3;
            entity_dof[enode] = element_dof[node1];
          }
        }
      } else if (index < 4) {
        for (index_t j3 = 0; j3 < order - 1; j3++) {
          for (index_t j1 = 0; j1 < order - 1; j1++) {
            index_t j2 = (index % 2) * (order - 1);

            index_t i1, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j1, j3, &i1,
                                               &i3);

            // Get the node on the face
            index_t node2 = (offset + order * (order - 1) * (order - 1)) +
                            (j1 + (order - 1) * j2 + order * (order - 1) * j3);
            index_t enode = i1 + (order - 1) * i3;
            entity_dof[enode] = element_dof[node2];
          }
        }
      } else {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++) {
            index_t j3 = (index % 2) * (order - 1);

            index_t i1, i2;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j1, j2, &i1,
                                               &i2);

            // Get the node on the face
            index_t node3 =
                (offset + 2 * order * (order - 1) * (order - 1)) +
                (j1 + (order - 1) * j2 + (order - 1) * (order - 1) * j3);
            index_t enode = i1 + (order - 1) * i2;
            entity_dof[enode] = element_dof[node3];
          }
        }
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
   * @param element_sign Sign indices for each degree of freedom
   */
  template <index_t offset>
  static void set_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const index_t entity_dof[],
                             index_t element_dof[], int element_sign[]) {
    if (entity == ET::FACE) {
      // Loop over the reference face and transform to the local face
      if (index < 2) {
        for (index_t j3 = 0; j3 < order - 1; j3++) {
          for (index_t j2 = 0; j2 < order - 1; j2++) {
            index_t j1 = (index % 2) * (order - 1);

            index_t i2, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j2, j3, &i2,
                                               &i3);

            // Get the node on the face
            index_t node1 = offset + j1 + order * j2 + (order - 1) * order * j3;
            index_t enode = i2 + (order - 1) * i3;
            element_dof[node1] = entity_dof[enode];

            // Check if the orientation is flipped relative to the reference
            if (orient >= 4) {
              element_sign[node1] = -1;
            } else {
              element_sign[node1] = 1;
            }
          }
        }
      } else if (index < 4) {
        for (index_t j3 = 0; j3 < order - 1; j3++) {
          for (index_t j1 = 0; j1 < order - 1; j1++) {
            index_t j2 = (index % 2) * (order - 1);

            index_t i1, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j1, j3, &i1,
                                               &i3);

            // Get the node on the face
            index_t node2 = (offset + order * (order - 1) * (order - 1)) +
                            (j1 + (order - 1) * j2 + order * (order - 1) * j3);
            index_t enode = i1 + (order - 1) * i3;
            element_dof[node2] = entity_dof[enode];

            // Check if the orientation is flipped relative to the reference
            if (orient >= 4) {
              element_sign[node2] = -1;
            } else {
              element_sign[node2] = 1;
            }
          }
        }
      } else {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++) {
            index_t j3 = (index % 2) * (order - 1);

            index_t i1, i2;
            ET::get_coords_on_quad_ref_element(orient, order - 2, j1, j2, &i1,
                                               &i2);

            // Get the node on the face
            index_t node3 =
                (offset + 2 * order * (order - 1) * (order - 1)) +
                (j1 + (order - 1) * j2 + (order - 1) * (order - 1) * j3);
            index_t enode = i1 + (order - 1) * i2;
            element_dof[node3] = entity_dof[enode];

            // Check if the orientation is flipped relative to the reference
            if (orient >= 4) {
              element_sign[node3] = -1;
            } else {
              element_sign[node3] = 1;
            }
          }
        }
      }
    } else if (entity == ET::VOLUME) {
      index_t enode = 0;

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 1; j1 < order - 1; j1++, enode++) {
            index_t node1 = offset + j1 + order * j2 + (order - 1) * order * j3;

            element_dof[node1] = entity_dof[enode];
            element_sign[node1] = 1;
          }
        }
      }

      for (index_t j3 = 0; j3 < order - 1; j3++) {
        for (index_t j2 = 1; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, enode++) {
            index_t node2 = (offset + order * (order - 1) * (order - 1)) +
                            (j1 + (order - 1) * j2 + order * (order - 1) * j3);

            element_dof[node2] = entity_dof[enode];
            element_sign[node2] = 1;
          }
        }
      }

      for (index_t j3 = 0; j3 < order; j3++) {
        for (index_t j2 = 0; j2 < order - 1; j2++) {
          for (index_t j1 = 0; j1 < order - 1; j1++, enode++) {
            index_t node3 =
                (offset + 2 * order * (order - 1) * (order - 1)) +
                (j1 + (order - 1) * j2 + (order - 1) * (order - 1) * j3);

            element_dof[node3] = entity_dof[enode];
            element_sign[node3] = 1;
          }
        }
      }
    }
  }

  template <class Quadrature, index_t offset, class SolnType>
  static void interp(index_t n, const SolnType sol, HdivSpace<T, dim>& out) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    Vec<T, dim>& u = out.get_value();
    T& div = out.get_div();

    u.zero();
    div = 0.0;

    const double* knots = get_qhdiv_knots<degree>();

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

  template <class Quadrature, index_t offset, class SolnType>
  static void add(index_t n, const HdivSpace<T, dim>& in, SolnType res) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    const Vec<T, dim>& u = in.get_value();
    const T& div = in.get_div();

    const double* knots = get_qhdiv_knots<degree>();

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

    const double* knots = get_qhdiv_knots<degree>();

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