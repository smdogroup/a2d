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
  template <index_t offset>
  static void get_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const index_t element_dof[],
                             index_t entity_dof[]) {
    if (entity == ET::VERTEX) {
      index_t node = (order - 1) * ET::HEX_VERTS_CART[index][0] +
                     (order - 1) * order * ET::HEX_VERTS_CART[index][1] +
                     (order - 1) * order * order * ET::HEX_VERTS_CART[index][2];

      for (index_t i = 0; i < C; i++) {
        entity_dof[i] = element_dof[offset + C * node + i];
      }
    } else if (entity == ET::EDGE) {
      // Get the start and end location - flip the orientation if orient = 1
      // (reversed edges)
      index_t v1 = ET::HEX_EDGE_VERTS[index][orient];
      index_t v2 = ET::HEX_EDGE_VERTS[index][(orient + 1) % 2];

      index_t start = (order - 1) * ET::HEX_VERTS_CART[v1][0] +
                      (order - 1) * order * ET::HEX_VERTS_CART[v1][1] +
                      (order - 1) * order * order * ET::HEX_VERTS_CART[v1][2];
      index_t incr =
          (ET::HEX_VERTS_CART[v2][0] - ET::HEX_VERTS_CART[v1][0]) +
          order * (ET::HEX_VERTS_CART[v2][1] - ET::HEX_VERTS_CART[v1][1]) +
          order * order *
              (ET::HEX_VERTS_CART[v2][2] - ET::HEX_VERTS_CART[v1][2]);

      for (index_t k = 1; k < order - 1; k++) {
        index_t node = start + k * incr;
        for (index_t i = 0; i < C; i++) {
          entity_dof[C * (k - 1) + i] = element_dof[offset + C * node + i];
        }
      }
    } else if (entity == ET::FACE) {
      // Loop over the reference face and transform to the local face
      if (index < 2) {
        for (index_t k3 = 0; k3 < order - 2; k3++) {
          for (index_t k2 = 0; k2 < order - 2; k2++) {
            index_t j1 = (index % 2) * (order - 1);
            index_t j2 = k2 + 1;
            index_t j3 = k3 + 1;

            index_t i2, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j2, j3, &i2,
                                               &i3);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i2 - 1) + (order - 2) * (i3 - 1);
            for (index_t i = 0; i < C; i++) {
              entity_dof[C * enode + i] = element_dof[offset + C * node + i];
            }
          }
        }
      } else if (index < 4) {
        for (index_t k3 = 0; k3 < order - 2; k3++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = (index % 2) * (order - 1);
            index_t j3 = k3 + 1;

            index_t i1, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j1, j3, &i1,
                                               &i3);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i1 - 1) + (order - 2) * (i3 - 1);
            for (index_t i = 0; i < C; i++) {
              entity_dof[C * enode + i] = element_dof[offset + C * node + i];
            }
          }
        }
      } else {
        for (index_t k2 = 0; k2 < order - 2; k2++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = k2 + 1;
            index_t j3 = (index % 2) * (order - 1);

            index_t i1, i2;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j1, j2, &i1,
                                               &i2);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i1 - 1) + (order - 2) * (i2 - 1);
            for (index_t i = 0; i < C; i++) {
              entity_dof[C * enode + i] = element_dof[offset + C * node + i];
            }
          }
        }
      }
    } else if (entity == ET::VOLUME) {
      for (index_t k3 = 0; k3 < order - 2; k3++) {
        for (index_t k2 = 0; k2 < order - 2; k2++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = k2 + 1;
            index_t j3 = k3 + 1;
            index_t e = k1 + k2 * (order - 2) + k3 * (order - 2) * (order - 2);
            index_t node = j1 + j2 * order + j3 * order * order;

            for (index_t i = 0; i < C; i++) {
              entity_dof[C * e + i] = element_dof[offset + C * node + i];
            }
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
    if (entity == ET::VERTEX) {
      index_t node = (order - 1) * ET::HEX_VERTS_CART[index][0] +
                     (order - 1) * order * ET::HEX_VERTS_CART[index][1] +
                     (order - 1) * order * order * ET::HEX_VERTS_CART[index][2];

      for (index_t i = 0; i < C; i++) {
        element_dof[offset + C * node + i] = entity_dof[i];
        element_sign[offset + C * node + i] = 1;
      }
    } else if (entity == ET::EDGE) {
      // Get the start and end location - flip the orientation if orient =
      // 1 (reversed edges)
      index_t v1 = ET::HEX_EDGE_VERTS[index][orient];
      index_t v2 = ET::HEX_EDGE_VERTS[index][(orient + 1) % 2];

      index_t start = (order - 1) * ET::HEX_VERTS_CART[v1][0] +
                      (order - 1) * order * ET::HEX_VERTS_CART[v1][1] +
                      (order - 1) * order * order * ET::HEX_VERTS_CART[v1][2];
      index_t incr =
          (ET::HEX_VERTS_CART[v2][0] - ET::HEX_VERTS_CART[v1][0]) +
          order * (ET::HEX_VERTS_CART[v2][1] - ET::HEX_VERTS_CART[v1][1]) +
          order * order *
              (ET::HEX_VERTS_CART[v2][2] - ET::HEX_VERTS_CART[v1][2]);

      for (index_t k = 1; k < order - 1; k++) {
        index_t node = start + k * incr;
        for (index_t i = 0; i < C; i++) {
          element_dof[offset + C * node + i] = entity_dof[C * (k - 1) + i];
          element_sign[offset + C * node + i] = 1;
        }
      }
    } else if (entity == ET::FACE) {
      // Loop over the reference face and transform to the local face
      if (index < 2) {
        for (index_t k3 = 0; k3 < order - 2; k3++) {
          for (index_t k2 = 0; k2 < order - 2; k2++) {
            index_t j1 = (index % 2) * (order - 1);
            index_t j2 = k2 + 1;
            index_t j3 = k3 + 1;

            index_t i2, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j2, j3, &i2,
                                               &i3);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i2 - 1) + (order - 2) * (i3 - 1);
            for (index_t i = 0; i < C; i++) {
              element_dof[offset + C * node + i] = entity_dof[C * enode + i];
              element_sign[offset + C * node + i] = 1;
            }
          }
        }
      } else if (index < 4) {
        for (index_t k3 = 0; k3 < order - 2; k3++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = (index % 2) * (order - 1);
            index_t j3 = k3 + 1;

            index_t i1, i3;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j1, j3, &i1,
                                               &i3);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i1 - 1) + (order - 2) * (i3 - 1);
            for (index_t i = 0; i < C; i++) {
              element_dof[offset + C * node + i] = entity_dof[C * enode + i];
              element_sign[offset + C * node + i] = 1;
            }
          }
        }
      } else {
        for (index_t k2 = 0; k2 < order - 2; k2++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = k2 + 1;
            index_t j3 = (index % 2) * (order - 1);

            index_t i1, i2;
            ET::get_coords_on_quad_ref_element(orient, order - 1, j1, j2, &i1,
                                               &i2);

            // Get the node on the face
            index_t node = j1 + order * j2 + order * order * j3;
            index_t enode = (i1 - 1) + (order - 2) * (i2 - 1);
            for (index_t i = 0; i < C; i++) {
              element_dof[offset + C * node + i] = entity_dof[C * enode + i];
              element_sign[offset + C * node + i] = 1;
            }
          }
        }
      }
    } else if (entity == ET::VOLUME) {
      for (index_t k3 = 0; k3 < order - 2; k3++) {
        for (index_t k2 = 0; k2 < order - 2; k2++) {
          for (index_t k1 = 0; k1 < order - 2; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = k2 + 1;
            index_t j3 = k3 + 1;
            index_t e = k1 + k2 * (order - 2) + k3 * (order - 2) * (order - 2);
            index_t node = j1 + j2 * order + j3 * order * order;

            for (index_t i = 0; i < C; i++) {
              element_dof[offset + C * node + i] = entity_dof[C * e + i];
              element_sign[offset + C * node + i] = 1;
            }
          }
        }
      }
    }
  }

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

            grad(i, 0) =
                d1[j1] * n2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];

            grad(i, 1) =
                n1[j1] * d2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];

            grad(i, 2) =
                n1[j1] * n2[j2] * d3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];
          }
        }
      }
    }
  }

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

  // Compute the full matrix of basis functions
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
  template <index_t offset>
  static void get_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const index_t element_dof[],
                             index_t entity_dof[]) {
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
  template <index_t offset>
  static void set_entity_dof(ET::ElementEntity entity, index_t index,
                             index_t orient, const index_t entity_dof[],
                             index_t element_dof[], int element_sign[]) {
    if (entity == ET::VOLUME) {
      for (index_t i = 0; i < ndof; i++) {
        element_dof[offset + i] = entity_dof[i];
        element_sign[offset + i] = 1;
      }
    }
  }

  template <class Quadrature, index_t offset, class SolnType>
  static void interp(index_t n, const SolnType sol, L2Space<T, C, dim>& out) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    Vec<T, C>& u = out.get_value();

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    lagrange_basis<order>(pt[0], n1);
    lagrange_basis<order>(pt[1], n2);
    lagrange_basis<order>(pt[2], n3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t i = 0; i < C; i++) {
            u(i) +=
                n1[j1] * n2[j2] * n3[j3] *
                sol[offset + C * (j1 + j2 * order + j3 * order * order) + i];
          }
        }
      }
    }
  }

  template <class Quadrature, index_t offset, class SolnType>
  static void add(index_t n, const L2Space<T, C, dim>& in, SolnType res) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    const Vec<T, C>& u = in.get_value();
    u.zero();

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
                n1[j1] * n2[j2] * n3[j3] * u(i);
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
    double d1[order], d2[order], d3[order];
    lagrange_basis<order>(pt[0], n1, d1);
    lagrange_basis<order>(pt[1], n2, d2);
    lagrange_basis<order>(pt[2], n3, d3);

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