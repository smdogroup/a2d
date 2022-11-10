#ifndef A2D_HEXA_BASIS_H
#define A2D_HEXA_BASIS_H

#include "multiphysics/febasis.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

template <typename T, index_t C, index_t degree>
class LagrangeHexBasis {
 public:
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
  index_t get_entity_ndof(TopoEntity entity, index_t index) {
    if (entity == VERTEX) {
      return C;
    } else if (entity == EDGE) {
      return C * (order - 2);
    } else if (entity == FACE) {
      return C * (order - 2) * (order - 2);
    } else if (entity == VOLUME) {
      return C * (order - 2) * (order - 2) * (order - 2);
    }
  }

  /**
   * @brief Get the DOF from the element dof
   *
   * @tparam offset The offset into the basis
   * @param entity The type of topological entity (vertex, edge, face or volume)
   * @param index The index of the topological entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */

  template <index_t offset>
  void get_dof(TopoEntity entity, index_t index, index_t orient,
               const index_t element_dof[], index_t entity_dof[]) {
    if (entity == VERTEX) {
      for (index i = 0; i < C; i++) {
        entity_dof[i] =
            element_dof[offset +
                        C * (j1 + j2 * order + j3 * order * order * order) + i];
      }
    } else if (entity == EDGE) {
      for (index i = 0; i < C; i++) {
        entity_dof[i] =
            element_dof[offset +
                        C * (j1 + j2 * order + j3 * order * order * order) + i];
      }
    } else if (entity == FACE) {
      entity_dof[i] =
          element_dof[offset +
                      C * (j1 + j2 * order + j3 * order * order * order) + i];

    } else if (entity == VOLUME) {
      for (index_t k3 = 0; k3 < order; k3++) {
        for (index_t k2 = 0; k2 < order; k2++) {
          for (index_t k1 = 0; k1 < order; k1++) {
            index_t j1 = k1 + 1;
            index_t j2 = k2 + 1;
            index_t j3 = k3 + 1;
            index_t e = k1 + k2 * (order - 2) + k3 * (order - 2) * (order - 2);

            entity_dof[e] =
                element_dof[offset +
                            C * (j1 + j2 * order + j3 * order * order * order) +
                            i];
          }
        }
      }
    }
  }

  template <index_t offset>
  void set_dof(TopoEntity domain, index_t index, index_t orient,
               index_t entity_dof[], index_t element_dof[],
               int element_sign[]) {}

  template <class Quadrature, index_t offset, class SolnType>
  static void interp(index_t n, const SolnType sol, H1Space<T, C, dim>& out) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    Vec<T, C>& u = out.get_value();
    Mat<T, C, dim>& grad = out.get_grad();

    // Evaluate the basis functions
    double n1[order], n2[order], n3[order];
    double d1[order], d2[order], d3[order];
    lagrange_basis<order>(pt[0], n1, d1);
    lagrange_basis<order>(pt[1], n2, d2);
    lagrange_basis<order>(pt[2], n3, d3);

    for (index_t j3 = 0; j3 < order; j3++) {
      for (index_t j2 = 0; j2 < order; j2++) {
        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t i = 0; i < C : i++) {
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
          for (index_t i = 0; i < C : i++) {
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

}  // namespace A2D

#endif  // A2D_HEXA_BASIS_H