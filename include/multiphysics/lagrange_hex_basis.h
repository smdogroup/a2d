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

  static constexpr BasisType get_basis_type() { return H1; }

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
    pt[1] = pts[(n % (order * order)) / order];
    pt[2] = pts[n / (order * order)];
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
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      const index_t q2dim = Quadrature::tensor_dim2;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 0-direction
        T u0[order * order * q0dim];
        T u0x[order * order * q0dim];
        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          lagrange_basis<order>(pt0, n0, d0);

          for (index_t j2 = 0; j2 < order; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val(0.0), derx(0.0);
              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                val += n0[j0] * sol[offset + C * node + i];
                derx += d0[j0] * sol[offset + C * node + i];
              }

              u0[j1 + order * (j2 + order * q0)] = val;
              u0x[j1 + order * (j2 + order * q0)] = derx;
            }
          }
        }

        // Interpolate along the 1-direction
        T u1[order * q0dim * q1dim];
        T u1x[order * q0dim * q1dim];
        T u1y[order * q0dim * q1dim];
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order], d1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1, d1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            for (index_t j2 = 0; j2 < order; j2++) {
              T val(0.0), derx(0.0), dery(0.0);
              for (index_t j1 = 0; j1 < order; j1++) {
                val += n1[j1] * u0[j1 + order * (j2 + order * q0)];
                derx += n1[j1] * u0x[j1 + order * (j2 + order * q0)];
                dery += d1[j1] * u0[j1 + order * (j2 + order * q0)];
              }

              u1[j2 + order * (q0 + q0dim * q1)] = val;
              u1x[j2 + order * (q0 + q0dim * q1)] = derx;
              u1y[j2 + order * (q0 + q0dim * q1)] = dery;
            }
          }
        }

        // Interpolate along the 2-direction
        for (index_t q2 = 0; q2 < q2dim; q2++) {
          double n2[order], d2[order];
          const double pt2 = Quadrature::get_tensor_point(2, q2);
          lagrange_basis<order>(pt2, n2, d2);

          for (index_t q1 = 0; q1 < q1dim; q1++) {
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              T val(0.0), derx(0.0), dery(0.0), derz(0.0);
              for (index_t j2 = 0; j2 < order; j2++) {
                val += n2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
                derx += n2[j2] * u1x[j2 + order * (q0 + q0dim * q1)];
                dery += n2[j2] * u1y[j2 + order * (q0 + q0dim * q1)];
                derz += d2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
              }

              const index_t qindex = Quadrature::get_tensor_index(q0, q1, q2);
              FiniteElementSpace& s = out.get(qindex);
              H1Space<T, C, dim>& h1 = s.template get<space>();
              typename H1Space<T, C, dim>::VarType& u = h1.get_value();
              typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

              if constexpr (C == 1) {
                u = val;
                grad(0) = derx;
                grad(1) = dery;
                grad(2) = derz;

              } else {
                u(i) = val;
                grad(i, 0) = derx;
                grad(i, 1) = dery;
                grad(i, 2) = derz;
              }
            }
          }
        }
      }
    } else {
      for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
        // Get the quadrature point
        double pt[dim];
        Quadrature::get_point(q, pt);

        // Evaluate the basis functions
        double n0[order], d0[order];
        double n1[order], d1[order];
        double n2[order], d2[order];
        lagrange_basis<order>(pt[0], n0, d0);
        lagrange_basis<order>(pt[1], n1, d1);
        lagrange_basis<order>(pt[2], n2, d2);

        FiniteElementSpace& s = out.get(q);
        H1Space<T, C, dim>& h1 = s.template get<space>();
        typename H1Space<T, C, dim>::VarType& u = h1.get_value();
        typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

        if constexpr (C == 1) {
          u = 0.0;
        } else {
          u.zero();
        }
        grad.zero();

        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order; j1++) {
            double n1n2 = n1[j1] * n2[j2];
            double d1n2 = d1[j1] * n2[j2];
            double n1d2 = n1[j1] * d2[j2];
            for (index_t j0 = 0; j0 < order; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double N = n0[j0] * n1n2;
              double dx = d0[j0] * n1n2;
              double dy = n0[j0] * d1n2;
              double dz = n0[j0] * n1d2;

              if constexpr (C == 1) {
                const T val = sol[offset + node];

                u += N * val;
                grad(0) += dx * val;
                grad(1) += dy * val;
                grad(2) += dz * val;
              } else {
                for (index_t i = 0; i < C; i++) {
                  const T val = sol[offset + C * node + i];

                  u(i) += N * val;
                  grad(i, 0) += dx * val;
                  grad(i, 1) += dy * val;
                  grad(i, 2) += dz * val;
                }
              }
            }
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
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      const index_t q2dim = Quadrature::tensor_dim2;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 2-direction
        T u1[order * q0dim * q1dim];
        T u1x[order * q0dim * q1dim];
        T u1y[order * q0dim * q1dim];
        std::fill(u1, u1 + order * q0dim * q1dim, T(0.0));
        std::fill(u1x, u1x + order * q0dim * q1dim, T(0.0));
        std::fill(u1y, u1y + order * q0dim * q1dim, T(0.0));

        for (index_t q2 = 0; q2 < q2dim; q2++) {
          double n2[order], d2[order];
          const double pt2 = Quadrature::get_tensor_point(2, q2);
          lagrange_basis<order>(pt2, n2, d2);

          for (index_t q1 = 0; q1 < q1dim; q1++) {
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              const index_t qindex = Quadrature::get_tensor_index(q0, q1, q2);

              const FiniteElementSpace& s = in.get(qindex);
              const H1Space<T, C, dim>& h1 = s.template get<space>();
              const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
              const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

              T val, derx, dery, derz;
              if constexpr (C == 1) {
                val = u;
                derx = grad(0);
                dery = grad(1);
                derz = grad(2);
              } else {
                val = u(i);
                derx = grad(i, 0);
                dery = grad(i, 1);
                derz = grad(i, 2);
              }
              for (index_t j2 = 0; j2 < order; j2++) {
                u1[j2 + order * (q0 + q0dim * q1)] +=
                    n2[j2] * val + d2[j2] * derz;
                u1x[j2 + order * (q0 + q0dim * q1)] += n2[j2] * derx;
                u1y[j2 + order * (q0 + q0dim * q1)] += n2[j2] * dery;
              }
            }
          }
        }

        T u0[order * order * q0dim];
        T u0x[order * order * q0dim];
        std::fill(u0, u0 + order * order * q0dim, T(0.0));
        std::fill(u0x, u0x + order * order * q0dim, T(0.0));
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order], d1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1, d1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            for (index_t j2 = 0; j2 < order; j2++) {
              T val = u1[j2 + order * (q0 + q0dim * q1)];
              T derx = u1x[j2 + order * (q0 + q0dim * q1)];
              T dery = u1y[j2 + order * (q0 + q0dim * q1)];

              for (index_t j1 = 0; j1 < order; j1++) {
                u0[j1 + order * (j2 + order * q0)] += n1[j1] * val;
                u0x[j1 + order * (j2 + order * q0)] += n1[j1] * derx;
                u0[j1 + order * (j2 + order * q0)] += d1[j1] * dery;
              }
            }
          }
        }

        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          lagrange_basis<order>(pt0, n0, d0);

          for (index_t j2 = 0; j2 < order; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val = u0[j1 + order * (j2 + order * q0)];
              T derx = u0x[j1 + order * (j2 + order * q0)];

              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                res[offset + C * node + i] += n0[j0] * val + d0[j0] * derx;
              }
            }
          }
        }
      }
    } else {
      for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
        // Get the quadrature point
        double pt[dim];
        Quadrature::get_point(q, pt);

        // Evaluate the basis functions
        double n0[order], d0[order];
        double n1[order], d1[order];
        double n2[order], d2[order];
        lagrange_basis<order>(pt[0], n0, d0);
        lagrange_basis<order>(pt[1], n1, d1);
        lagrange_basis<order>(pt[2], n2, d2);

        const FiniteElementSpace& s = in.get(q);
        const H1Space<T, C, dim>& h1 = s.template get<space>();
        const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
        const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order; j1++) {
            double n1n2 = n1[j1] * n2[j2];
            double d1n2 = d1[j1] * n2[j2];
            double n1d2 = n1[j1] * d2[j2];
            for (index_t j0 = 0; j0 < order; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double N = n0[j0] * n1n2;
              double dx = d0[j0] * n1n2;
              double dy = n0[j0] * d1n2;
              double dz = n0[j0] * n1d2;

              if constexpr (C == 1) {
                res[offset + node] +=
                    N * u + dx * grad(0) + dy * grad(1) + dz * grad(2);
              } else {
                for (index_t i = 0; i < C; i++) {
                  res[offset + C * node + i] += N * u(i) + dx * grad(i, 0) +
                                                dy * grad(i, 1) +
                                                dz * grad(i, 2);
                }
              }
            }
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
  static const index_t ncomp_per_stride = ncomp / stride;  // = 4

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
    double n0[order], n1[order], n2[order];
    double d0[order], d1[order], d2[order];
    lagrange_basis<order>(pt[0], n0, d0);
    lagrange_basis<order>(pt[1], n1, d1);
    lagrange_basis<order>(pt[2], n2, d2);

    for (index_t j2 = 0; j2 < order; j2++) {
      for (index_t j1 = 0; j1 < order; j1++) {
        double n1n2 = n1[j1] * n2[j2];
        double d1n2 = d1[j1] * n2[j2];
        double n1d2 = n1[j1] * d2[j2];

        for (index_t j0 = 0; j0 < order; j0++) {
          const index_t node = j0 + order * (j1 + order * j2);
          N[(dim + 1) * node] = n0[j0] * n1n2;
          N[(dim + 1) * node + 1] = d0[j0] * n1n2;
          N[(dim + 1) * node + 2] = n0[j0] * d1n2;
          N[(dim + 1) * node + 3] = n0[j0] * n1d2;
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

  static constexpr BasisType get_basis_type() { return L2; }

  /**
   * @brief Degree of freedom handling on the vertices, edges, faces and
   * volume
   *
   * @param entity The type of topological entity (vertex, edge, face or
   * volume)
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
   * @brief Get the parametric point location associated with the given
   * degree of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  static void get_dof_point(index_t index, double pt[]) {
    constexpr const double* pts = get_gauss_quadrature_pts<order>();

    index_t n = index / C;
    pt[0] = pts[n % order];
    pt[1] = pts[(n % (order * order)) / order];
    pt[2] = pts[n / (order * order)];
  }

  template <index_t space, class Quadrature, class FiniteElementSpace,
            index_t offset, class SolnType>
  static void interp(const SolnType& sol,
                     QptSpace<Quadrature, FiniteElementSpace>& out) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      const index_t q2dim = Quadrature::tensor_dim2;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 0-direction
        T u0[order * order * q0dim];
        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          lagrange_basis<order>(pt0, n0);

          for (index_t j2 = 0; j2 < order; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val(0.0);
              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                val += n0[j0] * sol[offset + C * node + i];
              }

              u0[j1 + order * (j2 + order * q0)] = val;
            }
          }
        }

        // Interpolate along the 1-direction
        T u1[order * q0dim * q1dim];
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            for (index_t j2 = 0; j2 < order; j2++) {
              T val(0.0), derx(0.0), dery(0.0);
              for (index_t j1 = 0; j1 < order; j1++) {
                val += n1[j1] * u0[j1 + order * (j2 + order * q0)];
              }

              u1[j2 + order * (q0 + q0dim * q1)] = val;
            }
          }
        }

        // Interpolate along the 2-direction
        for (index_t q2 = 0; q2 < q2dim; q2++) {
          double n2[order];
          const double pt2 = Quadrature::get_tensor_point(2, q2);
          lagrange_basis<order>(pt2, n2);

          for (index_t q1 = 0; q1 < q1dim; q1++) {
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              T val(0.0);
              for (index_t j2 = 0; j2 < order; j2++) {
                val += n2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
              }

              const index_t qindex = Quadrature::get_tensor_index(q0, q1, q2);
              FiniteElementSpace& s = out.get(qindex);
              L2Space<T, C, dim>& l2 = s.template get<space>();
              typename L2Space<T, C, dim>::VarType& u = l2.get_value();

              if constexpr (C == 1) {
                u = val;
              } else {
                u(i) = val;
              }
            }
          }
        }
      }
    } else {
      for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
        // Get the quadrature point
        double pt[dim];
        Quadrature::get_point(q, pt);

        FiniteElementSpace& s = out.get(q);
        L2Space<T, C, dim>& l2 = s.template get<space>();
        typename L2Space<T, C, dim>::VarType& u = l2.get_value();

        if constexpr (C == 1) {
          u = 0.0;
        } else {
          u.zero();
        }

        // Get the quadrature knot locations
        constexpr const double* knots = get_gauss_quadrature_pts<order>();

        // Evaluate the basis functions
        double n0[order], n1[order], n2[order];
        lagrange_basis<order>(knots, pt[0], n0);
        lagrange_basis<order>(knots, pt[1], n1);
        lagrange_basis<order>(knots, pt[2], n2);

        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order; j1++) {
            double n1n2 = n1[j1] * n2[j2];
            for (index_t j0 = 0; j0 < order; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double n0n1n2 = n0[j0] * n1n2;

              if constexpr (C == 1) {
                u += n0n1n2 * sol[offset + node];
              } else {
                for (index_t i = 0; i < C; i++) {
                  u(i) += n0n1n2 * sol[offset + C * node + i];
                }
              }
            }
          }
        }
      }
    }
  }

  template <index_t space, class Quadrature, class FiniteElementSpace,
            index_t offset, class SolnType>
  static void add(const QptSpace<Quadrature, FiniteElementSpace>& in,
                  SolnType& res) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      const index_t q2dim = Quadrature::tensor_dim2;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 2-direction
        T u1[order * q0dim * q1dim];
        std::fill(u1, u1 + order * q0dim * q1dim, T(0.0));

        for (index_t q2 = 0; q2 < q2dim; q2++) {
          double n2[order];
          const double pt2 = Quadrature::get_tensor_point(2, q2);
          lagrange_basis<order>(pt2, n2);

          for (index_t q1 = 0; q1 < q1dim; q1++) {
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              const index_t qindex = Quadrature::get_tensor_index(q0, q1, q2);

              const FiniteElementSpace& s = in.get(qindex);
              const L2Space<T, C, dim>& l2 = s.template get<space>();
              const typename L2Space<T, C, dim>::VarType& u = l2.get_value();

              T val;
              if constexpr (C == 1) {
                val = u;
              } else {
                val = u(i);
              }
              for (index_t j2 = 0; j2 < order; j2++) {
                u1[j2 + order * (q0 + q0dim * q1)] += n2[j2] * val;
              }
            }
          }
        }

        T u0[order * order * q0dim];
        std::fill(u0, u0 + order * order * q0dim, T(0.0));
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            for (index_t j2 = 0; j2 < order; j2++) {
              T val = u1[j2 + order * (q0 + q0dim * q1)];

              for (index_t j1 = 0; j1 < order; j1++) {
                u0[j1 + order * (j2 + order * q0)] += n1[j1] * val;
              }
            }
          }
        }

        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          lagrange_basis<order>(pt0, n0, d0);

          for (index_t j2 = 0; j2 < order; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val = u0[j1 + order * (j2 + order * q0)];

              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                res[offset + C * node + i] += n0[j0] * val;
              }
            }
          }
        }
      }
    } else {
      for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
        // Get the quadrature point
        double pt[dim];
        Quadrature::get_point(q, pt);

        const FiniteElementSpace& s = in.get(q);
        const L2Space<T, C, dim>& l2 = s.template get<space>();
        const typename L2Space<T, C, dim>::VarType& u = l2.get_value();

        // Get the quadrature knot locations
        constexpr const double* knots = get_gauss_quadrature_pts<order>();

        // Evaluate the basis functions
        double n0[order], n1[order], n2[order];
        lagrange_basis<order>(knots, pt[0], n0);
        lagrange_basis<order>(knots, pt[1], n1);
        lagrange_basis<order>(knots, pt[2], n2);

        for (index_t j2 = 0; j2 < order; j2++) {
          for (index_t j1 = 0; j1 < order; j1++) {
            double n1n2 = n1[j1] * n2[j2];
            for (index_t j0 = 0; j0 < order; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double n0n1n2 = n0[j0] * n1n2;

              if constexpr (C == 1) {
                res[offset + node] += n0n1n2 * u;
              } else {
                for (index_t i = 0; i < C; i++) {
                  res[offset + C * node + i] += n0n1n2 * u(i);
                }
              }
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

    // Get the quadrature knot locations
    constexpr const double* knots = get_gauss_quadrature_pts<order>();

    // Evaluate the basis functions
    double n0[order], n1[order], n2[order];
    lagrange_basis<order>(knots, pt[0], n0);
    lagrange_basis<order>(knots, pt[1], n1);
    lagrange_basis<order>(knots, pt[2], n2);

    for (index_t j2 = 0; j2 < order; j2++) {
      for (index_t j1 = 0; j1 < order; j1++) {
        double n1n2 = n1[j1] * n2[j2];
        for (index_t j0 = 0; j0 < order; j0++) {
          const index_t node = j0 + order * (j1 + order * j2);
          N[node] = n1[j0] * n1n2;
        }
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_LAGRANGE_HEX_BASIS_H