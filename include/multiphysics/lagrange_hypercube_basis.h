#ifndef A2D_LAGRANGE_HEX_BASIS_H
#define A2D_LAGRANGE_HEX_BASIS_H

#include "a2denum.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

// Compute x^y at compile time, where x and y are integer types
// Usage: index_t power = indexpow<x, y>::value
template <index_t x, index_t y>
struct indexpow;

template <index_t x>
struct indexpow<x, 0> {
  static const int value = 1;
};

template <index_t x, index_t y>
struct indexpow {
  static const int value = x * indexpow<x, y - 1>::value;
};

/**
 * Line (dim == 1), Quad (dim == 2) or Hex (dim == 3) high order Lagrange basis
 */
template <typename T, index_t Dim, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
class LagrangeH1HypercubeBasis {
 public:
  using ET = ElementTypes;

  static constexpr index_t dim = Dim;
  static_assert(dim == 1 or dim == 2 or dim == 3, "unsupported dim");

  static const index_t order = degree + 1;  // Number of nodes along each edge

  static const index_t ndof =
      C * indexpow<order, dim>::value;  // Total number of degrees of freedom

  // Number of components
  static const index_t ncomp = H1Space<T, C, dim>::ncomp;

  // Get the type of basis class implemented
  static constexpr BasisType get_basis_type() { return H1; }

  // Define the equivalent low-order basis class if any
  using LOrderBasis = LagrangeH1HypercubeBasis<T, dim, C, 1, interp_type>;

  /**
   * @brief Degree of freedom handling on the topological entities
   *
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the topological entity (e.g. edge index)
   * @return The number of degrees of freedom on a given entity
   */
  A2D_INLINE_FUNCTION static index_t get_entity_ndof(ET::ElementEntity entity,
                                                     index_t index) {
    switch (entity) {
      case ET::Vertex:
        return C;

      // Only 3D element has edge
      case ET::Edge:
        if constexpr (dim == 3) {
          return C * (order - 2);
        }

      // Bound for n-dimensional element has n-1 dimensions
      case ET::Bound:
        return C * indexpow<order - 2, dim - 1>::value;

      // Domain for n-dimensional element has n dimensions
      case ET::Domain:
        return C * indexpow<order - 2, dim>::value;
    }
    return 0;
  }

  /**
   * @brief Get the entity DOF from the element DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */
  template <index_t offset, class ElemDof, class EntityDof>
  A2D_INLINE_FUNCTION static void get_entity_dof(ET::ElementEntity entity,
                                                 index_t index,
                                                 const ElemDof& element_dof,
                                                 EntityDof& entity_dof) {
    switch (entity) {
      case ET::Vertex:
        if constexpr (dim == 2) {
          ET::get_quad_vert_dof<offset, C, order, order, ElemDof, EntityDof>(
              index, element_dof, entity_dof);
        } else if constexpr (dim == 3) {
          ET::get_hex_vert_dof<offset, C, order, order, order, ElemDof,
                               EntityDof>(index, element_dof, entity_dof);
        }
        break;

      case ET::Edge:
        if constexpr (dim == 3) {
          const bool endp = false;
          ET::get_hex_edge_dof<offset, endp, C, order, order, order, ElemDof,
                               EntityDof>(index, element_dof, entity_dof);
        }
        break;

      case ET::Bound:
        if constexpr (dim == 1) {
          ET::get_line_bound_dof<offset, C, order, ElemDof, EntityDof>(
              index, element_dof, entity_dof);
        } else if constexpr (dim == 2) {
          const bool endp = false;
          ET::get_quad_bound_dof<offset, endp, C, order, order, ElemDof,
                                 EntityDof>(index, element_dof, entity_dof);
        } else if constexpr (dim == 3) {
          const bool endp = false;
          ET::get_hex_bound_dof<offset, endp, C, order, order, order, ElemDof,
                                EntityDof>(index, element_dof, entity_dof);
        }
        break;

      case ET::Domain:
        if constexpr (dim == 1) {
          const bool endp = false;
          ET::get_line_domain_dof<offset, endp, C, order, ElemDof, EntityDof>(
              element_dof, entity_dof);
        } else if constexpr (dim == 2) {
          const bool endp = false;
          ET::get_quad_domain_dof<offset, endp, C, order, order, ElemDof,
                                  EntityDof>(element_dof, entity_dof);
        } else if constexpr (dim == 3) {
          const bool endp = false;
          ET::get_hex_domain_dof<offset, endp, C, order, order, order, ElemDof,
                                 EntityDof>(element_dof, entity_dof);
        }
        break;
    }
  }

  /**
   * @brief Set the element DOF and signs from the entity DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param entity_dof Entity DOF in the global orientation
   * @param element_dof Degrees of freedom for this element
   */
  template <index_t offset, class EntityDof, class ElemDof>
  A2D_INLINE_FUNCTION static void set_entity_dof(ET::ElementEntity entity,
                                                 index_t index, index_t orient,
                                                 const EntityDof& entity_dof,
                                                 ElemDof& element_dof) {
    switch (entity) {
      case ET::Vertex:
        if constexpr (dim == 2) {
          ET::set_quad_vert_dof<offset, C, order, order, EntityDof, ElemDof>(
              index, entity_dof, element_dof);
        } else if constexpr (dim == 3) {
          ET::set_hex_vert_dof<offset, C, order, order, order, EntityDof,
                               ElemDof>(index, entity_dof, element_dof);
        }
        break;

      case ET::Edge:
        if constexpr (dim == 3) {
          const bool endp = false;
          ET::set_hex_edge_dof<offset, endp, C, order, order, order, EntityDof,
                               ElemDof>(index, orient, entity_dof, element_dof);
        }
        break;

      case ET::Bound:
        if constexpr (dim == 1) {
          ET::set_line_bound_dof<offset, C, order, EntityDof, ElemDof>(
              index, orient, entity_dof, element_dof);
        } else if constexpr (dim == 2) {
          const bool endp = false;
          ET::set_quad_bound_dof<offset, endp, C, order, order, EntityDof,
                                 ElemDof>(index, orient, entity_dof,
                                          element_dof);
        } else if constexpr (dim == 3) {
          const bool endp = false;
          ET::set_hex_bound_dof<offset, endp, C, order, order, order, EntityDof,
                                ElemDof>(index, orient, entity_dof,
                                         element_dof);
        }
        break;

      case ET::Domain:
        if constexpr (dim == 1) {
          const bool endp = false;
          ET::set_line_domain_dof<offset, endp, C, order, EntityDof, ElemDof>(
              orient, entity_dof, element_dof);
        } else if constexpr (dim == 2) {
          const bool endp = false;
          ET::set_quad_domain_dof<offset, endp, C, order, order, EntityDof,
                                  ElemDof>(orient, entity_dof, element_dof);
        } else if constexpr (dim == 3) {
          const bool endp = false;
          ET::set_hex_domain_dof<offset, endp, C, order, order, order,
                                 EntityDof, ElemDof>(entity_dof, element_dof);
        }
        break;
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
  A2D_INLINE_FUNCTION static void set_entity_signs(ET::ElementEntity entity,
                                                   index_t index,
                                                   index_t orient,
                                                   int signs[]) {
    int sgns[ndof];
    const index_t entity_ndof = get_entity_ndof(entity, index);
    for (index_t i = 0; i < entity_ndof; i++) {
      sgns[i] = 1;
    }
    if (entity_ndof != 0) {
      set_entity_dof<offset>(entity, index, orient, sgns, signs);
    }
  }

  /**
   * @brief Get the number of low order elements defined by this basis
   */
  A2D_INLINE_FUNCTION constexpr static index_t get_num_lorder_elements() {
    return indexpow<degree, dim>::value;  // degree ** dim
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
  A2D_INLINE_FUNCTION static void get_lorder_dof(const index_t n,
                                                 const HOrderDof& hdof,
                                                 LOrderDof& ldof) {
    if constexpr (dim == 1) {
      // lorder element index:
      // n = i
      for (index_t ii = 0; ii < 2; ii++) {
        const index_t lnode = ii;
        const index_t hnode = n + ii;

        for (index_t p = 0; p < C; p++) {
          ldof[lorder_offset + C * lnode + p] =
              hdof[horder_offset + C * hnode + p];
        }
      }
    } else if constexpr (dim == 2) {
      // lorder element index:
      // n = i + j * degree, 0 <= i, j < degree
      const index_t i = n % degree;
      const index_t j = n / degree;

      for (index_t jj = 0; jj < 2; jj++) {
        for (index_t ii = 0; ii < 2; ii++) {
          const index_t lnode = ii + 2 * jj;
          const index_t hnode = (i + ii) + order * (j + jj);

          for (index_t p = 0; p < C; p++) {
            ldof[lorder_offset + C * lnode + p] =
                hdof[horder_offset + C * hnode + p];
          }
        }
      }
    } else {  // dim == 3
      // lorder element index:
      // n = i + j * degree + k * degree * degree, 0 <= i, j, k < degree
      const index_t i = n % degree;
      const index_t j = (n % (degree * degree)) / degree;
      const index_t k = n / (degree * degree);

      for (index_t kk = 0; kk < 2; kk++) {
        for (index_t jj = 0; jj < 2; jj++) {
          for (index_t ii = 0; ii < 2; ii++) {
            const index_t lnode = ii + 2 * (jj + 2 * kk);
            const index_t hnode =
                (i + ii) + order * ((j + jj) + order * (k + kk));

            for (index_t p = 0; p < C; p++) {
              ldof[lorder_offset + C * lnode + p] =
                  hdof[horder_offset + C * hnode + p];
            }
          }
        }
      }
    }
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
  A2D_INLINE_FUNCTION static void get_lorder_signs(const index_t n,
                                                   const int horder_signs[],
                                                   int signs[]) {
    for (index_t p = 0; p < indexpow<2, dim>::value * C; p++) {
      signs[lorder_offset + p] = 1;
    }
  }

  /**
   * @brief Get the parametric point location associated with the given
   * degree of freedom
   *
   * TODO: This code assumes we're using GLL interpolation. That is not
   * necessarily the case. Generalize this?
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  A2D_INLINE_FUNCTION static void get_dof_point(index_t index, double pt[]) {
    // Get the quadrature knot locations
    constexpr const double* pts = get_interpolation_pts<order, interp_type>();
    index_t n = index / C;

    if constexpr (dim == 1) {
      pt[0] = pts[n];
    } else if constexpr (dim == 2) {
      pt[0] = pts[n % order];
      pt[1] = pts[n / order];
    } else {  // dim == 3
      pt[0] = pts[n % order];
      pt[1] = pts[(n % (order * order)) / order];
      pt[2] = pts[n / (order * order)];
    }
  }

  /**
   * @brief Interpolate over all quadrature points using
   *
   * This function just calls the interpolation  the basis functions are
   * interpolated
   *
   * @tparam space The finite element space index
   * @tparam QptGeoSpace Geometry object (not used for this basis)
   * @tparam Quadrature The quadrature scheme
   * @tparam FiniteElementSpace The finite element space object
   * @tparam offset Offset index into the solution degrees of freedom
   * @tparam SolnType The solution vector type for the dof
   * @param geo The geometry object
   * @param sol The solution array
   * @param out The finite element space object at all quadrature points
   */
  template <index_t space, class QptGeoSpace, class Quadrature,
            class FiniteElementSpace, index_t offset, class SolnType>
  A2D_INLINE_FUNCTION static void interp(
      const QptGeoSpace& geo, const SolnType& sol,
      QptSpace<Quadrature, FiniteElementSpace>& out) {
    interp<space, Quadrature, FiniteElementSpace, offset, SolnType>(sol, out);
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
  A2D_INLINE_FUNCTION static void interp(
      const SolnType& sol, QptSpace<Quadrature, FiniteElementSpace>& out) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      index_t q2dim = 0;
      if constexpr (dim == 3) {
        q2dim = Quadrature::tensor_dim2;
      }

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 0-direction
        T u0[indexpow<order, dim - 1>::value * q0dim];
        T u0x[indexpow<order, dim - 1>::value * q0dim];
        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          interpolation_basis<order, interp_type>(pt0, n0, d0);

          const index_t outer = dim == 2 ? 1 : order;
          for (index_t j2 = 0; j2 < outer; j2++) {  // j2 == 0 for dim == 2
            for (index_t j1 = 0; j1 < order; j1++) {
              T val(0.0), derx(0.0);
              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                val += n0[j0] * sol[offset + C * node + i];
                derx += d0[j0] * sol[offset + C * node + i];
              }

              u0[j1 + order * (j2 + outer * q0)] = val;
              u0x[j1 + order * (j2 + outer * q0)] = derx;
            }
          }
        }

        // Interpolate along the 1-direction, for dim == 3 only
        T u1[order * q0dim * q1dim];
        T u1x[order * q0dim * q1dim];
        T u1y[order * q0dim * q1dim];
        if constexpr (dim == 3) {
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order], d1[order];
            const double pt1 = Quadrature::get_tensor_point(1, q1);
            interpolation_basis<order, interp_type>(pt1, n1, d1);

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
        }

        // For dim == 3: interpolate along the 2-direction, or
        // for dim == 2: interpolate along the 1-direction
        const index_t qouter = dim == 2 ? 1 : q2dim;
        double n2[order], d2[order];
        for (index_t q2 = 0; q2 < qouter; q2++) {  // q2 == 0 for dim == 2
          if constexpr (dim == 3) {
            const double pt2 = Quadrature::get_tensor_point(2, q2);
            interpolation_basis<order, interp_type>(pt2, n2, d2);
          }
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order], d1[order];
            if constexpr (dim == 2) {
              const double pt1 = Quadrature::get_tensor_point(1, q1);
              interpolation_basis<order, interp_type>(pt1, n1, d1);
            }
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              T val(0.0), derx(0.0), dery(0.0), derz(0.0);
              for (index_t j2 = 0; j2 < order; j2++) {
                if constexpr (dim == 3) {
                  val += n2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
                  derx += n2[j2] * u1x[j2 + order * (q0 + q0dim * q1)];
                  dery += n2[j2] * u1y[j2 + order * (q0 + q0dim * q1)];
                  derz += d2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
                } else {  // dim == 2
                  val += n1[j2] * u0[j2 + order * q0];
                  derx += n1[j2] * u0x[j2 + order * q0];
                  dery += d1[j2] * u0[j2 + order * q0];
                }
              }

              index_t qindex;
              if constexpr (dim == 2) {
                qindex = Quadrature::get_tensor_index(q0, q1);
              } else {  // dim == 3
                qindex = Quadrature::get_tensor_index(q0, q1, q2);
              }

              FiniteElementSpace& s = out.get(qindex);
              H1Space<T, C, dim>& h1 = s.template get<space>();
              typename H1Space<T, C, dim>::VarType& u = h1.get_value();
              typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

              if constexpr (C == 1) {
                u = val;
                grad(0) = derx;
                grad(1) = dery;
                if constexpr (dim == 3) {
                  grad(2) = derz;
                }

              } else {
                u(i) = val;
                grad(i, 0) = derx;
                grad(i, 1) = dery;
                if constexpr (dim == 3) {
                  grad(i, 2) = derz;
                }
              }
            }
          }
        }
      }
    } else {  // is not tensor_product
      for (index_t q = 0; q < Quadrature::get_num_points(); q++) {
        // Get the quadrature point
        double pt[dim];
        Quadrature::get_point(q, pt);

        // Evaluate the basis functions
        double n0[order], d0[order];
        double n1[order], d1[order];
        double n2[order], d2[order];

        if constexpr (dim >= 1) {
          interpolation_basis<order, interp_type>(pt[0], n0, d0);
        }
        if constexpr (dim >= 2) {
          interpolation_basis<order, interp_type>(pt[1], n1, d1);
        }
        if constexpr (dim >= 3) {
          interpolation_basis<order, interp_type>(pt[2], n2, d2);
        }

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

        // | dim | j2range | j1range | j0range |
        // -------------------------------------
        // | 3   | order   | order   | order   |
        // | 2   | 1       | order   | order   |
        // | 1   | 1       | 1       | order   |
        constexpr index_t j2range = dim >= 3 ? order : 1;
        constexpr index_t j1range = dim >= 2 ? order : 1;
        constexpr index_t j0range = dim >= 1 ? order : 1;

        for (index_t j2 = 0; j2 < j2range; j2++) {
          for (index_t j1 = 0; j1 < j1range; j1++) {
            double n1n2, d1n2, n1d2;
            if constexpr (dim == 3) {
              n1n2 = n1[j1] * n2[j2];
              d1n2 = d1[j1] * n2[j2];
              n1d2 = n1[j1] * d2[j2];
            }
            for (index_t j0 = 0; j0 < j0range; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double N, dx, dy, dz;
              if constexpr (dim == 3) {
                N = n0[j0] * n1n2;
                dx = d0[j0] * n1n2;
                dy = n0[j0] * d1n2;
                dz = n0[j0] * n1d2;
              } else if constexpr (dim == 2) {
                N = n0[j0] * n1[j1];
                dx = d0[j0] * n1[j1];
                dy = n0[j0] * d1[j1];
              } else {  // dim == 1
                N = n0[j0];
                dx = d0[j0];
              }

              if constexpr (C == 1) {
                const T val = sol[offset + node];

                u += N * val;
                if constexpr (dim >= 1) {
                  grad(0) += dx * val;
                }
                if constexpr (dim >= 2) {
                  grad(1) += dy * val;
                }
                if constexpr (dim >= 3) {
                  grad(2) += dz * val;
                }
              } else {
                for (index_t i = 0; i < C; i++) {
                  const T val = sol[offset + C * node + i];

                  u(i) += N * val;
                  if constexpr (dim >= 1) {
                    grad(i, 0) += dx * val;
                  }
                  if constexpr (dim >= 2) {
                    grad(i, 1) += dy * val;
                  }
                  if constexpr (dim >= 3) {
                    grad(i, 2) += dz * val;
                  }
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
   * @tparam QptGeoSpace Geometry object (not used for this basis)
   * @tparam Quadrature The quadrature object
   * @tparam FiniteElementSpace The finite element space object
   * @tparam offset Degree of freedom offset into the array
   * @tparam SolnType Solution array type
   * @param in The finite element space output object
   * @param res The residual array - same shape as the solution array
   */
  template <index_t space, class QptGeoSpace, class Quadrature,
            class FiniteElementSpace, index_t offset, class SolnType>
  A2D_INLINE_FUNCTION static void add(
      const QptSpace<Quadrature, FiniteElementSpace>& in, SolnType& res) {
    add<space, Quadrature, FiniteElementSpace, offset, SolnType>(in, res);
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
  A2D_INLINE_FUNCTION static void add(
      const QptSpace<Quadrature, FiniteElementSpace>& in, SolnType& res) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      index_t q2dim = 0;
      if constexpr (dim == 3) {
        q2dim = Quadrature::tensor_dim2;
      }

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 2-direction
        // index_t u0size, u1size;
        // if constexpr (dim == 2) {
        //   u0size = order * q0dim;
        //   u1size = order * q0dim;
        // } else {
        //   u0size = order * order * q0dim;
        //   u1size = order * q0dim * q1dim;
        // }

        constexpr index_t u0size = order * order * q0dim;
        constexpr index_t u1size = order * q0dim * q1dim;

        T u1[u1size];
        T u1x[u1size];
        T u1y[u1size];  // not used for dim == 2
        std::fill(u1, u1 + u1size, T(0.0));
        std::fill(u1x, u1x + u1size, T(0.0));
        std::fill(u1y, u1y + u1size, T(0.0));

        T u0[u0size];
        T u0x[u0size];
        std::fill(u0, u0 + u0size, T(0.0));
        std::fill(u0x, u0x + u0size, T(0.0));

        const index_t qouter = dim == 2 ? 1 : q2dim;
        for (index_t q2 = 0; q2 < qouter; q2++) {
          double n2[order], d2[order];
          if constexpr (dim == 3) {
            const double pt2 = Quadrature::get_tensor_point(2, q2);
            interpolation_basis<order, interp_type>(pt2, n2, d2);
          }
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order], d1[order];
            if constexpr (dim == 2) {
              const double pt1 = Quadrature::get_tensor_point(1, q1);
              interpolation_basis<order, interp_type>(pt1, n1, d1);
            }
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              index_t qindex;
              if constexpr (dim == 2) {
                qindex = Quadrature::get_tensor_index(q0, q1);
              } else {
                qindex = Quadrature::get_tensor_index(q0, q1, q2);
              }

              const FiniteElementSpace& s = in.get(qindex);
              const H1Space<T, C, dim>& h1 = s.template get<space>();
              const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
              const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

              T val, derx, dery, derz;
              if constexpr (C == 1) {
                val = u;
                derx = grad(0);
                dery = grad(1);
                if constexpr (dim == 3) {
                  derz = grad(2);
                }
              } else {
                val = u(i);
                derx = grad(i, 0);
                dery = grad(i, 1);
                if constexpr (dim == 3) {
                  derz = grad(i, 2);
                }
              }
              for (index_t j2 = 0; j2 < order; j2++) {
                if constexpr (dim == 2) {
                  u0[j2 + order * q0] += n1[j2] * val + d1[j2] * dery;
                  u0x[j2 + order * q0] += n1[j2] * derx;
                } else {
                  u1[j2 + order * (q0 + q0dim * q1)] +=
                      n2[j2] * val + d2[j2] * derz;
                  u1x[j2 + order * (q0 + q0dim * q1)] += n2[j2] * derx;
                  u1y[j2 + order * (q0 + q0dim * q1)] += n2[j2] * dery;
                }
              }
            }
          }
        }

        if constexpr (dim == 3) {
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order], d1[order];
            const double pt1 = Quadrature::get_tensor_point(1, q1);
            interpolation_basis<order, interp_type>(pt1, n1, d1);

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
        }

        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          interpolation_basis<order, interp_type>(pt0, n0, d0);

          const index_t outer = dim == 2 ? 1 : order;
          for (index_t j2 = 0; j2 < outer; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val = u0[j1 + order * (j2 + outer * q0)];
              T derx = u0x[j1 + order * (j2 + outer * q0)];

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

        if constexpr (dim >= 1) {
          interpolation_basis<order, interp_type>(pt[0], n0, d0);
        }
        if constexpr (dim >= 2) {
          interpolation_basis<order, interp_type>(pt[1], n1, d1);
        }
        if constexpr (dim >= 3) {
          interpolation_basis<order, interp_type>(pt[2], n2, d2);
        }

        const FiniteElementSpace& s = in.get(q);
        const H1Space<T, C, dim>& h1 = s.template get<space>();
        const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
        const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

        // | dim | j2range | j1range | j0range |
        // -------------------------------------
        // | 3   | order   | order   | order   |
        // | 2   | 1       | order   | order   |
        // | 1   | 1       | 1       | order   |
        constexpr index_t j2range = dim >= 3 ? order : 1;
        constexpr index_t j1range = dim >= 2 ? order : 1;
        constexpr index_t j0range = dim >= 1 ? order : 1;

        for (index_t j2 = 0; j2 < j2range; j2++) {
          for (index_t j1 = 0; j1 < j1range; j1++) {
            double n1n2, d1n2, n1d2;
            if constexpr (dim == 3) {
              n1n2 = n1[j1] * n2[j2];
              d1n2 = d1[j1] * n2[j2];
              n1d2 = n1[j1] * d2[j2];
            }
            for (index_t j0 = 0; j0 < j0range; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double N, dx, dy, dz;
              if constexpr (dim == 3) {
                N = n0[j0] * n1n2;
                dx = d0[j0] * n1n2;
                dy = n0[j0] * d1n2;
                dz = n0[j0] * n1d2;
              } else if constexpr (dim == 2) {
                N = n0[j0] * n1[j1];
                dx = d0[j0] * n1[j1];
                dy = n0[j0] * d1[j1];
              } else {  // dim == 1
                N = n0[j0];
                dx = d0[j0];
              }

              if constexpr (C == 1) {
                res[offset + node] += N * u + dx * grad(0);
                if constexpr (dim == 2) {
                  res[offset + node] += dy * grad(1);
                }
                if constexpr (dim == 3) {
                  res[offset + node] += dy * grad(1) + dz * grad(2);
                }
              } else {
                for (index_t i = 0; i < C; i++) {
                  res[offset + C * node + i] += N * u(i) + dx * grad(i, 0);
                  if constexpr (dim == 2) {
                    res[offset + C * node + i] += dy * grad(i, 1);
                  }
                  if constexpr (dim == 3) {
                    res[offset + C * node + i] +=
                        dy * grad(i, 1) + dz * grad(i, 2);
                  }
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
  static const index_t basis_size = (dim + 1) * indexpow<order, dim>::value;

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
  A2D_INLINE_FUNCTION static void basis(index_t n, BasisType N) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    // Evaluate the basis functions
    double n0[order], n1[order], n2[order];
    double d0[order], d1[order], d2[order];

    if constexpr (dim >= 1) {
      interpolation_basis<order, interp_type>(pt[0], n0, d0);
    }
    if constexpr (dim >= 2) {
      interpolation_basis<order, interp_type>(pt[1], n1, d1);
    }
    if constexpr (dim == 3) {
      interpolation_basis<order, interp_type>(pt[2], n2, d2);
    }

    // | dim | j2range | j1range | j0range |
    // -------------------------------------
    // | 3   | order   | order   | order   |
    // | 2   | 1       | order   | order   |
    // | 1   | 1       | 1       | order   |
    constexpr index_t j2range = dim >= 3 ? order : 1;
    constexpr index_t j1range = dim >= 2 ? order : 1;
    constexpr index_t j0range = dim >= 1 ? order : 1;

    for (index_t j2 = 0; j2 < j2range; j2++) {
      for (index_t j1 = 0; j1 < j1range; j1++) {
        double n1n2, d1n2, n1d2;

        if constexpr (dim == 3) {
          n1n2 = n1[j1] * n2[j2];
          d1n2 = d1[j1] * n2[j2];
          n1d2 = n1[j1] * d2[j2];
        }

        for (index_t j0 = 0; j0 < j0range; j0++) {
          const index_t node = j0 + order * (j1 + order * j2);

          if constexpr (dim == 1) {
            N[(dim + 1) * node] = n0[j0];
            N[(dim + 1) * node + 1] = d0[j0];
          } else if constexpr (dim == 2) {
            N[(dim + 1) * node] = n0[j0] * n1[j1];
            N[(dim + 1) * node + 1] = d0[j0] * n1[j1];
            N[(dim + 1) * node + 2] = n0[j0] * d1[j1];
          } else {  // dim == 3
            N[(dim + 1) * node] = n0[j0] * n1n2;
            N[(dim + 1) * node + 1] = d0[j0] * n1n2;
            N[(dim + 1) * node + 2] = n0[j0] * d1n2;
            N[(dim + 1) * node + 3] = n0[j0] * n1d2;
          }
        }
      }
    }
  }
};

template <typename T, index_t Dim, index_t C, index_t degree,
          InterpolationType interp_type = GAUSS_INTERPOLATION>
class LagrangeL2HypercubeBasis {
 public:
  using ET = ElementTypes;

  // The topological dimension (dimension of the reference element)
  static constexpr index_t dim = Dim;

  static const index_t order = degree + 1;  // Number of nodes along each edge

  static const index_t ndof =
      C * indexpow<order, dim>::value;  // Total number of degrees of freedom

  // Number of components
  static const index_t ncomp = L2Space<T, C, dim>::ncomp;

  // Get the type of basis class implemented
  A2D_INLINE_FUNCTION static constexpr BasisType get_basis_type() { return L2; }

  // Define the equivalent low-order basis class if any
  using LOrderBasis = LagrangeL2HypercubeBasis<T, dim, C, 0, interp_type>;

  /**
   * @brief Degree of freedom handling on the topological entities
   *
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the topological entity (e.g. edge index)
   * @return The number of degrees of freedom
   */
  A2D_INLINE_FUNCTION static index_t get_entity_ndof(ET::ElementEntity entity,
                                                     index_t index) {
    if (entity == ET::Domain) {
      return ndof;
    }
    return 0;
  }

  /**
   * @brief Get the entity DOF from the element DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param element_dof Degrees of freedom for this element
   * @param entity_dof Entity DOF in the global orientation
   */
  template <index_t offset, class ElemDof, class EntityDof>
  A2D_INLINE_FUNCTION static void get_entity_dof(ET::ElementEntity entity,
                                                 index_t index,
                                                 const ElemDof& element_dof,
                                                 EntityDof& entity_dof) {
    if (entity == ET::Domain) {
      for (index_t i = 0; i < ndof; i++) {
        entity_dof[i] = element_dof[offset + i];
      }
    }
  }

  /**
   * @brief Set the element DOF and signs from the entity DOF
   *
   * @tparam offset The offset into the basis
   * @param entity The type of topological entity (domain, bound, edge, vertex)
   * @param index The index of the entity (e.g. edge index)
   * @param orient Orientation flag indicating the relative orientation
   * @param entity_dof Entity DOF in the global orientation
   * @param element_dof Degrees of freedom for this element
   * @param element_sign Sign indices for each degree of freedom
   */
  template <index_t offset, class EntityDof, class ElemDof>
  A2D_INLINE_FUNCTION static void set_entity_dof(ET::ElementEntity entity,
                                                 index_t index, index_t orient,
                                                 const EntityDof& entity_dof,
                                                 ElemDof& element_dof) {
    if (entity == ET::Domain) {
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
  A2D_INLINE_FUNCTION static void set_entity_signs(ET::ElementEntity entity,
                                                   index_t index,
                                                   index_t orient,
                                                   int signs[]) {
    int sgns[ndof];
    const index_t entity_ndof = get_entity_ndof(entity, index);
    for (index_t i = 0; i < entity_ndof; i++) {
      sgns[i] = 1;
    }
    if (entity_ndof) {
      set_entity_dof<offset>(entity, index, orient, sgns, signs);
    }
  }

  /**
   * @brief Get the number of low order elements defined by this basis
   */
  A2D_INLINE_FUNCTION constexpr static index_t get_num_lorder_elements() {
    return indexpow<order, dim>::value;
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
  A2D_INLINE_FUNCTION static void get_lorder_dof(const index_t n,
                                                 const HOrderDof& hdof,
                                                 LOrderDof& ldof) {
    index_t hnode;

    if constexpr (dim == 1) {
      // lorder element index:
      // n = i
      hnode = n;
    } else if constexpr (dim == 2) {
      // lorder element index:
      // n = i + j * order, 0 <= i, j < order
      const index_t i = n % order;
      const index_t j = n / order;
      hnode = i + order * j;
    } else {  // dim == 3
      // lorder element index:
      // n = i + j * order + k * order * order, 0 <= i, j, k < order
      const index_t i = n % order;
      const index_t j = (n % (order * order)) / order;
      const index_t k = n / (order * order);
      hnode = i + order * (j + order * k);
    }

    for (index_t p = 0; p < C; p++) {
      ldof[lorder_offset + p] = hdof[horder_offset + C * hnode + p];
    }
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
  A2D_INLINE_FUNCTION static void get_lorder_signs(const index_t n,
                                                   const int horder_signs[],
                                                   int signs[]) {
    for (index_t p = 0; p < C; p++) {
      signs[lorder_offset + p] = 1;
    }
  }

  /**
   * @brief Get the parametric point location associated with the given
   * degree of freedom
   *
   * @param index The index for the dof
   * @param pt The parametric point location of dimension dim
   */
  A2D_INLINE_FUNCTION static void get_dof_point(index_t index, double pt[]) {
    // Get the quadrature knot locations
    constexpr const double* pts = get_interpolation_pts<order, interp_type>();

    index_t n = index / C;

    if constexpr (dim == 1) {
      pt[0] = pts[n];
    } else if constexpr (dim == 2) {
      pt[0] = pts[n % order];
      pt[1] = pts[n / order];
    } else {
      pt[0] = pts[n % order];
      pt[1] = pts[(n % (order * order)) / order];
      pt[2] = pts[n / (order * order)];
    }
  }

  /**
   * @brief Interpolate over all quadrature points using
   *
   * This function just calls the interpolation  the basis functions are
   * interpolated
   *
   * @tparam space The finite element space index
   * @tparam QptGeoSpace Geometry object (not used for this basis)
   * @tparam Quadrature The quadrature scheme
   * @tparam FiniteElementSpace The finite element space object
   * @tparam offset Offset index into the solution degrees of freedom
   * @tparam SolnType The solution vector type for the dof
   * @param geo The geometry object
   * @param sol The solution array
   * @param out The finite element space object at all quadrature points
   */
  template <index_t space, class QptGeoSpace, class Quadrature,
            class FiniteElementSpace, index_t offset, class SolnType>
  A2D_INLINE_FUNCTION static void interp(
      const QptGeoSpace& geo, const SolnType& sol,
      QptSpace<Quadrature, FiniteElementSpace>& out) {
    interp<space, Quadrature, FiniteElementSpace, offset, SolnType>(sol, out);
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
  A2D_INLINE_FUNCTION static void interp(
      const SolnType& sol, QptSpace<Quadrature, FiniteElementSpace>& out) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      index_t q2dim = 0;
      if constexpr (dim == 3) {
        q2dim = Quadrature::tensor_dim2;
      }

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 0-direction
        T u0[indexpow<order, dim - 1>::value * q0dim];
        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          interpolation_basis<order, interp_type>(pt0, n0);

          const index_t outer = dim == 2 ? 1 : order;
          for (index_t j2 = 0; j2 < outer; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val(0.0);
              for (index_t j0 = 0; j0 < order; j0++) {
                const index_t node = j0 + order * (j1 + order * j2);
                val += n0[j0] * sol[offset + C * node + i];
              }

              u0[j1 + order * (j2 + outer * q0)] = val;
            }
          }
        }

        // Interpolate along the 1-direction, for dim == 3 only
        T u1[order * q0dim * q1dim];
        if constexpr (dim == 3) {
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order];
            const double pt1 = Quadrature::get_tensor_point(1, q1);
            interpolation_basis<order, interp_type>(pt1, n1);

            for (index_t q0 = 0; q0 < q0dim; q0++) {
              for (index_t j2 = 0; j2 < order; j2++) {
                T val(0.0);
                for (index_t j1 = 0; j1 < order; j1++) {
                  val += n1[j1] * u0[j1 + order * (j2 + order * q0)];
                }

                u1[j2 + order * (q0 + q0dim * q1)] = val;
              }
            }
          }
        }

        // For dim == 3: interpolate along the 2-direction, or
        // for dim == 2: interpolate along the 1-direction
        const index_t qouter = dim == 2 ? 1 : q2dim;
        double n2[order];
        for (index_t q2 = 0; q2 < qouter; q2++) {  // q2 == 0 for dim == 2
          if constexpr (dim == 3) {
            const double pt2 = Quadrature::get_tensor_point(2, q2);
            interpolation_basis<order, interp_type>(pt2, n2);
          }
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order];
            if constexpr (dim == 2) {
              const double pt1 = Quadrature::get_tensor_point(1, q1);
              interpolation_basis<order, interp_type>(pt1, n1);
            }
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              T val(0.0);
              for (index_t j2 = 0; j2 < order; j2++) {
                if constexpr (dim == 3) {
                  val += n2[j2] * u1[j2 + order * (q0 + q0dim * q1)];
                } else {
                  val += n1[j2] * u0[j2 + order * q0];
                }
              }

              index_t qindex;
              if constexpr (dim == 2) {
                qindex = Quadrature::get_tensor_index(q0, q1);
              } else {
                qindex = Quadrature::get_tensor_index(q0, q1, q2);
              }

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

        // Evaluate the basis functions
        double n0[order], n1[order], n2[order];

        if constexpr (dim >= 1) {
          interpolation_basis<order, interp_type>(pt[0], n0);
        }
        if constexpr (dim >= 2) {
          interpolation_basis<order, interp_type>(pt[1], n1);
        }
        if constexpr (dim == 3) {
          interpolation_basis<order, interp_type>(pt[2], n2);
        }

        // | dim | j2range | j1range | j0range |
        // -------------------------------------
        // | 3   | order   | order   | order   |
        // | 2   | 1       | order   | order   |
        // | 1   | 1       | 1       | order   |
        constexpr index_t j2range = dim >= 3 ? order : 1;
        constexpr index_t j1range = dim >= 2 ? order : 1;
        constexpr index_t j0range = dim >= 1 ? order : 1;

        for (index_t j2 = 0; j2 < j2range; j2++) {
          for (index_t j1 = 0; j1 < j1range; j1++) {
            double n1n2;
            if constexpr (dim == 3) {
              n1n2 = n1[j1] * n2[j2];
            }
            for (index_t j0 = 0; j0 < j0range; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double coeff;
              if constexpr (dim == 3) {
                coeff = n0[j0] * n1n2;
              } else if constexpr (dim == 2) {
                coeff = n0[j0] * n1[j1];
              } else {  // dim == 1
                coeff = n0[j0];
              }
              if constexpr (C == 1) {
                u += coeff * sol[offset + node];
              } else {
                for (index_t i = 0; i < C; i++) {
                  u(i) += coeff * sol[offset + C * node + i];
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
   * @tparam QptGeoSpace Geometry object (not used for this basis)
   * @tparam Quadrature The quadrature object
   * @tparam FiniteElementSpace The finite element space object
   * @tparam offset Degree of freedom offset into the array
   * @tparam SolnType Solution array type
   * @param in The finite element space output object
   * @param res The residual array - same shape as the solution array
   */
  template <index_t space, class QptGeoSpace, class Quadrature,
            class FiniteElementSpace, index_t offset, class SolnType>
  A2D_INLINE_FUNCTION static void add(
      const QptSpace<Quadrature, FiniteElementSpace>& in, SolnType& res) {
    add<space, Quadrature, FiniteElementSpace, offset, SolnType>(in, res);
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
  A2D_INLINE_FUNCTION static void add(
      const QptSpace<Quadrature, FiniteElementSpace>& in, SolnType& res) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;
      index_t q2dim = 0;
      if constexpr (dim == 3) {
        q2dim = Quadrature::tensor_dim2;
      }

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 2-direction
        index_t u0size, u1size;

        if constexpr (dim == 2) {
          u0size = order * q0dim;
          u1size = order * q0dim;
        } else {  // dim == 3
          u0size = order * order * q0dim;
          u1size = order * q0dim * q1dim;
        }

        T u0[u0size];
        std::fill(u0, u0 + u0size, T(0.0));
        T u1[u1size];
        std::fill(u1, u1 + u1size, T(0.0));

        const index_t qouter = dim == 2 ? 1 : q2dim;
        for (index_t q2 = 0; q2 < qouter; q2++) {
          double n2[order];
          if constexpr (dim == 3) {
            const double pt2 = Quadrature::get_tensor_point(2, q2);
            interpolation_basis<order, interp_type>(pt2, n2);
          }

          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order];
            if constexpr (dim == 2) {
              const double pt1 = Quadrature::get_tensor_point(1, q1);
              interpolation_basis<order, interp_type>(pt1, n1);
            }
            for (index_t q0 = 0; q0 < q0dim; q0++) {
              index_t qindex;
              if constexpr (dim == 2) {
                qindex = Quadrature::get_tensor_index(q0, q1);
              } else {
                qindex = Quadrature::get_tensor_index(q0, q1, q2);
              }

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
                if constexpr (dim == 2) {
                  u0[j2 + order * q0] += n1[j2] * val;
                } else {
                  u1[j2 + order * (q0 + q0dim * q1)] += n2[j2] * val;
                }
              }
            }
          }
        }

        if constexpr (dim == 3) {
          for (index_t q1 = 0; q1 < q1dim; q1++) {
            double n1[order];
            const double pt1 = Quadrature::get_tensor_point(1, q1);
            interpolation_basis<order, interp_type>(pt1, n1);

            for (index_t q0 = 0; q0 < q0dim; q0++) {
              for (index_t j2 = 0; j2 < order; j2++) {
                T val = u1[j2 + order * (q0 + q0dim * q1)];

                for (index_t j1 = 0; j1 < order; j1++) {
                  u0[j1 + order * (j2 + order * q0)] += n1[j1] * val;
                }
              }
            }
          }
        }

        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          interpolation_basis<order, interp_type>(pt0, n0);

          const index_t outer = dim == 2 ? 1 : order;
          for (index_t j2 = 0; j2 < outer; j2++) {
            for (index_t j1 = 0; j1 < order; j1++) {
              T val = u0[j1 + order * (j2 + outer * q0)];

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

        // Evaluate the basis functions
        double n0[order], n1[order], n2[order];
        if constexpr (dim >= 1) {
          interpolation_basis<order, interp_type>(pt[0], n0);
        }
        if constexpr (dim >= 2) {
          interpolation_basis<order, interp_type>(pt[1], n1);
        }
        if constexpr (dim >= 3) {
          interpolation_basis<order, interp_type>(pt[2], n2);
        }

        // | dim | j2range | j1range | j0range |
        // -------------------------------------
        // | 3   | order   | order   | order   |
        // | 2   | 1       | order   | order   |
        // | 1   | 1       | 1       | order   |
        constexpr index_t j2range = dim >= 3 ? order : 1;
        constexpr index_t j1range = dim >= 2 ? order : 1;
        constexpr index_t j0range = dim >= 1 ? order : 1;

        for (index_t j2 = 0; j2 < j2range; j2++) {
          for (index_t j1 = 0; j1 < j1range; j1++) {
            double n1n2;
            if constexpr (dim == 3) {
              n1n2 = n1[j1] * n2[j2];
            }
            for (index_t j0 = 0; j0 < j0range; j0++) {
              const index_t node = j0 + order * (j1 + order * j2);
              double coeff;
              if constexpr (dim == 3) {
                coeff = n0[j0] * n1n2;
              } else if constexpr (dim == 2) {
                coeff = n0[j0] * n1[j1];
              } else {  // dim == 1
                coeff = n0[j0];
              }

              if constexpr (C == 1) {
                res[offset + node] += coeff * u;
              } else {
                for (index_t i = 0; i < C; i++) {
                  res[offset + C * node + i] += coeff * u(i);
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
  static const index_t basis_size = indexpow<order, dim>::value;

  // Set the derived quantities - number of dof for each stride
  static const index_t ndof_per_stride = ndof / stride;

  // Number of components per stride
  static const index_t ncomp_per_stride = ncomp / stride;

  // Compute the full matrix of basis functions
  template <class Quadrature, class BasisType>
  A2D_INLINE_FUNCTION static void basis(index_t n, BasisType N) {
    double pt[dim];
    Quadrature::get_point(n, pt);

    // Evaluate the basis functions
    double n0[order], n1[order], n2[order];

    if constexpr (dim >= 1) {
      interpolation_basis<order, interp_type>(pt[0], n0);
    }
    if constexpr (dim >= 2) {
      interpolation_basis<order, interp_type>(pt[1], n1);
    }
    if constexpr (dim >= 3) {
      interpolation_basis<order, interp_type>(pt[2], n2);
    }

    // | dim | j2range | j1range | j0range |
    // -------------------------------------
    // | 3   | order   | order   | order   |
    // | 2   | 1       | order   | order   |
    // | 1   | 1       | 1       | order   |
    constexpr index_t j2range = dim >= 3 ? order : 1;
    constexpr index_t j1range = dim >= 2 ? order : 1;
    constexpr index_t j0range = dim >= 1 ? order : 1;

    for (index_t j2 = 0; j2 < j2range; j2++) {
      for (index_t j1 = 0; j1 < j1range; j1++) {
        double n1n2;
        if constexpr (dim == 3) {
          n1n2 = n1[j1] * n2[j2];
        }
        for (index_t j0 = 0; j0 < j0range; j0++) {
          const index_t node = j0 + order * (j1 + order * j2);
          if constexpr (dim == 3) {
            N[node] = n0[j0] * n1n2;
          } else if constexpr (dim == 2) {
            N[node] = n0[j0] * n1[j1];
          } else {  // dim == 1
            N[node] = n0[j0];
          }
        }
      }
    }
  }
};

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeH1HexBasis =
    LagrangeH1HypercubeBasis<T, 3, C, degree, interp_type>;

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeL2HexBasis =
    LagrangeL2HypercubeBasis<T, 3, C, degree, interp_type>;

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeH1QuadBasis =
    LagrangeH1HypercubeBasis<T, 2, C, degree, interp_type>;

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeL2QuadBasis =
    LagrangeL2HypercubeBasis<T, 2, C, degree, interp_type>;

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeH1LineBasis =
    LagrangeH1HypercubeBasis<T, 1, C, degree, interp_type>;

template <typename T, index_t C, index_t degree,
          InterpolationType interp_type = GLL_INTERPOLATION>
using LagrangeL2LineBasis =
    LagrangeL2HypercubeBasis<T, 1, C, degree, interp_type>;

}  // namespace A2D

#endif  // A2D_LAGRANGE_HEX_BASIS_H