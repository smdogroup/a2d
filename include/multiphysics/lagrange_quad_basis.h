#ifndef A2D_LAGRANGE_QUAD_BASIS_H
#define A2D_LAGRANGE_QUAD_BASIS_H

#include "multiphysics/febasis.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/lagrange_tools.h"

namespace A2D {

template <typename T, index_t C, index_t degree>
class LagrangeH1QuadBasis {
 public:
  using ET = ElementTypes;

  static const index_t dim = 2;             // Spatial dimension
  static const index_t order = degree + 1;  // Number of nodes along each edge

  static const index_t ndof =
      C * order * order;  // Total number of degrees of freedom

  // Number of components
  static const index_t ncomp = H1Space<T, C, dim>::ncomp;

  // Get the type of basis class implemented
  static constexpr BasisType get_basis_type() { return H1; }

  // Define the equivalent low-order basis class if any
  using LOrderBasis = LagrangeH1HexBasis<T, C, 1>;

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
      ET::get_quad_vert_dof<offset, C, order, order, ElemDof, EntityDof>(
          index, element_dof, entity_dof);
    } else if (entity == ET::EDGE) {
      const bool endp = false;
      ET::get_quad_edge_dof<offset, endp, C, order, order, ElemDof, EntityDof>(
          index, element_dof, entity_dof);
    } else if (entity == ET::FACE) {
      const bool endp = false;
      ET::get_quad_face_dof<offset, endp, C, order, order, ElemDof, EntityDof>(
          index, element_dof, entity_dof);
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
      ET::set_quad_vert_dof<offset, C, order, order, EntityDof, ElemDof>(
          index, entity_dof, element_dof);
    } else if (entity == ET::EDGE) {
      const bool endp = false;
      ET::set_quad_edge_dof<offset, endp, C, order, order, EntityDof, ElemDof>(
          index, orient, entity_dof, element_dof);
    } else if (entity == ET::FACE) {
      const bool endp = false;
      ET::set_quad_face_dof<offset, endp, C, order, order, EntityDof, ElemDof>(
          index, orient, entity_dof, element_dof);
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
   * @brief Get the number of low order elements defined by this basis
   */
  constexpr static index_t get_num_lorder_elements() { return degree * degree; }

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
    const index_t j = n / (degree * degree);

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
    for (index_t p = 0; p < 4 * C; p++) {
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
  static void get_dof_point(index_t index, double pt[]) {
    constexpr const double* pts = get_gauss_lobatto_pts<order>();

    index_t n = index / C;
    pt[0] = pts[n % order];
    pt[1] = pts[n / order];
  }

  /**
   * @brief Interpolate over all quadrature points using
   *
   * This function just calls the interpolation the basis functions are
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
  static void interp(const QptGeoSpace& geo, const SolnType& sol,
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
  static void interp(const SolnType& sol,
                     QptSpace<Quadrature, FiniteElementSpace>& out) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 0-direction
        T u0[order * q0dim];
        T u0x[order * q0dim];
        for (index_t q0 = 0; q0 < q0dim; q0++) {
          double n0[order], d0[order];
          const double pt0 = Quadrature::get_tensor_point(0, q0);
          lagrange_basis<order>(pt0, n0, d0);

          for (index_t j1 = 0; j1 < order; j1++) {
            T val(0.0), derx(0.0);
            for (index_t j0 = 0; j0 < order; j0++) {
              const index_t node = j0 + order * j1;
              val += n0[j0] * sol[offset + C * node + i];
              derx += d0[j0] * sol[offset + C * node + i];
            }

            u0[j1 + order * q0] = val;
            u0x[j1 + order * q0] = derx;
          }
        }

        // Interpolate along the 1-direction
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order], d1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1, d1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            T val(0.0), derx(0.0), dery(0.0);
            for (index_t j1 = 0; j1 < order; j1++) {
              val += n1[j1] * u0[j1 + order * q0];
              derx += n1[j2] * u0x[j1 + order * q0];
              dery += d1[j2] * u0[j1 + order * q0];
            }

            const index_t qindex = Quadrature::get_tensor_index(q0, q1);
            FiniteElementSpace& s = out.get(qindex);
            H1Space<T, C, dim>& h1 = s.template get<space>();
            typename H1Space<T, C, dim>::VarType& u = h1.get_value();
            typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

            if constexpr (C == 1) {
              u = val;
              grad(0) = derx;
              grad(1) = dery;

            } else {
              u(i) = val;
              grad(i, 0) = derx;
              grad(i, 1) = dery;
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
        lagrange_basis<order>(pt[0], n0, d0);
        lagrange_basis<order>(pt[1], n1, d1);

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

        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t j0 = 0; j0 < order; j0++) {
            const index_t node = j0 + order * (j1 + order * j2);
            double N = n0[j0] * n1[j1];
            double dx = d0[j0] * n1[j1];
            double dy = n0[j0] * d1[j1];

            if constexpr (C == 1) {
              const T val = sol[offset + node];

              u += N * val;
              grad(0) += dx * val;
              grad(1) += dy * val;
            } else {
              for (index_t i = 0; i < C; i++) {
                const T val = sol[offset + C * node + i];

                u(i) += N * val;
                grad(i, 0) += dx * val;
                grad(i, 1) += dy * val;
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
  static void add(const QptSpace<Quadrature, FiniteElementSpace>& in,
                  SolnType& res) {
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
  static void add(const QptSpace<Quadrature, FiniteElementSpace>& in,
                  SolnType& res) {
    if constexpr (Quadrature::is_tensor_product) {
      const index_t q0dim = Quadrature::tensor_dim0;
      const index_t q1dim = Quadrature::tensor_dim1;

      for (index_t i = 0; i < C; i++) {
        // Interpolate along the 2-direction
        T u0[order * q0dim];
        T u0x[order * q0dim];
        std::fill(u0, u0 + order * q0dim, T(0.0));
        std::fill(u0x, u0x + order * q0dim, T(0.0));
        for (index_t q1 = 0; q1 < q1dim; q1++) {
          double n1[order], d1[order];
          const double pt1 = Quadrature::get_tensor_point(1, q1);
          lagrange_basis<order>(pt1, n1, d1);

          for (index_t q0 = 0; q0 < q0dim; q0++) {
            const index_t qindex = Quadrature::get_tensor_index(q0, q1);

            const FiniteElementSpace& s = in.get(qindex);
            const H1Space<T, C, dim>& h1 = s.template get<space>();
            const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
            const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

            T val, derx, dery;
            if constexpr (C == 1) {
              val = u;
              derx = grad(0);
              dery = grad(1);
            } else {
              val = u(i);
              derx = grad(i, 0);
              dery = grad(i, 1);
            }
            for (index_t j1 = 0; j1 < order; j1++) {
              u0[j1 + order * q0] += n1[j1] * val + d1[j1] * dery;
              u0x[j1 + order * q0] += n1[j1] * derx;
            }
          }
        }
      }

      for (index_t q0 = 0; q0 < q0dim; q0++) {
        double n0[order], d0[order];
        const double pt0 = Quadrature::get_tensor_point(0, q0);
        lagrange_basis<order>(pt0, n0, d0);

        for (index_t j1 = 0; j1 < order; j1++) {
          T val = u0[j1 + order * q0];
          T derx = u0x[j1 + order * q0];

          for (index_t j0 = 0; j0 < order; j0++) {
            const index_t node = j0 + order * j1;
            res[offset + C * node + i] += n0[j0] * val + d0[j0] * derx;
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
        lagrange_basis<order>(pt[0], n0, d0);
        lagrange_basis<order>(pt[1], n1, d1);

        const FiniteElementSpace& s = in.get(q);
        const H1Space<T, C, dim>& h1 = s.template get<space>();
        const typename H1Space<T, C, dim>::VarType& u = h1.get_value();
        const typename H1Space<T, C, dim>::GradType& grad = h1.get_grad();

        for (index_t j1 = 0; j1 < order; j1++) {
          for (index_t j0 = 0; j0 < order; j0++) {
            const index_t node = j0 + order * (j1 + order * j2);
            double N = n0[j0] * n1[j1];
            double dx = d0[j0] * n1[j1];
            double dy = n0[j0] * d1[j1];

            if constexpr (C == 1) {
              res[offset + node] += N * u + dx * grad(0) + dy * grad(1);
            } else {
              for (index_t i = 0; i < C; i++) {
                res[offset + C * node + i] +=
                    N * u(i) + dx * grad(i, 0) + dy * grad(i, 1);
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
  static const index_t basis_size = (dim + 1) * order * order;

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
    double n0[order], n1[order];
    double d0[order], d1[order];
    lagrange_basis<order>(pt[0], n0, d0);
    lagrange_basis<order>(pt[1], n1, d1);

    for (index_t j1 = 0; j1 < order; j1++) {
      for (index_t j0 = 0; j0 < order; j0++) {
        const index_t node = j0 + order * j1;
        N[(dim + 1) * node] = n0[j0] * n1[j1];
        N[(dim + 1) * node + 1] = d0[j0] * n1[j1];
        N[(dim + 1) * node + 2] = n0[j0] * d1[j1];
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_LAGRANGE_QUAD_BASIS_H