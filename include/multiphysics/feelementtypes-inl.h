#ifndef A2D_FE_ELEMENT_TYPES_INL_H
#define A2D_FE_ELEMENT_TYPES_INL_H

#include "feelementtypes.h"

namespace A2D {

/**
 * @brief Get the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the line
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param v Bound index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, class ElemDof,
          class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_line_bound_dof(index_t b,
                                                      const ElemDof& element,
                                                      EntityDof& entity) {
  const index_t start = offset + ndof * (nx - 1) * b;
  for (index_t i = 0; i < ndof; i++) {
    entity[i] = element[start + i];
  }
}

/**
 * @brief Set the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the line
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param b Bound index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, class EntityDof,
          class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_line_bound_dof(index_t b,
                                                      const index_t orient,
                                                      const EntityDof& entity,
                                                      ElemDof& element) {
  const index_t start = offset + ndof * (nx - 1) * ((b + orient) % 2);
  for (index_t i = 0; i < ndof; i++) {
    element[start + i] = entity[i];
  }
}

/**
 * @brief Get the degrees of freedom from the edge
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the line
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, class ElemDof,
          class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_line_domain_dof(const ElemDof& element,
                                                       EntityDof& entity) {
  if constexpr (ends) {
    for (index_t j = 0; j < nx; j++) {
      for (index_t i = 0; i < ndof; i++) {
        entity[ndof * j + i] = element[offset + ndof * j + i];
      }
    }
  } else {
    for (index_t j = 1; j < nx - 1; j++) {
      for (index_t i = 0; i < ndof; i++) {
        entity[ndof * (j - 1) + i] = element[offset + ndof * j + i];
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom from the edge
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the line
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, class EntityDof,
          class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_line_domain_dof(const index_t orient,
                                                       const EntityDof& entity,
                                                       ElemDof& element) {
  index_t index;
  if constexpr (ends) {
    for (index_t j = 0; j < nx; j++) {
      for (index_t i = 0; i < ndof; i++) {
        index = offset + ndof * j;
        if (orient == 1) {
          index = offset + ndof * (nx - 1 - j);
        }
        element[index + i] = entity[ndof * j + i];
      }
    }
  } else {
    for (index_t j = 1; j < nx - 1; j++) {
      for (index_t i = 0; i < ndof; i++) {
        index = offset + ndof * j;
        if (orient == 1) {
          index = offset + ndof * (nx - 1 - j);
        }
        element[index + i] = entity[ndof * (j - 1) + i];
      }
    }
  }
}

/**
 * @brief Given a reference domain, find the orientation
 *
 * @param ref_domain_verts
 * @param domain_verts
 * @return index_t
 */
KOKKOS_FUNCTION inline index_t ElementTypes::get_quad_domain_orientation(
    const index_t ref_domain_verts[], const index_t domain_verts[]) {
  static const index_t QUAD_DOMAIN_ORIENTATIONS[][4] = {
      {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
      {0, 3, 2, 1}, {3, 2, 1, 0}, {2, 1, 0, 3}, {1, 0, 3, 2}};

  index_t orient = 0;
  for (; orient < NUM_QUAD_DOMAIN_ORIENTATIONS; orient++) {
    if (ref_domain_verts[0] ==
            domain_verts[QUAD_DOMAIN_ORIENTATIONS[orient][0]] &&
        ref_domain_verts[1] ==
            domain_verts[QUAD_DOMAIN_ORIENTATIONS[orient][1]] &&
        ref_domain_verts[2] ==
            domain_verts[QUAD_DOMAIN_ORIENTATIONS[orient][2]] &&
        ref_domain_verts[3] ==
            domain_verts[QUAD_DOMAIN_ORIENTATIONS[orient][3]]) {
      break;
    }
  }
  return orient;
}

/**
 * @brief Get the coords on quad ref element object
 *
 * @param orient The orientation
 * @param hx The edge length along the x-direction
 * @param hy The edge length along the y-direction
 * @param x The 1st coordinate on the transformed face
 * @param y The 2nd coordinate on the transformed face
 * @param u The 1st coordinate on the reference face
 * @param v The 2nd coordinate on the reference face
 */
KOKKOS_FUNCTION inline void ElementTypes::get_coords_on_quad_ref_element(
    const index_t orient, const index_t hx, const index_t hy, const index_t x,
    const index_t y, index_t* u, index_t* v) {
  if (orient == 0) {
    *u = x;
    *v = y;
  } else if (orient == 1) {
    *u = hy - y;
    *v = x;
  } else if (orient == 2) {
    *u = hx - x;
    *v = hy - y;
  } else if (orient == 3) {
    *u = y;
    *v = hx - x;
  } else if (orient == 4) {
    *u = y;
    *v = x;
  } else if (orient == 5) {
    *u = x;
    *v = hy - y;
  } else if (orient == 6) {
    *u = hy - y;
    *v = hx - x;
  } else if (orient == 7) {
    *u = hx - x;
    *v = y;
  } else {
    *u = 0;
    *v = 0;
  }
}

KOKKOS_FUNCTION inline index_t ElementTypes::get_index_on_quad_ref_element(
    const index_t orient, const index_t hx, const index_t hy, const index_t x,
    const index_t y) {
  if (orient == 0) {
    // *u = x;
    // *v = y;
    return x + hx * y;
  } else if (orient == 1) {
    // *u = hy - y;
    // *v = x;
    return (hy - 1 - y) + hy * x;
  } else if (orient == 2) {
    // *u = hx - x;
    // *v = hy - y;
    return (hx - 1 - x) + hx * (hy - 1 - y);
  } else if (orient == 3) {
    // *u = y;
    // *v = hx - x;
    return y * hy * (hx - 1 - x);
  } else if (orient == 4) {
    // *u = y;
    // *v = x;
    return y + hy * x;
  } else if (orient == 5) {
    // *u = x;
    // *v = hy - y;
    return x + hx * (hy - 1 - y);
  } else if (orient == 6) {
    // *u = hy - y;
    // *v = hx - x;
    return (hy - 1 - y) + hy * (hx - 1 - x);
  } else if (orient == 7) {
    // *u = hx - x;
    // *v = y;
    return (hx - 1 - x) + hx * y;
  } else {
    return 0;
  }
}

/**
 * @brief Get the local node index (without offset) for a node on an element
 * with nx and ny nodes along the local directions
 *
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @param i Index along the x-direction
 * @param j Index along the y-direction
 */
template <index_t nx, index_t ny>
KOKKOS_FUNCTION int ElementTypes::get_quad_node(const int i, const int j) {
  return i + nx * j;
}

/**
 * @brief Get the edge length between the verties v1 and v2
 *
 * This only works if the verties are connected by a single edge
 *
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @param v0 First vertex
 * @param v1 Second vertex
 * @return The number of nodes along the edge
 */
template <index_t nx, index_t ny>
KOKKOS_FUNCTION index_t ElementTypes::get_quad_bound_length(const index_t v0,
                                                            const index_t v1) {
  if (get_quad_verts_cart(v0, 0) != get_quad_verts_cart(v1, 0)) {
    return nx;
  } else {
    return ny;
  }
}

/**
 * @brief Get the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param v Vertex index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, index_t ny, class ElemDof,
          class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_quad_vert_dof(index_t v,
                                                     const ElemDof& element,
                                                     EntityDof& entity) {
  const index_t node =
      get_quad_node<nx, ny>((nx - 1) * get_quad_verts_cart(v, 0),
                            (ny - 1) * get_quad_verts_cart(v, 1));

  const index_t start = offset + ndof * node;
  for (index_t i = 0; i < ndof; i++) {
    entity[i] = element[start + i];
  }
}

/**
 * @brief Get the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam EntityDof Entity degrees of freedom array type
 * @tparam ElemDof Element degrees of freedom array type
 * @param v Vertex index
 * @param entity Entity degree of freedom array
 * @param element Element degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, index_t ny, class EntityDof,
          class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_quad_vert_dof(index_t v,
                                                     const EntityDof& entity,
                                                     ElemDof& element) {
  const index_t node =
      get_quad_node<nx, ny>((nx - 1) * get_quad_verts_cart(v, 0),
                            (ny - 1) * get_quad_verts_cart(v, 1));

  const index_t start = offset + ndof * node;
  for (index_t i = 0; i < ndof; i++) {
    element[start + i] = entity[i];
  }
}

/**
 * @brief Get the degrees of freedom from the edge
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param b Bound index
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_quad_bound_dof(const index_t b,
                                                      const ElemDof& element,
                                                      EntityDof& entity) {
  // Get the first and last vertices on the edge
  const index_t v0 = QUAD_BOUND_VERTS[b][0];
  const index_t v1 = QUAD_BOUND_VERTS[b][1];

  // Get the starting index on the element
  const index_t start =
      get_quad_node<nx, ny>((nx - 1) * get_quad_verts_cart(v0, 0),
                            (ny - 1) * get_quad_verts_cart(v0, 1));

  // Find the number of nodes along the u-edge
  const index_t nu = get_quad_bound_length<nx, ny>(v0, v1);

  // Get the increment
  const int incr = get_quad_node<nx, ny>(
      get_quad_verts_cart(v1, 0) - get_quad_verts_cart(v0, 0),
      get_quad_verts_cart(v1, 1) - get_quad_verts_cart(v0, 1));

  if constexpr (ends) {
    index_t index = start;
    for (index_t u = 0; u < nu; u++, index += incr) {
      const index_t entity_index = ndof * u;
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        entity[entity_index + i] = element[elem_index + i];
      }
    }
  } else {
    index_t index = start + incr;
    for (index_t u = 1; u < nu - 1; u++, index += incr) {
      const index_t entity_index = ndof * (u - 1);
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        entity[entity_index + i] = element[elem_index + i];
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom from the edge
 *
 * The orient argument indicates whether the edges have the same orientation
 * (orient = 0) or opposite orientations (orient = 1)
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param b Bound index
 * @param orient Relative orientation between edges
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_quad_bound_dof(const index_t b,
                                                      const index_t orient,
                                                      const EntityDof& entity,
                                                      ElemDof& element) {
  // Get the first and last vertices on the edge
  const index_t v0 = QUAD_BOUND_VERTS[b][orient];
  const index_t v1 = QUAD_BOUND_VERTS[b][(orient + 1) % 2];

  // Get the starting index on the element
  const index_t start =
      get_quad_node<nx, ny>((nx - 1) * get_quad_verts_cart(v0, 0),
                            (ny - 1) * get_quad_verts_cart(v0, 1));

  // Find the number of nodes along the u-edge
  const index_t nu = get_quad_bound_length<nx, ny>(v0, v1);

  // Get the increment
  const int incr = get_quad_node<nx, ny>(
      get_quad_verts_cart(v1, 0) - get_quad_verts_cart(v0, 0),
      get_quad_verts_cart(v1, 1) - get_quad_verts_cart(v0, 1));

  if constexpr (ends) {
    index_t index = start;
    for (index_t u = 0; u < nu; u++, index += incr) {
      const index_t entity_index = ndof * u;
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        element[elem_index + i] = entity[entity_index + i];
      }
    }
  } else {
    index_t index = start + incr;
    for (index_t u = 1; u < nu - 1; u++, index += incr) {
      const index_t entity_index = ndof * (u - 1);
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        element[elem_index + i] = entity[entity_index + i];
      }
    }
  }
}

/**
 * @brief Get the degrees of freedom from the quad domain
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_quad_domain_dof(const ElemDof& element,
                                                       EntityDof& entity) {
  if constexpr (ends) {
    for (index_t v = 0; v < ny; v++) {
      for (index_t u = 0; u < nx; u++) {
        const index_t entity_index = ndof * get_quad_node<nx, ny>(u, v);
        const index_t elem_index = offset + entity_index;

        for (index_t i = 0; i < ndof; i++) {
          entity[entity_index + i] = element[elem_index + i];
        }
      }
    }
  } else {
    for (index_t v = 1; v < ny - 1; v++) {
      for (index_t u = 1; u < nx - 1; u++) {
        const index_t entity_index =
            ndof * get_quad_node<nx - 2, ny - 2>(u - 1, v - 1);
        const index_t elem_index = offset + ndof * get_quad_node<nx, ny>(u, v);

        for (index_t i = 0; i < ndof; i++) {
          entity[entity_index + i] = element[elem_index + i];
        }
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom from the domain into the element
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param orient Relative orientation between domains
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_quad_domain_dof(const index_t orient,
                                                       const EntityDof& entity,
                                                       ElemDof& element) {
  if constexpr (ends) {
    for (index_t v = 0; v < ny; v++) {
      for (index_t u = 0; u < nx; u++) {
        const index_t entity_index =
            ndof * get_index_on_quad_ref_element(orient, nx, ny, u, v);

        const index_t index = get_quad_node<nx, ny>(u, v);
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          element[elem_index + i] = entity[entity_index + i];
        }
      }
    }
  } else {
    for (index_t v = 1; v < ny - 1; v++) {
      for (index_t u = 1; u < nx - 1; u++) {
        const index_t entity_index =
            ndof *
            get_index_on_quad_ref_element(orient, nx - 2, ny - 2, u - 1, v - 1);
        const index_t index = get_quad_node<nx, ny>(u, v);
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          element[elem_index + i] = entity[entity_index + i];
        }
      }
    }
  }
}

/**
 * @brief Get the local node index (without offset) for a node on an element
 * with nx, ny and nz nodes along the local directions
 *
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @param i Index along the x-direction
 * @param j Index along the y-direction
 * @param k Index along the z-direction
 */
template <index_t nx, index_t ny, index_t nz>
KOKKOS_FUNCTION int ElementTypes::get_hex_node(const int i, const int j,
                                               const int k) {
  return i + nx * (j + ny * k);
}

/**
 * @brief Get the edge length between the verties v1 and v2
 *
 * This only works if the verties are connected by a single edge
 *
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @param v0 First vertex
 * @param v1 Second vertex
 * @return The number of nodes along the edge
 */
template <index_t nx, index_t ny, index_t nz>
KOKKOS_FUNCTION index_t ElementTypes::get_hex_edge_length(const index_t v0,
                                                          const index_t v1) {
  if (HEX_VERTS_CART[v0][0] != HEX_VERTS_CART[v1][0]) {
    return nx;
  } else if (HEX_VERTS_CART[v0][1] != HEX_VERTS_CART[v1][1]) {
    return ny;
  } else {
    return nz;
  }
}

/**
 * @brief Get the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param v Vertex index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
          class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_hex_vert_dof(index_t v,
                                                    const ElemDof& element,
                                                    EntityDof& entity) {
  const index_t node = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v][0], (ny - 1) * HEX_VERTS_CART[v][1],
      (nz - 1) * HEX_VERTS_CART[v][2]);

  const index_t start = offset + ndof * node;
  for (index_t i = 0; i < ndof; i++) {
    entity[i] = element[start + i];
  }
}

/**
 * @brief Get the degrees of freedom associated with the vertex
 *
 * @tparam offset Offset into the global degree of freedom array
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam EntityDof Entity degrees of freedom array type
 * @tparam ElemDof Element degrees of freedom array type
 * @param v Vertex index
 * @param entity Entity degree of freedom array
 * @param element Element degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
          class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_hex_vert_dof(index_t v,
                                                    const EntityDof& entity,
                                                    ElemDof& element) {
  const index_t node = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v][0], (ny - 1) * HEX_VERTS_CART[v][1],
      (nz - 1) * HEX_VERTS_CART[v][2]);

  const index_t start = offset + ndof * node;
  for (index_t i = 0; i < ndof; i++) {
    element[start + i] = entity[i];
  }
}

/**
 * @brief Get the degrees of freedom from the edge
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param e Edge index
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_hex_edge_dof(const index_t e,
                                                    const ElemDof& element,
                                                    EntityDof& entity) {
  // Get the first and last vertices on the edge
  const index_t v0 = HEX_EDGE_VERTS[e][0];
  const index_t v1 = HEX_EDGE_VERTS[e][1];

  // Get the starting index on the element
  const index_t start = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v0][0], (ny - 1) * HEX_VERTS_CART[v0][1],
      (nz - 1) * HEX_VERTS_CART[v0][2]);

  // Find the number of nodes along the u-edge
  const index_t nu = get_hex_edge_length<nx, ny, nz>(v0, v1);

  // Get the increment
  const int incr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v1][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v1][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v1][2] - HEX_VERTS_CART[v0][2]);

  if constexpr (ends) {
    index_t index = start;
    for (index_t u = 0; u < nu; u++, index += incr) {
      const index_t entity_index = ndof * u;
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        entity[entity_index + i] = element[elem_index + i];
      }
    }
  } else {
    index_t index = start + incr;
    for (index_t u = 1; u < nu - 1; u++, index += incr) {
      const index_t entity_index = ndof * (u - 1);
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        entity[entity_index + i] = element[elem_index + i];
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom from the edge
 *
 * The orient argument indicates whether the edges have the same orientation
 * (orient = 0) or opposite orientations (orient = 1)
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param e Edge index
 * @param orient Relative orientation between edges
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_hex_edge_dof(const index_t e,
                                                    const index_t orient,
                                                    const EntityDof& entity,
                                                    ElemDof& element) {
  // Get the first and last vertices on the edge
  const index_t v0 = HEX_EDGE_VERTS[e][orient];
  const index_t v1 = HEX_EDGE_VERTS[e][(orient + 1) % 2];

  // Get the starting index on the element
  const index_t start = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v0][0], (ny - 1) * HEX_VERTS_CART[v0][1],
      (nz - 1) * HEX_VERTS_CART[v0][2]);

  // Find the number of nodes along the u-edge
  const index_t nu = get_hex_edge_length<nx, ny, nz>(v0, v1);

  // Get the increment
  const int incr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v1][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v1][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v1][2] - HEX_VERTS_CART[v0][2]);

  if constexpr (ends) {
    index_t index = start;
    for (index_t u = 0; u < nu; u++, index += incr) {
      const index_t entity_index = ndof * u;
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        element[elem_index + i] = entity[entity_index + i];
      }
    }
  } else {
    index_t index = start + incr;
    for (index_t u = 1; u < nu - 1; u++, index += incr) {
      const index_t entity_index = ndof * (u - 1);
      const index_t elem_index = offset + ndof * index;

      for (index_t i = 0; i < ndof; i++) {
        element[elem_index + i] = entity[entity_index + i];
      }
    }
  }
}

/**
 * @brief Get the degrees of freedom from the hex bound
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param b Bound index
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_hex_bound_dof(const index_t b,
                                                     const ElemDof& element,
                                                     EntityDof& entity) {
  // Get the origin and bound vert directions
  const index_t v0 = HEX_BOUND_VERTS[b][0];  // Root vertex
  const index_t v1 = HEX_BOUND_VERTS[b][1];  // Vertex along the bound u dir
  const index_t v3 = HEX_BOUND_VERTS[b][3];  // Vertex along the bound v dir

  // Get the root node location
  const index_t start = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v0][0], (ny - 1) * HEX_VERTS_CART[v0][1],
      (nz - 1) * HEX_VERTS_CART[v0][2]);

  // Find the number of nodes along the u and v edges
  const index_t nu = get_hex_edge_length<nx, ny, nz>(v0, v1);
  const index_t nv = get_hex_edge_length<nx, ny, nz>(v0, v3);

  // Find the increment in u
  const int uincr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v1][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v1][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v1][2] - HEX_VERTS_CART[v0][2]);

  // Find the increment in v
  const int vincr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v3][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v3][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v3][2] - HEX_VERTS_CART[v0][2]);

  if constexpr (ends) {
    for (index_t v = 0; v < nv; v++) {
      for (index_t u = 0; u < nu; u++) {
        const index_t entity_index = ndof * (u + nu * v);

        const index_t index = start + u * uincr + v * vincr;
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          entity[entity_index + i] = element[elem_index + i];
        }
      }
    }
  } else {
    for (index_t v = 1; v < nv - 1; v++) {
      for (index_t u = 1; u < nu - 1; u++) {
        const index_t entity_index = ndof * (u - 1 + (nu - 2) * (v - 1));
        const index_t index = start + u * uincr + v * vincr;
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          entity[entity_index + i] = element[elem_index + i];
        }
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom from the bound into the element
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param b Bound index
 * @param orient Relative orientation between bounds
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_hex_bound_dof(const index_t b,
                                                     const index_t orient,
                                                     const EntityDof& entity,
                                                     ElemDof& element) {
  // Get the origin and bound vert directions
  const index_t v0 = HEX_BOUND_VERTS[b][0];  // Root vertex
  const index_t v1 = HEX_BOUND_VERTS[b][1];  // Vertex along the bound u dir
  const index_t v3 = HEX_BOUND_VERTS[b][3];  // Vertex along the bound v dir

  // Get the root node location
  const index_t start = get_hex_node<nx, ny, nz>(
      (nx - 1) * HEX_VERTS_CART[v0][0], (ny - 1) * HEX_VERTS_CART[v0][1],
      (nz - 1) * HEX_VERTS_CART[v0][2]);

  // Find the number of nodes along the u and v edges
  const index_t nu = get_hex_edge_length<nx, ny, nz>(v0, v1);
  const index_t nv = get_hex_edge_length<nx, ny, nz>(v0, v3);

  // Find the increment in u
  const int uincr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v1][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v1][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v1][2] - HEX_VERTS_CART[v0][2]);

  // Find the increment in v
  const int vincr =
      get_hex_node<nx, ny, nz>(HEX_VERTS_CART[v3][0] - HEX_VERTS_CART[v0][0],
                               HEX_VERTS_CART[v3][1] - HEX_VERTS_CART[v0][1],
                               HEX_VERTS_CART[v3][2] - HEX_VERTS_CART[v0][2]);

  if constexpr (ends) {
    for (index_t v = 0; v < nv; v++) {
      for (index_t u = 0; u < nu; u++) {
        const index_t entity_index =
            ndof * get_index_on_quad_ref_element(orient, nu, nv, u, v);

        const index_t index = start + u * uincr + v * vincr;
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          element[elem_index + i] = entity[entity_index + i];
        }
      }
    }
  } else {
    for (index_t v = 1; v < nv - 1; v++) {
      for (index_t u = 1; u < nu - 1; u++) {
        const index_t entity_index =
            ndof *
            get_index_on_quad_ref_element(orient, nu - 2, nv - 2, u - 1, v - 1);
        const index_t index = start + u * uincr + v * vincr;
        const index_t elem_index = offset + ndof * index;

        for (index_t i = 0; i < ndof; i++) {
          element[elem_index + i] = entity[entity_index + i];
        }
      }
    }
  }
}

/**
 * @brief Get the degrees of freedom from the volume
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class ElemDof, class EntityDof>
KOKKOS_FUNCTION void ElementTypes::get_hex_domain_dof(const ElemDof& element,
                                                      EntityDof& entity) {
  if constexpr (ends) {
    for (index_t w = 0; w < nz; w++) {
      for (index_t v = 0; v < ny; v++) {
        for (index_t u = 0; u < nx; u++) {
          const index_t entity_index = ndof * get_hex_node<nx, ny, nz>(u, v, w);
          const index_t elem_index = offset + entity_index;

          for (index_t i = 0; i < ndof; i++) {
            entity[entity_index + i] = element[elem_index + i];
          }
        }
      }
    }
  } else {
    for (index_t w = 1; w < nz - 1; w++) {
      for (index_t v = 1; v < ny - 1; v++) {
        for (index_t u = 1; u < nx - 1; u++) {
          const index_t entity_index =
              ndof * get_hex_node<nx - 2, ny - 2, nz - 2>(u - 1, v - 1, w - 1);
          const index_t elem_index =
              offset + ndof * get_hex_node<nx, ny, nz>(u, v, w);

          for (index_t i = 0; i < ndof; i++) {
            entity[entity_index + i] = element[elem_index + i];
          }
        }
      }
    }
  }
}

/**
 * @brief Set the degrees of freedom into the element array
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class EntityDof, class ElemDof>
KOKKOS_FUNCTION void ElementTypes::set_hex_domain_dof(const EntityDof& entity,
                                                      ElemDof& element) {
  if constexpr (ends) {
    for (index_t w = 0; w < nz; w++) {
      for (index_t v = 0; v < ny; v++) {
        for (index_t u = 0; u < nx; u++) {
          const index_t entity_index = ndof * get_hex_node<nx, ny, nz>(u, v, w);
          const index_t elem_index = offset + entity_index;

          for (index_t i = 0; i < ndof; i++) {
            element[elem_index + i] = entity[entity_index + i];
          }
        }
      }
    }
  } else {
    for (index_t w = 1; w < nz - 1; w++) {
      for (index_t v = 1; v < ny - 1; v++) {
        for (index_t u = 1; u < nx - 1; u++) {
          const index_t entity_index =
              ndof * get_hex_node<nx - 2, ny - 2, nz - 2>(u - 1, v - 1, w - 1);
          const index_t elem_index =
              offset + ndof * get_hex_node<nx, ny, nz>(u, v, w);

          for (index_t i = 0; i < ndof; i++) {
            element[elem_index + i] = entity[entity_index + i];
          }
        }
      }
    }
  }
}

};  // namespace A2D

#endif  // A2D_FE_ELEMENT_TYPES_INL_H