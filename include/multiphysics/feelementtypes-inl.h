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
 * @param v Vertex index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, class ElemDof,
          class EntityDof>
void ElementTypes::get_line_vert_dof(index_t v, const ElemDof& element,
                                     EntityDof& entity) {
  const index_t start = offset + ndof * (nx - 1) * v;
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
 * @param v Vertex index
 * @param element Element degree of freedom array
 * @param entity Entity degree of freedom array
 */
template <index_t offset, index_t ndof, index_t nx, class EntityDof,
          class ElemDof>
void ElementTypes::set_line_vert_dof(index_t v, const EntityDof& entity,
                                     ElemDof& element) {
  const index_t start = offset + ndof * (nx - 1) * v;
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
void ElementTypes::get_line_edge_dof(const ElemDof& element,
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
        entity[ndof * j + i] = element[offset + ndof * j + i];
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
void ElementTypes::set_line_edge_dof(const index_t orient,
                                     const EntityDof& entity,
                                     ElemDof& element) {
  if constexpr (ends) {
    for (index_t j = 0; j < nx; j++) {
      for (index_t i = 0; i < ndof; i++) {
        element[offset + ndof * j + i] = entity[ndof * j + i];
      }
    }
  } else {
    for (index_t j = 1; j < nx - 1; j++) {
      for (index_t i = 0; i < ndof; i++) {
        element[offset + ndof * j + i] = entity[ndof * j + i];
      }
    }
  }
}

/**
 * @brief Given a reference face, find the orientation
 *
 * @param ref
 * @param face
 * @return index_t
 */
inline index_t ElementTypes::get_quad_face_orientation(const index_t ref[],
                                                       const index_t face[]) {
  index_t orient = 0;
  for (; orient < NUM_QUAD_FACE_ORIENTATIONS; orient++) {
    if (ref[0] == face[QUAD_FACE_ORIENTATIONS[orient][0]] &&
        ref[1] == face[QUAD_FACE_ORIENTATIONS[orient][1]] &&
        ref[2] == face[QUAD_FACE_ORIENTATIONS[orient][2]] &&
        ref[3] == face[QUAD_FACE_ORIENTATIONS[orient][3]]) {
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
inline void ElementTypes::get_coords_on_quad_ref_element(
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

inline index_t ElementTypes::get_index_on_quad_ref_element(const index_t orient,
                                                           const index_t hx,
                                                           const index_t hy,
                                                           const index_t x,
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
int ElementTypes::get_quad_node(const int i, const int j) {
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
index_t ElementTypes::get_quad_edge_length(const index_t v0, const index_t v1) {
  if (QUAD_VERTS_CART[v0][0] != QUAD_VERTS_CART[v1][0]) {
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
void ElementTypes::get_quad_vert_dof(index_t v, const ElemDof& element,
                                     EntityDof& entity) {
  const index_t node = get_quad_node<nx, ny>((nx - 1) * QUAD_VERTS_CART[v][0],
                                             (ny - 1) * QUAD_VERTS_CART[v][1]);

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
void ElementTypes::set_quad_vert_dof(index_t v, const EntityDof& entity,
                                     ElemDof& element) {
  const index_t node = get_quad_node<nx, ny>((nx - 1) * QUAD_VERTS_CART[v][0],
                                             (ny - 1) * QUAD_VERTS_CART[v][1]);

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
 * @param e Edge index
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class ElemDof, class EntityDof>
void ElementTypes::get_quad_edge_dof(const index_t e, const ElemDof& element,
                                     EntityDof& entity) {
  // Get the first and last vertices on the edge
  const index_t v0 = QUAD_EDGE_VERTS[e][0];
  const index_t v1 = QUAD_EDGE_VERTS[e][1];

  // Get the starting index on the element
  const index_t start = get_quad_node<nx, ny>(
      (nx - 1) * QUAD_VERTS_CART[v0][0], (ny - 1) * QUAD_VERTS_CART[v0][1]);

  // Find the number of nodes along the u-edge
  const index_t nu = get_quad_edge_length<nx, ny>(v0, v1);

  // Get the increment
  const int incr =
      get_quad_node<nx, ny>(QUAD_VERTS_CART[v1][0] - QUAD_VERTS_CART[v0][0],
                            QUAD_VERTS_CART[v1][1] - QUAD_VERTS_CART[v0][1]);

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
 * @param e Edge index
 * @param orient Relative orientation between edges
 * @param element Element degrees of freedom
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class EntityDof, class ElemDof>
void ElementTypes::set_quad_edge_dof(const index_t e, const index_t orient,
                                     const EntityDof& entity,
                                     ElemDof& element) {
  // Get the first and last vertices on the edge
  const index_t v0 = QUAD_EDGE_VERTS[e][orient];
  const index_t v1 = QUAD_EDGE_VERTS[e][(orient + 1) % 2];

  // Get the starting index on the element
  const index_t start = get_quad_node<nx, ny>(
      (nx - 1) * QUAD_VERTS_CART[v0][0], (ny - 1) * QUAD_VERTS_CART[v0][1]);

  // Find the number of nodes along the u-edge
  const index_t nu = get_quad_edge_length<nx, ny>(v0, v1);

  // Get the increment
  const int incr =
      get_quad_node<nx, ny>(QUAD_VERTS_CART[v1][0] - QUAD_VERTS_CART[v0][0],
                            QUAD_VERTS_CART[v1][1] - QUAD_VERTS_CART[v0][1]);

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
 * @brief Get the degrees of freedom from the quad face
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param f Face index
 * @param element Element degrees of freedmo
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class ElemDof, class EntityDof>
void ElementTypes::get_quad_face_dof(const ElemDof& element,
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
 * @brief Set the degrees of freeom from the face into the element
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param f Face index
 * @param orient Relative orientation between faces
 * @param element Element degrees of freedmo
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          class EntityDof, class ElemDof>
void ElementTypes::set_quad_face_dof(const index_t orient,
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
int ElementTypes::get_hex_node(const int i, const int j, const int k) {
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
index_t ElementTypes::get_hex_edge_length(const index_t v0, const index_t v1) {
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
void ElementTypes::get_hex_vert_dof(index_t v, const ElemDof& element,
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
void ElementTypes::set_hex_vert_dof(index_t v, const EntityDof& entity,
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
void ElementTypes::get_hex_edge_dof(const index_t e, const ElemDof& element,
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
void ElementTypes::set_hex_edge_dof(const index_t e, const index_t orient,
                                    const EntityDof& entity, ElemDof& element) {
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
 * @brief Get the degrees of freedom from the hex face
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param f Face index
 * @param element Element degrees of freedmo
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class ElemDof, class EntityDof>
void ElementTypes::get_hex_face_dof(const index_t f, const ElemDof& element,
                                    EntityDof& entity) {
  // Get the origin and face vert directions
  const index_t v0 = HEX_FACE_VERTS[f][0];  // Root vertex
  const index_t v1 = HEX_FACE_VERTS[f][1];  // Vertex along the face u dir
  const index_t v3 = HEX_FACE_VERTS[f][3];  // Vertex along the face v dir

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
 * @brief Set the degrees of freeom from the face into the element
 *
 * @tparam offset Offset into the element dof array
 * @tparam ends Include the end points of the edge or not
 * @tparam ndof The nuber of degrees of freedom at each node
 * @tparam nx Number of nodes along the x-direction
 * @tparam ny Number of nodes along the y-direction
 * @tparam nz Number of nodes along the z-direction
 * @tparam ElemDof Element degrees of freedom array type
 * @tparam EntityDof Entity degrees of freedom array type
 * @param f Face index
 * @param orient Relative orientation between faces
 * @param element Element degrees of freedmo
 * @param entity Entity degrees of freedom
 */
template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
          index_t nz, class EntityDof, class ElemDof>
void ElementTypes::set_hex_face_dof(const index_t f, const index_t orient,
                                    const EntityDof& entity, ElemDof& element) {
  // Get the origin and face vert directions
  const index_t v0 = HEX_FACE_VERTS[f][0];  // Root vertex
  const index_t v1 = HEX_FACE_VERTS[f][1];  // Vertex along the face u dir
  const index_t v3 = HEX_FACE_VERTS[f][3];  // Vertex along the face v dir

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
 * @brief Get the degrees of freeom from the volume
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
void ElementTypes::get_hex_volume_dof(const ElemDof& element,
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
void ElementTypes::set_hex_volume_dof(const EntityDof& entity,
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