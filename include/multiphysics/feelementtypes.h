#ifndef A2D_FE_ELEMENT_TYPES_H
#define A2D_FE_ELEMENT_TYPES_H

#include "a2dobjs.h"

/**
 * @brief Topological entities
 *
 * We first define our naming conventions for the geometric entities for
 * elements with different dimension (0, 1, 2, 3, etc.):
 *
 * |                  | 3D element  | 2D element  | 1D element  | 0D element  |
 * | Entity dimension | (Hex etc.)  | (Quad etc.) | (Line etc.) |             |
 *  --------------------------------------------------------------------------
 * | 3D               | domain      |             |             |             |
 * | 2D               | bound       | domain      |             |             |
 * | 1D               | edge        | bound       | domain      |             |
 * | 0D               | vertex      | vertex      | bound       | domain      |
 *
 * Note: in our terminology, bound and boundary are different and need not to be
 * used interchangeably. Bound is the boundary of a single element, while
 * boundary is the boundary of the entire mesh that boundary conditions may be
 * applied to.
 */

namespace A2D {

class ElementTypes {
 public:
  /**
   * @brief Element types implemented by A2D
   */
  //   enum ElementReferenceDomain {
  //     NODE,
  //     LINE,
  //     TRIANGLE,
  //     QUADRILATERAL,
  //     TETRAHEDRAL,
  //     HEXAHEDRAL,
  //     WEDGE,
  //     PYRAMID
  //   };

  /**
   * @brief Elements implemented by A2D
   */
  enum class Element { Tri, Quad, Tet, Hex, Wedge, Pyrmd };

  /**
   * @brief Element topological entities for ordering
   */
  enum ElementEntity { Domain, Bound, Edge, Vertex };

  // Maximum number of vertices of a bound, used as the second array dimension
  static constexpr index_t MAX_BOUND_VERTS = 4;

  // Maximum number of edges of a bound, used as the second array dimension
  static constexpr index_t MAX_BOUND_EDGES = 4;

  /**
   * @brief Line element
   *
   * Line is a 1D element, its entities include:
   * - domain
   * - bound
   *
   * The bounds of a line is
   *
   *  0 ------------- 1
   */

  // Get the degrees of freedom associated with the bound
  template <index_t offset, index_t ndof, index_t nx, class ElemDof,
            class EntityDof>
  KOKKOS_FUNCTION static void get_line_bound_dof(index_t b,
                                                 const ElemDof& element,
                                                 EntityDof& entity);

  // Get the degrees of freedom associated with a bound
  template <index_t offset, index_t ndof, index_t nx, class EntityDof,
            class ElemDof>
  KOKKOS_FUNCTION static void set_line_bound_dof(index_t b,
                                                 const index_t orient,
                                                 const EntityDof& entity,
                                                 ElemDof& element);

  // Get the degrees of freedom from the domain
  template <index_t offset, bool ends, index_t ndof, index_t nx, class ElemDof,
            class EntityDof>
  KOKKOS_FUNCTION static void get_line_domain_dof(const ElemDof& element,
                                                  EntityDof& entity);

  // Set the degrees of freedom from the domain
  template <index_t offset, bool ends, index_t ndof, index_t nx,
            class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_line_domain_dof(const index_t orient,
                                                  const EntityDof& entity,
                                                  ElemDof& element);

  /**
   * @brief Triangle element
   *
   * Triangle is a 2D element, its entities include:
   * - domain
   * - bound
   * - vertex
   *
   * The vertices of the triangle are
   *
   *     2
   *     | .
   *     |    .
   * (1) |       . (0)
   *     |          .
   *     |             .
   *     0 ------------- 1
   *            (2)
   *
   * The bounds of the triangle are
   * Idx    Bound
   * (0)    1 -> 2
   * (1)    2 -> 0
   * (2)    0 -> 1
   */
  static const index_t TRI_NBOUNDS = 3;
  static const index_t TRI_NVERTS = 3;

  // Given bound index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_tri_bound_nverts() {
    static constexpr index_t TRI_BOUND_NVERTS[] = {2, 2, 2};
    return TRI_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_tri_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t TRI_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {1, 2, NO_INDEX, NO_INDEX},
        {2, 0, NO_INDEX, NO_INDEX},
        {0, 1, NO_INDEX, NO_INDEX}};
    return TRI_BOUND_VERTS;
  }

  /**
   * @brief Quadrilateral element
   *
   * Quadrilateral is a 2D element, its entities include:
   * - domain
   * - bound
   * - vertex
   *
   * The vertices of the quadrilateral element are
   *
   *           (1)
   *      3 ----------- 2
   *      |             |
   *      |             | (3)
   *  (2) |             |
   *      |             |
   *      0 ----------- 1
   *            (0)
   *
   * The bounds are
   * Idx    Bound
   * (0)    0 -> 1
   * (1)    3 -> 2
   * (2)    0 -> 3
   * (3)    1 -> 2
   */
  static const index_t QUAD_NBOUNDS = 4;
  static const index_t QUAD_NVERTS = 4;

  // Given bound index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_quad_bound_nverts() {
    static constexpr index_t QUAD_BOUND_NVERTS[] = {2, 2, 2, 2};
    return QUAD_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_quad_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t QUAD_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {0, 1, NO_INDEX, NO_INDEX},
        {3, 2, NO_INDEX, NO_INDEX},
        {0, 3, NO_INDEX, NO_INDEX},
        {1, 2, NO_INDEX, NO_INDEX}};
    return QUAD_BOUND_VERTS;
  }

  // Cartesian coordinates of the vertices in the reference element
  KOKKOS_FUNCTION static const index_t (*get_quad_verts_cart())[2] {
    static constexpr index_t QUAD_VERTS_CART[][2] = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1}};
    return QUAD_VERTS_CART;
  }

  static const index_t NUM_QUAD_DOMAIN_ORIENTATIONS = 8;

  // Given a reference domain, find the orientation
  KOKKOS_FUNCTION static index_t get_quad_domain_orientation(
      const index_t ref_domain_verts[], const index_t domain_verts[]);

  // Get the coords on quad ref element object
  KOKKOS_FUNCTION static void get_coords_on_quad_ref_element(
      const index_t orient, const index_t hx, const index_t hy, const index_t x,
      const index_t y, index_t* u, index_t* v);

  KOKKOS_FUNCTION static index_t get_index_on_quad_ref_element(
      const index_t orient, const index_t hx, const index_t hy, const index_t x,
      const index_t y);

  // Get the local node index (without offset) for a node on an element
  // with nx and ny nodes along the local directions
  template <index_t nx, index_t ny>
  KOKKOS_FUNCTION static inline int get_quad_node(const int i, const int j);

  // Get the bound length between the vertices v1 and v2
  template <index_t nx, index_t ny>
  KOKKOS_FUNCTION static inline index_t get_quad_bound_length(const index_t v0,
                                                              const index_t v1);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, class ElemDof,
            class EntityDof>
  KOKKOS_FUNCTION static void get_quad_vert_dof(index_t v,
                                                const ElemDof& element,
                                                EntityDof& entity);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_quad_vert_dof(index_t v,
                                                const EntityDof& entity,
                                                ElemDof& element);

  // Get the degrees of freedom from the bound
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_quad_bound_dof(const index_t b,
                                                 const ElemDof& element,
                                                 EntityDof& entity);

  // Set the degrees of freedom from the bound
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_quad_bound_dof(const index_t b,
                                                 const index_t orient,
                                                 const EntityDof& entity,
                                                 ElemDof& element);

  // Get the degrees of freedom from the quad domain
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_quad_domain_dof(const ElemDof& element,
                                                  EntityDof& entity);

  // Set the degrees of freedom from the domain into the element
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_quad_domain_dof(const index_t orient,
                                                  const EntityDof& entity,
                                                  ElemDof& element);

  /**
   * @brief Tetrahedral properties
   *
   * Tetrahedral is a 3D element, its entities include:
   * - domain
   * - bound
   * - edge
   * - vertex
   *
   * The vertices of the tetrahedral element are
   *
   *       3
   *     / . \
   *    /  .  \
   *   /   .   \
   *  /    .    \
   * 0 ----.---- 2
   *  \    .    /
   *   \   .   /
   *    \  .  /
   *     \ . /
   *       1
   *
   * The edges of the tetrahedral are
   * Idx    Edge
   * (0)    0 -> 1
   * (1)    1 -> 2
   * (2)    2 -> 0
   * (3)    0 -> 3
   * (4)    1 -> 3
   * (5)    2 -> 3
   *
   * The bounds of the element are
   * Idx    Bound          Edges
   * (0)    1 -> 2 -> 3    1, 5, -4
   * (1)    0 -> 3 -> 2    3, -5, 2
   * (2)    0 -> 1 -> 3    0, 4, -3
   * (3)    0 -> 2 -> 1    -2, -1, -0
   */
  // Number of vertices, edges and bounds
  static const index_t TET_NBOUNDS = 4;
  static const index_t TET_NEDGES = 6;
  static const index_t TET_NVERTS = 4;

  // Given edge index, return edge vertex indices
  KOKKOS_FUNCTION static const index_t (*get_tet_edge_verts())[2] {
    static constexpr index_t TET_EDGE_VERTS[][2] = {{0, 1}, {1, 2}, {2, 0},
                                                    {0, 3}, {1, 3}, {2, 3}};
    return TET_EDGE_VERTS;
  }

  // Given bounds index, return edge indices
  KOKKOS_FUNCTION static const index_t* get_tet_bound_nedges() {
    static constexpr index_t TET_BOUND_NEDGES[] = {3, 3, 3, 3};
    return TET_BOUND_NEDGES;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_tet_bound_edges())[MAX_BOUND_EDGES] {
    static constexpr index_t TET_BOUND_EDGES[][MAX_BOUND_EDGES] = {
        {1, 5, 4, NO_INDEX},
        {3, 5, 2, NO_INDEX},
        {0, 4, 3, NO_INDEX},
        {2, 1, 0, NO_INDEX}};
    return TET_BOUND_EDGES;
  }

  // Given bounds index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_tet_bound_nverts() {
    static constexpr index_t TET_BOUND_NVERTS[] = {3, 3, 3, 3};
    return TET_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_tet_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t TET_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {1, 2, 3, NO_INDEX},
        {0, 3, 2, NO_INDEX},
        {0, 1, 3, NO_INDEX},
        {0, 2, 1, NO_INDEX}};
    return TET_BOUND_VERTS;
  }

  /**
   * @brief Hexahedral properties
   *
   * Hexahedral is a 3D element, its entities include:
   * - domain
   * - bound
   * - edge
   * - vertex
   *
   *        7 --------------- 6
   *       / |              / |
   *      /  |             /  |
   *     /   |            /   |
   *    4 -------------- 5    |
   *    |    |           |    |
   *    |    3 ----------|--- 2
   *    |   /            |   /
   *    |  /             |  /
   *    | /              | /
   *    0 -------------- 1
   *
   * The edges are
   * Idx    Edge
   * (0)    0 -> 1
   * (1)    3 -> 2
   * (2)    4 -> 5
   * (3)    7 -> 6
   * (4)    0 -> 3
   * (5)    1 -> 2
   * (6)    4 -> 7
   * (7)    5 -> 6
   * (8)    0 -> 4
   * (9)    1 -> 5
   * (10)   3 -> 7
   * (11)   2 -> 6
   *
   * The bounds are
   * Idx    Bound               Edges
   * (0)    0 -> 4 -> 7 -> 3    8, 6, -10, -4
   * (1)    1 -> 2 -> 6 -> 5    5, 11, -7, -9
   * (2)    0 -> 1 -> 5 -> 4    0, 9, -2, -8
   * (3)    3 -> 7 -> 6 -> 2    10, 3, -11, -1
   * (4)    0 -> 3 -> 2 -> 1    4, 1, -5, -0
   * (5)    4 -> 5 -> 6 -> 7    2, 7, -3, -6
   */
  // Number of vertices, edges and bounds
  static const index_t HEX_NBOUNDS = 6;
  static const index_t HEX_NEDGES = 12;
  static const index_t HEX_NVERTS = 8;

  // Given edge index, return edge vertex indices
  KOKKOS_FUNCTION static const index_t (*get_hex_edge_verts())[2] {
    static constexpr index_t HEX_EDGE_VERTS[][2] = {
        {0, 1}, {3, 2}, {4, 5}, {7, 6}, {0, 3}, {1, 2},
        {4, 7}, {5, 6}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
    return HEX_EDGE_VERTS;
  }

  // Given bound index, return edge indices
  KOKKOS_FUNCTION static const index_t* get_hex_bound_nedges() {
    static constexpr index_t HEX_BOUND_NEDGES[] = {4, 4, 4, 4, 4, 4};
    return HEX_BOUND_NEDGES;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_hex_bound_edges())[MAX_BOUND_EDGES] {
    static constexpr index_t HEX_BOUND_EDGES[][MAX_BOUND_EDGES] = {
        {8, 10, 4, 6},  {5, 7, 9, 11}, {0, 2, 8, 9},
        {1, 3, 11, 10}, {4, 5, 0, 1},  {2, 3, 6, 7}};
    return HEX_BOUND_EDGES;
  }

  // Given bound index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_hex_bound_nverts() {
    static constexpr index_t HEX_BOUND_NVERTS[] = {4, 4, 4, 4, 4, 4};
    return HEX_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_hex_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t HEX_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
        {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};
    return HEX_BOUND_VERTS;
  }

  // Cartesian coordinates of the vertices in the reference element
  KOKKOS_FUNCTION static const index_t (*get_hex_verts_cart())[3] {
    static constexpr index_t HEX_VERTS_CART[][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    return HEX_VERTS_CART;
  }

  // Get the local node index (without offset) for a node on an element
  // with nx, ny and nz nodes along the local directions
  template <index_t nx, index_t ny, index_t nz>
  KOKKOS_FUNCTION static int get_hex_node(const int i, const int j,
                                          const int k);

  // Get the edge length between the verties v1 and v2
  template <index_t nx, index_t ny, index_t nz>
  KOKKOS_FUNCTION static index_t get_hex_edge_length(const index_t v0,
                                                     const index_t v1);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
            class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_hex_vert_dof(index_t v,
                                               const ElemDof& element,
                                               EntityDof& entity);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
            class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_hex_vert_dof(index_t v,
                                               const EntityDof& entity,
                                               ElemDof& element);

  // Get the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_hex_edge_dof(const index_t e,
                                               const ElemDof& element,
                                               EntityDof& entity);

  // Set the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_hex_edge_dof(const index_t e,
                                               const index_t orient,
                                               const EntityDof& entity,
                                               ElemDof& element);

  // Get the degrees of freedom from the hex bound
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_hex_bound_dof(const index_t f,
                                                const ElemDof& element,
                                                EntityDof& entity);

  // Set the degrees of freedom from the bound into the element
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_hex_bound_dof(const index_t f,
                                                const index_t orient,
                                                const EntityDof& entity,
                                                ElemDof& element);

  // Get the degrees of freedom from the domain
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  KOKKOS_FUNCTION static void get_hex_domain_dof(const ElemDof& element,
                                                 EntityDof& entity);

  // Set the degrees of freedom into the element array
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  KOKKOS_FUNCTION static void set_hex_domain_dof(const EntityDof& entity,
                                                 ElemDof& element);

  /**
   * @brief Wedge properties
   *
   * Wedge is a 3D element, its entities include:
   * - domain
   * - bound
   * - edge
   * - vertex
   *
   *        2 ------------- 5
   *      / |             / |
   *     /  |            /  |
   *    /   |           /   |
   *   0 ---|----------3    |
   *    \   |           \   |
   *     \  |            \  |
   *      \ |             \ |
   *        1 ------------- 4
   *
   * The edges are
   * Idx    Edge
   * (0)    0 -> 1
   * (1)    1 -> 2
   * (2)    2 -> 0
   * (3)    3 -> 4
   * (4)    4 -> 5
   * (5)    5 -> 3
   * (6)    0 -> 3
   * (7)    1 -> 4
   * (8)    2 -> 5
   *
   * The bounds are
   * Idx    Bound               Edges
   * (0)    0 -> 1 -> 2         0, 1, 2
   * (1)    3 -> 4 -> 5         3, 4, 5
   * (2)    0 -> 3 -> 4 -> 1    6, 3, -7, -0
   * (3)    1 -> 4 -> 5 -> 2    7, 4, -8, -1
   * (4)    0 -> 2 -> 5 -> 3    -2, 8, 5, -6
   */
  static const index_t WEDGE_NBOUNDS = 5;
  static const index_t WEDGE_NEDGES = 9;
  static const index_t WEDGE_NVERTS = 6;

  // Given edge index, return edge vertex indices
  KOKKOS_FUNCTION static const index_t (*get_wedge_edge_verts())[2] {
    static constexpr index_t WEDGE_EDGE_VERTS[][2] = {
        {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    return WEDGE_EDGE_VERTS;
  }

  // Given bound index, return edge indices
  KOKKOS_FUNCTION static const index_t* get_wedge_bound_nedges() {
    static constexpr index_t WEDGE_BOUND_NEDGES[] = {3, 3, 4, 4, 4};
    return WEDGE_BOUND_NEDGES;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_wedge_bound_edges())[MAX_BOUND_EDGES] {
    static constexpr index_t WEDGE_BOUND_EDGES[][MAX_BOUND_EDGES] = {
        {0, 1, 2, NO_INDEX},
        {3, 4, 5, NO_INDEX},
        {6, 3, 7, 0},
        {7, 4, 8, 1},
        {2, 8, 5, 6}};
    return WEDGE_BOUND_EDGES;
  }

  // Given bound index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_wedge_bound_nverts() {
    static constexpr index_t WEDGE_BOUND_NVERTS[] = {3, 3, 4, 4, 4};
    return WEDGE_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_wedge_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t WEDGE_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {0, 1, 2, NO_INDEX},
        {3, 4, 5, NO_INDEX},
        {0, 3, 4, 1},
        {1, 4, 5, 2},
        {0, 2, 5, 3}};
    return WEDGE_BOUND_VERTS;
  }

  /**
   * @brief Pyramid properties
   *
   * Pyramid is a 3D element, its entities include:
   * - domain
   * - bound
   * - edge
   * - vertex
   *
   *           4
   *        . / \ .
   *     .   /   \   .
   *   0 ---/-----\--- 3
   *   |   /       \   |
   *   |  /         \  |
   *   | /           \ |
   *   |/             \|
   *   1 --------------2
   *
   * The edges are
   * Idx    Edge
   * (0)    0 -> 1
   * (1)    1 -> 2
   * (2)    2 -> 3
   * (3)    3 -> 0
   * (4)    0 -> 4
   * (5)    1 -> 4
   * (6)    2 -> 4
   * (7)    3 -> 4
   *
   * The bounds are
   * Idx    Bound               Edges
   * (0)    0 -> 1 -> 4         0, 5, -4
   * (1)    1 -> 2 -> 4         1, 6, -5
   * (2)    2 -> 3 -> 4         2, 7, -6
   * (3)    0 -> 4 -> 3         4, -7, 3
   * (4)    0 -> 3 -> 2 -> 1    -3, -2, -1, -0
   */
  static const index_t PYRMD_NBOUNDS = 5;
  static const index_t PYRMD_NEDGES = 8;
  static const index_t PYRMD_NVERTS = 5;

  // Given edge index, return edge vertex indices
  KOKKOS_FUNCTION static const index_t (*get_pyrmd_edge_verts())[2] {
    static constexpr index_t PYRMD_EDGE_VERTS[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    return PYRMD_EDGE_VERTS;
  }

  // Given bound index, return edge indices
  KOKKOS_FUNCTION static const index_t* get_pyrmd_bound_nedges() {
    static constexpr index_t PYRMD_BOUND_NEDGES[] = {3, 3, 3, 3, 4};
    return PYRMD_BOUND_NEDGES;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_pyrmd_bound_edges())[MAX_BOUND_EDGES] {
    static constexpr index_t PYRMD_BOUND_EDGES[][MAX_BOUND_EDGES] = {
        {0, 5, 4, NO_INDEX},
        {1, 6, 5, NO_INDEX},
        {2, 7, 6, NO_INDEX},
        {4, 7, 3, NO_INDEX},
        {3, 2, 1, 0}};
    return PYRMD_BOUND_EDGES;
  }

  // Given bound index, return number of vertices/vertex indices
  KOKKOS_FUNCTION static const index_t* get_pyrmd_bound_nverts() {
    static constexpr index_t PYRMD_BOUND_NVERTS[] = {3, 3, 3, 3, 4};
    return PYRMD_BOUND_NVERTS;
  }
  KOKKOS_FUNCTION static const index_t (
      *get_pyrmd_bound_verts())[MAX_BOUND_VERTS] {
    static constexpr index_t PYRMD_BOUND_VERTS[][MAX_BOUND_VERTS] = {
        {0, 1, 4, NO_INDEX},
        {1, 2, 4, NO_INDEX},
        {2, 3, 4, NO_INDEX},
        {0, 4, 3, NO_INDEX},
        {0, 3, 2, 1}};
    return PYRMD_BOUND_VERTS;
  }
};  // namespace A2D

}  // namespace A2D

#include "feelementtypes-inl.h"

#endif  //  A2D_FE_ELEMENT_TYPES_H