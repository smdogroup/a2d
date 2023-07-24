#ifndef A2D_FE_ELEMENT_TYPES_H
#define A2D_FE_ELEMENT_TYPES_H

#include "a2dobjs.h"

namespace A2D {

class ElementTypes {
 public:
  static const index_t MAX_ELEMENT_EDGES = 12;

  /**
   * @brief Element types implemented by A2D
   */
  enum ElementReferenceDomain {
    NODE,
    LINE,
    TRIANGLE,
    QUADRILATERAL,
    TETRAHEDRAL,
    HEXAHEDRAL,
    WEDGE,
    PYRAMID
  };

  /**
   * @brief Element entity for ordering
   */
  enum ElementEntity { VERTEX, EDGE, FACE, VOLUME };

  // Maximum number of vertices of a face, for alignment purpose
  static constexpr index_t MAX_FACE_VERTS = 4;

  // Maximum number of edges of a face, for alignment purpose
  static constexpr index_t MAX_FACE_EDGES = 4;

  // Maximum number of faces of an element
  static constexpr index_t MAX_FACES = 6;
  static constexpr index_t NONE_FACE_NQUANTS[MAX_FACES] = {0, 0, 0, 0, 0, 0};

  /**
   * @brief Triangle element
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
   * The edges of the triangle are
   * Idx    Edge
   * (0)    1 -> 2
   * (1)    2 -> 0
   * (2)    0 -> 1
   *
   * The faces are (note that for 2d elements face and edge are the same)
   * Idx    Face      Edges
   * (0)    1 -> 2    0,
   * (1)    2 -> 0    1,
   * (2)    0 -> 1    2,
   */
  static const index_t TRI_VERTS = 3;
  static const index_t TRI_EDGES = 3;
  static const index_t TRI_FACES = 3;

  // Given edge index, return edge vertex indices
  static constexpr index_t TRI_EDGE_VERTS[][2] = {{1, 2}, {2, 0}, {0, 1}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t TRI_FACE_NVERTS[] = {2, 2, 2};
  static constexpr index_t TRI_FACE_VERTS[][MAX_FACE_VERTS] = {
      {1, 2, NO_INDEX, NO_INDEX},
      {2, 0, NO_INDEX, NO_INDEX},
      {0, 1, NO_INDEX, NO_INDEX}};

  // Given face index, return edge indices
  static constexpr index_t TRI_FACE_NEDGES[] = {1, 1, 1};
  static constexpr index_t TRI_FACE_EDGES[][MAX_FACE_EDGES] = {
      {0, NO_INDEX, NO_INDEX, NO_INDEX},
      {1, NO_INDEX, NO_INDEX, NO_INDEX},
      {2, NO_INDEX, NO_INDEX, NO_INDEX}};

  /**
   * @brief Quadrilateral element
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
   * The edges of the quadrilateral are
   * Idx    Edge
   * (0)    0 -> 1
   * (1)    3 -> 2
   * (2)    0 -> 3
   * (3)    1 -> 2
   *
   * The faces are (note that for 2d elements face and edge are the same)
   * Idx    Face      Edges
   * (0)    0 -> 1    0,
   * (1)    3 -> 2    1,
   * (2)    0 -> 3    2,
   * (3)    1 -> 2    3,
   */
  static const index_t QUAD_VERTS = 4;
  static const index_t QUAD_EDGES = 4;
  static const index_t QUAD_FACES = 4;

  // Given edge index, return edge vertex indices
  static constexpr index_t QUAD_EDGE_VERTS[][2] = {
      {0, 1}, {3, 2}, {0, 3}, {1, 2}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t QUAD_FACE_NVERTS[] = {2, 2, 2, 2};
  static constexpr index_t QUAD_FACE_VERTS[][MAX_FACE_VERTS] = {
      {0, 1, NO_INDEX, NO_INDEX},
      {3, 2, NO_INDEX, NO_INDEX},
      {0, 3, NO_INDEX, NO_INDEX},
      {1, 2, NO_INDEX, NO_INDEX}};

  // Given face index, return edge indices
  static constexpr index_t QUAD_FACE_NEDGES[] = {1, 1, 1, 1};
  static constexpr index_t QUAD_FACE_EDGES[][MAX_FACE_EDGES] = {
      {0, NO_INDEX, NO_INDEX, NO_INDEX},
      {1, NO_INDEX, NO_INDEX, NO_INDEX},
      {2, NO_INDEX, NO_INDEX, NO_INDEX},
      {3, NO_INDEX, NO_INDEX, NO_INDEX}};

  // Cartesian coordinates of the vertices in the reference element
  static constexpr index_t QUAD_VERTS_CART[][2] = {
      {0, 0}, {1, 0}, {1, 1}, {0, 1}};

  static const index_t NUM_QUAD_FACE_ORIENTATIONS = 8;
  static constexpr index_t QUAD_FACE_ORIENTATIONS[][4] = {
      {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
      {0, 3, 2, 1}, {3, 2, 1, 0}, {2, 1, 0, 3}, {1, 0, 3, 2}};

  // Given a reference face, find the orientation
  static index_t get_quad_face_orientation(const index_t ref[],
                                           const index_t face[]);

  // Get the coords on quad ref element object
  static void get_coords_on_quad_ref_element(const index_t orient,
                                             const index_t hx, const index_t hy,
                                             const index_t x, const index_t y,
                                             index_t* u, index_t* v);

  static index_t get_index_on_quad_ref_element(const index_t orient,
                                               const index_t hx,
                                               const index_t hy,
                                               const index_t x,
                                               const index_t y);

  // Get the local node index (without offset) for a node on an element
  // with nx and ny nodes along the local directions
  template <index_t nx, index_t ny>
  static inline int get_quad_node(const int i, const int j);

  // Get the edge length between the verties v1 and v2
  template <index_t nx, index_t ny>
  static inline index_t get_quad_edge_length(const index_t v0,
                                             const index_t v1);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, class ElemDof,
            class EntityDof>
  static void get_quad_vert_dof(index_t v, const ElemDof& element,
                                EntityDof& entity);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  static void set_quad_vert_dof(index_t v, const EntityDof& entity,
                                ElemDof& element);

  // Get the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class ElemDof, class EntityDof>
  static void get_quad_edge_dof(const index_t e, const ElemDof& element,
                                EntityDof& entity);

  // Set the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  static void set_quad_edge_dof(const index_t e, const index_t orient,
                                const EntityDof& entity, ElemDof& element);

  // Get the degrees of freedom from the quad face
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  static void get_quad_face_dof(const ElemDof& element, EntityDof& entity);

  // Set the degrees of freeom from the face into the element
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            class EntityDof, class ElemDof>
  static void set_quad_face_dof(const index_t orient, const EntityDof& entity,
                                ElemDof& element);

  /**
   * @brief Tetrahedral properties
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
   * The faces of the element are
   * Idx    Face
   * (0)    1 -> 2 -> 3    1, 5, -4
   * (1)    0 -> 3 -> 2    3, -5, 2
   * (2)    0 -> 1 -> 3    0, 4, -3
   * (3)    0 -> 2 -> 1    -2, -1, -0
   */
  // Number of vertices, edges and faces
  static const index_t TET_VERTS = 4;
  static const index_t TET_EDGES = 6;
  static const index_t TET_FACES = 4;

  // Given edge index, return edge vertex indices
  static constexpr index_t TET_EDGE_VERTS[][2] = {{0, 1}, {1, 2}, {2, 0},
                                                  {0, 3}, {1, 3}, {2, 3}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t TET_FACE_NVERTS[] = {3, 3, 3, 3};
  static constexpr index_t TET_FACE_VERTS[][MAX_FACE_VERTS] = {
      {1, 2, 3, NO_INDEX},
      {0, 3, 2, NO_INDEX},
      {0, 1, 3, NO_INDEX},
      {0, 2, 1, NO_INDEX}};

  // Given face index, return edge indices
  static constexpr index_t TET_FACE_NEDGES[] = {3, 3, 3, 3};
  static constexpr index_t TET_FACE_EDGES[][MAX_FACE_EDGES] = {
      {1, 5, 4, NO_INDEX},
      {3, 5, 2, NO_INDEX},
      {0, 4, 3, NO_INDEX},
      {2, 1, 0, NO_INDEX}};

  /**
   * @brief Hexahedral properties
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
   * The faces are
   * Idx    Face                Edges
   * (0)    0 -> 4 -> 7 -> 3    8, 6, -10, -4
   * (1)    1 -> 2 -> 6 -> 5    5, 11, -7, -9
   * (2)    0 -> 1 -> 5 -> 4    0, 9, -2, -8
   * (3)    3 -> 7 -> 6 -> 2    10, 3, -11, -1
   * (4)    0 -> 3 -> 2 -> 1    4, 1, -5, -0
   * (5)    4 -> 5 -> 6 -> 7    2, 7, -3, -6
   */
  // Number of vertices, edges and faces
  static const index_t HEX_VERTS = 8;
  static const index_t HEX_EDGES = 12;
  static const index_t HEX_FACES = 6;

  // Given edge index, return edge vertex indices
  static constexpr index_t HEX_EDGE_VERTS[][2] = {
      {0, 1}, {3, 2}, {4, 5}, {7, 6}, {0, 3}, {1, 2},
      {4, 7}, {5, 6}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t HEX_FACE_NVERTS[] = {4, 4, 4, 4, 4, 4};
  static constexpr index_t HEX_FACE_VERTS[][MAX_FACE_VERTS] = {
      {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
      {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

  // Given face index, return edge indices
  static constexpr index_t HEX_FACE_NEDGES[] = {4, 4, 4, 4, 4, 4};
  static constexpr index_t HEX_FACE_EDGES[][MAX_FACE_EDGES] = {
      {8, 10, 4, 6},  {5, 7, 9, 11}, {0, 2, 8, 9},
      {1, 3, 11, 10}, {4, 5, 0, 1},  {2, 3, 6, 7}};

  // Cartesian coordinates of the vertices in the reference element
  static constexpr index_t HEX_VERTS_CART[][3] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

  // Get the local node index (without offset) for a node on an element
  // with nx, ny and nz nodes along the local directions
  template <index_t nx, index_t ny, index_t nz>
  static int get_hex_node(const int i, const int j, const int k);

  // Get the edge length between the verties v1 and v2
  template <index_t nx, index_t ny, index_t nz>
  static index_t get_hex_edge_length(const index_t v0, const index_t v1);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
            class ElemDof, class EntityDof>
  static void get_hex_vert_dof(index_t v, const ElemDof& element,
                               EntityDof& entity);

  // Get the degrees of freedom associated with the vertex
  template <index_t offset, index_t ndof, index_t nx, index_t ny, index_t nz,
            class EntityDof, class ElemDof>
  static void set_hex_vert_dof(index_t v, const EntityDof& entity,
                               ElemDof& element);

  // Get the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  static void get_hex_edge_dof(const index_t e, const ElemDof& element,
                               EntityDof& entity);

  // Set the degrees of freedom from the edge
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  static void set_hex_edge_dof(const index_t e, const index_t orient,
                               const EntityDof& entity, ElemDof& element);

  // Get the degrees of freedom from the hex face
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  static void get_hex_face_dof(const index_t f, const ElemDof& element,
                               EntityDof& entity);

  // Set the degrees of freeom from the face into the element
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  static void set_hex_face_dof(const index_t f, const index_t orient,
                               const EntityDof& entity, ElemDof& element);

  // Get the degrees of freeom from the volume
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class ElemDof, class EntityDof>
  static void get_hex_volume_dof(const ElemDof& element, EntityDof& entity);

  // Set the degrees of freedom into the element array
  template <index_t offset, bool ends, index_t ndof, index_t nx, index_t ny,
            index_t nz, class EntityDof, class ElemDof>
  static void set_hex_volume_dof(const EntityDof& entity, ElemDof& element);

  /**
   * @brief Wedge properties
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
   * The faces are
   * Idx    Face                Edges
   * (0)    0 -> 1 -> 2         0, 1, 2
   * (1)    3 -> 4 -> 5         3, 4, 5
   * (2)    0 -> 3 -> 4 -> 1    6, 3, -7, -0
   * (3)    1 -> 4 -> 5 -> 2    7, 4, -8, -1
   * (4)    0 -> 2 -> 5 -> 3    -2, 8, 5, -6
   */
  static const index_t WEDGE_VERTS = 6;
  static const index_t WEDGE_EDGES = 9;
  static const index_t WEDGE_FACES = 5;

  // Given edge index, return edge vertex indices
  static constexpr index_t WEDGE_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t WEDGE_FACE_NVERTS[] = {3, 3, 4, 4, 4};
  static constexpr index_t WEDGE_FACE_VERTS[][MAX_FACE_VERTS] = {
      {0, 1, 2, NO_INDEX},
      {3, 4, 5, NO_INDEX},
      {0, 3, 4, 1},
      {1, 4, 5, 2},
      {0, 2, 5, 3}};

  // Given face index, return edge indices
  static constexpr index_t WEDGE_FACE_NEDGES[] = {3, 3, 4, 4, 4};
  static constexpr index_t WEDGE_FACE_EDGES[][MAX_FACE_EDGES] = {
      {0, 1, 2, NO_INDEX},
      {3, 4, 5, NO_INDEX},
      {6, 3, 7, 0},
      {7, 4, 8, 1},
      {2, 8, 5, 6}};

  /**
   * @brief Pyramid properties
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
   * The faces are
   * Idx    Face                Edges
   * (0)    0 -> 1 -> 4         0, 5, -4
   * (1)    1 -> 2 -> 4         1, 6, -5
   * (2)    2 -> 3 -> 4         2, 7, -6
   * (3)    0 -> 4 -> 3         4, -7, 3
   * (4)    0 -> 3 -> 2 -> 1    -3, -2, -1, -0
   */
  static const index_t PYRMD_VERTS = 5;
  static const index_t PYRMD_EDGES = 8;
  static const index_t PYRMD_FACES = 5;

  // Given edge index, return edge vertex indices
  static constexpr index_t PYRMD_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};

  // Given face index, return number of vertices/vertex indices
  static constexpr index_t PYRMD_FACE_NVERTS[] = {3, 3, 3, 3, 4};
  static constexpr index_t PYRMD_FACE_VERTS[][MAX_FACE_VERTS] = {
      {0, 1, 4, NO_INDEX},
      {1, 2, 4, NO_INDEX},
      {2, 3, 4, NO_INDEX},
      {0, 4, 3, NO_INDEX},
      {0, 3, 2, 1}};

  // Given face index, return edge indices
  static constexpr index_t PYRMD_FACE_NEDGES[] = {3, 3, 3, 3, 4};
  static constexpr index_t PYRMD_FACE_EDGES[][MAX_FACE_EDGES] = {
      {0, 5, 4, NO_INDEX},
      {1, 6, 5, NO_INDEX},
      {2, 7, 6, NO_INDEX},
      {4, 7, 3, NO_INDEX},
      {3, 2, 1, 0}};
};

}  // namespace A2D

#include "feelementtypes-inl.h"

#endif  //  A2D_FE_ELEMENT_TYPES_H