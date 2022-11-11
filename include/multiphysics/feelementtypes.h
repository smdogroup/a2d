#ifndef A2D_FE_ELEMENT_TYPES_H
#define A2D_FE_ELEMENT_TYPES_H

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

  /**
   * @brief Triangle element
   *
   *  The vertices of the triangle are
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
   *  The edges of the triangle are
   *  (0)  1 -> 2
   *  (1)  2 -> 1
   *  (2)  0 -> 1
   */
  static const index_t TRI_VERTS = 3;
  static const index_t TRI_EDGES = 3;
  static constexpr index_t TRI_EDGE_VERTS[][2] = {{1, 2}, {2, 0}, {0, 1}};

  /**
   * @brief Quadrilateral element
   *
   * The vertices of the quadrilateral element are
   *
   *           (2)
   *      3 ----------- 2
   *      |             |
   *      |             | (1)
   *  (3) |             |
   *      |             |
   *      0 ----------- 1
   *            (0)
   *
   *  The edges of the quadrilateral are
   *  (0)  0 -> 1
   *  (1)  1 -> 2
   *  (2)  2 -> 3
   *  (3)  3 -> 0
   */
  static const index_t QUAD_VERTS = 4;
  static const index_t QUAD_EDGES = 4;
  static constexpr index_t QUAD_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}};

  static constexpr index_t QUAD_VERTS_CART[][2] = {
      {0, 0}, {1, 0}, {1, 1}, {1, 0}};

  static const index_t NUM_QUAD_FACE_ORIENTATIONS = 8;
  // static constexpr index_t QUAD_FACE_ORIENTATIONS[][4] = {
  //     {0, 1, 2, 3}, {2, 0, 3, 1}, {3, 2, 1, 0}, {1, 3, 0, 2},
  //     {0, 2, 1, 3}, {2, 3, 0, 1}, {3, 1, 2, 0}, {1, 0, 3, 2}};
  static constexpr index_t QUAD_FACE_ORIENTATIONS[][4] = {
      {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
      {0, 3, 2, 1}, {3, 2, 1, 0}, {2, 1, 0, 3}, {1, 0, 3, 2}};

  /**
   * @brief Given a reference face, find the orientation
   *
   * @param ref
   * @param face
   * @return index_t
   */
  static index_t get_quad_face_orientation(const index_t ref[],
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
   * @param h The edge length
   * @param x The 1st coordinate on the transformed face
   * @param y The 2nd coordinate on the transformed face
   * @param u The 1st coordinate on the reference face
   * @param v The 2nd coordinate on the reference face
   */
  static void get_coords_on_quad_ref_element(const index_t orient,
                                             const index_t h, const index_t x,
                                             const index_t y, index_t *u,
                                             index_t *v) {
    if (orient == 0) {
      *u = x;
      *v = y;
    } else if (orient == 1) {
      *u = h - y;
      *v = x;
    } else if (orient == 2) {
      *u = h - x;
      *v = h - y;
    } else if (orient == 3) {
      *u = y;
      *v = h - x;
    } else if (orient == 4) {
      *u = y;
      *v = x;
    } else if (orient == 5) {
      *u = x;
      *v = h - y;
    } else if (orient == 6) {
      *u = h - y;
      *v = h - x;
    } else if (orient == 7) {
      *u = h - x;
      *v = y;
    } else {
      *u = 0;
      *v = 0;
    }
  }

  /**
   * @brief Tetrahedral properties
   * The vertices of the tetrahedral element are
   *
   *     3
   *    /. \
   *   / .  \
   *  /  .   \
   * 0 --.---- 2
   *  \  .   /
   *   \ . /
   *     1
   *
   * The edges of the tetrahedral are
   * 0 -> 1
   * 1 -> 2
   * 2 -> 0
   * 0 -> 3
   * 1 -> 3
   * 2 -> 3
   *
   * The faces of the element are
   * 1 -> 2 -> 3
   * 0 -> 3 -> 2
   * 0 -> 1 -> 3
   * 0 -> 2 -> 1
   */
  static const index_t TET_VERTS = 4;
  static const index_t TET_EDGES = 6;
  static const index_t TET_FACES = 4;
  static constexpr index_t TET_EDGE_VERTS[][2] = {{0, 1}, {1, 2}, {2, 0},
                                                  {0, 3}, {1, 3}, {2, 3}};
  static constexpr index_t TET_FACE_VERTS[][3] = {
      {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

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
   */
  static const index_t HEX_VERTS = 8;
  static const index_t HEX_EDGES = 12;
  static const index_t HEX_FACES = 6;
  static constexpr index_t HEX_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
      {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  static constexpr index_t HEX_FACE_VERTS[][4] = {{0, 4, 7, 3}, {1, 2, 6, 5},
                                                  {0, 1, 5, 4}, {2, 3, 7, 6},
                                                  {0, 3, 2, 1}, {4, 5, 6, 7}};

  static constexpr index_t HEX_VERTS_CART[][3] = {
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

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
   */
  static const index_t WEDGE_VERTS = 6;
  static const index_t WEDGE_EDGES = 9;
  static const index_t WEDGE_FACES = 5;
  static constexpr index_t WEDGE_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};

  static const index_t WEDGE_TRI_FACES = 2;
  static constexpr index_t WEDGE_TRI_FACE_VERTS[][3] = {{0, 1, 2}, {3, 4, 5}};
  static const index_t WEDGE_QUAD_FACES = 2;
  static constexpr index_t WEDGE_QUAD_FACE_VERTS[][4] = {
      {0, 3, 4, 1}, {1, 4, 5, 2}, {0, 2, 5, 3}};

  static constexpr index_t WEDGE_EDGE_TO_ADJ_FACES[][2] = {
      {2, 4}, {1, 4}, {3, 4}, {0, 4}, {2, 5}, {1, 5},
      {3, 5}, {0, 5}, {0, 2}, {1, 2}, {1, 3}, {0, 3}};
  static constexpr index_t WEDGE_EDGE_TO_ADJ_FACE_EDGE[][2] = {
      {0, 3}, {0, 2}, {0, 1}, {3, 0}, {2, 0}, {2, 1},
      {2, 2}, {1, 3}, {0, 3}, {3, 1}, {1, 3}, {2, 1}};

  // static constexpr index_t WEDGE_EDGE_TO_ADJ_FACES[][2] = {
  //     {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 4}, {2, 3}, {3,
  //     4}};
  // static constexpr index_t WEDGE_EDGE_TO_ADJ_FACE_EDGE[][2] = {
  //     {0, 3}, {1, 3}, {2, 0}, {0, 1}, {1, 1}, {2, 2}, {0, 3}, {2, 0}, {2,
  //     1}};

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
   */
  static const index_t PYRMD_VERTS = 5;
  static const index_t PYRMD_EDGES = 8;
  static const index_t PYRMD_FACES = 5;
  static constexpr index_t PYRMD_EDGE_VERTS[][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
  static const index_t PYRMD_TRI_FACES = 4;
  static constexpr index_t PYRMD_TRI_FACE_VERTS[][3] = {
      {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {0, 4, 3}};
  static const index_t PYRMD_QUAD_FACES = 4;
  static constexpr index_t PYRMD_QUAD_FACE_VERTS[] = {0, 3, 2, 1};

  static constexpr index_t PYRMD_EDGE_TO_ADJ_FACES[][2] = {
      {0, 4}, {1, 4}, {2, 4}, {3, 4}, {0, 3}, {0, 1}, {1, 2}, {2, 3}};
  static constexpr index_t PYRMD_EDGE_TO_ADJ_FACE_EDGE[][2] = {
      {0, 3}, {0, 2}, {0, 1}, {2, 0}, {2, 0}, {1, 2}, {1, 2}, {1, 1}};
};

}  // namespace A2D

#endif  //  A2D_FE_ELEMENT_TYPES_H