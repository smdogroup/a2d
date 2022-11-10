#ifndef A2D_FE_ELEMENT_TYPES_H
#define A2D_FE_ELEMENT_TYPES_H

namespace A2D {

/**
 * @brief Connecitivity types for all the elements implemented by A2D
 *
 */
enum A2D_ELEMENT_TYPES {
  POINT,
  EDGE,
  TRIANGLE,
  QUADRILATERAL,
  TETRAHEDRAL,
  HEXAHEDRAL,
  WEDGE,
  PYRAMID
};

class ElementTypes {
 public:
  static const index_t MAX_ELEMENT_EDGES = 12;

  /**
   * @brief Triangle element
   *
   *  The vertices of the triangle are
   *
   *         2
   *        / \
   *       /   \
   * (1)  /     \ (0)
   *     /       \
   *    /         \
   *   0 --------- 1
   *        (2)
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

  static constexpr index_t TET_EDGE_TO_ADJ_FACES[][2] = {
      {2, 3}, {0, 3}, {1, 3}, {1, 2}, {0, 2}, {0, 1}};
  static constexpr index_t TET_EDGE_TO_ADJ_FACE_EDGE[][2] = {
      {0, 2}, {0, 1}, {2, 0}, {0, 2}, {2, 1}, {1, 1}};

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

  // Faces adjacent to the given edge index
  static constexpr index_t HEX_EDGE_TO_ADJ_FACES[][2] = {
      {2, 4}, {1, 4}, {3, 4}, {0, 4}, {2, 5}, {1, 5},
      {3, 5}, {0, 5}, {0, 2}, {1, 2}, {1, 3}, {0, 3}};

  // Index of the edge on each adjacent face
  static constexpr index_t HEX_EDGE_TO_ADJ_FACE_EDGE[][2] = {
      {0, 3}, {0, 2}, {0, 1}, {3, 0}, {2, 0}, {2, 1},
      {2, 2}, {1, 3}, {0, 3}, {3, 1}, {1, 3}, {2, 1}};

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
      {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 4}, {2, 3}, {3, 4}};
  static constexpr index_t WEDGE_EDGE_TO_ADJ_FACE_EDGE[][2] = {
      {0, 3}, {1, 3}, {2, 0}, {0, 1}, {1, 1}, {2, 2}, {0, 3}, {2, 0}, {2, 1}};

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