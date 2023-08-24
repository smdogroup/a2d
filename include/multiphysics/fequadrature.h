#ifndef A2D_FE_QUADRATURE_H
#define A2D_FE_QUADRATURE_H

#include "multiphysics/lagrange_tools.h"

namespace A2D {

/*
  Quadrature class for triangles
*/
const double TriangleWts3[] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
const double TrianglePts3[] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0,
                               1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

class TriQuadrature3 {
 public:
  static index_t get_num_points() { return 3; }
  static void get_point(const index_t n, double pt[]) {
    pt[0] = TrianglePts3[2 * n];
    pt[1] = TrianglePts3[2 * n + 1];
  }
  static double get_weight(const index_t n) { return TriangleWts3[n]; }
};

template <index_t order>
class LineGaussQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = false;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order;

  static index_t get_num_points() { return order; }

  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n];
  }
  static double get_weight(const index_t n) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Wts[n];
  }
};

template <index_t order>
class QuadGaussQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = true;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order * order;

  /// @brief Number of points along each direction
  static const index_t tensor_dim0 = order;
  static const index_t tensor_dim1 = order;

  /**
   * @brief Get the quadrature point along the given dimension
   *
   * @param dim Dimension 0 or 1
   * @param pt The point along the direction
   * @return The quadrature point along the direction
   */
  static double get_tensor_point(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Pts[pt];
  }

  /**
   * @brief Get the quadrature weight along the given dimension
   *
   * @param dim Dimension 0, 1 or 2
   * @param pt The point along the direction
   * @return The weight factor along the direction
   */
  static double get_tensor_weight(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Wts[pt];
  }

  /**
   * @brief Get the quadrature point index
   *
   * @param q0 Quadrature index along 0-direction
   * @param q1 Quadrature index along 1-direction
   * @return The quadrature index
   */
  static index_t get_tensor_index(const index_t q0, const index_t q1) {
    return q0 + order * q1;
  }

  static index_t get_num_points() { return order * order; }
  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n % order];
    pt[1] = QuadratureData::Pts[n / order];
  }
  static double get_weight(const index_t n) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Wts[n % order] * QuadratureData::Wts[n / order];
  }
};

template <index_t order>
class HexGaussQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = true;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order * order * order;

  /// @brief Number of points along each direction
  static const index_t tensor_dim0 = order;
  static const index_t tensor_dim1 = order;
  static const index_t tensor_dim2 = order;

  /**
   * @brief Get the quadrature point along the given dimension
   *
   * @param dim Dimension 0, 1 or 2
   * @param pt The point along the direction
   * @return The quadrature point along the direction
   */
  static double get_tensor_point(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Pts[pt];
  }

  /**
   * @brief Get the quadrature weight along the given dimension
   *
   * @param dim Dimension 0, 1 or 2
   * @param pt The point along the direction
   * @return The weight factor along the direction
   */
  static double get_tensor_weight(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Wts[pt];
  }

  /**
   * @brief Get the quadrature point index
   *
   * @param q0 Quadrature index along 0-direction
   * @param q1 Quadrature index along 1-direction
   * @param q2 Quadrature index along 2-direction
   * @return The quadrature index
   */
  static index_t get_tensor_index(const index_t q0, const index_t q1,
                                  const index_t q2) {
    return q0 + order * (q1 + order * q2);
  }

  static index_t get_num_points() { return order * order * order; }
  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n % order];
    pt[1] = QuadratureData::Pts[(n % (order * order)) / order];
    pt[2] = QuadratureData::Pts[n / (order * order)];
  }
  static double get_weight(const index_t n) {
    using QuadratureData = get_interpolation<order, GAUSS_INTERPOLATION>;
    return QuadratureData::Wts[n % order] *
           QuadratureData::Wts[(n % (order * order)) / order] *
           QuadratureData::Wts[n / (order * order)];
  }
};

template <index_t order>
class LineGaussLobattoQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = false;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order;

  static index_t get_num_points() { return order; }
  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GLL_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n];
  }
  static double get_weight(const index_t n) {
    return GaussLobattoQuadData<order>::Wts[n];
  }
};

template <index_t order>
class QuadGaussLobattoQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = true;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order * order;

  /// @brief Number of points along each direction
  static const index_t tensor_dim0 = order;
  static const index_t tensor_dim1 = order;

  /**
   * @brief Get the quadrature point along the given dimension
   *
   * @param dim Dimension 0 or 1
   * @param pt The point along the direction
   * @return The quadrature point along the direction
   */
  static double get_tensor_point(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GLL_INTERPOLATION>;
    return QuadratureData::Pts[pt];
  }

  /**
   * @brief Get the quadrature weight along the given dimension
   *
   * @param dim Dimension 0 or 1
   * @param pt The point along the direction
   * @return The weight factor along the direction
   */
  static double get_tensor_weight(const index_t dim, const index_t pt) {
    return GaussLobattoQuadData<order>::Wts[pt];
  }

  /**
   * @brief Get the quadrature point index
   *
   * @param q0 Quadrature index along 0-direction
   * @param q1 Quadrature index along 1-direction
   * @return The quadrature index
   */
  static index_t get_tensor_index(const index_t q0, const index_t q1) {
    return q0 + order * q1;
  }

  static index_t get_num_points() { return order * order; }
  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GLL_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n % order];
    pt[1] = QuadratureData::Pts[n / order];
  }
  static double get_weight(const index_t n) {
    return GaussLobattoQuadData<order>::Wts[n % order] *
           GaussLobattoQuadData<order>::Wts[n / order];
  }
};

template <index_t order>
class HexGaussLobattoQuadrature {
 public:
  /// @brief  Is this a tensor product implementation
  static const bool is_tensor_product = true;

  /// @brief The total number of quadrature points
  static const index_t num_quad_points = order * order * order;

  /// @brief Number of points along each direction
  static const index_t tensor_dim0 = order;
  static const index_t tensor_dim1 = order;
  static const index_t tensor_dim2 = order;

  /**
   * @brief Get the quadrature point along the given dimension
   *
   * @param dim Dimension 0, 1 or 2
   * @param pt The point along the direction
   * @return The quadrature point along the direction
   */
  static double get_tensor_point(const index_t dim, const index_t pt) {
    using QuadratureData = get_interpolation<order, GLL_INTERPOLATION>;
    return QuadratureData::Pts[pt];
  }

  /**
   * @brief Get the quadrature weight along the given dimension
   *
   * @param dim Dimension 0, 1 or 2
   * @param pt The point along the direction
   * @return The weight factor along the direction
   */
  static double get_tensor_weight(const index_t dim, const index_t pt) {
    return GaussLobattoQuadData<order>::Wts[pt];
  }

  /**
   * @brief Get the quadrature point index
   *
   * @param q0 Quadrature index along 0-direction
   * @param q1 Quadrature index along 1-direction
   * @param q2 Quadrature index along 2-direction
   * @return The quadrature index
   */
  static index_t get_tensor_index(const index_t q0, const index_t q1,
                                  const index_t q2) {
    return q0 + order * (q1 + order * q2);
  }

  static index_t get_num_points() { return order * order * order; }
  static void get_point(const index_t n, double pt[]) {
    using QuadratureData = get_interpolation<order, GLL_INTERPOLATION>;
    pt[0] = QuadratureData::Pts[n % order];
    pt[1] = QuadratureData::Pts[(n % (order * order)) / order];
    pt[2] = QuadratureData::Pts[n / (order * order)];
  }
  static double get_weight(const index_t n) {
    return GaussLobattoQuadData<order>::Wts[n % order] *
           GaussLobattoQuadData<order>::Wts[(n % (order * order)) / order] *
           GaussLobattoQuadData<order>::Wts[n / (order * order)];
    return 0.0;
  }
};

}  // namespace A2D

#endif  // A2D_FE_QUADRATURE_H