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
class HexQuadrature {
 public:
  static index_t get_num_points() { return order * order * order; }
  static void get_point(const index_t n, double pt[]) {
    constexpr const double* pts = get_gauss_quadrature_pts<order>();
    pt[0] = pts[n % order];
    pt[1] = pts[(n % order * order) / order];
    pt[2] = pts[n / (order * order)];
  }
  static double get_weight(const index_t n) {
    constexpr const double* wts = get_gauss_quadrature_wts<order>();
    return wts[n % order] * wts[(n % (order * order)) / order] *
           wts[n / (order * order)];
  }
};

}  // namespace A2D

#endif  // A2D_FE_QUADRATURE_H