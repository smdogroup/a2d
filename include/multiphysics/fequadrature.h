#ifndef A2D_FE_QUADRATURE_H
#define A2D_FE_QUADRATURE_H

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

/*
  Quadrature class for quadrilaterals
*/
const double GaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double GaussQuadWts2[] = {1.0, 1.0};

class Quad4ptQuadrature {
 public:
  static index_t get_num_points() { return 4; };
  static void get_point(const index_t index, double pt[]) {
    pt[0] = GaussQuadPts2[index % 2];
    pt[1] = GaussQuadPts2[(index % 4) / 2];
  }

  static double get_weight(const index_t index) {
    return (GaussQuadWts2[index % 2] * GaussQuadWts2[(index % 4) / 2]);
  }
};

}  // namespace A2D

#endif  // A2D_FE_QUADRATURE_H