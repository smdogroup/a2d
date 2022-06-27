#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <cstddef>

#include "a2dobjs.h"
#include "block_numeric.h"
#include "parallel.h"

namespace A2D {

const double GaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double GaussQuadWts2[] = {1.0, 1.0};

class Quad4ptQuadrature {
 public:
  static const index_t NUM_QUAD_PTS = 4;

  static void getQuadPoint(const index_t index, double pt[]) {
    pt[0] = GaussQuadPts2[index % 2];
    pt[1] = GaussQuadPts2[(index % 4) / 2];
  }

  static double getQuadWeight(const index_t index) {
    return (GaussQuadWts2[index % 2] * GaussQuadWts2[(index % 4) / 2]);
  }
};

class QuadBiLinear {
 public:
  static const index_t NUM_NODES = 4;

  static void evalBasis(const double pt[], double N[]) {
    N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
    N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
    N[2] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);
    N[3] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  }

  static void evalBasisDeriv(const double pt[], double Nxy[]) {
    // Nx
    Nxy[0] = -0.25 * (1.0 - pt[1]);
    Nxy[1] = 0.25 * (1.0 - pt[1]);
    Nxy[2] = 0.25 * (1.0 + pt[1]);
    Nxy[3] = -0.25 * (1.0 + pt[1]);

    // Ny
    Nxy[4] = -0.25 * (1.0 - pt[0]);
    Nxy[5] = -0.25 * (1.0 + pt[0]);
    Nxy[6] = 0.25 * (1.0 + pt[0]);
    Nxy[7] = 0.25 * (1.0 - pt[0]);
  }
};

class Hex8ptQuadrature {
 public:
  static const index_t NUM_QUAD_PTS = 8;

  static void getQuadPoint(const index_t index, double pt[]) {
    pt[0] = GaussQuadPts2[index % 2];
    pt[1] = GaussQuadPts2[(index % 4) / 2];
    pt[2] = GaussQuadPts2[index / 4];
  }

  static double getQuadWeight(const index_t index) {
    return (GaussQuadWts2[index % 2] * GaussQuadWts2[(index % 4) / 2] *
            GaussQuadWts2[index / 4]);
  }
};

class HexTriLinear {
 public:
  static const index_t NUM_NODES = 8;

  static void evalBasis(const double pt[], double N[]) {
    N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
    N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
    N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
    N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
    N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
    N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
    N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
    N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
  }

  static void evalBasisDeriv(const double pt[], double Nx[], double Ny[],
                             double Nz[]) {
    Nx[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
    Nx[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
    Nx[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
    Nx[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
    Nx[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
    Nx[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
    Nx[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
    Nx[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

    Ny[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
    Ny[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
    Ny[2] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
    Ny[3] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
    Ny[4] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
    Ny[5] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
    Ny[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
    Ny[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

    Nz[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
    Nz[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
    Nz[2] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
    Nz[3] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
    Nz[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
    Nz[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
    Nz[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
    Nz[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
  }

  static void evalBasisDeriv(const double pt[], double Nxyz[]) {
    // Nx
    Nxyz[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
    Nxyz[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
    Nxyz[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
    Nxyz[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
    Nxyz[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
    Nxyz[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
    Nxyz[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
    Nxyz[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

    // Ny
    Nxyz[8] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
    Nxyz[9] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
    Nxyz[10] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
    Nxyz[11] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
    Nxyz[12] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
    Nxyz[13] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
    Nxyz[14] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
    Nxyz[15] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

    // Nz
    Nxyz[16] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
    Nxyz[17] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
    Nxyz[18] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
    Nxyz[19] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
    Nxyz[20] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
    Nxyz[21] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
    Nxyz[22] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
    Nxyz[23] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
  }
};

const double TetrahedronWts5[] = {-4.0 / 5.0, 9.0 / 20.0};
const double TetrahedronPts5[] = {
    0.25,      0.25,      0.25, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.5, 1.0 / 6.0,
    1.0 / 6.0, 1.0 / 6.0, 0.5,  1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 0.5};

class Tetra5ptQuadrature {
 public:
  static const index_t NUM_QUAD_PTS = 5;

  static void getQuadPoint(index_t n, double pt[]) {
    if (n == 0) {
      pt[0] = TetrahedronPts5[0];
      pt[1] = TetrahedronPts5[1];
      pt[2] = TetrahedronPts5[2];
    } else if (n == 1) {
      pt[0] = TetrahedronPts5[3];
      pt[1] = TetrahedronPts5[4];
      pt[2] = TetrahedronPts5[5];
    } else if (n == 2) {
      pt[0] = TetrahedronPts5[6];
      pt[1] = TetrahedronPts5[7];
      pt[2] = TetrahedronPts5[8];
    } else if (n == 3) {
      pt[0] = TetrahedronPts5[9];
      pt[1] = TetrahedronPts5[10];
      pt[2] = TetrahedronPts5[11];
    } else if (n == 4) {
      pt[0] = TetrahedronPts5[12];
      pt[1] = TetrahedronPts5[13];
      pt[2] = TetrahedronPts5[14];
    }
  }

  static double getQuadWeight(index_t n) {
    if (n == 0) {
      return TetrahedronWts5[0];
    } else if (n == 1 || n == 2 || n == 3 || n == 4) {
      return TetrahedronWts5[1];
    }
    return 0.0;
  }
};

class TetraQuadraticBasis {
 public:
  static const index_t NUM_NODES = 10;
  static void evalBasis(const double pt[], double N[]) {
    double l0 = 1.0 - pt[0] - pt[1] - pt[2];
    double l1 = pt[0];
    double l2 = pt[1];
    double l3 = pt[2];

    // Corner nodes
    N[0] = l0 * (2.0 * l0 - 1.0);
    N[1] = l1 * (2.0 * l1 - 1.0);
    N[2] = l2 * (2.0 * l2 - 1.0);
    N[3] = l3 * (2.0 * l3 - 1.0);

    // Mid-side nodes
    N[4] = 4.0 * l1 * l0;
    N[5] = 4.0 * l1 * l2;
    N[6] = 4.0 * l2 * l0;
    N[7] = 4.0 * l3 * l0;
    N[8] = 4.0 * l1 * l3;
    N[9] = 4.0 * l3 * l2;
  }
  static void evalBasisDeriv(const double pt[], double Nxyz[]) {
    // Nx
    Nxyz[0] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Nxyz[1] = 4.0 * pt[0] - 1.0;
    Nxyz[2] = 0.0;
    Nxyz[3] = 0.0;
    Nxyz[4] = -4.0 * (2.0 * pt[0] + pt[1] + pt[2] - 1.0);
    Nxyz[5] = 4.0 * pt[1];
    Nxyz[6] = -4.0 * pt[1];
    Nxyz[7] = -4.0 * pt[2];
    Nxyz[8] = 4.0 * pt[2];
    Nxyz[9] = 0.0;

    // Ny
    Nxyz[10] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Nxyz[11] = 0.0;
    Nxyz[12] = 4.0 * pt[1] - 1.0;
    Nxyz[13] = 0.0;
    Nxyz[14] = -4.0 * pt[0];
    Nxyz[15] = 4.0 * pt[0];
    Nxyz[16] = -4.0 * (pt[0] + 2.0 * pt[1] + pt[2] - 1.0);
    Nxyz[17] = -4.0 * pt[2];
    Nxyz[18] = 0.0;
    Nxyz[19] = 4.0 * pt[2];

    // Nz
    Nxyz[20] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Nxyz[21] = 0.0;
    Nxyz[22] = 0.0;
    Nxyz[23] = 4.0 * pt[2] - 1.0;
    Nxyz[24] = -4.0 * pt[0];
    Nxyz[25] = 0.0;
    Nxyz[26] = -4.0 * pt[1];
    Nxyz[27] = -4.0 * (pt[0] + pt[1] + 2.0 * pt[2] - 1.0);
    Nxyz[28] = 4.0 * pt[0];
    Nxyz[29] = 4.0 * pt[1];
  }
};

}  // namespace A2D

#endif  // QUADRATURE_H
