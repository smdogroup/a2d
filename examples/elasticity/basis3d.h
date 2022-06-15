#ifndef BRICK_BASIS_3D_H
#define BRICK_BASIS_3D_H

#include <cstddef>

#include "parallel.h"

namespace A2D {

const double GaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double GaussQuadWts2[] = {1.0, 1.0};

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
  static void evalBasisDeriv(const double pt[], double Nx[], double Ny[],
                             double Nz[]) {
    Nx[0] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Nx[1] = 4.0 * pt[0] - 1.0;
    Nx[2] = 0.0;
    Nx[3] = 0.0;
    Nx[4] = -4.0 * (2.0 * pt[0] + pt[1] + pt[2] - 1.0);
    Nx[5] = 4.0 * pt[1];
    Nx[6] = -4.0 * pt[1];
    Nx[7] = -4.0 * pt[2];
    Nx[8] = 4.0 * pt[2];
    Nx[9] = 0.0;

    Ny[0] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Ny[1] = 0.0;
    Ny[2] = 4.0 * pt[1] - 1.0;
    Ny[3] = 0.0;
    Ny[4] = -4.0 * pt[0];
    Ny[5] = 4.0 * pt[0];
    Ny[6] = -4.0 * (pt[0] + 2.0 * pt[1] + pt[2] - 1.0);
    Ny[7] = -4.0 * pt[2];
    Ny[8] = 0.0;
    Ny[9] = 4.0 * pt[2];

    Nz[0] = 4.0 * pt[0] + 4.0 * pt[1] + 4.0 * pt[2] - 3.0;
    Nz[1] = 0.0;
    Nz[2] = 0.0;
    Nz[3] = 4.0 * pt[2] - 1.0;
    Nz[4] = -4.0 * pt[0];
    Nz[5] = 0.0;
    Nz[6] = -4.0 * pt[1];
    Nz[7] = -4.0 * (pt[0] + pt[1] + 2.0 * pt[2] - 1.0);
    Nz[8] = 4.0 * pt[0];
    Nz[9] = 4.0 * pt[1];
  }
};

template <class Basis, class Quadrature>
class Basis3D {
 public:
  static const index_t NUM_NODES = Basis::NUM_NODES;
  static const index_t SPATIAL_DIM = 3;
  typedef Quadrature quadrature;

  /*
    Interpolate from element-oriented data to quadrature point data
  */
  template <const index_t M, class ElementArray, class QuadPointArray>
  static void interp(ElementArray& input, QuadPointArray& output) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double N[Basis::NUM_NODES];
      Basis::evalBasis(pt, N);

      const A2D::index_t npts = input.extent(0);
      A2D::parallel_for(npts, [&, N](A2D::index_t i) -> void {
        for (index_t ii = 0; ii < M; ii++) {
          output(i, j, ii) = 0.0;
          for (index_t kk = 0; kk < NUM_NODES; kk++) {
            output(i, j, ii) += N[kk] * input(i, kk, ii);
          }
        }
      });
    }
  }

  /*
    Compute the Jacobian transformation at each quadrature point
  */
  template <typename T, class ElementNodeArray, class QuadPointDetJArray,
            class QuadPointJacobianArray>
  static void compute_jtrans(ElementNodeArray& X, QuadPointDetJArray& detJ,
                             QuadPointJacobianArray& Jinv) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double Nx[Basis::NUM_NODES], Ny[Basis::NUM_NODES], Nz[Basis::NUM_NODES];
      Basis::evalBasisDeriv(pt, Nx, Ny, Nz);

      const A2D::index_t npts = X.extent(0);
      A2D::parallel_for(npts, [&, Nx, Ny, Nz](A2D::index_t i) -> void {
        // Compute the Jacobian transformation
        A2D::Mat<T, 3, 3> J;
        for (index_t ii = 0; ii < 3; ii++) {
          J(ii, 0u) = 0.0;
          J(ii, 1u) = 0.0;
          J(ii, 2u) = 0.0;
          for (index_t kk = 0; kk < NUM_NODES; kk++) {
            J(ii, 0u) += Nx[kk] * X(i, kk, ii);
            J(ii, 1u) += Ny[kk] * X(i, kk, ii);
            J(ii, 2u) += Nz[kk] * X(i, kk, ii);
          }
        }

        // Compute the 3x3 matrix inverse
        A2D::Mat<T, 3, 3> jinv;
        A2D::Mat3x3Inverse(J, jinv);
        A2D::Mat3x3Det(J, detJ(i, j));

        // Copy values of the inverse of the Jacobian
        for (index_t ii = 0; ii < 3; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Jinv(i, j, ii, jj) = jinv(ii, jj);
          }
        }
      });
    }
  }

  template <typename T, const index_t num_vars, class ElementSolutionArray,
            class QuadPointGradientArray>
  static void gradient(ElementSolutionArray& U, QuadPointGradientArray& Uxi) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double Nx[Basis::NUM_NODES], Ny[Basis::NUM_NODES], Nz[Basis::NUM_NODES];
      Basis::evalBasisDeriv(pt, Nx, Ny, Nz);

      const A2D::index_t npts = U.extent(0);
      A2D::parallel_for(npts, [&, Nx, Ny, Nz](A2D::index_t i) -> void {
        for (index_t ii = 0; ii < num_vars; ii++) {
          Uxi(i, j, ii, 0u) = 0.0;
          Uxi(i, j, ii, 1u) = 0.0;
          Uxi(i, j, ii, 2u) = 0.0;

          for (index_t kk = 0; kk < NUM_NODES; kk++) {
            Uxi(i, j, ii, 0u) += Nx[kk] * U(i, kk, ii);
            Uxi(i, j, ii, 1u) += Ny[kk] * U(i, kk, ii);
            Uxi(i, j, ii, 2u) += Nz[kk] * U(i, kk, ii);
          }
        }
      });
    }
  }

  template <typename T, index_t M, class FunctorType, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray>
  static T integrate(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                     QuadPointGradientArray& Uxi,
                     const FunctorType& integrand) {
    T value = 0.0;
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);

      const A2D::index_t npts = detJ.extent(0);
      value += A2D::parallel_reduce<T>(npts, [&](A2D::index_t i) -> T {
        // Extract Jinv
        A2D::Mat<T, 3, 3> Jinv0;
        for (index_t ii = 0; ii < 3; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        A2D::Mat<T, M, 3> Uxi0;
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        T wdetJ = weight * detJ(i, j);
        return integrand(i, j, wdetJ, Jinv0, Uxi0);
      });
    }

    return value;
  }

  /*
    Residuals that depend on U
  */
  template <typename T, index_t M, class FunctorType, class QuadPointDetJArray,
            class QuadPointSolutionArray, class ElementResidualArray>
  static void residuals(QuadPointDetJArray& detJ, QuadPointSolutionArray& Uq,
                        const FunctorType& resfunc, ElementResidualArray& res) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[Basis::NUM_NODES];
      Basis::evalBasis(pt, N);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(npts, [&, N](A2D::index_t i) -> void {
        A2D::Vec<T, M> U0, Ub;
        for (index_t ii = 0; ii < M; ii++) {
          U0(ii) = Uq(i, j, ii);
        }

        T wdetJ = weight * detJ(i, j);
        resfunc(i, j, wdetJ, U0, Ub);

        auto resi = MakeSlice(res, i);
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t k = 0; k < NUM_NODES; k++) {
            resi(k, ii) += N[k] * Ub(ii);
          }
        }
      });
    }
  }

  /*
    Residuals that depend on U,xi
  */
  template <typename T, index_t M, class FunctorType, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void residuals(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi, const FunctorType& resfunc,
                        ElementResidualArray& res) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nx[Basis::NUM_NODES], Ny[Basis::NUM_NODES], Nz[Basis::NUM_NODES];
      Basis::evalBasisDeriv(pt, Nx, Ny, Nz);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(npts, [&, Nx, Ny, Nz](A2D::index_t i) -> void {
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, M, 3> Uxi0, Uxib;

        // Extract Jinv
        for (index_t ii = 0; ii < 3; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        T wdetJ = weight * detJ(i, j);
        resfunc(i, j, wdetJ, Jinv0, Uxi0, Uxib);

        auto resi = MakeSlice(res, i);
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t k = 0; k < NUM_NODES; k++) {
            resi(k, ii) += Nx[k] * Uxib(ii, 0u) + Ny[k] * Uxib(ii, 1u) +
                           Nz[k] * Uxib(ii, 2u);
          }
        }
      });
    }
  }

  template <typename T, index_t M, class FunctorType, class QuadPointDetJArray,
            class QuadPointSolutionArray, class ElementJacArray>
  static void jacobians(QuadPointDetJArray& detJ, QuadPointSolutionArray& Uq,
                        const FunctorType& jacfunc, ElementJacArray& jac) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[Basis::NUM_NODES];
      Basis::evalBasis(pt, N);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(npts, [&, N](A2D::index_t i) -> void {
        A2D::Vec<T, M> U0, Ub;
        for (index_t ii = 0; ii < M; ii++) {
          U0(ii) = Uq(i, j, ii);
        }

        // The Jacobian of the energy
        A2D::SymmMat<T, M> ja;
        T wdetJ = weight * detJ(i, j);
        jacfunc(i, j, wdetJ, U0, ja);

        auto jaci = MakeSlice(jac, i);
        for (index_t ky = 0; ky < NUM_NODES; ky++) {
          for (index_t iy = 0; iy < M; iy++) {
            for (index_t ix = 0; ix < M; ix++) {
              T n = N[ky] * ja(iy, ix);
              for (index_t kx = 0; kx < NUM_NODES; kx++) {
                jac(i, ky, kx, iy, ix) += N[kx] * n;
              }
            }
          }
        }
      });
    }
  }

  template <typename T, index_t M, class FunctorType, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointGradientArray,
            class ElementResidualArray>
  static void jacobians(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi, const FunctorType& jacfunc,
                        ElementResidualArray& jac) {
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nx[Basis::NUM_NODES], Ny[Basis::NUM_NODES], Nz[Basis::NUM_NODES];
      Basis::evalBasisDeriv(pt, Nx, Ny, Nz);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(npts, [&, Nx, Ny, Nz](A2D::index_t i) -> void {
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, M, 3> Uxi0, Uxib;

        // Extract Jinv
        for (index_t ii = 0; ii < 3; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        // The Jacobian of the energy
        A2D::SymmTensor<T, M, 3> ja;

        T wdetJ = weight * detJ(i, j);
        jacfunc(i, j, wdetJ, Jinv0, Uxi0, Uxib, ja);

        auto jaci = MakeSlice(jac, i);
        for (index_t ky = 0; ky < NUM_NODES; ky++) {
          for (index_t iy = 0; iy < M; iy++) {
            for (index_t ix = 0; ix < M; ix++) {
              T nx = Nx[ky] * ja(iy, 0u, ix, 0u) + Ny[ky] * ja(iy, 1u, ix, 0u) +
                     Nz[ky] * ja(iy, 2u, ix, 0u);
              T ny = Nx[ky] * ja(iy, 0u, ix, 1u) + Ny[ky] * ja(iy, 1u, ix, 1u) +
                     Nz[ky] * ja(iy, 2u, ix, 1u);
              T nz = Nx[ky] * ja(iy, 0u, ix, 2u) + Ny[ky] * ja(iy, 1u, ix, 2u) +
                     Nz[ky] * ja(iy, 2u, ix, 2u);

              for (index_t kx = 0; kx < NUM_NODES; kx++) {
                jaci(ky, kx, iy, ix) += Nx[kx] * nx + Ny[kx] * ny + Nz[kx] * nz;
              }
            }
          }
        }
      });
    }
  }
};

}  // namespace A2D

#endif  // BRICK_BASIS_3D_H
