#ifndef BRICK_BASIS_3D_H
#define BRICK_BASIS_3D_H

#include <cstddef>

const double GaussQuadPts2[] = {-0.577350269189626, 0.577350269189626};
const double GaussQuadWts2[] = {1.0, 1.0};

class HexQuadrature {
 public:
  static const int NUM_QUAD_PTS = 8;

  static unsigned int getNumQuadPoints() { return NUM_QUAD_PTS; }
  static void getQuadPoint(const int index, double pt[]) {
    pt[0] = GaussQuadPts2[index % 2];
    pt[1] = GaussQuadPts2[(index % 4) / 2];
    pt[2] = GaussQuadPts2[index / 4];
  }
  static double getQuadWeight(const int index) {
    return (GaussQuadWts2[index % 2] * GaussQuadWts2[(index % 4) / 2] *
            GaussQuadWts2[index / 4]);
  }
};

template <class Quadrature>
class HexBasis {
 public:
  static const int NUM_NODES = 8;
  static const int SPATIAL_DIM = 3;
  typedef Quadrature quadrature;

  /*
    Interpolate from element-oriented data to quadrature point data
  */
  template <const int num_vars, class ElementArray, class QuadPointArray>
  static void interp(ElementArray& input, QuadPointArray& output) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double N[8];
      N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
      N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

      for (std::size_t i = 0; i < input.extent(0); i++) {
        for (int ii = 0; ii < num_vars; ii++) {
          output(i, j, ii) = N[0] * input(i, 0, ii) + N[1] * input(i, 1, ii) +
                             N[2] * input(i, 2, ii) + N[3] * input(i, 3, ii) +
                             N[4] * input(i, 4, ii) + N[5] * input(i, 5, ii) +
                             N[6] * input(i, 6, ii) + N[7] * input(i, 7, ii);
        }
      }
    }
  }

  /*
    Compute the Jacobian transformation at each quadrature point
  */
  template <typename T, class ElementNodeArray, class QuadPointDetJArray,
            class QuadPointJacobianArray>
  static void compute_jtrans(ElementNodeArray& X, QuadPointDetJArray& detJ,
                             QuadPointJacobianArray& Jinv) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double Nx[8];
      Nx[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
      Nx[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

      double Ny[8];
      Ny[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[2] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[3] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[4] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
      Ny[5] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

      double Nz[8];
      Nz[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[2] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[3] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
      Nz[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

      for (std::size_t i = 0; i < X.extent(0); i++) {
        // Compute the Jacobian transformation
        A2D::Mat<T, 3, 3> J;
        for (int ii = 0; ii < 3; ii++) {
          J(ii, 0) = Nx[0] * X(i, 0, ii) + Nx[1] * X(i, 1, ii) +
                     Nx[2] * X(i, 2, ii) + Nx[3] * X(i, 3, ii) +
                     Nx[4] * X(i, 4, ii) + Nx[5] * X(i, 5, ii) +
                     Nx[6] * X(i, 6, ii) + Nx[7] * X(i, 7, ii);

          J(ii, 1) = Ny[0] * X(i, 0, ii) + Ny[1] * X(i, 1, ii) +
                     Ny[2] * X(i, 2, ii) + Ny[3] * X(i, 3, ii) +
                     Ny[4] * X(i, 4, ii) + Ny[5] * X(i, 5, ii) +
                     Ny[6] * X(i, 6, ii) + Ny[7] * X(i, 7, ii);

          J(ii, 2) = Nz[0] * X(i, 0, ii) + Nz[1] * X(i, 1, ii) +
                     Nz[2] * X(i, 2, ii) + Nz[3] * X(i, 3, ii) +
                     Nz[4] * X(i, 4, ii) + Nz[5] * X(i, 5, ii) +
                     Nz[6] * X(i, 6, ii) + Nz[7] * X(i, 7, ii);
        }

        // Compute the 3x3 matrix inverse
        A2D::Mat<T, 3, 3> jinv;
        A2D::Mat3x3Inverse(J, jinv);
        A2D::Mat3x3Det(J, detJ(i, j));

        // Copy values of the inverse of the Jacobian
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Jinv(i, j, ii, jj) = jinv(ii, jj);
          }
        }
      }
    }
  }

  template <typename T, const int num_vars, class ElementSolutionArray,
            class QuadPointGradientArray>
  static void gradient(ElementSolutionArray& U, QuadPointGradientArray& Uxi) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);

      double Nx[8];
      Nx[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
      Nx[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

      double Ny[8];
      Ny[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[2] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[3] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[4] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
      Ny[5] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

      double Nz[8];
      Nz[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[2] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[3] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
      Nz[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

      for (std::size_t i = 0; i < U.extent(0); i++) {
        for (int ii = 0; ii < num_vars; ii++) {
          Uxi(i, j, ii, 0) = Nx[0] * U(i, 0, ii) + Nx[1] * U(i, 1, ii) +
                             Nx[2] * U(i, 2, ii) + Nx[3] * U(i, 3, ii) +
                             Nx[4] * U(i, 4, ii) + Nx[5] * U(i, 5, ii) +
                             Nx[6] * U(i, 6, ii) + Nx[7] * U(i, 7, ii);

          Uxi(i, j, ii, 1) = Ny[0] * U(i, 0, ii) + Ny[1] * U(i, 1, ii) +
                             Ny[2] * U(i, 2, ii) + Ny[3] * U(i, 3, ii) +
                             Ny[4] * U(i, 4, ii) + Ny[5] * U(i, 5, ii) +
                             Ny[6] * U(i, 6, ii) + Ny[7] * U(i, 7, ii);

          Uxi(i, j, ii, 2) = Nz[0] * U(i, 0, ii) + Nz[1] * U(i, 1, ii) +
                             Nz[2] * U(i, 2, ii) + Nz[3] * U(i, 3, ii) +
                             Nz[4] * U(i, 4, ii) + Nz[5] * U(i, 5, ii) +
                             Nz[6] * U(i, 6, ii) + Nz[7] * U(i, 7, ii);
        }
      }
    }
  }

  template <typename T, class Model, class QuadPointModelDataArray,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray>
  static void energy(QuadPointModelDataArray& Edata, QuadPointDetJArray& detJ,
                     QuadPointJacobianArray& Jinv, QuadPointGradientArray& Uxi,
                     T& energy) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double weight = Quadrature::getQuadWeight(j);

      for (std::size_t i = 0; i < Uxi.extent(0); i++) {
        // Extract Jinv
        A2D::Mat<T, 3, 3> Jinv0;
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0;
        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        T wdetJ = weight * detJ(i, j);
        energy += Model::compute_energy(i, j, Edata, wdetJ, Jinv0, Uxi0);
      }
    }
  }

  /*
    Residuals that depend on U
  */
  template <typename T, class Model, class QuadPointModelDataArray,
            class QuadPointDetJArray, class QuadPointSolutionArray,
            class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointSolutionArray& Uq, ElementResidualArray& res) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[8];
      N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
      N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

      for (std::size_t i = 0; i < detJ.extent(0); i++) {
        A2D::Vec<T, Model::NUM_VARS> U0, Ub;
        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          U0(ii) = Uq(i, j, ii);
        }

        Model::compute_residual(i, j, data, weight * detJ(i, j), U0, Ub);

        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          for (int k = 0; k < 8; k++) {
            res(i, k, ii) += N[k] * Ub(ii);
          }
        }
      }
    }
  }

  /*
    Residuals that depend on U,xi
  */
  template <typename T, class Model, class QuadPointModelDataArray,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& res) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nx[8];
      Nx[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
      Nx[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

      double Ny[8];
      Ny[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[2] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[3] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[4] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
      Ny[5] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

      double Nz[8];
      Nz[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[2] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[3] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
      Nz[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

      for (std::size_t i = 0; i < detJ.extent(0); i++) {
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0, Uxib;

        // Extract Jinv
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        Model::compute_residual(i, j, data, weight * detJ(i, j), Jinv0, Uxi0,
                                Uxib);

        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          for (int k = 0; k < 8; k++) {
            res(i, k, ii) +=
                Nx[k] * Uxib(ii, 0) + Ny[k] * Uxib(ii, 1) + Nz[k] * Uxib(ii, 2);
          }
        }
      }
    }
  }

  template <typename T, class Model, class QuadPointModelDataArray,
            class QuadPointDetJArray, class QuadPointSolutionArray,
            class ElementResidualArray>
  static void jacobians(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointSolutionArray& Uq, ElementResidualArray& jac) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[8];
      N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
      N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

      for (std::size_t i = 0; i < detJ.extent(0); i++) {
        A2D::Vec<T, Model::NUM_VARS> U0, Ub;
        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          U0(ii) = Uq(i, j, ii);
        }

        // The Jacobian of the energy
        A2D::Mat<T, Model::NUM_VARS, Model::NUM_VARS> ja;

        Model::compute_jacobian(i, j, data, weight * detJ(i, j), U0, Ub, ja);

        for (int ky = 0; ky < 8; ky++) {
          for (int iy = 0; iy < Model::NUM_VARS; iy++) {
            for (int ix = 0; ix < Model::NUM_VARS; ix++) {
              T n = N[ky] * ja(iy, ix);
              for (int kx = 0; kx < 8; kx++) {
                jac(i, ky, kx, iy, ix) += N[kx] * n;
              }
            }
          }
        }
      }
    }
  }

  template <typename T, class Model, class QuadPointModelDataArray,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void jacobians(QuadPointModelDataArray& Edata,
                        QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& jac) {
    for (std::size_t j = 0; j < Quadrature::getNumQuadPoints(); j++) {
      double pt[3];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nx[8];
      Nx[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[1] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
      Nx[2] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[3] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
      Nx[4] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[5] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
      Nx[6] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
      Nx[7] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);

      double Ny[8];
      Ny[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[2] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
      Ny[3] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
      Ny[4] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
      Ny[5] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
      Ny[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);

      double Nz[8];
      Nz[0] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[1] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[2] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[3] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);
      Nz[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nz[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nz[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
      Nz[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

      for (std::size_t i = 0; i < detJ.extent(0); i++) {
        A2D::Mat<T, 3, 3> Jinv0;
        A2D::Mat<T, Model::NUM_VARS, 3> Uxi0, Uxib;

        // Extract Jinv
        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        for (int ii = 0; ii < Model::NUM_VARS; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        // The Jacobian of the energy
        A2D::SymmTensor<T, Model::NUM_VARS, 3> ja;

        Model::compute_jacobian(i, j, Edata, weight * detJ(i, j), Jinv0, Uxi0,
                                Uxib, ja);

        for (int ky = 0; ky < 8; ky++) {
          for (int iy = 0; iy < Model::NUM_VARS; iy++) {
            for (int ix = 0; ix < Model::NUM_VARS; ix++) {
              T nx = Nx[ky] * ja(iy, 0, ix, 0) + Ny[ky] * ja(iy, 1, ix, 0) +
                     Nz[ky] * ja(iy, 2, ix, 0);
              T ny = Nx[ky] * ja(iy, 0, ix, 1) + Ny[ky] * ja(iy, 1, ix, 1) +
                     Nz[ky] * ja(iy, 2, ix, 1);
              T nz = Nx[ky] * ja(iy, 0, ix, 2) + Ny[ky] * ja(iy, 1, ix, 2) +
                     Nz[ky] * ja(iy, 2, ix, 2);

              for (int kx = 0; kx < 8; kx++) {
                jac(i, ky, kx, iy, ix) +=
                    Nx[kx] * nx + Ny[kx] * ny + Nz[kx] * nz;
              }
            }
          }
        }
      }
    }
  }
};

#endif  // BRICK_BASIS_3D_H