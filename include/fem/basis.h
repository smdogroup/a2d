#ifndef A2D_BASIS_H
#define A2D_BASIS_H

#include <cstddef>

#include "a2dtmp2d.h"
#include "a2dtmp3d.h"
#include "block_numeric.h"
#include "parallel.h"
#include "quadrature.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/**
 * @brief A collection of basis operations.
 *
 * @tparam spatial_dim 2 or 3
 * @tparam BasisFunc basis function class that defines N and Nx evaluation
 *                   given quadrature
 * @tparam Quadrature quadrature class that defines the quadrature scheme
 */
template <index_t spatial_dim, class BasisFunc, class Quadrature>
class BasisOps {
 public:
  static_assert(spatial_dim == 3 or spatial_dim == 2,
                "spatial_dim must be 2 or 3");
  static const index_t NUM_NODES = BasisFunc::NUM_NODES;
  static const index_t SPATIAL_DIM = spatial_dim;
  typedef Quadrature quadrature;

  /**
   * @brief Interpolate from element-oriented data to quadrature point data
   *
   * @tparam ElementArray array of element data, shape: (nelems, nodes_per_elem,
   *                      ndata_per_nodes)
   * @tparam QuadPointArray array of quadrature data, shape: (nelems,
   *                        quad_pts_per_elem, ndata_per_nodes)
   */
  template <index_t ndata_per_nodes, class ElementArray, class QuadPointArray>
  static void interp(ElementArray& input, QuadPointArray& output) {
    Timer timer("BasisOps::interp()");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);

      double N[BasisFunc::NUM_NODES];
      BasisFunc::evalBasis(pt, N);

      const A2D::index_t npts = input.extent(0);
      auto lam = A2D_LAMBDA(A2D::index_t i)->void {
        for (index_t ii = 0; ii < ndata_per_nodes; ii++) {
          output(i, j, ii) = 0.0;
          for (index_t kk = 0; kk < NUM_NODES; kk++) {
            output(i, j, ii) += N[kk] * input(i, kk, ii);
          }
        }
      };
      A2D::parallel_for(npts, lam);
    }
  }

  /**
   * @brief Add the contributions to the element-oriented data from quad point
   *        data
   *
   * @tparam QuadPointArray array of quadrature data, shape: (nelems,
   *                        quad_pts_per_elem, ndata_per_nodes)
   * @tparam ElementArray array of element data, shape: (nelems, nodes_per_elem,
   *                      ndata_per_nodes)
   */
  template <const index_t ndata_per_nodes, class QuadPointArray,
            class ElementArray>
  static void interpReverseAdd(QuadPointArray& input, ElementArray& output) {
    Timer timer("BasisOps::interpReverseAdd()");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);

      double N[BasisFunc::NUM_NODES];
      BasisFunc::evalBasis(pt, N);

      const A2D::index_t npts = input.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            for (index_t ii = 0; ii < ndata_per_nodes; ii++) {
              for (index_t kk = 0; kk < NUM_NODES; kk++) {
                output(i, kk, ii) += N[kk] * input(i, j, ii);
              }
            }
          });
    }
  }

  /**
   * @brief Compute the Jacobian transformation at each quadrature point
   *        data
   *
   * @tparam ElementNodeArray array of element nodal data, shape: (nelems,
   *                          nodes_per_elem, spatial_dim)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   */
  template <typename T, class ElementNodeArray, class QuadPointDetJArray,
            class QuadPointJacobianArray>
  static void compute_jtrans(ElementNodeArray& X, QuadPointDetJArray& detJ,
                             QuadPointJacobianArray& Jinv) {
    Timer timer("BasisOps::compute_jtrans()");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);

      double Nxyz[BasisFunc::NUM_NODES * SPATIAL_DIM];
      BasisFunc::evalBasisDeriv(pt, Nxyz);

      const A2D::index_t npts = X.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            // Compute the Jacobian transformation
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> J;
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
                J(ii, idim) = 0.0;
              }
              for (index_t kk = 0; kk < NUM_NODES; kk++) {
                for (index_t ll = 0; ll < SPATIAL_DIM; ll++) {
                  J(ii, ll) +=
                      Nxyz[kk + ll * BasisFunc::NUM_NODES] * X(i, kk, ii);
                }
              }
            }

            // Compute the matrix inverse
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> jinv;
            MatInverse(J, jinv);
            MatDet(J, detJ(i, j));

            // Copy values of the inverse of the Jacobian
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Jinv(i, j, ii, jj) = jinv(ii, jj);
              }
            }
          });
    }
  }

  /**
   * @brief Compute gradient of nodal solution with respect to local coordinates
   *        data
   *
   * @tparam vars_per_node number of components in each nodal solution
   * @tparam ElementSolutionArray solution array organized in element-wise,
   *                              shape: (nelems, nodes_per_elem, vars_per_node)
   * @tparam QuadPointGradientArray array of the solution gradient (w.r.t. local
   *                                coordinates) at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node,
   *                                spatial_tim)
   */
  template <typename T, const index_t vars_per_node, class ElementSolutionArray,
            class QuadPointGradientArray>
  static void gradient(ElementSolutionArray& U, QuadPointGradientArray& Uxi) {
    Timer timer("BasisOps::gradient()");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);

      double Nxyz[BasisFunc::NUM_NODES * SPATIAL_DIM];
      BasisFunc::evalBasisDeriv(pt, Nxyz);

      const A2D::index_t npts = U.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
                Uxi(i, j, ii, idim) = 0.0;
              }

              for (index_t kk = 0; kk < NUM_NODES; kk++) {
                for (index_t ll = 0; ll < SPATIAL_DIM; ll++) {
                  Uxi(i, j, ii, ll) +=
                      Nxyz[kk + ll * BasisFunc::NUM_NODES] * U(i, kk, ii);
                }
              }
            }
          });
    }
  }

  /**
   * @brief Integrate a value over all the elements in the domain
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ) -> val,
   *                     where i is nelems index, j is quad_pt index, wdetJ is
   *                     weight * detJ
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   */
  template <typename T, class FunctorType, class QuadPointDetJArray>
  static T integrate(QuadPointDetJArray& detJ, const FunctorType& integrand) {
    Timer timer("BasisOps::integrate(1)");
    T value = 0.0;
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);

      const A2D::index_t npts = detJ.extent(0);
      value += A2D::parallel_reduce<T>(
          npts, [j, weight, detJ, integrand](A2D::index_t i) -> T {
            T wdetJ = weight * detJ(i, j);
            return integrand(i, j, wdetJ);
          });
    }

    return value;
  }

  /**
   * @brief Integrate a value over all the elements in the domain
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, Jinv0,
   *                     Uxi0) -> val, where i is nelems index, j is quad_pt
   *                     index, wdetJ is weight * detJ, Jinv0 is 3x3 or 2x2
   *                     Jacobian inverse matrix, Uxi0 is solution gradient
   *                     matrix
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   * @tparam QuadPointGradientArray array of the solution gradient (w.r.t. local
   *                                coordinates) at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node,
   *                                spatial_tim)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray>
  static T integrate(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                     QuadPointGradientArray& Uxi,
                     const FunctorType& integrand) {
    Timer timer("BasisOps::integrate(2)");
    T value = 0.0;
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);

      const A2D::index_t npts = detJ.extent(0);
      value += A2D::parallel_reduce<T>(
          npts, [j, Jinv, Uxi, weight, detJ, integrand](A2D::index_t i) -> T {
            // Extract Jinv
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Jinv0;
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Jinv0(ii, jj) = Jinv(i, j, ii, jj);
              }
            }

            // Extract Uxi0
            A2D::Mat<T, vars_per_node, SPATIAL_DIM> Uxi0;
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Uxi0(ii, jj) = Uxi(i, j, ii, jj);
              }
            }

            T wdetJ = weight * detJ(i, j);
            return integrand(i, j, wdetJ, Jinv0, Uxi0);
          });
    }

    return value;
  }

  /**
   * @brief compute maximum value over all elements in the domain
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, Jinv0,
   *                     Uxi0) -> val, where i is nelems index, j is quad_pt
   *                     index, wdetJ is weight * detJ, Jinv0 is 3x3 or 2x2
   *                     Jacobian inverse matrix, Uxi0 is solution gradient
   *                     matrix of shape (vars_per_node, spatial_dim)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   * @tparam QuadPointGradientArray array of the solution gradient (w.r.t. local
   *                                coordinates) at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node,
   *                                spatial_tim)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray>
  static T maximum(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                   QuadPointGradientArray& Uxi, const FunctorType& func) {
    Timer timer("BasisOps::maximum()");
    T value = -1e20;
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);

      const A2D::index_t npts = detJ.extent(0);
      // value = A2D::parallel_reduce_max<T>(npts, [&](A2D::index_t i) -> T {
      for (A2D::index_t i = 0; i < npts; i++) {
        // Extract Jinv
        A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Jinv0;
        for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
          for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
            Jinv0(ii, jj) = Jinv(i, j, ii, jj);
          }
        }

        // Extract Uxi0
        A2D::Mat<T, vars_per_node, SPATIAL_DIM> Uxi0;
        for (index_t ii = 0; ii < vars_per_node; ii++) {
          for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
            Uxi0(ii, jj) = Uxi(i, j, ii, jj);
          }
        }

        T wdetJ = weight * detJ(i, j);
        T val = func(i, j, wdetJ, Jinv0, Uxi0);
        if (A2D::RealPart(val) > A2D::RealPart(value)) {
          value = val;
        }
      }
      // );
    }

    return value;
  }

  /**
   * @brief Compute element-wise residuals that depend on U
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, U0,
   *                     Ub) -> val, where i is nelems index, j is quad_pt
   *                     index, wdetJ is weight * detJ, U0 and Ub are vectors
   *                     of size vars_per_node
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointSolutionArray array of solution at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node)
   * @tparam ElementResidualArray array of residual value at each node, shape:
   *                              (nelems, nodes_per_elem, vars_per_node)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointSolutionArray,
            class ElementResidualArray>
  static void residuals(QuadPointDetJArray& detJ, QuadPointSolutionArray& Uq,
                        const FunctorType& resfunc, ElementResidualArray& res) {
    Timer timer("BasisOps::residuals(1)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[BasisFunc::NUM_NODES];
      BasisFunc::evalBasis(pt, N);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            A2D::Vec<T, vars_per_node> U0, Ub;
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              U0(ii) = Uq(i, j, ii);
            }

            T wdetJ = weight * detJ(i, j);
            resfunc(i, j, wdetJ, U0, Ub);

            auto resi = Kokkos::subview(res, i, Kokkos::ALL, Kokkos::ALL);
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t k = 0; k < NUM_NODES; k++) {
                resi(k, ii) += N[k] * Ub(ii);
              }
            }
          });
    }
  }

  /**
   * @brief Compute element-wise residuals that depend on U,xi
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, Jinv0,
   *                     Uxi0, Uxib) -> val, where i is nelems index, j is
   *                     quad_pt index, wdetJ is weight * detJ, Jinv0 is 3x3 or
   *                     2x2 Jacobian inverse matrix, Uxi0 and Uxib are solution
   *                     gradient and seed matrix of shape (vars_per_node,
   *                     spatial_dim)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   * @tparam QuadPointGradientArray array of solution gradient at quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                vars_per_node, spatial_dim)
   * @tparam ElementResidualArray array of residual value at each node, shape:
   *                              (nelems, nodes_per_elem, vars_per_node)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void residuals(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi, const FunctorType& resfunc,
                        ElementResidualArray& res) {
    Timer timer("BasisOps::residuals(2)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nxyz[BasisFunc::NUM_NODES * SPATIAL_DIM];
      BasisFunc::evalBasisDeriv(pt, Nxyz);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Jinv0;
            A2D::Mat<T, vars_per_node, SPATIAL_DIM> Uxi0, Uxib;

            // Extract Jinv
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Jinv0(ii, jj) = Jinv(i, j, ii, jj);
              }
            }

            // Extract Uxi0
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Uxi0(ii, jj) = Uxi(i, j, ii, jj);
              }
            }

            T wdetJ = weight * detJ(i, j);
            resfunc(i, j, wdetJ, Jinv0, Uxi0, Uxib);

            auto resi = Kokkos::subview(res, i, Kokkos::ALL, Kokkos::ALL);
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t k = 0; k < NUM_NODES; k++) {
                for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
                  resi(k, ii) +=
                      Nxyz[k + idim * BasisFunc::NUM_NODES] * Uxib(ii, idim);
                }
              }
            }
          });
    }
  }

  /**
   * @brief Compute element-wise jacobians that depend on U
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, U0,
   *                     ja) -> val, where i is nelems index, j is quad_pt
   *                     index, wdetJ is weight * detJ, U0 is vector of
   *                     size vars_per_node, ja is matrix of shape
   *                     (vars_per_node, vars_per_node)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointSolutionArray array of solution at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node)
   * @tparam ElementJacArray array of jacobian for each element, shape: (nelems,
   *                         nodes_per_elem, nodes_per_elem, vars_per_node,
   *                         vars_per_node)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointSolutionArray,
            class ElementJacArray>
  static void jacobians(QuadPointDetJArray& detJ, QuadPointSolutionArray& Uq,
                        const FunctorType& jacfunc, ElementJacArray& jac) {
    Timer timer("BasisOps::jacobians(1)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double N[BasisFunc::NUM_NODES];
      BasisFunc::evalBasis(pt, N);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            A2D::Vec<T, vars_per_node> U0, Ub;
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              U0(ii) = Uq(i, j, ii);
            }

            // The Jacobian of the energy
            A2D::SymmMat<T, vars_per_node> ja;
            T wdetJ = weight * detJ(i, j);
            jacfunc(i, j, wdetJ, U0, ja);

            for (index_t ky = 0; ky < NUM_NODES; ky++) {
              for (index_t iy = 0; iy < vars_per_node; iy++) {
                for (index_t ix = 0; ix < vars_per_node; ix++) {
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

  /**
   * @brief Compute element-wise jacobians that depend on U,xi
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, Jinv0,
   *                     Uxi0, Uxib, ja) -> val, where i is nelems index, j is
   *                     quad_pt index, wdetJ is weight * detJ, Jinv0 is 3x3 or
   *                     2x2 Jacobian inverse matrix, Uxi0 and Uxib are solution
   *                     gradient and seed matrix of shape (vars_per_node,
   *                     spatial_dim), ja is matrix of shape (vars_per_node,
   *                     vars_per_node)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   * @tparam QuadPointGradientArray array of solution gradient at quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                vars_per_node, spatial_dim)
   * @tparam ElementJacArray array of jacobian for each element, shape: (nelems,
   *                         nodes_per_elem, nodes_per_elem, vars_per_node,
   *                         vars_per_node)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray, class ElementJacArray>
  static void jacobians(QuadPointDetJArray& detJ, QuadPointJacobianArray& Jinv,
                        QuadPointGradientArray& Uxi, const FunctorType& jacfunc,
                        ElementJacArray& jac) {
    Timer timer("BasisOps::jacobians(2)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double pt[SPATIAL_DIM];
      Quadrature::getQuadPoint(j, pt);
      double weight = Quadrature::getQuadWeight(j);

      double Nxyz[BasisFunc::NUM_NODES * SPATIAL_DIM];
      BasisFunc::evalBasisDeriv(pt, Nxyz);

      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i) {
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Jinv0;
            A2D::Mat<T, vars_per_node, SPATIAL_DIM> Uxi0, Uxib;

            // Extract Jinv
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Jinv0(ii, jj) = Jinv(i, j, ii, jj);
              }
            }

            // Extract Uxi0
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Uxi0(ii, jj) = Uxi(i, j, ii, jj);
              }
            }

            // The Jacobian of the energy
            A2D::SymmTensor<T, vars_per_node, SPATIAL_DIM> ja;

            T wdetJ = weight * detJ(i, j);
            jacfunc(i, j, wdetJ, Jinv0, Uxi0, Uxib, ja);

            T nxyz[SPATIAL_DIM];

            auto jaci = Kokkos::subview(jac, i, Kokkos::ALL, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL);
            for (index_t ky = 0; ky < NUM_NODES; ky++) {
              for (index_t iy = 0; iy < vars_per_node; iy++) {
                for (index_t ix = 0; ix < vars_per_node; ix++) {
                  for (index_t idim = 0; idim < SPATIAL_DIM; idim++) {
                    nxyz[idim] = 0.0;
                    for (index_t jdim = 0; jdim < SPATIAL_DIM; jdim++) {
                      nxyz[idim] += Nxyz[ky + jdim * BasisFunc::NUM_NODES] *
                                    ja(iy, jdim, ix, idim);
                    }
                  }

                  for (index_t kx = 0; kx < NUM_NODES; kx++) {
                    for (index_t iidim = 0; iidim < SPATIAL_DIM; iidim++) {
                      jaci(ky, kx, iy, ix) +=
                          Nxyz[kx + iidim * BasisFunc::NUM_NODES] * nxyz[iidim];
                    }
                  }
                }
              }
            }
          });
    }
  }

  /**
   * @brief Adjoint-residual products for residuals that depend on U
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, U0,
   *                     Psi0) -> val, where i is nelems index, j is quad_pt
   *                     index, wdetJ is weight * detJ, U0 and Psi0 are vectors
   *                     of size vars_per_node
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointSolutionArray array of solution at each quadrature, shape:
   *                                (nelems, quad_pts_per_elem, vars_per_node)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointSolutionArray>
  static void adjoint_product(QuadPointDetJArray& detJ,
                              QuadPointSolutionArray& Uq,
                              QuadPointSolutionArray& Psiq,
                              const FunctorType& func) {
    Timer timer("BasisOps::adjoint_product(1)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);
      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            A2D::Vec<T, vars_per_node> U0, Psi0;
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              U0(ii) = Uq(i, j, ii);
              Psi0(ii) = Psiq(i, j, ii);
            }

            T wdetJ = weight * detJ(i, j);
            func(i, j, wdetJ, U0, Psi0);
          });
    }
  }

  /*
    Adjoint-residual products for residuals that depend on U,xi
  */

  /**
   * @brief Adjoint-residual products for residuals that depend on U,xi
   *
   * @tparam FunctorType type of a callable that takes in (i, j, wdetJ, Jinv0,
   *                     Uxi0, Psi0) -> val, where i is nelems index, j is
   *                     quad_pt index, wdetJ is weight * detJ, Jinv0 is 3x3 or
   *                     2x2 Jacobian inverse matrix, Uxi0 and Psi0 are solution
   *                     and adjoint gradient of shape (vars_per_node,
   *                     spatial_dim)
   * @tparam QuadPointDetJArray array of Jacobian determinant at each quadrature
   *                            , shape: (nelems, quad_pts_per_elem)
   * @tparam QuadPointJacobianArray array of Jacobian matrix at each quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                spatial_dim, spatial_dim)
   * @tparam QuadPointGradientArray array of solution gradient at quadrature,
   *                                shape: (nelems, quad_pts_per_elem,
   *                                vars_per_node, spatial_dim)
   */
  template <typename T, index_t vars_per_node, class FunctorType,
            class QuadPointDetJArray, class QuadPointJacobianArray,
            class QuadPointGradientArray>
  static void adjoint_product(QuadPointDetJArray& detJ,
                              QuadPointJacobianArray& Jinv,
                              QuadPointGradientArray& Uxi,
                              QuadPointGradientArray& Psixi,
                              const FunctorType& func) {
    Timer timer("BasisOps::adjoint_product(2)");
    for (A2D::index_t j = 0; j < Quadrature::NUM_QUAD_PTS; j++) {
      double weight = Quadrature::getQuadWeight(j);
      const A2D::index_t npts = detJ.extent(0);
      A2D::parallel_for(
          npts, A2D_LAMBDA(A2D::index_t i)->void {
            A2D::Mat<T, SPATIAL_DIM, SPATIAL_DIM> Jinv0;
            A2D::Mat<T, vars_per_node, SPATIAL_DIM> Uxi0, Psi0;

            // Extract Jinv
            for (index_t ii = 0; ii < SPATIAL_DIM; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Jinv0(ii, jj) = Jinv(i, j, ii, jj);
              }
            }

            // Extract Uxi0
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Uxi0(ii, jj) = Uxi(i, j, ii, jj);
              }
            }

            // Extract Psi0
            for (index_t ii = 0; ii < vars_per_node; ii++) {
              for (index_t jj = 0; jj < SPATIAL_DIM; jj++) {
                Psi0(ii, jj) = Psixi(i, j, ii, jj);
              }
            }

            T wdetJ = weight * detJ(i, j);
            func(i, j, wdetJ, Jinv0, Uxi0, Psi0);
          });
    }
  }
};

}  // namespace A2D

#endif  // A2D_BASIS_H
