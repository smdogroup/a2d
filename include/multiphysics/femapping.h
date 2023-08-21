#ifndef A2D_FE_MAPPING_H
#define A2D_FE_MAPPING_H

#include "a2dobjs.h"
// #include "a2dvecops3d.h"
#include "multiphysics/febasis.h"

namespace A2D {

/**
 * @brief 3D or 2D volume transform transform
 *
 */
template <typename T, index_t dim>
class InteriorMapping {
 public:
  template <class FiniteElementGeometry>
  InteriorMapping(const FiniteElementGeometry& geo, T& detJ)
      : J(geo.template get<0>().get_grad()), detJ(detJ) {
    // Compute the inverse of the transformation
    A2D::MatInverse(J, Jinv);

    // Compute the determinant of the Jacobian matrix
    A2D::MatDet(J, detJ);
  }

  template <class FiniteElementSpace>
  void transform(const FiniteElementSpace& in, FiniteElementSpace& out) {
    in.transform(detJ, J, Jinv, out);
  }

  template <class FiniteElementSpace>
  void rtransform(const FiniteElementSpace& in, FiniteElementSpace& out) {
    in.rtransform(detJ, J, Jinv, out);
  }

  // TODO: maybe the best way is to put jtransform in FESpace?
  template <class PDEIntegrand>
  void jtransform(const typename PDEIntegrand::QMatType& mat_in,
                  typename PDEIntegrand::QMatType& mat_out) {
    constexpr index_t ncomp = PDEIntegrand::FiniteElementSpace::ncomp;

    typename PDEIntegrand::FiniteElementSpace pref, p, Jp;

    for (index_t k = 0; k < ncomp; k++) {
      pref.zero();
      pref[k] = T(1.0);
      transform(pref, p);

      // Compute Jp
      Jp.zero();
      // TODO: use MatVec operation instead
      for (index_t j = 0; j < ncomp; j++) {
        for (index_t i = 0; i < ncomp; i++) {
          Jp[i] += mat_in(i, j) * p[j];
        }
      }
      rtransform(Jp, pref);
      for (index_t i = 0; i < ncomp; i++) {
        mat_out(i, k) = pref[i];
      }
    }
  }

 private:
  const A2D::Mat<T, dim, dim>& J;
  A2D::Mat<T, dim, dim> Jinv;
  T& detJ;
};

/**
 * @brief surface transformation - TODO: this needs to be fixed
 */
template <typename T, index_t dim>
class SurfaceMapping {
 public:
  static const index_t dim_surf = dim - 1;

  template <class FiniteElementGeometry>
  SurfaceMapping(const FiniteElementGeometry& geo, T& detJ) : detJ(detJ) {
    const A2D::Mat<T, dim, dim>& Jxi = geo.template get<0>().get_grad();
    A2D::Vec<T, dim> x, y, nA;
    if constexpr (dim == 2) {
      // Find the nA = vector of the distorted bound
      nA(0) = Jxi(0, 0);
      nA(1) = Jxi(1, 0);
      detJ = std::sqrt(nA(0) * nA(0) + nA(1) * nA(1));
    } else if constexpr (dim == 3) {
      // Find the nA = (Area) * normal direction
      x(0) = Jxi(0, 0);  // d(x, y, z)/dxi
      x(1) = Jxi(1, 0);
      x(2) = Jxi(2, 0);

      y(0) = Jxi(0, 1);  // d(x, y, z)/deta
      y(1) = Jxi(1, 1);
      y(2) = Jxi(2, 1);

      nA(0) = x(1) * y(2) - x(2) * y(1);
      nA(1) = x(2) * y(0) - x(0) * y(2);
      nA(2) = x(0) * y(1) - x(1) * y(0);

      detJ = std::sqrt(nA(0) * nA(0) + nA(1) * nA(1) + nA(2) * nA(2));
    }
  }

  template <class FiniteElementSpace>
  void transform(const FiniteElementSpace& in, FiniteElementSpace& out) {
    in.transform(detJ, J, Jinv, out);
  }

  template <class FiniteElementSpace>
  void rtransform(const FiniteElementSpace& in, FiniteElementSpace& out) {
    in.rtransform(detJ, J, Jinv, out);
  }

 private:
  // Determinant of the Jacobian transformation
  T& detJ;

  // J with the normal direction added
  A2D::Mat<T, dim_surf, dim_surf> J, Jinv;  // TODO
};

}  // namespace A2D

#endif  //  A2D_FE_MAPPING_H
