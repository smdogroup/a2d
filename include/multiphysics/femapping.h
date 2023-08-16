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

 private:
  const A2D::Mat<T, dim, dim>& J;
  A2D::Mat<T, dim, dim> Jinv;
  T& detJ;
};

/**
 * @brief surface transformation - this needs to be fixed TODO: how?
 */
template <typename T, index_t D>
class SurfaceMapping {
 public:
  static const index_t dim = D;
  static const index_t dim_surf = D - 1;

  template <class FiniteElementGeometry>
  SurfaceMapping(const FiniteElementGeometry& geo, T& detJ) : detJ(detJ) {
    if constexpr (dim == 2) {
      detJ = 1.0;
    } else if constexpr (dim == 3) {
      const A2D::Mat<T, dim, dim>& Jxi = geo.template get<0>().get_grad();

      // Find the nA = (Area) * normal direction
      A2D::Vec<T, dim> x, y, nA;
      x(0) = Jxi(0, 0);
      x(1) = Jxi(1, 0);
      x(2) = Jxi(2, 0);

      y(0) = Jxi(0, 1);
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
  A2D::Mat<T, dim_surf, dim_surf> J, Jinv;  // TODO: where are these computed?
};

}  // namespace A2D

#endif  //  A2D_FE_MAPPING_H
