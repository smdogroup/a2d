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
class VolumeMapping {
 public:
  template <class FiniteElementGeometry>
  VolumeMapping(const FiniteElementGeometry& geo, T& detJ)
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
 * @brief 3D surface transformation
 *
 */
template <typename T>
class SurfaceTransform3D {
 public:
  static const index_t dim3 = 3;
  static const index_t dim2 = 2;

  template <class FiniteElementGeometry>
  SurfaceTransform3D(const FiniteElementGeometry& geo, T& detJ)
      : Jxi(geo.template get<0>().get_grad()), detJ(detJ) {
    // // Find the nA = (Area) * normal direction
    // // A2D::Vec<T, dim3> x0, x1, nA;
    // x0(0) = Jxi(0, 0);
    // x0(1) = Jxi(1, 0);
    // x0(2) = Jxi(2, 0);

    // x1(0) = Jxi(0, 1);
    // x1(1) = Jxi(1, 1);
    // x1(2) = Jxi(2, 1);
    // A2D::Vec3Cross(x1, x2, nA);

    // // Normalize the vector so we just have the normal
    // A2D::Vec<T, dim3> n;
    // A2D::Vec3Normalize(nA, n);

    // // Now initialize the Jacobian transformation
    // J(0, 0) = Jxi(0, 0);
    // J(1, 0) = Jxi(1, 0);
    // J(2, 0) = Jxi(2, 0);

    // J(0, 1) = Jxi(0, 1);
    // J(1, 1) = Jxi(1, 1);
    // J(2, 1) = Jxi(2, 1);

    // J(0, 2) = n(0);
    // J(1, 2) = n(1);
    // J(2, 2) = n(2);

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
  // The Jacobian transform is a 3 x 2
  const A2D::Mat<T, dim3, dim2>& Jxi;

  // Determinant of the Jacobian transformation
  T& detJ;

  // J with the normal direction added
  A2D::Mat<T, dim3, dim3> J, Jinv;
};

}  // namespace A2D

#endif  //  A2D_FE_MAPPING_H
