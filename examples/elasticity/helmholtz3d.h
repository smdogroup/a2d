#ifndef HELMHOLTZ_3D_H
#define HELMHOLTZ_3D_H

#include "a2dtmp.h"
#include "basis3d.h"

template <class Basis>
class HelmholtzPDE {
 public:
  static const int NUM_VARS = 1;

  template <typename T, class QuadPointModelDataArray, class QuadPointDetJArray,
            class QuadPointJacobianArray, class QuadPointSolutionArray,
            class QuadPointGradientArray, class ElementResidualArray>
  static void residuals(QuadPointModelDataArray& data, QuadPointDetJArray& detJ,
                        QuadPointJacobianArray& Jinv, QuadPointSolutionArray& U,
                        QuadPointGradientArray& Uxi,
                        ElementResidualArray& res) {
    Basis::template residuals<T, HelmholtzPDE<Basis>::Impl>(data, detJ, Jinv,
                                                            Uxi, res);
    Basis::template residuals<T, HelmholtzPDE<Basis>::Impl>(data, detJ, U, res);
  }

  class Impl {
   public:
    template <typename T, class IdxType, class QuadPointData>
    static T compute_energy(IdxType i, IdxType j, QuadPointData& data,
                            A2D::Mat<T, 3, 3>& Jinv, A2D::Vec<T, 3>& U,
                            A2D::Mat<T, 3, 3>& Uxi) {
      // A2D::Mat3x3VecMult(Uxi, Jinv, Ux);

      return 0.5
    }
  };
};

#endif  // HELMHOLTZ_3D_H