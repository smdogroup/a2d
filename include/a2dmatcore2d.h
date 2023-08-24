#ifndef A2D_MAT_CORE_2D_H
#define A2D_MAT_CORE_2D_H

#include "a2dobjs.h"

namespace A2D {

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2SymmMultCore(const AMatType& A,
                                             const BMatType& B, CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2SymmMultScaleCore(T scale, const AMatType& A,
                                                  const BMatType& B,
                                                  CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2SymmMultAddCore(const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2SymmMultSubCore(const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2SymmMultAddScaleCore(T scale, const AMatType& A,
                                                     const BMatType& B,
                                                     CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatMultCore(const AMatType& A,
                                            const BMatType& B, CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatMultScaleCore(T scale, const AMatType& A,
                                                 const BMatType& B,
                                                 CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatMultAddCore(const AMatType& A,
                                               const BMatType& B, CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatMultSubCore(const AMatType& A,
                                               const BMatType& B, CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatMultAddScaleCore(T scale, const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatTransMultCore(const AMatType& A,
                                                 const BMatType& B,
                                                 CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatTransMultScaleCore(T scale,
                                                      const AMatType& A,
                                                      const BMatType& B,
                                                      CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatTransMultAddCore(const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatTransMultSubCore(const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) -= (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Symm2x2MatTransMultAddScaleCore(T scale,
                                                         const AMatType& A,
                                                         const BMatType& B,
                                                         CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) += scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2SymmMultCore(const AMatType& A,
                                            const BMatType& B, CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2SymmMultScaleCore(T scale, const AMatType& A,
                                                 const BMatType& B,
                                                 CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2SymmMultAddCore(const AMatType& A,
                                               const BMatType& B, CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2SymmMultSubCore(const AMatType& A,
                                               const BMatType& B, CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2SymmMultAddScaleCore(T scale, const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2SymmMultCore(const AMatType& A,
                                                 const BMatType& B,
                                                 CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2SymmMultScaleCore(T scale,
                                                      const AMatType& A,
                                                      const BMatType& B,
                                                      CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2SymmMultAddCore(const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2SymmMultSubCore(const AMatType& A,
                                                    const BMatType& B,
                                                    CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2SymmMultAddScaleCore(T scale,
                                                         const AMatType& A,
                                                         const BMatType& B,
                                                         CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) += scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatMultCore(const AMatType& A, const BMatType& B,
                                           CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatMultScaleCore(T scale, const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatMultAddCore(const AMatType& A,
                                              const BMatType& B, CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatMultSubCore(const AMatType& A,
                                              const BMatType& B, CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatMultAddScaleCore(T scale, const AMatType& A,
                                                   const BMatType& B,
                                                   CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatTransMultCore(const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) = (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatTransMultScaleCore(T scale, const AMatType& A,
                                                     const BMatType& B,
                                                     CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) = scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatTransMultAddCore(const AMatType& A,
                                                   const BMatType& B,
                                                   CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) += (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatTransMultSubCore(const AMatType& A,
                                                   const BMatType& B,
                                                   CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) -= (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) -= (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void Mat2x2MatTransMultAddScaleCore(T scale,
                                                        const AMatType& A,
                                                        const BMatType& B,
                                                        CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(0, 1) * B(0, 1));
  C(0, 1) += scale * (A(0, 0) * B(1, 0) + A(0, 1) * B(1, 1));
  C(1, 0) += scale * (A(1, 0) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += scale * (A(1, 0) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatMultCore(const AMatType& A,
                                                const BMatType& B,
                                                CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) = (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatMultScaleCore(T scale, const AMatType& A,
                                                     const BMatType& B,
                                                     CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) = scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) = scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatMultAddCore(const AMatType& A,
                                                   const BMatType& B,
                                                   CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) += (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatMultSubCore(const AMatType& A,
                                                   const BMatType& B,
                                                   CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) -= (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) -= (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatMultAddScaleCore(T scale,
                                                        const AMatType& A,
                                                        const BMatType& B,
                                                        CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(1, 0));
  C(0, 1) += scale * (A(0, 0) * B(0, 1) + A(1, 0) * B(1, 1));
  C(1, 0) += scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(1, 0));
  C(1, 1) += scale * (A(0, 1) * B(0, 1) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatTransMultCore(const AMatType& A,
                                                     const BMatType& B,
                                                     CMatType& C) {
  C(0, 0) = (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1));
  C(0, 1) = (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1));
  C(1, 0) = (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatTransMultScaleCore(T scale,
                                                          const AMatType& A,
                                                          const BMatType& B,
                                                          CMatType& C) {
  C(0, 0) = scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1));
  C(0, 1) = scale * (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1));
  C(1, 0) = scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) = scale * (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatTransMultAddCore(const AMatType& A,
                                                        const BMatType& B,
                                                        CMatType& C) {
  C(0, 0) += (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1));
  C(0, 1) += (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1));
  C(1, 0) += (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatTransMultSubCore(const AMatType& A,
                                                        const BMatType& B,
                                                        CMatType& C) {
  C(0, 0) -= (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1));
  C(0, 1) -= (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1));
  C(1, 0) -= (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) -= (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1));
}

template <typename T, class AMatType, class BMatType, class CMatType>
KOKKOS_FUNCTION void MatTrans2x2MatTransMultAddScaleCore(T scale,
                                                             const AMatType& A,
                                                             const BMatType& B,
                                                             CMatType& C) {
  C(0, 0) += scale * (A(0, 0) * B(0, 0) + A(1, 0) * B(0, 1));
  C(0, 1) += scale * (A(0, 0) * B(1, 0) + A(1, 0) * B(1, 1));
  C(1, 0) += scale * (A(0, 1) * B(0, 0) + A(1, 1) * B(0, 1));
  C(1, 1) += scale * (A(0, 1) * B(1, 0) + A(1, 1) * B(1, 1));
}

}  // namespace A2D

#endif  // A2D_MAT_CORE_2D_H
