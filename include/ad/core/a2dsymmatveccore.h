#ifndef A2D_SYMMAT_VEC_CORE_H
#define A2D_SYMMAT_VEC_CORE_H

#include "../../a2ddefs.h"

namespace A2D {

    template <typename T, bool additive = false>
    A2D_FUNCTION void SymMatVecCore3x3(const T A[], const T x[], T b[]) {
        // compute b = A*x
        if constexpr (additive) {
            b[0] += A[0] * x[0] + A[1] * x[1] + A[3] * x[2];
            b[1] += A[1] * x[0] + A[2] * x[1] + A[4] * x[2];
            b[2] += A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
        } else {
            b[0] = A[0] * x[0] + A[1] * x[1] + A[3] * x[2];
            b[1] = A[1] * x[0] + A[2] * x[1] + A[4] * x[2];
            b[2] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
        }
        
    }

    template <typename T, bool additive = false>
    A2D_FUNCTION void SymMatVecCoreScale3x3(const T& scale, const T A[], const T x[], T b[]) {
        // compute b = A * x
       if constexpr (additive) {
            b[0] += scale * (A[0] * x[0] + A[1] * x[1] + A[3] * x[2]);
            b[1] += scale * (A[1] * x[0] + A[2] * x[1] + A[4] * x[2]);
            b[2] += scale * (A[3] * x[0] + A[4] * x[1] + A[5] * x[2]);
        } else {
            b[0] = scale * (A[0] * x[0] + A[1] * x[1] + A[3] * x[2]);
            b[1] = scale * (A[1] * x[0] + A[2] * x[1] + A[4] * x[2]);
            b[2] = scale * (A[3] * x[0] + A[4] * x[1] + A[5] * x[2]);
        }
    }
}

#endif // A2D_SYMMAT_VEC_CORE_H