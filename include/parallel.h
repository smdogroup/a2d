#ifndef A2D_PARALLEL_H
#define A2D_PARALLEL_H

#include <omp.h>

#include <string>

#include "a2ddefs.h"

namespace A2D {

template <typename IdxType, class FunctorType>
void parallel_for(const IdxType N, const FunctorType& func) {
  Kokkos::parallel_for(N, func);
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence();
#endif
}

template <typename T, class FunctorType>
T parallel_reduce(const index_t N, const FunctorType& func) {
  T sum = 0.0;

#pragma omp parallel
  {
    T part_sum = 0.0;

#pragma omp for
    for (index_t i = 0; i < N; i++) {
      part_sum += func(i);
    }

#pragma omp critical
    { sum += part_sum; }
  }

  return sum;
}

}  // namespace A2D

#endif  // A2D_PARALLEL_H
