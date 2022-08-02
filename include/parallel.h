#ifndef A2D_PARALLEL_H
#define A2D_PARALLEL_H

#include <omp.h>

#include <string>

#include "a2dobjs.h"

#ifdef A2D_USE_KOKKOS
#include "a2dkokkos.h"
#endif

namespace A2D {

template <typename IdxType, class FunctorType>
void parallel_for(const IdxType N, const FunctorType& func) {
#ifdef A2D_USE_KOKKOS
  Kokkos::parallel_for(N, func);
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence();
#endif
#else
#pragma omp parallel for
  for (index_t i = 0; i < N; ++i) {
    func(i);
  }
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
