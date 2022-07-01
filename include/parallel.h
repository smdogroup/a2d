#ifndef A2D_PARALLEL_H
#define A2D_PARALLEL_H

#include <omp.h>

namespace A2D {

template <class FunctorType>
void parallel_for(const A2D::index_t N, const FunctorType& func) {
#pragma omp parallel for
  for (index_t i = 0; i < N; ++i) {
    func(i);
  }
}

template <typename T, class FunctorType>
T parallel_reduce(const A2D::index_t N, const FunctorType& func) {
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
