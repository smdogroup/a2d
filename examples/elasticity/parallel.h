#include <omp.h>

namespace A2D {

template <class FunctorType>
void parallel_for(const A2D::index_t N, const FunctorType& func) {
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    func(i);
  }
}

template <typename T, class FunctorType>
T parallel_reduce(const A2D::index_t N, const FunctorType& func) {
  T sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < N; ++i) {
    sum += func(i);
  }
  return sum;
}

}  // namespace A2D
