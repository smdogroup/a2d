#include <omp.h>

namespace A2D {

template <class FunctorType>
void parallel_for(const A2D::index_t N, const FunctorType& func) {
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    func(i);
  }
}

}  // namespace A2D
