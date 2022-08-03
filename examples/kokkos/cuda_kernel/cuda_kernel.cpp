#include "Kokkos_Core.hpp"
#include "a2dobjs.h"
#include "multiarray.h"
#include "parallel.h"

class Ops {
 public:
  template <typename ArrayType, typename FunctorType>
  static void launch_kernel_with_functor(int N, ArrayType& input,
                                         ArrayType& output,
                                         const FunctorType& fun) {
    A2D::parallel_for(
        N, A2D_LAMBDA(int i) { fun(input(i, 0), output(i, 0)); });
  }
};

template <typename Ops, typename ArrayType>
void wrapper_function(int N, ArrayType& input, ArrayType& output) {
  Ops::template launch_kernel_with_functor(
      N, input, output, A2D_LAMBDA(double& x, double& y) { y = 4.2 * x; });
}

void use_a2d_multiarray() {
  Kokkos::initialize();
  {
#ifdef KOKKOS_ENABLE_CUDA
    printf("KOKKOS_ENABLE_CUDA is defined\n");
#endif
    // Vector size
    int N = 1 << 20;  // 1M

    // Allocate vectors on CUDA unified memory
    using VecLayout = A2D::CLayout<1>;
    using Vec = A2D::MultiArray<double, VecLayout>;

    VecLayout x_layout(N), y_layout(N);
    Vec x(x_layout), y(y_layout);

    // Initialize vector x on host
    for (int i = 0; i < N; i++) {
      x(i, 0) = 4.2;
    }

    // Initialize vector y on device: y(i) = x(i)
    // using range_policy = Kokkos::RangePolicy<Kokkos::Cuda>;
    wrapper_function<Ops>(N, x, y);

    // Sync
    Kokkos::fence();

    // Access vector y from host
    for (int i = 0; i < 5; i++) {
      printf("y(%d, 0) = %.2f\n", i, y(i, 0));  // Expect 4.2, get 0.0
    }
  }
  Kokkos::finalize();
}

int main() { use_a2d_multiarray(); }