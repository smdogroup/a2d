#include <Kokkos_Core.hpp>

void launch_cuda_kernel() {
  Kokkos::initialize();
  {
#ifdef KOKKOS_ENABLE_CUDA
    printf("KOKKOS_ENABLE_CUDA is defined\n");
#endif
    // Vector size
    int N = 1 << 20;  // 1M

    // Allocate vectors on CUDA unified memory
    using ViewType =
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::CudaUVMSpace>;
    ViewType x("x", N);
    ViewType y("y", N);

    // Initialize vector x on host
    for (int i = 0; i < N; i++) {
      x(i) = 4.2;
    }

    // Initialize vector y on device: y(i) = x(i)
    using range_policy = Kokkos::RangePolicy<Kokkos::Cuda>;
    Kokkos::parallel_for(
        "assign_value", range_policy(0, N),
        KOKKOS_LAMBDA(int i) { y(i) = x(i); });

    // Sync
    Kokkos::fence();

    // Access vector y from host
    for (int i = 0; i < 5; i++) {
      printf("y(%d) = %.2f\n", i, y(i));  // Expect 4.2, get 0.0
    }
  }
  Kokkos::finalize();
}

class Foo {
 public:
  using ViewType =
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::CudaUVMSpace>;
  KOKKOS_INLINE_FUNCTION Foo(int dim0) { view = ViewType("Foo", dim0); }
  ViewType view;
};

void inline_constructor() {
  Kokkos::initialize();
  { Foo f(1000); }
  Kokkos::finalize();
}

void random_on_kernel() {
  Kokkos::initialize();
  {
#ifdef KOKKOS_ENABLE_CUDA
    printf("KOKKOS_ENABLE_CUDA is defined\n");
#endif
    // Vector size
    int N = 1 << 20;  // 1M

    // Allocate vectors on CUDA unified memory
    using ViewType =
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::CudaUVMSpace>;
    ViewType x("x", N);
    ViewType y("y", N);

    // Initialize vector y on device: y(i) = rand()
    using range_policy = Kokkos::RangePolicy<Kokkos::Cuda>;
    Kokkos::parallel_for(
        "assign_value", range_policy(0, N),
        KOKKOS_LAMBDA(int i) { y(i) = x(i); });

    // Sync
    Kokkos::fence();

    // Access vector y from host
    for (int i = 0; i < 20; i++) {
      printf("y(%d) = %.2f\n", i, y(i));  // Expect 4.2, get 0.0
    }
  }
  Kokkos::finalize();
}

int main() {
  // launch_cuda_kernel();
  // inline_constructor();
}