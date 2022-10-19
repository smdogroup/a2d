#include <iostream>

#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include "Kokkos_UnorderedMap.hpp"

using namespace std;
using T = double;
using I = long int;

#if 0
void test_axpy(int argc, char* argv[]) {
  using ViewDevice_t = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::CudaSpace>;
  Kokkos::initialize(argc, argv);
  {
    if (argc == 1) {
      std::printf("axpy\nusage: ./scratchpad N\n");
      exit(-1);
    }

    // Allocate x and y on host
    I N = pow(2, atoi(argv[1]));
    I bytes = N * sizeof(T);
    double Mbytes = (double)bytes / 1024 / 1024;

    std::printf("Allocating x on host, size: %.2f MB\n", Mbytes);
    std::printf("Allocating y on host, size: %.2f MB\n", Mbytes);

    // Allocate x and y on device
    ViewDevice_t x_device("x_device", N);
    ViewDevice_t y_device("y_device", N);

    // Initialize x and y on host
    ViewDevice_t::HostMirror x = Kokkos::create_mirror_view(x_device);
    ViewDevice_t::HostMirror y = Kokkos::create_mirror_view(y_device);

    for (I i = 0; i < N; i++) {
      x(i) = 2.4;
      y(i) = 0.1;
    }

    // Copy value to device
    Kokkos::deep_copy(x_device, x);
    Kokkos::deep_copy(y_device, y);

    // Perform 100 axpy operations on device and time
    Kokkos::Timer timer;
    int repeat = 100;
    T alpha = 0.01;
    for (int i = 0; i < repeat; i++) {
      Kokkos::parallel_for(
          "axpy", N,
          KOKKOS_LAMBDA(I index) { y_device(index) += alpha * x_device(i); });
      Kokkos::fence();
    }

    double elapse = timer.seconds();
    printf("averaged time: %.8f ms\n", elapse * 1e3);

    // Compute bandwidth:
    // x is read once, y is read once and written once
    double bandwidth = 3 * Mbytes / 1024 * repeat / elapse;
    printf("averaged bandwidth: %.8f GB/s\n", bandwidth);

    // Copy results back to host
    Kokkos::deep_copy(y, y_device);

    // Check maximum error
    T max_err = T(0);
    T val;
    for (I i = 0; i < N; i++) {
      val = fabs(repeat * alpha * 2.4 + 0.1 - y(i));
      if (val > max_err) {
        max_err = val;
      }
    }

    printf("Maximum error: %20.10e\n", max_err);
  }
  Kokkos::finalize();
}

template <class ExecSpace, class AType, class xType, class yType>
void _launch_matvec(AType A, xType x, yType y, int repeat = 100) {
  // Get extents
  I M = A.extent(0);
  I N = A.extent(1);

  // Perform 100 axpy operations on device and time
  Kokkos::Timer timer;
  for (int i = 0; i < repeat; i++) {
    Kokkos::parallel_for(
        "matvec", Kokkos::RangePolicy<ExecSpace>(0, M), KOKKOS_LAMBDA(I index) {
          double row_sum = 0.0;
          for (I j = 0; j < N; j++) {
            row_sum += A(index, j) * x(j);
          }
          y(index) += row_sum;
        });
    Kokkos::fence();
  }

  double elapse = timer.seconds();
  printf("averaged time: %.8f ms\n", elapse * 1e3);

  // Calculate bandwidth.
  // Each matrix A row (each of length M) is read once.
  // The x vector (of length N) is read M times.
  // The y vector (of length N) is read once and write once
  double Mbytes = double(sizeof(T) * (2 * M * N + M)) / 1024 / 1024;
  double bandwidth = Mbytes / 1024 * repeat / elapse;
  printf("averaged bandwidth: %.8f GB/s\n", bandwidth);
}

void test_matvec(int argc, char* argv[]) {
  using vec_device_t = Kokkos::View<T*, Kokkos::LayoutLeft, Kokkos::CudaSpace>;
  using mat_device_t = Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::CudaSpace>;
  Kokkos::initialize();
  {
    if (argc < 2) {
      std::printf("compute y <- Ax + y, where A is M-by-N matrix\n");
      std::printf("usage:  ./scratchpad M [N]\n");
      exit(-1);
    }

    // Allocate A, x and y on host
    I M = pow(2, atoi(argv[1]));
    I N;
    if (argc >= 3) {
      N = pow(2, atoi(argv[2]));
    } else {
      N = M;
    }

    mat_device_t A_device("A_device", M, N);
    vec_device_t x_device("x", N);
    vec_device_t y_device("y", M);

    double A_mbytes = (double)M * N * sizeof(T) / 1024 / 1024;
    std::printf("\nAllocating A on host, size: %.2f MB\n", A_mbytes);

    mat_device_t::HostMirror A = Kokkos::create_mirror_view(A_device);
    vec_device_t::HostMirror x = Kokkos::create_mirror_view(x_device);
    vec_device_t::HostMirror y = Kokkos::create_mirror_view(y_device);

    for (I i = 0; i < M; i++) {
      for (I j = 0; j < N; j++) {
        A(i, j) = 1.0;
      }
    }
    for (I i = 0; i < N; i++) {
      x(i) = 1.0;
    }

    for (I i = 0; i < M; i++) {
      y(i) = 0.0;
    }

    // Copy value to device
    Kokkos::deep_copy(A_device, A);
    Kokkos::deep_copy(x_device, x);
    Kokkos::deep_copy(y_device, y);

    int repeat = 100;

    // Run matvec on device
    printf("\n====== CUDA ======\n\n");
    _launch_matvec<Kokkos::CudaSpace::execution_space>(A_device, x_device,
                                                       y_device, repeat);

    // Copy results back to host
    Kokkos::deep_copy(y, y_device);

    // Check maximum error
    T max_err = T(0);
    T val;
    for (I i = 0; i < M; i++) {
      val = fabs((double)repeat * (double)N - y(i));
      if (val > max_err) {
        max_err = val;
      }
    }
    printf("Maximum error: %20.10e\n", max_err);

    // Run matvec on host
    printf("\n====== OpenMP ======\n\n");
    _launch_matvec<Kokkos::HostSpace::execution_space>(A, x, y, repeat);

    // Check maximum error
    max_err = T(0);
    for (I i = 0; i < M; i++) {
      val = fabs(2 * (double)repeat * (double)N - y(i));
      if (val > max_err) {
        max_err = val;
      }
    }
    printf("Maximum error: %20.10e\n", max_err);
  }
  Kokkos::finalize();
}

template <typename T>
struct COO {
  T x;
  T y;
};

void test_unordered_set() {
  Kokkos::initialize();
  {
    using MemSpace = Kokkos::HostSpace;
    using ExecSpace = Kokkos::Cuda;
    using RangePolicy = Kokkos::RangePolicy<ExecSpace>;

    int repeat = 20;
    int set_capacity = 1;
    int rand_max = 10;
    Kokkos::UnorderedMap<COO<int>, void, ExecSpace> node_set;
    node_set.rehash(set_capacity);

    // Add value in parallel
    int fail = 0;
    Kokkos::parallel_reduce(
        RangePolicy(0, repeat),
        KOKKOS_LAMBDA(int i, int& error) {
          int x = i % rand_max;
          int y = i % rand_max;
          auto result = node_set.insert(COO<int>{x, y});
          error += result.failed();
        },
        fail);

    printf("total fail: %d\n", fail);

    // Print value in parallel
    int total = 0;
    Kokkos::parallel_reduce(
        RangePolicy(0, node_set.capacity()),
        KOKKOS_LAMBDA(I i, int& _total) {
          _total += 1;
          if (node_set.valid_at(i)) {
            auto key = node_set.key_at(i);
            printf("[%2d] (%d, %d)\n", (int)i, key.x, key.y);
          }
        },
        total);

    printf("total number of threads: %d\n", total);
  }
  Kokkos::finalize();
}
#endif

#define N1 3
#define N2 4
using SomeViewType = Kokkos::View<double* [N1][N2]>;

SomeViewType create_view(int n) {
  SomeViewType view("view", n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 4; k++) {
        view(i, j, k) = (i + 1) * (j + 1) * (k + 1);
      }
    }
  }
  return view;
}

void test_subview() {
  Kokkos::initialize();
  {
    int n = 4;
    SomeViewType array = create_view(n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < N1; j++) {
        for (int k = 0; k < N2; k++) {
          printf("array(%d, %d, %d) = %.2f\n", i, j, k, array(i, j, k));
        }
      }
    }

    printf("rank = %d\n", decltype(array)::rank);
    auto slice = Kokkos::subview(array, 2, Kokkos::ALL, Kokkos::ALL);
    printf("is layout of slice LayoutRight? %d\n",
           std::is_same<decltype(slice)::array_layout,
                        Kokkos::LayoutRight>::value);
    printf(
        "is layout of slice LayoutLeft? %d\n",
        std::is_same<decltype(slice)::array_layout, Kokkos::LayoutLeft>::value);

    for (int j = 0; j < N1; j++) {
      for (int k = 0; k < N2; k++) {
        printf("slice(%d, %d) = %.2f\n", j, k, slice(j, k));
      }
    }
  }
  Kokkos::finalize();
}

void test_sort() {
  Kokkos::initialize();
  {
    int n = 10;
    Kokkos::View<double*> array;
    printf("bool(array.data()): %d\n", bool(array.data()));
    printf("array.is_allocated(): %d\n", array.is_allocated());

    array = Kokkos::View<double*>("array", n);
    printf("bool(array.data()): %d\n", bool(array.data()));
    printf("array.is_allocated(): %d\n", array.is_allocated());
    for (int i = 0; i != n; i++) {
      array(i) = (double)std::rand() / (double)RAND_MAX;
    }

    for (int i = 0; i != n; i++) {
      printf("array(%d) = %.5f\n", i, array(i));
    }

    // std::sort(array.data(), array.data() + n);
    Kokkos::sort(array, 0, 5);

    printf("sorted:\n");
    for (int i = 0; i != n; i++) {
      printf("array(%d) = %.5f\n", i, array(i));
    }
  }
  Kokkos::finalize();
}

int main(int argc, char* argv[]) {
  // test_axpy(argc, argv);
  // test_matvec(argc, argv);
  // test_unordered_set();
  // test_subview();
  test_sort();
}