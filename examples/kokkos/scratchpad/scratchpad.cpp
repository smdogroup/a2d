#include <iostream>

#include "Kokkos_Core.hpp"

using namespace std;

void test_axpy(int argc, char* argv[]) {
  using T = double;
  using I = long int;
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
          "axpy", N, KOKKOS_LAMBDA(I index) {
            y_device(i) = y_device(index) + alpha * x_device(i);
          });
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

void test_matvec(int argc, char* argv[]) {
  using T = double;
  using I = long int;
  using vec_device_t = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::CudaSpace>;
  using mat_device_t =
      Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::CudaSpace>;
  Kokkos::initialize();
  {
    if (argc < 3) {
      std::printf("compute y <- Ax + y, where A is M-by-N matrix\n");
      std::printf("usage:  ./scratchpad M N\n");
      exit(-1);
    }

    // Allocate A, x and y on host
    I M = pow(2, atoi(argv[1]));
    I N = pow(2, atoi(argv[2]));

    mat_device_t A_device("A_device", M, N);
    vec_device_t x_device("x", N);
    vec_device_t y_device("y", M);

    double A_mbytes = (double)M * N * sizeof(T) / 1024 / 1024;
    std::printf("Allocating A on host, size: %.2f MB\n", A_mbytes);

    mat_device_t::HostMirror A = Kokkos::create_mirror_view(A_device);
    vec_device_t::HostMirror x = Kokkos::create_mirror_view(x_device);
    vec_device_t::HostMirror y = Kokkos::create_mirror_view(y_device);

    for (I i = 0; i < M; i++) {
      for (I j = 0; j < N; j++) {
        A(i, j) = 0.01;
      }
    }
    for (I i = 0; i < N; i++) {
      x(i) = 2.3;
    }

    for (I i = 0; i < M; i++) {
      y(i) = 0.1;
    }

    // Copy value to device
    Kokkos::deep_copy(A_device, A);
    Kokkos::deep_copy(x_device, x);
    Kokkos::deep_copy(y_device, y);

    // Perform 100 axpy operations on device and time
    Kokkos::Timer timer;
    int repeat = 100;
    for (int i = 0; i < repeat; i++) {
      Kokkos::parallel_for(
          "matvec", M, KOKKOS_LAMBDA(I index) {
            double row_sum = 0.0;
            for (I j = 0; j < N; j++) {
              row_sum += A_device(index, j) * x_device(j);
            }
            y_device(index) = row_sum;
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

    // Copy results back to host
    Kokkos::deep_copy(y, y_device);

    // Check maximum error
    T max_err = T(0);
    T val;
    for (I i = 0; i < M; i++) {
      val = fabs(repeat * 0.01 * 2.3 * N + 0.1 - y(i));
      if (val > max_err) {
        max_err = val;
      }
    }

    printf("Maximum error: %20.10e\n", max_err);
  }
  Kokkos::finalize();
}

int main(int argc, char* argv[]) {
  // test_axpy(argc, argv);
  test_matvec(argc, argv);
}