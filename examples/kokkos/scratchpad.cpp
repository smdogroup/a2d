#include <iostream>

#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"
#include "Kokkos_StdAlgorithms.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "a2dobjs.h"
#include "array.h"

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

void test_is_same_layout() {
  Kokkos::initialize();
  {
    using Array_t = Kokkos::View<T* [5], Kokkos::LayoutLeft>;
    Array_t array = Array_t("arr", 20);
    auto sub = Kokkos::subview(array, 0, Kokkos::ALL);
    bool is_flayout =
        std::is_same<decltype(sub)::array_layout, Kokkos::LayoutStride>::value;
    printf("is_flayout: %d\n", is_flayout);
  }
  Kokkos::finalize();
}

void test_complex() {
  using cT = A2D_complex_t<double>;
  cT x = cT(2.3, 4.5);
  cT s = A2D::sqrt(x);
  cT e = A2D::exp(x);
  cT l = A2D::log(x);
  std::printf("x      : %.5f + %.5fj\n", x.real(), x.imag());
  std::printf("sqrt(x): %.5f + %.5fj\n", s.real(), s.imag());
  std::printf("exp(x):  %.5f + %.5fj\n", e.real(), e.imag());
  std::printf("log(x):  %.5f + %.5fj\n", l.real(), l.imag());
}

template <class T>
struct is_complex : public std::false_type {};
template <class T>
struct is_complex<std::complex<T>> : public std::true_type {};

void test_is_complex() {
  using T1 = std::complex<double>;
  using T2 = std::complex<int>;
  using T3 = double;
  using T4 = std::complex<std::complex<float>>;

  std::cout << is_complex<T1>::value << std::endl;  // should be 1
  std::cout << is_complex<T2>::value << std::endl;  // should be 1
  std::cout << is_complex<T3>::value << std::endl;  // should be 0
  std::cout << is_complex<T4>::value << std::endl;  // should be 1
}

class BSRMatToy {
 public:
  BSRMatToy() {}
  void initialize() { perm = A2D::IdxArray1D_t("perm", 20); }
  void yell() {
    if (!perm.is_allocated()) {
      std::cout << "perm is not allocated\n";
    } else {
      std::cout << "perm is allocated\n";
    }
  }

 private:
  A2D::IdxArray1D_t perm;
};

void test_view_is_allocated() {
  BSRMatToy b;
  b.yell();
  b.initialize();
  b.yell();
  return;
}

void allocate_and_populate(A2D::IdxArray1D_t& cols) {
  cols = A2D::IdxArray1D_t("cols", 20);
  cols(0) = 233;
  return;
}

void test_smart_pointer_behavior() {
  A2D::IdxArray1D_t cols;
  std::cout << cols(0) << "\n";
  std::cout << "is allocated? " << cols.is_allocated() << "\n";
  allocate_and_populate(cols);
  std::cout << cols(0) << "\n";
  std::cout << "is allocated? " << cols.is_allocated() << "\n";
  return;
}

void test_copy() {
  using T = double;
  A2D::MultiArrayNew<T* [5]> array1("array1", 20);
  array1(0, 0) = 4.2;

  A2D::MultiArrayNew<T* [5]> array2;

  array2 = array1;
  std::cout << array2(0, 0) << "\n";
}

void test_parallel_for() {
  constexpr int N = 10;
  using T = double;
  A2D::MultiArrayNew<T*> array("array", N);
  A2D::MultiArrayNew<T*> array2("array", N);

  Kokkos::parallel_for(
      "Loop", N, KOKKOS_LAMBDA(const int i) { array(i) = 4.2; });

  for (int i = 0; i < N; i++) {
    std::cout << array(i) << "\n";
  }
}

class ElemVec {
 public:
  using T = double;
  ElemVec(A2D::index_t n) : array("array", n) {}

  T get(A2D::index_t i) const { return array(i); }
  void set(A2D::index_t i, const T val) const { array(i) = val; }

 private:
  A2D::MultiArrayNew<T*> array;
};

// class ElemVec2 {
//  public:
//   using T = double;
//   ElemVec2(A2D::index_t n) : array(n) {}

//   T get(A2D::index_t i) const { return array[i]; }
//   void set(A2D::index_t i, const T val) const { array[i] = val; }

//  private:
//   std::vector<T> array;
// };

class ElemVec3 {
 public:
  using T = double;
  using array_t = A2D::MultiArrayNew<T*>;

  ElemVec3(array_t& array) : array(array) {}

  T get(A2D::index_t i) const { return array(i); }
  void set(A2D::index_t i, const T val) const { array(i) = val; }

 private:
  array_t& array;
};

void test_modify_view_from_const_lambda() {
  const ElemVec elem_vec(42);
  elem_vec.set(0, 2.3);
  std::cout << elem_vec.get(0) << "\n";
}

struct Alpha {
  constexpr static T value = 0.01;
  constexpr static T value_array[] = {0.01, 0.02};

  KOKKOS_INLINE_FUNCTION static const T* get_alpha_array() {
    return Alpha::value_array;
  };
};

KOKKOS_INLINE_FUNCTION T get_alpha() { return Alpha::value; };

void test_cuda_axpy(int argc, char* argv[]) {
#ifdef KOKKOS_ENABLE_CUDA
  using ViewDevice_t = Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::CudaSpace>;
  if (argc == 1) {
    std::printf("axpy\nusage: ./scratchpad N, where mat size = 2^N\n");
    return;
  }

  // Allocate x and y on host
  I N = pow(2, atoi(argv[1]));
  I bytes = N * sizeof(T);
  double Mbytes = (double)bytes / 1024 / 1024;

  std::printf("Allocating x[%d] on host, size: %.2f MB\n", N, Mbytes);
  std::printf("Allocating y[%d] on host, size: %.2f MB\n", N, Mbytes);

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
  Kokkos::fence();

  // Perform 100 axpy operations on device and time
  Kokkos::Timer timer;
  int repeat = 100;
  T alpha = 0.01;
  for (int i = 0; i < repeat; i++) {
    Kokkos::parallel_for(
        "axpy", N, KOKKOS_LAMBDA(const I index) {
          // T alpha = Alpha::get_alpha_array()[0];  // This doesn't work
          T alpha = Alpha::value_array[0];  // This works
          y_device(index) += alpha * x_device(index);
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
  Kokkos::fence();

  // Print values
  int len = 10;
  if (N < len) {
    len = N;
  }

  for (I i = 0; i < len; i++) {
    printf("x[%2d]: %8.2f  y[%2d]: %8.2f\n", i, x(i), i, y(i));
  }

  // Check maximum error
  T max_err = T(0);
  T val = T(0);
  for (I i = 0; i < N; i++) {
    val = fabs(repeat * alpha * 2.4 + 0.1 - y(i));
    if (val > max_err) {
      max_err = val;
    }
  }

  printf("Maximum error: %20.10e\n", max_err);
#endif
}

void subview() {
  // constexpr A2D::index_t ndof_per_element = 0;  // won't work
  constexpr A2D::index_t ndof_per_element = 1;  // ok
  using ElementDofArray = A2D::MultiArrayNew<A2D::index_t* [ndof_per_element]>;
  A2D::index_t nelems = 10;
  ElementDofArray element_dof("element_dof", nelems);
  A2D::index_t elem = 1;
  auto elem_dof = Kokkos::subview(element_dof, elem, Kokkos::ALL);
}

void test_cuda_axpy_with_UVM(int argc, char* argv[]) {
#ifdef KOKKOS_ENABLE_CUDA
  using ViewDevice_t =
      Kokkos::View<T*, Kokkos::LayoutRight, Kokkos::CudaUVMSpace>;
  if (argc == 1) {
    std::printf("axpy\nusage: ./scratchpad N, where mat size = 2^N\n");
    return;
  }

  // Allocate x and y on host
  I N = pow(2, atoi(argv[1]));
  I bytes = N * sizeof(T);
  double Mbytes = (double)bytes / 1024 / 1024;

  std::printf("Allocating x[%d] on host, size: %.2f MB\n", N, Mbytes);
  std::printf("Allocating y[%d] on host, size: %.2f MB\n", N, Mbytes);

  // Allocate x and y on device
  ViewDevice_t x_device("x_device", N);
  ViewDevice_t y_device("y_device", N);

  for (I i = 0; i < N; i++) {
    x_device(i) = 2.4;
    y_device(i) = 0.1;
  }

  // Perform 100 axpy operations on device and time
  Kokkos::Timer timer;
  int repeat = 100;
  T alpha = 0.01;
  for (int i = 0; i < repeat; i++) {
    auto loop_body = KOKKOS_LAMBDA(const I index) {
      // T alpha = Alpha::get_alpha_array()[0];  // This doesn't work
      T alpha = Alpha::value_array[0];  // This works
      auto x_slice = Kokkos::subview(x_device, Kokkos::ALL);
      y_device(index) += alpha * x_slice(index);
    };
    Kokkos::parallel_for("axpy", N, loop_body);
    Kokkos::fence();
  }

  double elapse = timer.seconds();
  printf("averaged time: %.8f ms\n", elapse * 1e3);

  // Compute bandwidth:
  // x is read once, y is read once and written once
  double bandwidth = 3 * Mbytes / 1024 * repeat / elapse;
  printf("averaged bandwidth: %.8f GB/s\n", bandwidth);

  // Print values
  int len = 10;
  if (N < len) {
    len = N;
  }

  for (I i = 0; i < len; i++) {
    printf("x[%2d]: %8.2f  y[%2d]: %8.2f\n", i, x_device(i), i, y_device(i));
  }

  // Check maximum error
  T max_err = T(0);
  T val = T(0);
  for (I i = 0; i < N; i++) {
    val = fabs(repeat * alpha * 2.4 + 0.1 - y_device(i));
    if (val > max_err) {
      max_err = val;
    }
  }

  printf("Maximum error: %20.10e\n", max_err);
#endif
}

class ParallelVector {
#ifdef KOKKOS_ENABLE_CUDA
  using MemSpace = Kokkos::CudaUVMSpace;
#else
  using MemSpace = Kokkos::HostSpace;
#endif
  using data_t = Kokkos::View<double*, Kokkos::LayoutRight, MemSpace>;

  KOKKOS_FUNCTION void set_one_value(int i, double val) const { data(i) = val; }

 public:
  ParallelVector(int len) : len(len), data("data", len) {}

  void set_values(double val) {
    auto loop_body = KOKKOS_CLASS_LAMBDA(int i) { set_one_value(i, val); };

    Kokkos::parallel_for("loop", len, loop_body);
    Kokkos::fence();
  }

 private:
  int len;
  data_t data;
};

struct Head {
  Head(double val) : val(val) {}
  KOKKOS_FUNCTION double get_val() { return val; };

 private:
  double val;
};
struct Data {
  Data(Head& head, double val, int id = 0) : head(head), val(val), id(id) {}
  double val;
  int id;
  Head& head;
};

// This can build, won't work -  ``it is generally not valid to have any pointer
// or reference members in the functor''
void test_cuda_functor_pass_by_ref() {
#ifdef KOKKOS_ENABLE_CUDA
  Head head(5.6);
  Data data(head, 4.2);
  Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::CudaUVMSpace> view("view",
                                                                        10);
  auto loop_body = [=] __device__(int i) { view(i) = data.head.get_val(); };
  Kokkos::parallel_for("loop", 10, loop_body);
  Kokkos::fence();

  for (int i = 0; i < 10; i++) {
    std::printf("view[%2d]: %.10f\n", i, view[i]);
  }
#endif
}

int main(int argc, char* argv[]) {
  Kokkos::initialize();
  {  // test_axpy(argc, argv);
     // test_matvec(argc, argv);
     // test_unordered_set();
     // test_subview();
     // test_sort();
     // test_is_same_layout();
     // test_complex();
     // test_is_complex();
     // test_view_is_allocated();
     // test_smart_pointer_behavior();
     // test_copy();
     // test_parallel_for();
     // test_modify_view_from_const_lambda();
     // test_cuda_axpy(argc, argv);
     // subview();
     // test_cuda_axpy_with_UVM(argc, argv);
     // ParallelVector pv(10);
     // pv.set_values(4.2);
     // test_cuda_functor_pass_by_ref();
  }
  Kokkos::finalize();
}