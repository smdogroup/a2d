#include <iostream>
#include <memory>

#ifdef A2D_USE_KOKKOS
#include "Kokkos_Core.hpp"
#endif
#include "multiarray.h"
#include "parallel.h"

using namespace A2D;
using namespace std;

template <class Array>
void print_array(Array array, int dim0, int dim1) {
  for (int i = 0; i < dim0; i++) {
    for (int j = 0; j < dim1; j++) {
      std::printf("array(%d, %d) = %.2f\n", i, j, array(i, j));
    }
  }
}

template <class Array>
void print_array(Array array, int dim0, int dim1, int dim2) {
  for (int i = 0; i < dim0; i++) {
    for (int j = 0; j < dim1; j++) {
      for (int k = 0; k < dim2; k++) {
        std::printf("array(%d, %d, %d) = %.2f\n", i, j, k, array(i, j, k));
      }
    }
  }
}

void test_lambda() {
  using Layout = CLayout<3>;
  Layout layout(2);
  MultiArray<double, Layout> array1(layout);
  MultiArray<double, Layout> array2(layout);

  array1.random();
  array2.zero();

  printf("array1:\n");
  print_array(array1, 2, 3);
  printf("array2(before):\n");
  print_array(array2, 2, 3);

  auto lam = A2D_LAMBDA(int i) {
    array2(0, i) = array1(0, i);
    array2(1, i) = array1(1, i);
  };

  parallel_for(3, lam);

  printf("array2(after):\n");
  print_array(array2, 2, 3);
}

// void test_multiarray_kokkos() {
//   Kokkos::initialize();
//   {
//     // Create layout
//     using Layout = CLayout<3, 4>;
//     using T = double;
//     Layout layout(5);

//     // Create a2d and kokkos multiarray with the same layout
//     MultiArray<T, Layout> karray(layout);

//     // Proves that our dimension parsing is correct
//     using Ktype = MultiArray<T, Layout>::ViewType;
// #ifdef KOKKOS_ENABLE_CUDA
//     using Ktype_ref =
//         Kokkos::View<T* [3][4], Kokkos::LayoutRight, Kokkos::CudaUVMSpace>;
// #else
//     using Ktype_ref =
//         Kokkos::View<T* [3][4], Kokkos::LayoutRight, Kokkos::HostSpace>;
// #endif
//     static_assert(std::is_same_v<Ktype, Ktype_ref>);
//     assert(karray.view.extent(0) == 5);
//     assert(karray.view.extent(1) == 3);
//     assert(karray.view.extent(2) == 4);

//     Ktype* karray_ptr = new Ktype("label", 5);
//   }
//   Kokkos::finalize();
// }

void test_clayout() {
  using C_1dim_t = CLayout<>;
  using C_2dim_t = CLayout<3>;
  using C_3dim_t = CLayout<3, 4>;

  int N = 10;

  C_1dim_t layout1;  // shape: (N)
  C_2dim_t layout2;  // shape: (N, 3)
  C_3dim_t layout3;  // shape: (N, 3, 4)

  layout1 = C_1dim_t(N);
  layout2 = C_2dim_t(N);
  layout3 = C_3dim_t(N);

  printf("=========================================\n");
  printf("test CLayout\n");
  printf("=========================================\n");

  // Test: get_rank()
  printf("layout1.get_rank() = %d (expect 1)\n", layout1.get_rank());
  printf("layout2.get_rank() = %d (expect 2)\n", layout2.get_rank());
  printf("layout3.get_rank() = %d (expect 3)\n", layout3.get_rank());
  printf("\n");

  // Test: get_size(N)
  printf("layout1.get_size(%d) = %d (expect 10)\n", N, layout1.get_size(N));
  printf("layout2.get_size(%d) = %d (expect 30)\n", N, layout2.get_size(N));
  printf("layout3.get_size(%d) = %d (expect 120)\n", N, layout3.get_size(N));
  printf("\n");

  // Test: get_extent(index)
  printf("layout1.get_extent(0) = %d (expect 10)\n", layout1.get_extent(0));
  printf("\n");

  printf("layout2.get_extent(0) = %d (expect 10)\n", layout2.get_extent(0));
  printf("layout2.get_extent(1) = %d (expect 3)\n", layout2.get_extent(1));
  printf("\n");

  printf("layout3.get_extent(0) = %d (expect 10)\n", layout3.get_extent(0));
  printf("layout3.get_extent(1) = %d (expect 3)\n", layout3.get_extent(1));
  printf("layout3.get_extent(2) = %d (expect 4)\n", layout3.get_extent(2));
  printf("\n");

  // Test: compute_index(i1, ...idx)
  printf("layout1.compute_index(0) = %d (expect 0)\n",
         layout1.compute_index(0));
  printf("layout1.compute_index(1) = %d (expect 1)\n",
         layout1.compute_index(1));
  printf("layout1.compute_index(2) = %d (expect 2)\n",
         layout1.compute_index(2));
  printf("\n");

  printf("layout2.compute_index(0, 0) = %d (expect 0)\n",
         layout2.compute_index(0, 0));
  printf("layout2.compute_index(3, 2) = %d (expect 11)\n",
         layout2.compute_index(3, 2));
  printf("\n");

  printf("layout3.compute_index(0, 0, 0) = %d (expect 0)\n",
         layout3.compute_index(0, 0, 0));
  printf("layout3.compute_index(3, 2, 1) = %d (expect 45)\n",
         layout3.compute_index(3, 2, 1));
  printf("\n");

  // get_size()
  printf("layout1.get_size() = %d (expect 10)\n", layout1.get_size());
  printf("layout2.get_size() = %d (expect 30)\n", layout2.get_size());
  printf("layout3.get_size() = %d (expect 120)\n", layout3.get_size());
}

void test_flayout() {
  using F_1dim_t = FLayout<>;
  using F_2dim_t = FLayout<3>;
  using F_3dim_t = FLayout<3, 4>;

  int N = 10;

  F_1dim_t layout1(N);  // shape: (N)
  F_2dim_t layout2(N);  // shape: (N, 3)
  F_3dim_t layout3(N);  // shape: (N, 3, 4)

  printf("=========================================\n");
  printf("test FLayout\n");
  printf("=========================================\n");

  // Test: get_rank()
  printf("layout1.get_rank() = %d (expect 1)\n", layout1.get_rank());
  printf("layout2.get_rank() = %d (expect 2)\n", layout2.get_rank());
  printf("layout3.get_rank() = %d (expect 3)\n", layout3.get_rank());
  printf("\n");

  // Test: get_size(N)
  printf("layout1.get_size(%d) = %d (expect 10)\n", N, layout1.get_size(N));
  printf("layout2.get_size(%d) = %d (expect 30)\n", N, layout2.get_size(N));
  printf("layout3.get_size(%d) = %d (expect 120)\n", N, layout3.get_size(N));
  printf("\n");

  // Test: get_extent(index)
  printf("layout1.get_extent(0) = %d (expect 10)\n", layout1.get_extent(0));
  printf("\n");

  printf("layout2.get_extent(0) = %d (expect 10)\n", layout2.get_extent(0));
  printf("layout2.get_extent(1) = %d (expect 3)\n", layout2.get_extent(1));
  printf("\n");

  printf("layout3.get_extent(0) = %d (expect 10)\n", layout3.get_extent(0));
  printf("layout3.get_extent(1) = %d (expect 3)\n", layout3.get_extent(1));
  printf("layout3.get_extent(2) = %d (expect 4)\n", layout3.get_extent(2));
  printf("\n");

  // Test: compute_index(i1, ...idx)
  printf("layout1.compute_index(0) = %d (expect 0)\n",
         layout1.compute_index(0));
  printf("layout1.compute_index(1) = %d (expect 1)\n",
         layout1.compute_index(1));
  printf("layout1.compute_index(2) = %d (expect 2)\n",
         layout1.compute_index(2));
  printf("\n");

  printf("layout2.compute_index(0, 0) = %d (expect 0)\n",
         layout2.compute_index(0, 0));
  printf("layout2.compute_index(3, 2) = %d (expect 23)\n",
         layout2.compute_index(3, 2));
  printf("\n");

  printf("layout3.compute_index(0, 0, 0) = %d (expect 0)\n",
         layout3.compute_index(0, 0, 0));
  printf("layout3.compute_index(3, 2, 1) = %d (expect 53)\n",
         layout3.compute_index(3, 2, 1));
  printf("\n");

  // get_size()
  printf("layout1.get_size() = %d (expect 10)\n", layout1.get_size());
  printf("layout2.get_size() = %d (expect 30)\n", layout2.get_size());
  printf("layout3.get_size() = %d (expect 120)\n", layout3.get_size());
}

void test_uninit_multiarray() {
  using C_1dim_t = CLayout<>;
  using C_2dim_t = CLayout<3>;
  using C_3dim_t = CLayout<3, 4>;

  using T = double;

  using Array_1dim_t = MultiArray<T, C_1dim_t>;
  using Array_2dim_t = MultiArray<T, C_2dim_t>;
  using Array_3dim_t = MultiArray<T, C_3dim_t>;

  int N = 10;

  C_1dim_t layout1;  // shape: (N)
  C_2dim_t layout2;  // shape: (N, 3)
  C_3dim_t layout3;  // shape: (N, 3, 4)

  layout1 = C_1dim_t(N);
  layout2 = C_2dim_t(N);
  layout3 = C_3dim_t(N);

  Array_1dim_t array1;
  Array_2dim_t array2;
  Array_3dim_t array3;

  array1 = Array_1dim_t(layout1);
  array2 = Array_2dim_t(layout2);
  array3 = Array_3dim_t(layout3);
}

int main(int argc, char* argv[]) {
  //   test_clayout();
  //   test_flayout();
  test_uninit_multiarray();
  return 0;
}