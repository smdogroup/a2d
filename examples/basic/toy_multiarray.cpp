#include <iostream>
#include <memory>

#include "Kokkos_Core.hpp"
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

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
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

    auto lam = KOKKOS_LAMBDA(int i) {
      array2(0, i) = array1(0, i);
      array2(1, i) = array1(1, i);
    };

    parallel_for(3, lam);

    printf("array2(after):\n");
    print_array(array2, 2, 3);
  }
  Kokkos::finalize();
}
