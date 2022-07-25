#include <Kokkos_Core.hpp>

void body(void) {
  int dim0 = 2;
  constexpr int dim1 = 3;

  Kokkos::View<double* [dim1]> data("label", dim0);
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  body();

  Kokkos::finalize();
  return 0;
}
