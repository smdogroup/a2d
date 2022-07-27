#include "Kokkos_Core.hpp"
#include "multiarray.h"

using namespace A2D;
using namespace std;
using ViewType = Kokkos::View<double* [2][3]>;

struct ViewWrapper {
  using ViewType = Kokkos::View<double* [2][3]>;
  ViewWrapper(int leading_dim) {
    // Populate view
    // ...
  }
  ViewType view;
};

void test_view() {
  ViewType view1("label", 4);
  ViewType view2("label", 4);
  cout << "view1.data(): " << view1.data() << endl;
  cout << "view2.data(): " << view2.data() << endl;
  cout << "view1 = view2" << endl;
  view1 = view2;
  cout << "view1.data(): " << view1.data() << endl;
  cout << "view2.data(): " << view2.data() << endl;

  ViewType* view3 = new ViewType("label", 4);
  ViewType* view4 = new ViewType("label", 4);

  cout << "view3->data(): " << view3->data() << endl;
  cout << "view4->data(): " << view4->data() << endl;
  cout << "(*view3) = (*view4)" << endl;
  (*view3) = (*view4);
  cout << "view3->data(): " << view3->data() << endl;
  cout << "view4->data(): " << view4->data() << endl;
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  ViewType a, b;

  Kokkos::finalize();
  return 0;
}
