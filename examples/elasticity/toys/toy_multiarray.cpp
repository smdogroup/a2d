#include <iostream>
#include <memory>

#include "multiarray.h"

using namespace A2D;

int main(int argc, char* argv[]) {
  int nnodes = 20;

  CLayout<2, 3, 4> dv_layout(nnodes);
  auto dv = std::make_shared<MultiArray<double, CLayout<2, 3, 4>>>(dv_layout);

  for (int i = 0; i < 10; i++) {
    std::printf("index(%d,1,1,2) = %d\n", i,
                dv_layout.compute_index(i, 1, 1, 2));
  }
  std::cout << "done\n";
}