#include <iostream>
#include <string>

#include "utils/a2dprofiler.h"
#include "utils/a2dvtk.h"

using namespace A2D;
using namespace std;

int main() {
  Timer t("main()");
  { ReadVTK<4, 3, double, int> vtk_reader("tetra_3d_refine.vtk"); }
  return 0;
}