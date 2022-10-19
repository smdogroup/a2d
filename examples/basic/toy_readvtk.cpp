#include <iostream>
#include <string>

#include "array.h"
#include "utils/a2dprofiler.h"
#include "utils/a2dvtk.h"

using namespace A2D;
using namespace std;

int main() {
#ifdef A2D_USE_KOKKOS
  Kokkos::initialize();
#endif
  Timer t("main()");
  {
    // Read vtk
    const static int node_per_elem = 4;
    const static int spatial_dim = 3;
    using T = double;
    using I = int;

    ReadVTK<node_per_elem, spatial_dim, T, I> vtk_reader("tetra_3d.vtk");

    I nnodes = vtk_reader.nnodes;
    I nelems = vtk_reader.nelems;

    // Write a new vtk
    using ConnArray = A2D::MultiArrayNew<I* [node_per_elem]>;
    using NodeArray = A2D::MultiArrayNew<T* [spatial_dim]>;

    ConnArray conn = ConnArray("conn", vtk_reader.nelems);
    NodeArray X = NodeArray("X", vtk_reader.nnodes);

    vtk_reader.set_conn(conn);
    vtk_reader.set_X(X);

    ToVTK<ConnArray, NodeArray> vtk_writer(conn, X, 10, "tetra_3d_out.vtk");
    vtk_writer.write_mesh();
  }

#ifdef A2D_USE_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}