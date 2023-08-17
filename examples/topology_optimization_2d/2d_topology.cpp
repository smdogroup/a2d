#include "2d_topology.h"

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    using T = double;

    constexpr int degree = 3;

    // Number of elements in each dimension
    const index_t nx = 16, ny = 16;
    const double lx = 1.5, ly = 2.0;

    // Set up mesh
    const index_t nverts = (nx + 1) * (ny + 1);
    const index_t ntri = 0, nquad = nx * ny;
    index_t *tri = nullptr;
    std::vector<index_t> quad(4 * nquad);
    std::vector<double> Xloc(2 * nverts);
    MesherRect2D mesher(nx, ny, lx, ly);
    mesher.set_X_conn<index_t, double>(Xloc.data(), quad.data());
    MeshConnectivity2D conn(nverts, ntri, tri, nquad, quad.data());

    // Set up bcs and traction
    auto node_num = [](index_t i, index_t j) { return i + j * (nx + 1); };

    const index_t num_boundary_verts = (ny + 1);
    const index_t num_traction_verts = (ny + 1);

    index_t boundary_verts[num_boundary_verts];
    index_t traction_verts[num_traction_verts];

    for (index_t j = 0; j < ny + 1; j++) {
      boundary_verts[j] = node_num(0, j);
    }
    for (index_t j = 0; j < (ny + 1); j++) {
      traction_verts[j] = node_num(nx, j);
    }
    T tx_traction[2] = {0.0, -1.0};

    A2D::index_t bc_label =
        conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);
    A2D::index_t traction_label =
        conn.add_boundary_label_from_verts(num_traction_verts, traction_verts);

    A2D::DirichletBCInfo bcinfo;
    bcinfo.add_boundary_condition(bc_label);

    StopWatch watch;
    double t1 = watch.lap();
    TopoElasticityAnalysis2D<T, degree> prob(conn, bcinfo, nquad, quad.data(),
                                             Xloc.data(), traction_label,
                                             tx_traction);
    double t2 = watch.lap();

    CSCMat<T> &csc_mat = prob.get_csc_matrix();

    double t3 = watch.lap();
    prob.solve();
    double t4 = watch.lap();

    prob.tovtk("solution.vtk");

    printf("number of vertices: %d\n", conn.get_num_verts());
    printf("number of edges:    %d\n", conn.get_num_edges());
    printf("number of faces:    %d\n", conn.get_num_bounds());
    printf("number of elements: %d\n", conn.get_num_elements());
    printf("number of dofs:     %d\n", prob.get_mesh().get_num_dof());
    printf("matrix dimension:  (%d, %d)\n", csc_mat.nrows, csc_mat.ncols);

    printf("Setup time:                %12.5e s\n", t2 - t1);
    printf("Solve time:                %12.5e s\n", t4 - t3);
  }

  Kokkos::finalize();

  return 0;
}
