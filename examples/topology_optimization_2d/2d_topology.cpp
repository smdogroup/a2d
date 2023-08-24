#include "2d_topology.h"

#include <complex>
#include <cstdlib>
#include <string>

#include "../topology_optimization/topology_paropt_prob.h"

template <typename T>
class MeshCreator {
 public:
  MeshCreator(index_t nx, index_t ny)
      : nx(nx),
        ny(ny),
        nz(0),
        nelems(nx * ny),
        nverts((nx + 1) * (ny + 1)),
        elem(4 * nelems),
        Xloc(2 * nverts),
        tx_traction(2) {}
  MeshCreator(index_t nx, index_t ny, index_t nz)
      : nx(nx),
        ny(ny),
        nz(nz),
        nelems(nx * ny * nz),
        nverts((nx + 1) * (ny + 1) * (nz + 1)),
        elem(8 * nelems),
        Xloc(3 * nverts),
        tx_traction(3) {}

  index_t create_2d_mesh(double lx = 1.0, double ly = 1.0) {
    // Set up mesh
    MesherRect2D mesher(nx, ny, lx, ly);
    mesher.set_X_conn<index_t, double>(Xloc.data(), elem.data(), true, 0, 0.2);
    conn = new MeshConnectivity2D(nverts, index_t(0), (index_t *)nullptr,
                                  nelems, elem.data());

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
    tx_traction[0] = T(0.0);
    tx_traction[1] = T(-1.0);

    index_t bc_label =
        conn->add_boundary_label_from_verts(num_boundary_verts, boundary_verts);
    index_t traction_label =
        conn->add_boundary_label_from_verts(num_traction_verts, traction_verts);
    bcinfo.add_boundary_condition(bc_label);

    return traction_label;
  }

  index_t create_3d_mesh(double lx = 1.0, double ly = 1.0, double lz = 1.0) {
    // Set up mesh
    MesherBrick3D mesher(nx, ny, nz, lx, ly, lz);
    mesher.set_X_conn<index_t, double>(Xloc.data(), elem.data());
    conn = new MeshConnectivity3D(
        nverts, (index_t)0, (index_t *)nullptr, nelems, elem.data(), (index_t)0,
        (index_t *)nullptr, index_t(0), (index_t *)nullptr);

    const index_t num_boundary_verts = 2 * (ny + 1) * (nz + 1);
    const index_t num_traction_verts = 2 * (ny + 1) * (nz + 1);

    index_t boundary_verts[num_boundary_verts];
    index_t traction_verts[num_traction_verts];

    for (index_t k = 0, index = 0; k < nz + 1; k++) {
      for (index_t j = 0; j < ny + 1; j++, index++) {
        boundary_verts[index] = node_num(0, j, k);
      }
    }

    for (index_t k = 0, index = 0; k < nz + 1; k++) {
      for (index_t j = 0; j < ny + 1; j++, index++) {
        traction_verts[index] = node_num(nx, j, k);
      }
    }

    tx_traction[0] = T(0.0);
    tx_traction[1] = T(0.0);
    tx_traction[2] = T(-1.0);

    A2D::index_t bc_label =
        conn->add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

    A2D::index_t traction_label =
        conn->add_boundary_label_from_verts(num_traction_verts, traction_verts);
    bcinfo.add_boundary_condition(bc_label);

    return traction_label;
  }

  index_t nx, ny, nz;
  index_t nelems, nverts;
  std::vector<index_t> elem;
  std::vector<double> Xloc;
  std::vector<T> tx_traction;
  MeshConnectivityBase *conn;
  DirichletBCInfo bcinfo;

 private:
  index_t node_num(index_t i, index_t j) { return i + j * (nx + 1); }
  index_t node_num(index_t i, index_t j, index_t k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  }
};

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    using T = double;

    constexpr int degree = 2;

    constexpr int spatial_dim = 2;
    MeshCreator<T> creator(128, 32);
    index_t traction_label = creator.create_2d_mesh(4.0, 1.0);

    // constexpr int spatial_dim = 3;
    // MeshCreator<T> creator(16, 4, 4);
    // index_t traction_label = creator.create_3d_mesh(4.0, 1.0, 1.0);

    StopWatch watch;
    double t1 = watch.lap();
    TopoElasticityAnalysis<T, degree, spatial_dim> prob(
        *creator.conn, creator.bcinfo, creator.nelems, creator.elem.data(),
        creator.Xloc.data(), traction_label, creator.tx_traction.data());
    double t2 = watch.lap();
    printf("Setup time: %12.5e s\n", t2 - t1);

    prob.set_design_var(1.0);
    prob.solve();
    prob.tovtk("solution.vtk");

    CSCMat<T> &csc_mat = prob.get_csc_matrix();
    printf("number of vertices: %d\n", creator.conn->get_num_verts());
    printf("number of edges:    %d\n", creator.conn->get_num_edges());
    printf("number of faces:    %d\n", creator.conn->get_num_bounds());
    printf("number of elements: %d\n", creator.conn->get_num_elements());
    printf("number of dofs:     %d\n", prob.get_mesh().get_num_dof());
    printf("matrix dimension:  (%d, %d)\n", csc_mat.nrows, csc_mat.ncols);

    const std::string prefix = "results";
    int nvars = prob.get_num_design_vars();
    int ncon = 1;
    int nineq = 1;
    TopOptProb<decltype(prob)> topo_prob(prefix, MPI_COMM_WORLD, nvars, ncon,
                                         nineq, 1.0, 1.0, 0.4, prob, true, 10);
    topo_prob.checkGradients(1e-6);
  }

  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}
