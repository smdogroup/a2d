#include "topology_elasticity.h"
#include "utils/a2dmesh.h"
#include "utils/a2dprofiler.h"

int main(int argc, char *argv[]) {
  A2D::Timer timer("main()");
  A2D::Timer::set_threshold_ms(5.0);
  Kokkos::initialize();

  const A2D::index_t dim = 3;
  const A2D::index_t degree = 2;
  const A2D::index_t filter_degree = 1;

  using T = double;
  // using T = A2D_complex_t<double>;

  // Number of elements in each dimension
  const int nx = 20, ny = 10, nz = 10;
  const double lx = 2.0, ly = 1.0, lz = 1.0;

  // Set up mesh
  const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = nx * ny * nz;
  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8 * nhex];
  double Xloc[3 * nverts];
  A2D::MesherBrick3D mesher(nx, ny, nz, lx, ly, lz);
  mesher.set_X_conn<int, double>(Xloc, hex);
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Set up bcs
  auto node_num = [](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  const int num_boundary_verts = 2 * (ny + 1) * (nz + 1);
  const int num_traction_verts = 2 * (ny + 1) * (nz + 1);
  int boundary_verts[num_boundary_verts];
  int traction_verts[num_traction_verts];

  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(0, j, k);
    }
  }

  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      traction_verts[index] = node_num(nx, j, k);
    }
  }

  A2D::index_t bc_label =
      conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

  A2D::index_t traction_label =
      conn.add_boundary_label_from_verts(num_traction_verts, traction_verts);

  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label);

  // Set the traction components
  T t[3] = {1.0, 0.0, 0.0};

  // Set the body force components
  T tb[3] = {0.0, 0.0, 0.0};

  // Set the torque components
  T tt[3] = {0.0, 0.0, 0.0};
  T x0[3] = {0.0, 0.0, 0.0};

  // Create the finite-element model
  T E = 70.0e3, nu = 0.3, q = 5.0;
  T design_stress = 200.0, ks_penalty = 50.0;
  bool verbose = true;
  int amg_nlevels = 3, cg_it = 100;
  double cg_rtol = 1e-8, cg_atol = 1e-30;
  TopoElasticityAnalysis<T, degree, filter_degree> topo(
      conn, bcinfo, E, nu, q, tb, traction_label, t, x0, tt, verbose,
      amg_nlevels, cg_it, cg_rtol, cg_atol);

  // Print info
  int nvars = topo.get_num_design_vars();
  int ndof = topo.get_num_dofs();
  int nelems = topo.get_num_elements();
  std::printf("basis degree:                 %d\n", degree);
  std::printf("number of mesh elements:      %d\n", nelems);
  std::printf("number of design variables:   %d\n", nvars);
  std::printf("number of degrees of freedom: %d\n", ndof);

  // Set the geometry from the node locations
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<
      TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  topo.reset_geometry();
  std::vector<T> x(topo.get_num_design_vars(), T(1.0));

  // Set the design variables
  topo.set_design_vars(x);

  // Solve the problem
  topo.solve();

  // Sensitivity analysis
  auto dvdx = std::make_shared<A2D::MultiArrayNew<T *>>(
      "dfdx", topo.get_num_design_vars());
  T vol = topo.eval_volume();
  topo.add_volume_gradient(*dvdx);

  printf("v:      %20.10f\n", vol);
  printf("|dvdx|: %20.10f\n", A2D::BLAS::norm(*dvdx));

  // Write the solution to a vtk file
  topo.tovtk("toy_elasticity.vtk");

  return 0;
}
