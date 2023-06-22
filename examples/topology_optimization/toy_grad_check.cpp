#include "topology_elasticity.h"
#include "utils/a2dprofiler.h"

int main(int argc, char *argv[]) {
  A2D::Timer timer("main()");
  A2D::Timer::set_threshold_ms(5.0);
  Kokkos::initialize();

  const A2D::index_t dim = 3;
  const A2D::index_t degree = 4;
  const A2D::index_t filter_degree = 2;

  std::cout << "Topology linear elasticity\n";
  A2D::TopoLinearElasticity<std::complex<double>, dim> elasticity(70e3, 0.3,
                                                                  5.0);
  A2D::TestPDEImplementation<std::complex<double>>(elasticity);

  using T = std::complex<double>;

  // Number of elements in each dimension
  const int nx = 8, ny = 4, nz = 4;
  auto node_num = [](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  // Number of edges
  const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = nx * ny * nz;

  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8 * nhex];

  using ET = A2D::ElementTypes;

  for (int k = 0, e = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, e++) {
        for (int ii = 0; ii < ET::HEX_VERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  double Xloc[3 * nverts];
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (2.0 * i) / nx;
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j) / ny;
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k) / nz;
      }
    }
  }

  // Constrain the nodes at either end of the block
  const int num_boundary_verts = 2 * (ny + 1) * (nz + 1);
  int boundary_verts[num_boundary_verts];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(0, j, k);
    }
  }

  for (int k = 0, index = (ny + 1) * (nz + 1); k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(nx, j, k);
    }
  }

  const int num_surface_verts = (nx + 1) * (ny + 1);
  int surface_verts[num_surface_verts];
  for (int j = 0, index = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++, index++) {
      surface_verts[index] = node_num(i, j, nz);
    }
  }

  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  A2D::index_t bc_label =
      conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

  A2D::index_t traction_label =
      conn.add_boundary_label_from_verts(num_surface_verts, surface_verts);

  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label, basis);

  // Set the traction components
  T t[3] = {0.2, -0.3, 0.5};

  // Set the body force components
  T tb[3] = {0.1, 0.1, -0.2};

  // Set the torque components
  T tt[3] = {0.2, -0.1, 0.5};
  T x0[3] = {0.0, 0.1, 0.2};

  // Create the finite-element model
  T E = 70.0e3, nu = 0.3, q = 5.0;
  T design_stress = 200.0, ks_penalty = 50.0;
  bool verbose = true;
  int amg_nlevels = 3, cg_it = 100;
  double cg_rtol = 1e-8, cg_atol = 1e-30;
  TopoElasticityAnalysis<T, degree, filter_degree> topo(
      conn, bcinfo, E, nu, q, tb, traction_label, t, tt, x0, verbose,
      amg_nlevels, cg_it, cg_rtol, cg_atol);

  // Set the geometry from the node locations
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<
      TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  topo.reset_geometry();

  int ndvs = topo.get_num_design_vars();
  std::vector<T> x(ndvs, T(1.0));
  std::vector<T> px(ndvs);

  double dh = 1e-30;
  if constexpr (std::is_same<T, std::complex<double>>::value == true) {
    // Create a perturbed set of design variables
    for (A2D::index_t i = 0; i < ndvs; i++) {
      px[i] = -0.25 + 1.0 * (i % 3);
      x[i] = 1.0 + dh * px[i] * T(0.0, 1.0);
    }
  }

  // Set the design variables
  topo.set_design_vars(x);

  // Solve the problem
  topo.solve();

  std::vector<T> dVdx(ndvs, 0.0);
  T volume = topo.eval_volume();
  topo.add_volume_gradient(dVdx);
  std::cout << "Volume value = " << volume << std::endl;

  std::vector<T> dcdx(ndvs, 0.0);
  T compliance = topo.eval_compliance();
  topo.add_compliance_gradient(dcdx);
  std::cout << "Compliance value = " << compliance << std::endl;

  std::vector<T> dfdx(ndvs, 0.0);
  T aggregation = topo.eval_aggregation(design_stress, ks_penalty);
  topo.add_aggregation_gradient(design_stress, ks_penalty, dfdx);
  std::cout << "Aggregation value = " << aggregation << std::endl;

  if constexpr (std::is_same<T, std::complex<double>>::value == true) {
    double fd;
    T ans;

    fd = imag(volume) / dh;
    ans = 0.0;
    for (A2D::index_t i = 0; i < ndvs; i++) {
      ans += dVdx[i] * px[i];
    }

    std::cout << "Volume gradient check" << std::endl;
    std::cout << std::setw(15) << "Complex-step" << std::setw(15) << "Adjoint"
              << std::setw(15) << "Relative error" << std::endl;
    std::cout << std::setw(15) << fd << std::setw(15) << std::real(ans)
              << std::setw(15) << (fd - std::real(ans)) / fd << std::endl;

    fd = imag(compliance) / dh;
    ans = 0.0;
    for (A2D::index_t i = 0; i < ndvs; i++) {
      ans += dcdx[i] * px[i];
    }

    std::cout << "Compliance gradient check" << std::endl;
    std::cout << std::setw(15) << "Complex-step" << std::setw(15) << "Adjoint"
              << std::setw(15) << "Relative error" << std::endl;
    std::cout << std::setw(15) << fd << std::setw(15) << std::real(ans)
              << std::setw(15) << (fd - std::real(ans)) / fd << std::endl;

    fd = imag(aggregation) / dh;
    ans = 0.0;
    for (A2D::index_t i = 0; i < ndvs; i++) {
      ans += dfdx[i] * px[i];
    }

    std::cout << "Aggregation gradient check" << std::endl;
    std::cout << std::setw(15) << "Complex-step" << std::setw(15) << "Adjoint"
              << std::setw(15) << "Relative error" << std::endl;
    std::cout << std::setw(15) << fd << std::setw(15) << std::real(ans)
              << std::setw(15) << (fd - std::real(ans)) / fd << std::endl;
  }

  // Write the problem to a vtk file
  topo.tovtk("filename.vtk");

  return 0;
}
