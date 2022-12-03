#include "topology_heat.h"

// #include "topology_paropt_prob.h"

#if 0
void main_body(int argc, char *argv[]) {
  constexpr int DEGREE = 6;
  using I = A2D::index_t;
  using T = double;

  // Load connectivity and vertex coordinates from vtk
  std::string vtk_name = "3d_hex.vtk";
  if (argc > 1) {
    vtk_name = argv[1];
  }
  A2D::ReadVTK3D<I, T> readvtk(vtk_name);
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;

  // Construct connectivity
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Extract boundary vertices
  std::vector<int> ids{100, 101};
  std::vector<I> verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Set up boundary condition information
  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(verts.size(), verts.data()), basis);

  // Initialize the analysis instance
  T kappa = 1.0, q = 5.0, heat_source = 1.0;
  HeatTopoOpt<T, DEGREE> topo(conn, bcinfo, kappa, q, heat_source);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<HeatTopoOpt<T, DEGREE>::GeoBasis>(nhex, hex, Xloc,
                                                                elem_geo);
  topo.reset_geometry();

  // Initialize the optimization object
  int nvars = topo.get_num_design_vars();
  int ncon = 1;
  int nineq = 1;

  topo.set_design_vars(std::vector<T>(nvars, T(1.0)));
  topo.solve();
  topo.tovtk("heat_topo.vtk");
  T ref_hcomp = topo.eval_compliance();
  T domain_vol = topo.eval_volume();

  T volume_frac = 0.4;
  HeatTopOptProb *prob =
      new HeatTopOptProb(MPI_COMM_WORLD, nvars, ncon, nineq, ref_hcomp,
                         domain_vol, volume_frac, topo);

  prob->incref();

  // Sanity check
  prob->checkGradients(1e-6);
  exit(0);
}
#endif

int main(int argc, char *argv[]) {
  A2D::Timer timer("main()");
  // MPI_Init(&argc, &argv);
  Kokkos::initialize();
  {
    double kappa = 1.0, q = 1.2, heat_source = 3.4;
    A2D::HeatConduction<double, 3> heat(kappa, q, heat_source);
    A2D::TestPDEImplementation<double>(heat);
    test_heat_analysis(argc, argv);
    // main_body(argc, argv);
  }
  Kokkos::finalize();
  // MPI_Finalize();
  return 0;
}