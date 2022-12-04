#include "topology_heat.h"

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;

void main_body(int argc, char *argv[]) {
  constexpr int degree = 2;
  constexpr int filter_degree = 1;

  // Set up cmd arguments and defaults
  ArgumentParser parser(argc, argv);
  std::string vtk_name =
      parser.parse_option("--vtk", std::string("3d_hex_1000.vtk"));
  std::string prefix = parser.parse_option("--prefix", std::string("results"));
  int maxit = parser.parse_option("--maxit", 100);
  double bc_temp = parser.parse_option("--bc_temp", 0.0);
  double heat_source = parser.parse_option("--heat_source", 1.0);
  double ramp_q = parser.parse_option("--ramp_q", 5.0);
  bool check_grad_and_exit = parser.parse_option("--check_grad_and_exit");
  bool verbose = parser.parse_option("--verbose");
  parser.help_info();

  // Set up result directory
  if (!std::filesystem::is_directory(prefix)) {
    std::filesystem::create_directory(prefix);
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
  std::vector<int> ids{100};
  std::vector<I> verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Set up boundary condition information
  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(verts.size(), verts.data()), basis);

  // Initialize the analysis instance
  T kappa = 1.0;
  TopoHeatAnalysis<T, degree, filter_degree> topo(
      conn, bcinfo, kappa, heat_source, bc_temp, ramp_q, verbose);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<
      TopoHeatAnalysis<T, degree, filter_degree>::GeoBasis>(nhex, hex, Xloc,
                                                            elem_geo);
  topo.reset_geometry();

  // Get problem size
  int nvars = topo.get_num_design_vars();
  int ndof = topo.get_num_dofs();

  // Print info
  std::printf("number of design variables:   %d\n", nvars);
  std::printf("number of degrees of freedom: %d\n", ndof);

  // Initialize the optimization object
  int ncon = 1;
  int nineq = 1;
  topo.solve();
  T comp = topo.eval_compliance();
  T ref_hcomp = topo.eval_compliance();
  T domain_vol = topo.eval_volume();
  T volume_frac = 0.4;
  TopOptProb<TopoHeatAnalysis<T, degree, filter_degree>> *prob =
      new TopOptProb(prefix, MPI_COMM_WORLD, nvars, ncon, nineq, ref_hcomp,
                     domain_vol, volume_frac, topo, verbose);
  prob->incref();

  // Sanity check
  if (check_grad_and_exit) {
    prob->checkGradients(1e-6);
    return;
  }

  // Define optimizer and options
  ParOptOptions *options = new ParOptOptions;
  options->incref();
  ParOptOptimizer::addDefaultOptions(options);

  options->setOption("algorithm", "mma");
  options->setOption("mma_asymptote_contract", 0.7);
  options->setOption("mma_asymptote_relax", 1.2);
  options->setOption("mma_bound_relax", 0);
  options->setOption("mma_delta_regularization", 1e-05);
  options->setOption("mma_eps_regularization", 0.001);
  options->setOption("mma_infeas_tol", 1e-05);
  options->setOption("mma_init_asymptote_offset", 0.25);
  options->setOption("mma_l1_tol", 1e-06);
  options->setOption("mma_linfty_tol", 1e-06);
  options->setOption("mma_max_asymptote_offset", 10);
  options->setOption("mma_max_iterations", maxit);
  options->setOption("mma_min_asymptote_offset", 0.01);
  options->setOption("mma_use_constraint_linearization", true);

  ParOptOptimizer *opt = new ParOptOptimizer(prob, options);
  opt->incref();
  opt->optimize();

  // Delete objects
  prob->decref();
  options->incref();
  opt->decref();
  return;
}

int main(int argc, char *argv[]) {
  A2D::Timer timer("main()");
  MPI_Init(&argc, &argv);
  Kokkos::initialize();
  {
    // double kappa = 1.0, q = 1.2, heat_source = 3.4;
    // A2D::HeatConduction<double, 3> heat(kappa, q, heat_source);
    // A2D::TestPDEImplementation<double>(heat);
    // test_heat_analysis(argc, argv);
    main_body(argc, argv);
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}