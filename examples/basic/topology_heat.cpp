#include "topology_heat.h"

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;

void main_body(int argc, char *argv[]) {
  constexpr int degree = 5;
  constexpr int filter_degree = 4;

  // Set up cmd arguments and defaults
  ArgumentParser parser(argc, argv);
  std::string vtk_name =
      parser.parse_option("--vtk", std::string("3d_hex_1000.vtk"));
  std::string prefix = parser.parse_option("--prefix", std::string("results"));
  int maxit = parser.parse_option("--maxit", 200);
  int vtk_freq = parser.parse_option("--vtk_freq", 10);
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

  // Set up profiler
  A2D::TIMER_OUTPUT_FILE =
      std::filesystem::path(prefix) / std::filesystem::path("profile.log");
  A2D::Timer timer("main_body()");

  // Read in vtk
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

  // // Extract boundary vertices
  // std::vector<int> ids{100};
  // std::vector<I> verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Set up bc region
  T lower[3], upper[3];
  readvtk.get_bounds(lower, upper);
  T ratio = 0.4;
  T tol = 1e-6;
  T xmin = lower[0];
  T xmax = lower[0];
  T ymin = lower[1] + (1.0 - ratio) / 2 * (upper[1] - lower[1]);
  T ymax = lower[1] + (1.0 + ratio) / 2 * (upper[1] - lower[1]);
  T zmin = lower[2] + (1.0 - ratio) / 2 * (upper[2] - lower[2]);
  T zmax = lower[2] + (1.0 + ratio) / 2 * (upper[2] - lower[2]);
  std::vector<I> verts =
      readvtk.get_verts_within_box(xmin, xmax, ymin, ymax, zmin, zmax);

  // Write bc verts to vtk for debugging purpose
  std::vector<T> labels(3 * nverts, 0.0);
  for (auto it = verts.begin(); it != verts.end(); it++) {
    labels[3 * (*it)] = 1.0;
  }
  {
    A2D::VectorFieldToVTK fieldtovtk(
        nverts, Xloc, labels.data(),
        std::filesystem::path(prefix) / std::filesystem::path("bc_verts.vtk"));
  }

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
                     domain_vol, volume_frac, topo, verbose, vtk_freq);
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

  std::string out_path =
      std::filesystem::path(prefix) / std::filesystem::path("paropt.out");
  std::string tr_path =
      std::filesystem::path(prefix) / std::filesystem::path("paropt.tr");
  std::string mma_path =
      std::filesystem::path(prefix) / std::filesystem::path("paropt.mma");
  options->setOption("algorithm", "mma");
  options->setOption("mma_max_iterations", maxit);
  options->setOption("output_file", out_path.c_str());
  options->setOption("tr_output_file", tr_path.c_str());
  options->setOption("mma_output_file", mma_path.c_str());

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
  MPI_Init(&argc, &argv);
  Kokkos::initialize();
  {
    // double kappa = 1.0, q = 1.2, heat_source = 3.4;
    // A2D::HeatConduction<double, 3> heat(kappa, q, heat_source);
    // A2D::TestPDEImplementation<double>(heat);
    // test_heat_analysis(argc, argv);
    // test_find_bc_verts(argc, argv);
    main_body(argc, argv);
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}