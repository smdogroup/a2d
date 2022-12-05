#include "topology_elasticity.h"

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;

void main_body(int argc, char *argv[]) {
  constexpr int degree = 4;
  constexpr int filter_degree = 2;

  // Set up cmd arguments and defaults
  ArgumentParser parser(argc, argv);
  std::string vtk_name =
      parser.parse_option("--vtk", std::string("3d_hex.vtk"));
  std::string prefix = parser.parse_option("--prefix", std::string("results"));
  int maxit = parser.parse_option("--maxit", 200);
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
  T E = 70.0e3, nu = 0.3;
  TopoElasticityAnalysis<T, degree, filter_degree> topo(conn, bcinfo, E, nu,
                                                        ramp_q, verbose);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<
      TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  topo.reset_geometry();

  // Initialize the optimization object
  int nvars = topo.get_num_design_vars();
  int ncon = 1;
  int nineq = 1;
  topo.set_design_vars(std::vector<T>(nvars, T(1.0)));
  topo.solve();
  T ref_comp = topo.eval_compliance();
  T domain_vol = topo.eval_volume();
  T volume_frac = 0.4;
  TopOptProb<TopoElasticityAnalysis<T, degree, filter_degree>> *prob =
      new TopOptProb(prefix, MPI_COMM_WORLD, nvars, ncon, nineq, ref_comp,
                     domain_vol, volume_frac, topo, verbose);
  prob->incref();

  // Sanity check
  if (check_grad_and_exit) {
    prob->checkGradients(1e-6);
    return;
  }

  // Define paropt optimizer
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

  // options->setOption("algorithm", "tr");
  // options->setOption("output_level", 0);
  // options->setOption("norm_type", "l1");
  // options->setOption("tr_init_size", 0.05);
  // options->setOption("tr_min_size", 1e-3);
  // options->setOption("tr_max_size", 1.0);
  // options->setOption("tr_eta", 0.25);
  // options->setOption("tr_infeas_tol", 1e-6);
  // options->setOption("tr_l1_tol", 0.0);
  // options->setOption("tr_linfty_tol", 0.0);
  // options->setOption("tr_adaptive_gamma_update", true);
  // options->setOption("tr_accept_step_strategy", "penalty_method");
  // options->setOption("filter_sufficient_reduction", true);
  // options->setOption("filter_has_feas_restore_phase", true);
  // options->setOption("tr_use_soc", false);
  // options->setOption("tr_max_iterations", 100);
  // options->setOption("penalty_gamma", 50.0);
  // options->setOption("qn_subspace_size", 5);
  // options->setOption("qn_type", "bfgs");
  // options->setOption("qn_diag_type", "yty_over_yts");
  // options->setOption("abs_res_tol", 1e-8);
  // options->setOption("starting_point_strategy", "affine_step");
  // options->setOption("barrier_strategy", "mehrotra_predictor_corrector");
  // options->setOption("tr_steering_barrier_strategy",
  //                    "mehrotra_predictor_corrector");
  // options->setOption("tr_steering_starting_point_strategy", "affine_step");
  // options->setOption("use_line_search", false);
  // options->setOption("max_major_iters", 20);

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
  { main_body(argc, argv); }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}