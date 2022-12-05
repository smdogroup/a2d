#include "topology_elasticity.h"

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;
using fspath = fspath;

void main_body(int argc, char *argv[]) {
  constexpr int degree = 2;
  constexpr int filter_degree = 1;

  // Set up cmd arguments and defaults
  A2D::ArgumentParser parser(argc, argv);
  std::string vtk_name =
      parser.parse_option("--vtk", std::string("3d_hex.vtk"));
  std::string prefix = parser.parse_option("--prefix", std::string("results"));
  int maxit = parser.parse_option("--maxit", 200);
  int vtk_freq = parser.parse_option("--vtk_freq", 10);
  double ramp_q = parser.parse_option("--ramp_q", 5.0);
  bool check_grad_and_exit = parser.parse_option("--check_grad_and_exit");
  bool verbose = parser.parse_option("--verbose");
  int amg_nlevels = parser.parse_option("--amg_nlevels", 3);
  int cg_it = parser.parse_option("--cg_it", 100);
  double cg_rtol = parser.parse_option("--cg_rtol", 1e-8);
  double cg_atol = parser.parse_option("--cg_atol", 1e-30);
  parser.help_info();

  // Set up result directory
  if (!std::filesystem::is_directory(prefix)) {
    std::filesystem::create_directory(prefix);
  }

  // Save cmd arguments to txt
  A2D::save_cmd(argc, argv, (fspath(prefix) / fspath("cmd.txt")));

  // Set up profiler
  A2D::TIMER_OUTPUT_FILE = fspath(prefix) / fspath("profile.log");
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

  // Find bc vertices
  std::vector<int> ids{100,
                       101};  // fixed both left and right faces of the domain
  std::vector<I> verts = readvtk.get_verts_given_cell_entity_id(ids);
  I bc_label = conn.add_boundary_label_from_verts(verts.size(), verts.data());

  // Find traction vertices
  T lower[3], upper[3];
  readvtk.get_bounds(lower, upper);
  std::vector<I> traction_verts = readvtk.get_verts_within_box(
      lower[0], upper[0], upper[1], upper[1], lower[2], upper[2]);
  I traction_label = conn.add_boundary_label_from_verts(traction_verts.size(),
                                                        traction_verts.data());

  // Set traction components and body force components
  T t[3] = {0.0, -1.0, 0.0};
  T tb[3] = {0.0, 0.0, 0.0};

  // Write traction verts to vtk for debugging purpose
  std::vector<T> labels(3 * nverts, 0.0);
  for (auto it = verts.begin(); it != verts.end(); it++) {
    labels[3 * (*it)] = 1.0;  // bc verts
  }
  for (auto it = traction_verts.begin(); it != traction_verts.end(); it++) {
    labels[3 * (*it) + 1] = 1.0;  // traction verts
  }
  {
    A2D::VectorFieldToVTK fieldtovtk(
        nverts, Xloc, labels.data(),
        fspath(prefix) / fspath("bc_traction_verts.vtk"));
  }

  // Set up boundary condition information
  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label, basis);

  // Initialize the analysis instance
  T E = 70.0e3, nu = 0.3;
  TopoElasticityAnalysis<T, degree, filter_degree> topo(
      conn, bcinfo, E, nu, ramp_q, tb, traction_label, t, verbose, amg_nlevels,
      cg_it, cg_rtol, cg_atol);
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
                     domain_vol, volume_frac, topo, verbose, vtk_freq);
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

  std::string out_path = fspath(prefix) / fspath("paropt.out");
  std::string tr_path = fspath(prefix) / fspath("paropt.tr");
  std::string mma_path = fspath(prefix) / fspath("paropt.mma");
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