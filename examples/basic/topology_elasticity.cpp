#include "topology_elasticity.h"

#include "topology_paropt_prob.h"
#include "utils/a2dvtk.h"

using I = A2D::index_t;
using T = ParOptScalar;

void main_body(int argc, char *argv[]) {
  constexpr int degree = 3;
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
  T E = 70.0e3, nu = 0.3, q = 5.0;
  TopoAnalysis<T, degree> topo(conn, bcinfo, E, nu, q);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<TopoAnalysis<T, degree>::GeoBasis>(
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
  TopOptProb<degree> *prob =
      new TopOptProb<degree>(MPI_COMM_WORLD, nvars, ncon, nineq, ref_comp,
                             domain_vol, volume_frac, topo);
  prob->incref();

  // Sanity check
  prob->checkGradients(1e-6);
  exit(0);

  // Define paropt optimizer
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
  options->setOption("mma_max_iterations", 200);
  options->setOption("mma_min_asymptote_offset", 0.01);
  options->setOption("mma_use_constraint_linearization", true);

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