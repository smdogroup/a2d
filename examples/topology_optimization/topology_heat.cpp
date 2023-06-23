#include "topology_heat.h"

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;
using fspath = std::filesystem::path;

/**
 * @brief Helper function, save bc verts to vtk to verify visually.
 */
void verts_to_vtk(I nverts, std::vector<I> &bc_verts, T *Xloc,
                  std::string prefix) {
  std::vector<T> labels(3 * nverts, 0.0);
  for (auto it = bc_verts.begin(); it != bc_verts.end(); it++) {
    labels[3 * (*it)] = 1.0;  // bc verts
  }
  {
    A2D::VectorFieldToVTK fieldtovtk(nverts, Xloc, labels.data(),
                                     fspath(prefix) / fspath("bc_verts.vtk"));
  }
}

template <int degree>
void main_body(std::string prefix, double bc_temp, double heat_source,
               double ratio, int nx, int ny, int nz, double lx, double ly,
               double lz, int amg_nlevels, int cg_it, double cg_rtol,
               double cg_atol, bool verbose, double vol_frac, int maxit,
               int vtk_freq, double ramp_q, bool check_grad_and_exit) {
  // Set the lower order degree for the Bernstein polynomial
  constexpr int filter_degree = degree - 1;

  // Set up profiler
  A2D::Timer::set_log_path(fspath(prefix) / fspath("profile.log"));
  A2D::Timer timer("main_body()");

  // Create the mesh by code
  I nverts = (nx + 1) * (ny + 1) * (nz + 1);
  I nhex = nx * ny * nz;

  // helper functor
  auto node_num = [&](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  // Allocate temporary arrays
  I hex[8 * nhex];
  T Xloc[3 * nverts];

  // Populate hex
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

  // Populate Xloc
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (lx * i) / nx;
        Xloc[3 * node_num(i, j, k) + 1] = (ly * j) / ny;
        Xloc[3 * node_num(i, j, k) + 2] = (lz * k) / nz;
      }
    }
  }

  // Create A2D's connectivity object
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Find bc verts
  std::vector<I> bc_verts = A2D::get_verts_within_box(
      nverts, Xloc, 0.0, 0.0, (1.0 - ratio) / 2.0 * ly,
      (1.0 + ratio) / 2.0 * ly, (1.0 - ratio) / 2.0 * lz,
      (1.0 + ratio) / 2.0 * lz);
  I bc_label =
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data());

  // Save bcverts to vtk
  verts_to_vtk(nverts, bc_verts, Xloc, prefix);

  // Set up boundary condition information
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label);

  // Initialize the analysis instance
  T kappa = 1.0;
  TopoHeatAnalysis<T, degree, filter_degree> analysis(
      conn, bcinfo, kappa, heat_source, bc_temp, ramp_q, verbose, amg_nlevels,
      cg_it, cg_rtol, cg_atol);

  auto elem_geo = analysis.get_geometry();
  A2D::set_geo_from_hex_nodes<
      typename TopoHeatAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  analysis.reset_geometry();

  // Get problem size
  int nvars = analysis.get_num_design_vars();
  int ndof = analysis.get_num_dofs();
  int nelems = analysis.get_num_elements();

  // Print info
  std::printf("basis degree:                 %d\n", degree);
  std::printf("number of mesh elements:      %d\n", nelems);
  std::printf("number of design variables:   %d\n", nvars);
  std::printf("number of degrees of freedom: %d\n", ndof);

  // Initialize the optimization object
  int ncon = 1;
  int nineq = 1;
  analysis.set_design_vars(std::vector<T>(nvars, T(1.0)));
  analysis.solve();
  T comp = analysis.eval_compliance();
  T ref_hcomp = analysis.eval_compliance();
  T domain_vol = analysis.eval_volume();
  TopOptProb<TopoHeatAnalysis<T, degree, filter_degree>> *prob =
      new TopOptProb<TopoHeatAnalysis<T, degree, filter_degree>>(
          prefix, MPI_COMM_WORLD, nvars, ncon, nineq, ref_hcomp, domain_vol,
          vol_frac, analysis, verbose, vtk_freq);
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

  std::string out_path = fspath(prefix) / fspath("paropt.out");
  std::string tr_path = fspath(prefix) / fspath("paropt.tr");
  std::string mma_path = fspath(prefix) / fspath("paropt.mma");
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
    // Initialize argument parser
    A2D::ArgumentParser parser(argc, argv);

    // Set up cmd arguments and defaults

    // - Basics: basis degree and result prefix
    int degree = parser.parse_option("--degree", 2);
    std::string prefix =
        parser.parse_option("--prefix", std::string("results"));

    // - Problem specific settings
    double bc_temp = parser.parse_option("--bc_temp", 0.0);
    double heat_source = parser.parse_option("--heat_source", 1.0);
    double ratio = parser.parse_option("--ratio", 0.2);

    // - Domain specific settings
    int nx = parser.parse_option("--nx", 10);
    int ny = parser.parse_option("--ny", 10);
    int nz = parser.parse_option("--nz", 10);
    double lx = parser.parse_option("--lx", 1.0);
    double ly = parser.parse_option("--ly", 1.0);
    double lz = parser.parse_option("--lz", 1.0);

    // - Linear solver settings
    int amg_nlevels = parser.parse_option("--amg_nlevels", 3);
    int cg_it = parser.parse_option("--cg_it", 400);
    double cg_rtol = parser.parse_option("--cg_rtol", 1e-8);
    double cg_atol = parser.parse_option("--cg_atol", 1e-30);
    bool verbose = parser.parse_option("--verbose");

    // -- Optimization settings
    double vol_frac = parser.parse_option("--vol_frac", 0.4);
    int maxit = parser.parse_option("--maxit", 400);
    int vtk_freq = parser.parse_option("--vtk_freq", 10);
    double ramp_q = parser.parse_option("--ramp_q", 5.0);
    bool check_grad_and_exit = parser.parse_option("--check_grad_and_exit");

    // - Print info with invoked with -h or --help
    parser.help_info();

    // Set up result directory
    if (!std::filesystem::is_directory(prefix)) {
      std::filesystem::create_directory(prefix);
    }

    // Save cmd arguments to txt
    A2D::save_cmd(argc, argv, (fspath(prefix) / fspath("cmd.txt")));

    // Execute
    switch (degree) {
      case 2:
        main_body<2>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 3:
        main_body<3>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 4:
        main_body<4>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 5:
        main_body<5>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 6:
        main_body<6>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 7:
        main_body<7>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 8:
        main_body<8>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 9:
        main_body<9>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                     lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                     vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 10:
        main_body<10>(prefix, bc_temp, heat_source, ratio, nx, ny, nz, lx, ly,
                      lz, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose,
                      vol_frac, maxit, vtk_freq, ramp_q, check_grad_and_exit);
        break;

      default:
        std::printf("Specified degree (%d) is not precompiled.\n", degree);
        break;
    }
  }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}