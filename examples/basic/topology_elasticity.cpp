#include "topology_elasticity.h"

#include <cmath>

#include "topology_paropt_prob.h"
#include "utils/a2dparser.h"

using I = A2D::index_t;
using T = ParOptScalar;
using fspath = std::filesystem::path;

/**
 * @brief Helper function, save bc and traction verts to vtk to verify visually.
 */
void verts_to_vtk(I nverts, std::vector<I> &bc_verts,
                  std::vector<I> &traction_verts, T *Xloc, std::string prefix) {
  // Write traction verts to vtk for debugging purpose
  std::vector<T> labels(3 * nverts, 0.0);
  for (auto it = bc_verts.begin(); it != bc_verts.end(); it++) {
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
}

/**
 * @brief Helper function, set high order mesh for cylinder problem instance
 */
template <class GeoBasis, class GeoElemVec>
void set_geo_cylindrical(int nhex, double rout, double rin, double height,
                         int nelems_c, int nelems_r, int nelems_h,
                         GeoElemVec &elem_geo) {
  constexpr double pi = 3.141592653589793238462643383279502884;

  for (int e = 0; e < nhex; e++) {
    typename GeoElemVec::FEDof geo_dof(e, elem_geo);
    // Get the reference vertex of the element
    I ic = e % nelems_c;
    I ir = (e / nelems_c) / nelems_h;
    I ih = (e / nelems_c) % nelems_h;

    double theta = (T)ic / nelems_c * 2.0 * pi;
    double r = rout - T(ir + 1) / nelems_r * (rout - rin);
    double arc = r * 2.0 * pi / nelems_c;
    double width = (rout - rin) / nelems_r;
    double h = T(ih) / nelems_h * height;

    for (int ii = 0; ii < GeoBasis::ndof; ii++) {
      double pt[3];
      GeoBasis::get_dof_point(ii, pt);  // Get the ref coordinates for each dof

      double xi = pt[0], eta = pt[1], zeta = pt[2];

      double dr = (xi + 1.0) / 2.0 * width;
      double dtheta = (eta + 1.0) / 2.0 * arc / r;
      double dh = (zeta + 1.0) / 2.0 * height / nelems_h;

      // Set high order node coordinates
      switch (ii % 3) {
        case 0:  // x coordinate
          geo_dof[ii] = (r + dr) * std::cos(theta + dtheta);
          break;

        case 1:  // y coordinate
          geo_dof[ii] = (r + dr) * std::sin(theta + dtheta);
          break;

        case 2:  // z coordinate
          geo_dof[ii] = h + dh;
          break;
      }
    }
    elem_geo.set_element_values(e, geo_dof);
  }
}

/**
 * @brief Generate the analysis object based on vtk mesh input
 */
template <int degree, int filter_degree>
std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>>
create_analysis_vtk(std::string prefix, std::string vtk_name,
                    double vb_traction_frac, int amg_nlevels, int cg_it,
                    double cg_rtol, double cg_atol, bool verbose, int maxit,
                    int vtk_freq, double ramp_q, bool check_grad_and_exit) {
  // Load vtk
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
  std::vector<int> ids{100};  // Fix the left surface of domain
  std::vector<I> bc_verts = readvtk.get_verts_given_cell_entity_id(ids);
  I bc_label =
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data());

  // Find traction vertices
  T lower[3], upper[3];
  readvtk.get_bounds(lower, upper);
  T xmin = upper[0];
  T xmax = upper[0];
  T ymin = lower[1];
  T ymax = lower[1] + vb_traction_frac * (upper[1] - lower[1]);
  T zmin = lower[2];
  T zmax = upper[2];
  std::vector<I> traction_verts = A2D::get_verts_within_box(
      nverts, Xloc, xmin, xmax, ymin, ymax, zmin, zmax);
  I traction_label = conn.add_boundary_label_from_verts(traction_verts.size(),
                                                        traction_verts.data());

  // Save bc/traction verts to vtk
  verts_to_vtk(nverts, bc_verts, traction_verts, Xloc, prefix);

  // Set traction components and body force components
  T tx_traction[3] = {0.0, -1.0, 0.0};
  T tx_body[3] = {0.0, 0.0, 0.0};

  // Set up boundary condition information
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label);

  // Initialize the analysis instance
  T E = 70.0e3, nu = 0.3;
  std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>> analysis =
      std::make_shared<TopoElasticityAnalysis<T, degree, filter_degree>>(
          conn, bcinfo, E, nu, ramp_q, tx_body, traction_label, tx_traction,
          nullptr, nullptr, verbose, amg_nlevels, cg_it, cg_rtol, cg_atol);

  auto elem_geo = analysis->get_geometry();
  A2D::set_geo_from_hex_nodes<
      typename TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  analysis->reset_geometry();

  return analysis;
}

/**
 * @brief Create an analysis object with a cuboid domain
 */
template <int degree, int filter_degree>
std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>>
create_analysis_box(std::string prefix, double vb_traction_frac, int b_nx,
                    int b_ny, int b_nz, double b_lx, double b_ly, double b_lz,
                    int amg_nlevels, int cg_it, double cg_rtol, double cg_atol,
                    bool verbose, int maxit, int vtk_freq, double ramp_q,
                    bool check_grad_and_exit) {
  I nverts = (b_nx + 1) * (b_ny + 1) * (b_nz + 1);
  I nhex = b_nx * b_ny * b_nz;

  // helper functor
  auto node_num = [&](int i, int j, int k) {
    return i + j * (b_nx + 1) + k * (b_nx + 1) * (b_ny + 1);
  };

  // Compute number of vertices and hex elements
  nverts = (b_nx + 1) * (b_ny + 1) * (b_nz + 1);
  nhex = b_nx * b_ny * b_nz;

  // Allocate temporary arrays
  I hex[8 * nhex];
  T Xloc[3 * nverts];

  // Populate hex
  using ET = A2D::ElementTypes;
  for (int k = 0, e = 0; k < b_nz; k++) {
    for (int j = 0; j < b_ny; j++) {
      for (int i = 0; i < b_nx; i++, e++) {
        for (int ii = 0; ii < ET::HEX_VERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  // Populate Xloc
  for (int k = 0; k < b_nz + 1; k++) {
    for (int j = 0; j < b_ny + 1; j++) {
      for (int i = 0; i < b_nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (b_lx * i) / b_nx;
        Xloc[3 * node_num(i, j, k) + 1] = (b_ly * j) / b_ny;
        Xloc[3 * node_num(i, j, k) + 2] = (b_lz * k) / b_nz;
      }
    }
  }

  // Create A2D's connectivity object
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Find bc verts
  std::vector<I> bc_verts =
      A2D::get_verts_within_box(nverts, Xloc, 0.0, 0.0, 0.0, b_ly, 0.0, b_lz);
  I bc_label =
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data());

  // Find traction verts
  std::vector<I> traction_verts = A2D::get_verts_within_box(
      nverts, Xloc, b_lx, b_lx, 0.0, vb_traction_frac * b_ly, 0.0, b_lz);
  I traction_label = conn.add_boundary_label_from_verts(traction_verts.size(),
                                                        traction_verts.data());

  // Set traction components and body force components
  T tx_traction[3] = {0.0, -1.0, 0.0};
  T tx_body[3] = {0.0, 0.0, 0.0};

  // Save bc/traction verts to vtk
  verts_to_vtk(nverts, bc_verts, traction_verts, Xloc, prefix);

  // Set up boundary condition information
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label);

  // Initialize the analysis instance
  T E = 70.0e3, nu = 0.3;
  std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>> analysis =
      std::make_shared<TopoElasticityAnalysis<T, degree, filter_degree>>(
          conn, bcinfo, E, nu, ramp_q, tx_body, traction_label, tx_traction,
          nullptr, nullptr, verbose, amg_nlevels, cg_it, cg_rtol, cg_atol);

  auto elem_geo = analysis->get_geometry();
  A2D::set_geo_from_hex_nodes<
      typename TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, hex, Xloc, elem_geo);
  analysis->reset_geometry();

  // Save mesh
  A2D::ToVTK3D(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd, pyrmd,
               Xloc, fspath(prefix) / fspath("box_low_order_mesh.vtk"));
  analysis->tovtk(fspath(prefix) / fspath("box_high_order_mesh.vtk"));

  return analysis;
}

/**
 * @brief Create an analysis object with a hollow cylinder domain
 *
 * @param rout the outer radius
 * @param rin the inner radius, 0 < rin < rout
 * @param height the height
 * @param nelems_c number of elements around circumference
 * @param nelems_e number of elements along radial direction
 * @param nelems_h number of elements along height direction
 */
template <int degree, int filter_degree>
std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>>
create_analysis_cylinder(std::string prefix, double rout, double rin,
                         double height, int nelems_c, int nelems_r,
                         int nelems_h, int amg_nlevels, int cg_it,
                         double cg_rtol, double cg_atol, bool verbose,
                         int maxit, int vtk_freq, double ramp_q,
                         bool check_grad_and_exit) {
  // constants
  constexpr double pi = 3.141592653589793238462643383279502884;
  using ET = A2D::ElementTypes;

  // Dimension data
  I nhex = nelems_c * nelems_r * nelems_h;
  I nverts = nelems_c * (nelems_r + 1) * (nelems_h + 1);

  // Construct element connectivity
  std::vector<I> hex(nhex * ET::HEX_VERTS);
  I diff = nelems_c * (nelems_h + 1);
  I elem_idx = 0;
  for (I r = 0; r < nelems_r; r++) {
    for (I i = 0; i < nelems_h; i++) {
      for (I j = 0, jc = 1; j < nelems_c; j++, elem_idx++, jc++) {
        I ref_vert = r * nelems_c * (nelems_h + 1) + i * nelems_c;
        hex[ET::HEX_VERTS * elem_idx] = ref_vert + j + diff;
        hex[ET::HEX_VERTS * elem_idx + 1] = ref_vert + j;
        hex[ET::HEX_VERTS * elem_idx + 2] = ref_vert + jc % nelems_c;
        hex[ET::HEX_VERTS * elem_idx + 3] = ref_vert + jc % nelems_c + diff;

        ref_vert = r * nelems_c * (nelems_h + 1) + (i + 1) * nelems_c;
        hex[ET::HEX_VERTS * elem_idx + 4] = ref_vert + j + diff;
        hex[ET::HEX_VERTS * elem_idx + 5] = ref_vert + j;
        hex[ET::HEX_VERTS * elem_idx + 6] = ref_vert + jc % nelems_c;
        hex[ET::HEX_VERTS * elem_idx + 7] = ref_vert + jc % nelems_c + diff;
      }
    }
  }

  // Construct nodal location
  std::vector<T> Xloc(nverts * 3);
  I vert_idx = 0;
  for (I r = 0; r < nelems_r + 1; r++) {  // from outer to inner
    for (I i = 0; i < nelems_h + 1; i++) {
      for (I j = 0; j < nelems_c; j++, vert_idx++) {
        T deg = 2 * pi * (T)j / nelems_c;
        T rad = rout - (T)r / nelems_r * (rout - rin);
        Xloc[3 * vert_idx] = std::cos(deg) * rad;
        Xloc[3 * vert_idx + 1] = std::sin(deg) * rad;
        Xloc[3 * vert_idx + 2] = T(i) / nelems_h * height;
      }
    }
  }

  // Construct connectivity
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex.data(), nwedge,
                               wedge, npyrmd, pyrmd);
  A2D::ToVTK3D(nverts, ntets, tets, nhex, hex.data(), nwedge, wedge, npyrmd,
               pyrmd, Xloc.data(),
               fspath(prefix) / fspath("cylinder_low_order_mesh.vtk"));

  // Set points to apply bc and traction
  double dtheta =
      pi / 15.0;  // angle within witch to apply torque around the circumference
  std::vector<double> thetas = {0.0, 0.5 * pi, pi, 1.5 * pi};

  // Find traction vertices
  std::vector<I> traction_verts;
  for (double theta : thetas) {
    std::vector<I> tor_v = A2D::get_verts_cylindrical_coords(
        nverts, Xloc.data(), theta, theta + dtheta, rin, rout, height, height,
        1e-3);
    traction_verts.insert(traction_verts.end(), tor_v.begin(), tor_v.end());
  }
  I traction_label = conn.add_boundary_label_from_verts(traction_verts.size(),
                                                        traction_verts.data());

  // Find bc vertices - use same pattern as traction vertices
  std::vector<I> bc_verts;
  for (double theta : thetas) {
    std::vector<I> bc_v = A2D::get_verts_cylindrical_coords(
        nverts, Xloc.data(), theta, theta + dtheta, rin, rout, 0.0, 0.0, 1e-3);
    bc_verts.insert(bc_verts.end(), bc_v.begin(), bc_v.end());
  }

  I bc_label =
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data());

  // Save bc/traction verts to vtk
  verts_to_vtk(nverts, bc_verts, traction_verts, Xloc.data(), prefix);

  // Set traction components and body force components
  T tx_traction[3] = {0.0, 0.0, 0.0};
  T tx_body[3] = {0.0, 0.0, 0.0};
  T tx_torque[3] = {0.0, 0.0, 1.0};
  T x0_torque[3] = {0.0, 0.0, height};

  // Set up boundary condition information
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label);

  // Initialize the analysis instance
  T E = 70.0e3, nu = 0.3;
  std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>> analysis =
      std::make_shared<TopoElasticityAnalysis<T, degree, filter_degree>>(
          conn, bcinfo, E, nu, ramp_q, tx_body, traction_label, tx_traction,
          x0_torque, tx_torque, verbose, amg_nlevels, cg_it, cg_rtol, cg_atol);

  // Set high order element geometry and save mesh to vtk
  auto elem_geo = analysis->get_geometry();
  set_geo_cylindrical<
      typename TopoElasticityAnalysis<T, degree, filter_degree>::GeoBasis>(
      nhex, rout, rin, height, nelems_c, nelems_r, nelems_h, elem_geo);
  analysis->reset_geometry();
  analysis->tovtk(fspath(prefix) / fspath("cylinder_high_order_mesh.vtk"));

  return analysis;
}

template <int degree>
void main_body(std::string prefix, std::string domain, std::string vtk_name,
               double vb_traction_frac, int b_nx, int b_ny, int b_nz,
               double b_lx, double b_ly, double b_lz, double cy_rout,
               double cy_rin, double cy_height, int cy_nelems_c,
               int cy_nelems_r, int cy_nelems_h, int amg_nlevels, int cg_it,
               double cg_rtol, double cg_atol, bool verbose,
               std::string optimizer, double vol_frac, int maxit, int vtk_freq,
               double ramp_q, bool check_grad_and_exit) {
  // Set the lower order degree for the Bernstein polynomial
  constexpr int filter_degree = degree - 1;

  // Set up profiler
  A2D::TIMER_OUTPUT_FILE = fspath(prefix) / fspath("profile.log");
  A2D::Timer timer("main_body()");

  // Create mesh: either load vtk or generate by code
  std::shared_ptr<TopoElasticityAnalysis<T, degree, filter_degree>> analysis;

  if (domain == "vtk") {
    analysis = create_analysis_vtk<degree, filter_degree>(
        prefix, vtk_name, vb_traction_frac, amg_nlevels, cg_it, cg_rtol,
        cg_atol, verbose, maxit, vtk_freq, ramp_q, check_grad_and_exit);
  } else if (domain == "box") {
    analysis = create_analysis_box<degree, filter_degree>(
        prefix, vb_traction_frac, b_nx, b_ny, b_nz, b_lx, b_ly, b_lz,
        amg_nlevels, cg_it, cg_rtol, cg_atol, verbose, maxit, vtk_freq, ramp_q,
        check_grad_and_exit);
  } else if (domain == "cylinder") {
    analysis = create_analysis_cylinder<degree, filter_degree>(
        prefix, cy_rout, cy_rin, cy_height, cy_nelems_c, cy_nelems_r,
        cy_nelems_h, amg_nlevels, cg_it, cg_rtol, cg_atol, verbose, maxit,
        vtk_freq, ramp_q, check_grad_and_exit);
  } else {
    std::printf("Invalid domain: %s!\n", domain.c_str());
    exit(-1);
  }

  // Get problem size
  int nvars = analysis->get_num_design_vars();
  int ndof = analysis->get_num_dofs();
  int nelems = analysis->get_num_elements();

  // Print info
  std::printf("basis degree:                 %d\n", degree);
  std::printf("number of mesh elements:      %d\n", nelems);
  std::printf("number of design variables:   %d\n", nvars);
  std::printf("number of degrees of freedom: %d\n", ndof);

  // Initialize the optimization object
  int ncon = 1;
  int nineq = 1;
  analysis->set_design_vars(std::vector<T>(nvars, T(1.0)));
  analysis->solve();
  T ref_comp = analysis->eval_compliance();
  T domain_vol = analysis->eval_volume();
  TopOptProb<TopoElasticityAnalysis<T, degree, filter_degree>> *prob =
      new TopOptProb<TopoElasticityAnalysis<T, degree, filter_degree>>(
          prefix, MPI_COMM_WORLD, nvars, ncon, nineq, ref_comp, domain_vol,
          vol_frac, *analysis, verbose, vtk_freq);
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
  options->setOption("output_file", out_path.c_str());
  options->setOption("tr_output_file", tr_path.c_str());
  options->setOption("mma_output_file", mma_path.c_str());

  if (optimizer == "mma") {
    options->setOption("algorithm", "mma");
    options->setOption("mma_max_iterations", maxit);
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
    options->setOption("mma_use_constraint_linearization", true);
  } else if (optimizer == "tr") {
    options->setOption("algorithm", "tr");
    options->setOption("output_level", 0);
    options->setOption("norm_type", "l1");
    options->setOption("tr_init_size", 0.05);
    options->setOption("tr_min_size", 1e-3);
    options->setOption("tr_max_size", 1.0);
    options->setOption("tr_eta", 0.25);
    options->setOption("tr_infeas_tol", 1e-6);
    options->setOption("tr_l1_tol", 0.0);
    options->setOption("tr_linfty_tol", 0.0);
    options->setOption("tr_adaptive_gamma_update", true);
    options->setOption("tr_accept_step_strategy", "penalty_method");
    options->setOption("filter_sufficient_reduction", true);
    options->setOption("filter_has_feas_restore_phase", true);
    options->setOption("tr_use_soc", false);
    options->setOption("tr_max_iterations", maxit);
    options->setOption("penalty_gamma", 50.0);
    options->setOption("qn_subspace_size", 5);
    options->setOption("qn_type", "bfgs");
    options->setOption("qn_diag_type", "yty_over_yts");
    options->setOption("abs_res_tol", 1e-8);
    options->setOption("starting_point_strategy", "affine_step");
    options->setOption("barrier_strategy", "mehrotra_predictor_corrector");
    options->setOption("tr_steering_barrier_strategy",
                       "mehrotra_predictor_corrector");
    options->setOption("tr_steering_starting_point_strategy", "affine_step");
    options->setOption("use_line_search", false);
    options->setOption("max_major_iters", 20);
  }

  ParOptOptimizer *opt = new ParOptOptimizer(prob, options);
  opt->incref();
  opt->optimize();

  // Delete objects
  prob->decref();
  options->incref();
  opt->decref();
  return;
}

/**
 * @brief Check if target equals one of valid_vals
 */
template <typename EntryType>
void assert_option_in(std::string target, std::vector<EntryType> valid_vals) {
  auto domain_it = std::find(valid_vals.begin(), valid_vals.end(), target);
  if (domain_it == valid_vals.end()) {
    std::printf("Agrument value %s is invalid! Valid options are: ",
                target.c_str());
    for (auto it = valid_vals.begin(); it != valid_vals.end(); it++) {
      std::printf("%s, ", it->c_str());
    }
    std::printf("\b\b.\n");
    exit(-1);
  }
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

    // - Problem specific settings: mesh, domain, bc, parameter, etc.
    std::string domain = parser.parse_option("--domain", std::string("vtk"));
    std::string vtk_name =
        parser.parse_option("--vtk", std::string("3d_hex.vtk"));

    // - Domain specific settings - vtk or box
    double vb_traction_frac = parser.parse_option("--vb_traction_frac", 0.2);
    int b_nx = parser.parse_option("--b_nx", 10);
    int b_ny = parser.parse_option("--b_ny", 5);
    int b_nz = parser.parse_option("--b_nz", 3);
    double b_lx = parser.parse_option("--b_lx", 1.0);
    double b_ly = parser.parse_option("--b_ly", 0.5);
    double b_lz = parser.parse_option("--b_lz", 0.3);

    // - Domain specific settings - cylinder
    double cy_rout = parser.parse_option("--cy_rout", 1.0);
    double cy_rin = parser.parse_option("--cy_rin", 0.8);
    double cy_height = parser.parse_option("--cy_height", 4.0);
    int cy_nelems_c = parser.parse_option("--cy_nelems_c", 32);
    int cy_nelems_r = parser.parse_option("--cy_nelems_r", 1);
    int cy_nelems_h = parser.parse_option("--cy_nelems_h", 20);

    // - Linear solver settings
    int amg_nlevels = parser.parse_option("--amg_nlevels", 3);
    int cg_it = parser.parse_option("--cg_it", 400);
    double cg_rtol = parser.parse_option("--cg_rtol", 1e-8);
    double cg_atol = parser.parse_option("--cg_atol", 1e-30);
    bool verbose = parser.parse_option("--verbose");

    // - Optimization settings
    std::string optimizer =
        parser.parse_option("--optimizer", std::string("mma"));
    double vol_frac = parser.parse_option("--vol_frac", 0.3);
    int maxit = parser.parse_option("--maxit", 400);
    int vtk_freq = parser.parse_option("--vtk_freq", 10);
    double ramp_q = parser.parse_option("--ramp_q", 5.0);
    bool check_grad_and_exit = parser.parse_option("--check_grad_and_exit");

    // - Print info with invoked with -h or --help
    parser.help_info();

    // Check values
    std::vector<std::string> valid_domains = {"box", "cylinder", "vtk"};
    assert_option_in(domain, valid_domains);
    std::vector<std::string> valid_opt = {"mma", "tr"};
    assert_option_in(optimizer, valid_opt);

    // Set up result directory
    if (!std::filesystem::is_directory(prefix)) {
      std::filesystem::create_directory(prefix);
    }

    // Save cmd arguments to txt
    A2D::save_cmd(argc, argv, (fspath(prefix) / fspath("cmd.txt")));

    // Execute
    switch (degree) {
      case 2:
        main_body<2>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 3:
        main_body<3>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 4:
        main_body<4>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 5:
        main_body<5>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 6:
        main_body<6>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 7:
        main_body<7>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 8:
        main_body<8>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 9:
        main_body<9>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                     b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                     cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                     cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                     vtk_freq, ramp_q, check_grad_and_exit);
        break;

      case 10:
        main_body<10>(prefix, domain, vtk_name, vb_traction_frac, b_nx, b_ny,
                      b_nz, b_lx, b_ly, b_lz, cy_rout, cy_rin, cy_height,
                      cy_nelems_c, cy_nelems_r, cy_nelems_h, amg_nlevels, cg_it,
                      cg_rtol, cg_atol, verbose, optimizer, vol_frac, maxit,
                      vtk_freq, ramp_q, check_grad_and_exit);
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