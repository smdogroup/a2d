#include <cassert>

#include "a2d.h"
#include "amgx_c.h"
#include "cuda_runtime.h"
#include "toolkit.h"
#include "ticktock.h"

using namespace A2D;

#define MAX_MSG_LEN 4096

/* print error message and exit */
void errAndExit(const char* err) {
  printf("%s\n", err);
  fflush(stdout);
  exit(1);
}

/* print callback (could be customized) */
void print_callback(const char* msg, int length) {
    printf("%s", msg);
}

/* parse parameters */
int findParamIndex(char** argv, int argc, const char* parm) {
  int count = 0;
  int index = -1;

  for (int i = 0; i < argc; i++) {
    if (strncmp(argv[i], parm, 100) == 0) {
      index = i;
      count++;
    }
  }

  if (count == 0 || count == 1) {
    return index;
  } else {
    char msg[MAX_MSG_LEN];
    sprintf(msg,
            "ERROR: parameter %s has been specified more than once, exiting\n",
            parm);
    print_callback(msg, MAX_MSG_LEN);
    exit(1);
  }

  return -1;
}

typedef index_t I;
typedef double T;
static const int SPATIAL_DIM = 3;
typedef BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature> Basis;
typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> ElasticityPDE;

template <class ConnArray, class XArray>
void populate_X_conn_3d(const int nx, const int ny, const int nz, XArray& X,
                        ConnArray& conn) {
  // Set X
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        X(node, 0) = 1.0 * i / nx;
        X(node, 1) = 1.0 * j / ny;
        X(node, 2) = 1.0 * k / nz;
      }
    }
  }

  // Set connectivity
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int elem = i + nx * (j + ny * k);

        int conn_coord[8];
        for (int kk = 0, index = 0; kk < 2; kk++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int ii = 0; ii < 2; ii++, index++) {
              conn_coord[index] =
                  (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
            }
          }
        }

        // Convert to the correct connectivity
        conn(elem, 0) = conn_coord[0];
        conn(elem, 1) = conn_coord[1];
        conn(elem, 2) = conn_coord[3];
        conn(elem, 3) = conn_coord[2];

        conn(elem, 4) = conn_coord[4];
        conn(elem, 5) = conn_coord[5];
        conn(elem, 6) = conn_coord[7];
        conn(elem, 7) = conn_coord[6];
      }
    }
  }
}

template <class BcsArray>
void populate_bcs_3d(const int nx, const int ny, const int nz, BcsArray& bcs) {
  index_t index = 0;
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      int i = 0;
      int node = i + (nx + 1) * (j + (ny + 1) * k);

      // Set the boundary conditions
      bcs(index, 0) = node;
      for (int ii = 0; ii < 3; ii++) {
        bcs(index, 1) |= 1U << ii;
      }
      index++;
    }
  }
}

template <class DvArray>
void set_dv_3d(const int nx, const int ny, const int nz, DvArray& x) {
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);
        x(node, 0) = 1.0;
        // if (i < nx / 2 and j < ny / 2 and k < nz / 2) {
        //   x(node, 0) = 1e-3;
        // } else {
        //   x(node, 0) = 1.0;
        // }
      }
    }
  }
}

template <class Model, class RhsArray>
void set_force_3d(const int nx, const int ny, const int nz, Model model,
                  RhsArray& residual) {
  residual->zero();
  for (int k = nz / 4; k < 3 * nz / 4; k++) {
    int node = nx + (nx + 1) * (0 + (ny + 1) * k);
    (*residual)(node, 1) = -1e2;
  }
  model->zero_bcs(residual);
}

// input contians the argv and argc from main
template <class Matrix, class Vector>
int amgx_solve(Matrix& mat, Vector& rhs, int argc, char** argv) {
  // Init
  // parameters parsing and initialization
  int pidx = 0;
  int pidy = 0;

  AMGX_Mode mode;
  AMGX_config_handle cfg;
  AMGX_resources_handle rsrc;
  AMGX_matrix_handle A;
  AMGX_vector_handle b, x;
  AMGX_solver_handle solver;

  // status handling
  AMGX_SOLVE_STATUS status;

  AMGX_SAFE_CALL(AMGX_initialize());
  /* system */
  AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
  AMGX_SAFE_CALL(AMGX_install_signal_handler());

  /* get mode */
  if ((pidx = findParamIndex(argv, argc, "-mode")) != -1) {
    if (strncmp(argv[pidx + 1], "hDDI", 100) == 0) {
      mode = AMGX_mode_hDDI;
    } else if (strncmp(argv[pidx + 1], "hDFI", 100) == 0) {
      mode = AMGX_mode_hDFI;
    } else if (strncmp(argv[pidx + 1], "hFFI", 100) == 0) {
      mode = AMGX_mode_hFFI;
    } else if (strncmp(argv[pidx + 1], "dDDI", 100) == 0) {
      mode = AMGX_mode_dDDI;
    } else if (strncmp(argv[pidx + 1], "dDFI", 100) == 0) {
      mode = AMGX_mode_dDFI;
    } else if (strncmp(argv[pidx + 1], "dFFI", 100) == 0) {
      mode = AMGX_mode_dFFI;
    } else if (strncmp(argv[pidx + 1], "hCCI", 100) == 0) {
      mode = AMGX_mode_hZZI;
    } else if (strncmp(argv[pidx + 1], "hZCI", 100) == 0) {
      mode = AMGX_mode_hZCI;
    } else if (strncmp(argv[pidx + 1], "hZZI", 100) == 0) {
      mode = AMGX_mode_hZZI;
    } else if (strncmp(argv[pidx + 1], "dCCI", 100) == 0) {
      mode = AMGX_mode_dCCI;
    } else if (strncmp(argv[pidx + 1], "dZCI", 100) == 0) {
      mode = AMGX_mode_dZCI;
    } else if (strncmp(argv[pidx + 1], "dZZI", 100) == 0) {
      mode = AMGX_mode_dZZI;
    } else {
      errAndExit("ERROR: invalid mode");
    }
  } else {
    printf("Warning: No mode specified, using dDDI by default.\n");
    mode = AMGX_mode_dDDI;
  }

  /* get max_it number of outer (non-linear) iteration */

  if ((pidx = findParamIndex(argv, argc, "-it")) != -1) {
    if ((pidy = findParamIndex(argv, argc, "-gpu")) != -1) {
      errAndExit(
          "ERROR: -gpu and -it options are not compatible, you must choose one "
          "or the other option");
    }
    srand(0);
  }

  /* create config */
  pidx = findParamIndex(argv, argc, "-amg");
  pidy = findParamIndex(argv, argc, "-c");

  if ((pidx != -1) && (pidy != -1)) {
    printf("%s\n", argv[pidx + 1]);
    AMGX_SAFE_CALL(AMGX_config_create_from_file_and_string(&cfg, argv[pidy + 1],
                                                           argv[pidx + 1]));
  } else if (pidy != -1) {
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, argv[pidy + 1]));
  } else if (pidx != -1) {
    printf("%s\n", argv[pidx + 1]);
    AMGX_SAFE_CALL(AMGX_config_create(&cfg, argv[pidx + 1]));
  } else {
    errAndExit("ERROR: no config was specified");
    return 1;
  }

  // /* create config */
  // pidy = findParamIndex(argv, argc, "-c");
  // if (pidy == -1) {
  //   std::cout << "Missing config file" << std::endl;
  //   return 1;
  // }else {
  //   std::string config_file = argv[pidy + 1];
  //   AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, config_file.c_str()));
  // }

  // // Read configuration from file
  // const char cfg_file[] = "config.json";
  // AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, cfg_file));

  // switch on internal error handling (no need to use AMGX_SAFE_CALL after
  // this point)
  AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));

  // Create resources, matrix, vector and solver
  AMGX_resources_create_simple(&rsrc, cfg);
  AMGX_matrix_create(&A, rsrc, mode);
  AMGX_vector_create(&x, rsrc, mode);
  AMGX_vector_create(&b, rsrc, mode);
  AMGX_solver_create(&solver, rsrc, mode, cfg);

  // Set A, b and x
  assert(mat->nbrows == mat->nbcols);                    // Required by AMGX
  assert(mat->Avals.extent(1) == mat->Avals.extent(2));  // Required by AMGX
  int n = mat->nbrows;
  int nnz = mat->nnz;
  int block_dimx = mat->Avals.extent(1);
  int block_dimy = mat->Avals.extent(2);
  // Warning: reinterpret_cast is never a good idea here, use with caution!!
  AMGX_matrix_upload_all(
      A, n, nnz, block_dimx, block_dimy, reinterpret_cast<int*>(mat->rowp),
      reinterpret_cast<int*>(mat->cols), mat->Avals.data, mat->diag);

  AMGX_vector_upload(b, n, block_dimx, rhs->data);
  AMGX_vector_set_zero(x, n, block_dimx);

  // Set up solver and solve
  AMGX_solver_setup(solver, A);
  AMGX_solver_solve(solver, b, x);

  // Get solve status
  AMGX_solver_get_status(solver, &status);

  // Destroy resources, matrix, vectors and solver
  AMGX_solver_destroy(solver);
  AMGX_vector_destroy(x);
  AMGX_vector_destroy(b);
  AMGX_matrix_destroy(A);
  AMGX_resources_destroy(rsrc);
  AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
  AMGX_SAFE_CALL(AMGX_finalize());

  return status;
}

int main(int argc, char **argv) {
  const index_t nx = 8;
  const index_t ny = 8;
  const index_t nz = 8;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);

  std::shared_ptr<FEModel<I, T, ElasticityPDE>> model =
      std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  populate_bcs_3d(nx, ny, nz, bcs);

  // Set the connectivity
  auto conn = element->get_conn();

  // Set the node locations
  auto X = model->get_nodes();

  populate_X_conn_3d(nx, ny, nz, X, conn);

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();

  // Set the element
  T q = 5.0, E = 70e3, nu = 0.3;
  T density = 1.0, design_stress = 1e3;
  auto constitutive = std::make_shared<TopoIsoConstitutive<I, T, Basis>>(
      element, q, E, nu, density, design_stress);
  model->add_constitutive(constitutive);

  // Create the design vector
  A2D::CLayout<1> design_layout(model->nnodes);
  auto x = std::make_shared<A2D::MultiArray<T, A2D::CLayout<1>>>(design_layout);

  // Set the design variable values
  set_dv_3d(nx, ny, nz, *x);

  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, ElasticityPDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // ================================================================
  // Set up input parameters for the solver
  // Compute the Jacobian matrix
  auto J = model->new_matrix();
  model->jacobian(J);
  // Set the residuals and apply the boundary conditions
  auto residual = model->new_solution();
  set_force_3d(nx, ny, nz, model, residual);
  // Initialize the solution vector
  auto solution = model->new_solution();
  solution->zero();
  // ================================================================

  // Conjugate gradient solution
  std::cout << "----------------------------------------------------"
            << std::endl;
  std::cout << "Conjugate gradient methond start" << std::endl;
  TICK("Conjugate");
  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.0;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);
  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 80;
  amg->cg(*residual, *solution, monitor, max_iters);
  TOCK("Conjugate");
  std::cout << "----------------------------------------------------"
            << std::endl;

  std::cout << "----------------------------------------------------"
            << std::endl;
  std::cout << "AMGX solver start" << std::endl;
  // Compute the solution using AMGX
  TICK("AMGX:");
  amgx_solve(J, residual, argc, argv);
  TOCK("AMGX:");
  std::cout << "----------------------------------------------------"
            << std::endl;

  // ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(conn,
  //                                                                        X);
  // vtk.write_mesh();
  // vtk.write_sol("x", *x, 0);
  // vtk.write_sol("ux", *solution, 0);
  // vtk.write_sol("uy", *solution, 1);
  // vtk.write_sol("uz", *solution, 2);

  return 0;
}
