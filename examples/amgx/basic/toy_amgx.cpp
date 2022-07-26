#include "amgx_c.h"

int main(int argc, char* argv[]) {
  // Init
  AMGX_SAFE_CALL(AMGX_initialize());
  // System
  // AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
  // AMGX_SAFE_CALL(AMGX_install_signal_handler());
  // library handles

  AMGX_Mode mode = AMGX_mode_dDDI;
  AMGX_config_handle cfg;
  AMGX_resources_handle rsrc;
  AMGX_matrix_handle A;
  AMGX_vector_handle b, x;
  AMGX_solver_handle solver;

  // status handling
  AMGX_SOLVE_STATUS status;

  // Read configuration from file
  const char cfg_file[] = "config.json";
  AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, cfg_file));

  // switch on internal error handling (no need to use AMGX_SAFE_CALL after
  // this point)
  AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));

  // Create resources, matrix, vector and solver
  AMGX_resources_create_simple(&rsrc, cfg);
  AMGX_matrix_create(&A, rsrc, mode);
  AMGX_vector_create(&x, rsrc, mode);
  AMGX_vector_create(&b, rsrc, mode);
  AMGX_solver_create(&solver, rsrc, mode, cfg);

  // read A, set x to zero
  const char mat_file[] = "matrix.mtx";
  AMGX_read_system(A, b, x, mat_file);
  int n, block_dimx, block_dimy;
  AMGX_matrix_get_size(A, &n, &block_dimx, &block_dimy);
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