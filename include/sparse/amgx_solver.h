#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amgx_c.h"
#include "cuda_runtime.h"

/*
This is a solver for a sparse matrix A and a vector b using the AMGX library.
The matrix A is assumed to be in BSR format.
usage: ./out [-mode dDDI | dDFI | dDDII | dDFII ...] [-gpu]
-mode: the mode of the solver. default is dDDI.
-gpu: load the matrix from the device memory, default is from the host memory.
*/

using namespace A2D;

#define MAX_MSG_LEN 4096

/* CUDA error macro */
#define CUDA_SAFE_CALL(call) \
  { cuda_safe_call(call); }

void cuda_safe_call(cudaError_t err) {
  if (cudaSuccess != err) {
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__,
            __LINE__, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  while (0)
    ;
}

/* print error message and exit */
void errAndExit(const char* err) {
  printf("%s\n", err);
  fflush(stdout);
  exit(1);
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
    exit(1);
  }
  return -1;
}

/* input contians the argv and argc from main */
template <class Matrix, class Vector>
int amgxSolver(Matrix& mat, Vector& rhs, Vector& sol, int argc, char** argv) {
  /* parameter parsing */
  int pidx = 0;
  /* input matrix and rhs/solution vectors */
  int n, nnz;
  int block_dimx, block_dimy, block_size;
  int* row_ptrs = NULL;
  int* col_indices = NULL;
  int* h_row_ptrs = NULL;
  int* h_col_indices = NULL;
  int* d_row_ptrs = NULL;
  int* d_col_indices = NULL;
  void *values = NULL, *diag = NULL, *dh_x = NULL, *dh_b = NULL;
  void *h_values = NULL, *h_diag = NULL, *h_x = NULL, *h_b = NULL;
  void *d_values = NULL, *d_diag = NULL, *d_x = NULL, *d_b = NULL;
  int sizeof_m_val, sizeof_v_val;

  /* Set up A, b and x */
  assert(mat->nbrows == mat->nbcols);
  assert(mat->Avals.extent(1) == mat->Avals.extent(2));
  n = mat->nbrows;
  nnz = mat->nnz;
  block_dimx = mat->Avals.extent(1);
  block_dimy = mat->Avals.extent(2);
  block_size = block_dimx * block_dimy;
  h_row_ptrs = (int*)(mat->rowp);
  h_col_indices = (int*)(mat->cols);
  h_diag = mat->diag;
  h_values = mat->Avals.data;
  h_b = rhs->data;
  h_x = sol->data;

  /* library handles */
  AMGX_Mode mode;
  AMGX_config_handle cfg;
  AMGX_resources_handle rsrc;
  AMGX_matrix_handle A;
  AMGX_vector_handle b, x;
  AMGX_solver_handle solver;
  /* status handling */
  AMGX_SOLVE_STATUS status;

  /* init */
  AMGX_SAFE_CALL(AMGX_initialize());
  /* system */
  // AMGX_SAFE_CALL(AMGX_install_signal_handler());

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
      mode = AMGX_mode_hZZI;
    } else {
      errAndExit("ERROR: invalid mode");
    }
  } else {
    printf("No mode specified, using dDDI by default.\n");
    mode = AMGX_mode_dDDI;
  }

  sizeof_m_val =
      ((AMGX_GET_MODE_VAL(AMGX_MatPrecision, mode) == AMGX_matDouble))
          ? sizeof(double)
          : sizeof(float);
  sizeof_v_val =
      ((AMGX_GET_MODE_VAL(AMGX_VecPrecision, mode) == AMGX_vecDouble))
          ? sizeof(double)
          : sizeof(float);

  /* create config */
  AMGX_SAFE_CALL(
      AMGX_config_create(&cfg,
                         "config_version=2,"

                         /* cg accelerator config */
                         "solver(main)=CG,"
                         //  "preconditioner(amg)=NOSOLVER,"
                         "preconditioner(amg)=AMG,"

                         /* outer solver setup */
                         "main:max_iters=10000,"
                         "main:tolerance=1e-8,"
                         "main:norm=L2,"
                         "main:use_scalar_norm=1,"
                         "main:convergence=RELATIVE_INI_CORE, "

                         /* printing obtions */
                         "main:monitor_residual=1,"
                         "main:obtain_timings=1,"
                         "main:store_res_history=1,"
                         /* print out the residuals at each iteration */
                         //  "main:print_solve_stats=1,"
                         //  "main:print_grid_stats=1,"

                         /* preconditioner setup */
                         "amg:algorithm=CLASSICAL,"
                         "amg:max_iters=2,"
                         "amg:presweeps=1,"
                         "amg:postsweeps=1,"
                         "amg:cycle=V"));

  /* or read configuration from file */
  // const char cfg_file[] = "config.json";
  // AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, cfg_file));

  /* create resources, matrix, vector and solver */
  AMGX_resources_create_simple(&rsrc, cfg);
  AMGX_matrix_create(&A, rsrc, mode);
  AMGX_vector_create(&x, rsrc, mode);
  AMGX_vector_create(&b, rsrc, mode);
  AMGX_solver_create(&solver, rsrc, mode, cfg);

  /* set memeory pointers */
  if ((pidx = findParamIndex(argv, argc, "-gpu")) != -1) {
    /* allocate memory and copy the data to the GPU */
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_x, n * block_dimx * sizeof_v_val));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_b, n * block_dimy * sizeof_v_val));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_col_indices, nnz * sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_row_ptrs, (n + 1) * sizeof(int)));
    CUDA_SAFE_CALL(
        cudaMalloc((void**)&d_values, nnz * block_size * sizeof_m_val));
    CUDA_SAFE_CALL(
        cudaMemcpy(d_x, h_x, n * block_dimx * sizeof_v_val, cudaMemcpyDefault));
    CUDA_SAFE_CALL(
        cudaMemcpy(d_b, h_b, n * block_dimy * sizeof_v_val, cudaMemcpyDefault));
    CUDA_SAFE_CALL(cudaMemcpy(d_col_indices, h_col_indices, nnz * sizeof(int),
                              cudaMemcpyDefault));
    CUDA_SAFE_CALL(cudaMemcpy(d_row_ptrs, h_row_ptrs, (n + 1) * sizeof(int),
                              cudaMemcpyDefault));
    CUDA_SAFE_CALL(cudaMemcpy(d_values, h_values,
                              nnz * block_size * sizeof_m_val,
                              cudaMemcpyDefault));

    if (h_diag != NULL) {
      CUDA_SAFE_CALL(cudaMalloc(&d_diag, n * block_size * sizeof_m_val));
      CUDA_SAFE_CALL(cudaMemcpy(d_diag, h_diag, n * block_size * sizeof_m_val,
                                cudaMemcpyDefault));
    }

    /* set pointers to point to GPU (device) memory */
    row_ptrs = d_row_ptrs;
    col_indices = d_col_indices;
    values = d_values;
    diag = d_diag;
    dh_x = d_x;
    dh_b = d_b;

  } else {
    /* pin the memory to improve performance */
    AMGX_SAFE_CALL(AMGX_pin_memory(h_x, n * block_dimx * sizeof_v_val));
    AMGX_SAFE_CALL(AMGX_pin_memory(h_b, n * block_dimx * sizeof_v_val));
    AMGX_SAFE_CALL(AMGX_pin_memory(h_col_indices, nnz * sizeof(int)));
    AMGX_SAFE_CALL(AMGX_pin_memory(h_row_ptrs, (n + 1) * sizeof(int)));
    AMGX_SAFE_CALL(AMGX_pin_memory(h_values, nnz * block_size * sizeof_m_val));

    if (h_diag != NULL) {
      AMGX_SAFE_CALL(AMGX_pin_memory(h_diag, n * block_size * sizeof_m_val));
    }

    /* set pointers to point to CPU (host) memory */
    row_ptrs = h_row_ptrs;
    col_indices = h_col_indices;
    values = h_values;
    diag = h_diag;
    dh_x = h_x;
    dh_b = h_b;
  }

  /* upload the matrix to the solver from the host to the device */
  AMGX_matrix_upload_all(A, n, nnz, block_dimx, block_dimy, row_ptrs,
                         col_indices, values, diag);
  AMGX_vector_bind(x, A);
  AMGX_vector_bind(b, A);
  /* upload the vector (and the connectivity information) */
  AMGX_vector_upload(x, n, block_dimx, dh_x);
  AMGX_vector_upload(b, n, block_dimx, dh_b);
  /* solver setup */
  AMGX_solver_setup(solver, A);
  std::cout << "------------------------------------------------------------"
            << std::endl;
  /* solver solve */
  AMGX_solver_solve(solver, b, x);
  /* download the solution from the device to the host */
  AMGX_vector_download(x, sol->data);
  /* Get solve status */
  AMGX_solver_get_status(solver, &status);

  /* print the residual history */
  // int nit;
  // double res;
  // AMGX_solver_get_iterations_number(solver, &nit);
  // // printf("Iterations: %d, Residual: %f\n", nit, res);
  // printf("Residual history: ");
  // for (int i = 0; i < nit+1; i++) {
  //   printf("residual from iteration %d=", i);
  //   AMGX_solver_get_iteration_residual(solver, i, 0, &res);
  //   printf("%f \n", (float)(res));
  // }

  /* print the solution */
  double eps = 1e-8;
  int nit = 0;
  double res_init = 0.0, res_final = 0.0;
  AMGX_solver_get_iterations_number(solver, &nit);
  AMGX_solver_get_iteration_residual(solver, 0, 0, &res_init);
  AMGX_solver_get_iteration_residual(solver, nit, 0, &res_final);
  std::cout << "------------------------------------------------------------"
            << std::endl;
  std::cout << ((status == 0) ? "AMGX: converged!" : "AMGX: not converged!")
            << std::endl;
  std::cout << "    Total Iterations: " << nit << "/10000" << std::endl;
  std::cout << "    Initial Residual: \t\t\t" << std::setprecision(6)
            << std::scientific << std::setw(15) << res_init << std::fixed
            << std::endl;
  std::cout << "    Final Residual: \t\t\t" << std::setprecision(6)
            << std::scientific << std::setw(15) << res_final << std::fixed
            << std::endl;
  std::cout << "    Total Reduction in Residual: \t" << std::setprecision(6)
            << std::scientific << std::setw(15)
            << ((res_init > eps) ? res_final / res_init : res_init)
            << std::fixed << std::endl;
  if (status != 0) {
    std::cout << "    Solver failed with status " << status << "\n"
              << "------------------------------------------------------------"
              << std::endl;
    return 1;
  }
  std::cout << "------------------------------------------------------------"
            << std::endl;

  if ((pidx = findParamIndex(argv, argc, "-gpu")) != -1) {
    /* deallocate GPU (device) memory */
    CUDA_SAFE_CALL(cudaFree(d_x));
    CUDA_SAFE_CALL(cudaFree(d_b));
    CUDA_SAFE_CALL(cudaFree(d_row_ptrs));
    CUDA_SAFE_CALL(cudaFree(d_col_indices));
    CUDA_SAFE_CALL(cudaFree(d_values));

    if (d_diag != NULL) {
      CUDA_SAFE_CALL(cudaFree(d_diag));
    }
  } else {
    /* unpin the memory */
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_x));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_b));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_values));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_row_ptrs));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_col_indices));

    if (h_diag != NULL) {
      AMGX_SAFE_CALL(AMGX_unpin_memory(h_diag));
    }
  }

  /* reconstruct the global matrix and write it to a file */
  // AMGX_write_system_distributed(A, b, x, "output_system.mtx", nrings, nranks,
  // partition_sizes, partition_vector_size, partition_vector);

  /* Destroy resources, matrix, vectors and solver */
  AMGX_solver_destroy(solver);
  AMGX_vector_destroy(x);
  AMGX_vector_destroy(b);
  AMGX_matrix_destroy(A);
  AMGX_resources_destroy(rsrc);
  AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
  AMGX_SAFE_CALL(AMGX_finalize());
  CUDA_SAFE_CALL(cudaDeviceReset());
  return status;
}
