#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>

#include "ticktock.h"
#if (defined(_MSC_VER) && (_MSC_VER < 1600))

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

#else
#include <stdint.h>
#endif
#include <string.h>

#include "a2d.h"
#include "amgx_c.h"
#include "cuda_runtime.h"
#include "toolkit.h"

using namespace A2D;

#define MAX_MSG_LEN 4096

/* print error message and exit */
void errAndExit(const char* err) {
  printf("%s\n", err);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}

/* print callback (could be customized) */
void print_callback(const char* msg, int length) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    printf("%s", msg);
  }
}

/* print usage and exit */
void printUsageAndExit() {
  char msg[MAX_MSG_LEN] =
      "Usage: mpirun [-n nranks] ./amgx_mpi_capi_agg [-mode [dDDI | dDFI | "
      "dFFI]] [-m file] [-c config_file] [-amg \"variable1=value1 ... "
      "variable3=value3\"] [-gpu] [-it k]\n";
  strcat(msg, "     -mode:   select the solver mode\n");
  strcat(msg,
         "     -c:      set the amg solver options from the config file\n");
  strcat(msg,
         "     -amg:    set the amg solver options from the command line\n");
  strcat(msg, "     -gpu:    load the matrix from the device memory\n");
  strcat(msg,
         "     -it k:   set the number k of outer (non-linear) iteration\n");
  print_callback(msg, MAX_MSG_LEN);
  exit(0);
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

void block_element_random_update(int block_dimx, int block_dimy, int i, void* v,
                                 AMGX_Mode mode) {
  int j1, j2, block_size;
  block_size = block_dimy * block_dimx;

  for (j1 = 0; j1 < block_dimy; j1++) {
    for (j2 = 0; j2 < block_dimx; j2++) {
      if ((AMGX_GET_MODE_VAL(AMGX_MatPrecision, mode) == AMGX_matFloat)) {
        float* t = (float*)v;
        t[i * block_size + j1 * block_dimx + j2] *=
            (1.0f + (rand() % 10) * (1e-6f));
      } else {
        double* t = (double*)v;
        t[i * block_size + j1 * block_dimx + j2] *=
            (1.0 + (rand() % 10) * (1e-12));
      }
    }
  }
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
int amgx_solve(Matrix& mat, Vector& rhs, Vector& sol, int argc, char** argv) {
  // parameter parsing
  int pidx = 0;
  int pidy = 0;
  // MPI (with CUDA GPUs)
  int rank = 0;
  int lrank = 0;
  int nranks = 0;
  int gpu_count = 0;
  MPI_Comm amgx_mpi_comm = MPI_COMM_WORLD;
  // number of outer (non-linear) iterations
  int i = 0;
  int k = 0;
  int max_it = 0;
  // versions
  int major, minor;
  char *ver, *date, *time;
  // input matrix and rhs/solution
  int n, nnz, block_dimx, block_dimy, block_size, num_neighbors, ncol;
  int*row_ptrs = NULL, *neighbors = NULL;
  int* col_indices = NULL;
  int *h_row_ptrs = NULL, *h_col_indices = NULL;
  int *d_row_ptrs = NULL, *d_col_indices = NULL;
  void *values = NULL, *diag = NULL, *dh_x = NULL, *dh_b = NULL;
  void *h_values = NULL, *h_diag = NULL, *h_x = NULL, *h_b = NULL;
  void *d_values = NULL, *d_diag = NULL, *d_x = NULL, *d_b = NULL;
  int sizeof_m_val;
  int sizeof_v_val;
  int* partition_sizes = NULL;
  int* partition_vector = NULL;
  int partition_vector_size = 0;

  /* Set up A, b and x */
  assert(mat->nbrows == mat->nbcols);                    // Required by AMGX
  assert(mat->Avals.extent(1) == mat->Avals.extent(2));  // Required by AMGX
  n = mat->nbrows;
  ncol = mat->nbcols;
  nnz = mat->nnz;
  block_dimx = mat->Avals.extent(1);
  block_dimy = mat->Avals.extent(2);
  block_size = block_dimx * block_dimy;
  h_row_ptrs = reinterpret_cast<int*>(mat->rowp);
  h_col_indices = reinterpret_cast<int*>(mat->cols);
  // h_row_ptrs = mat->rowp;
  // h_col_indices = mat->cols;
  // print out the row pointers
  // for (int i = 0; i < n + 1; i++) {
  //   printf("%d\n", ((double*)h_row_ptrs)[i]);
  // }
  
  h_diag = mat->diag;
  h_values = mat->Avals.data;
  h_b = rhs->data;
  h_x = sol->data;
  __CHECK__

  for (int i = 0; i < block_size; i++) {
    if (((double*)h_x)[i] < 1e-6) {
      ((double*)h_x)[i] = 0.0;
    }
    std::cout << " " << ((double*)h_x)[i];
  }
  std::cout << std::endl;

  // AMGX_SAFE_CALL(AMGX_pin_memory(h_x, n * block_dimx * sizeof_v_val));
  // AMGX_SAFE_CALL(AMGX_pin_memory(h_b, n * block_dimx * sizeof_v_val));
  // AMGX_SAFE_CALL(AMGX_pin_memory(h_col_indices, nnz * sizeof(int64_t)));
  // AMGX_SAFE_CALL(AMGX_pin_memory(h_row_ptrs, (n + 1) * sizeof(int)));
  // AMGX_SAFE_CALL(AMGX_pin_memory(h_values, nnz * block_size * sizeof_m_val));

  // AMGX_matrix_upload_all(
  //     A, n, nnz, block_dimx, block_dimy, reinterpret_cast<int*>(mat->rowp),
  //     reinterpret_cast<int*>(mat->cols), mat->Avals.data, mat->diag);
  // AMGX_matrix_upload_all_global(A, nglobal, n, nnz, block_dimx, block_dimy,
  //                               row_ptrs, col_indices, values, diag, nrings,
  //                               nrings, partition_vector);

  // library handles
  AMGX_Mode mode;
  AMGX_config_handle cfg;
  AMGX_resources_handle rsrc;
  AMGX_matrix_handle A;
  AMGX_vector_handle b, x;
  AMGX_solver_handle solver;
  // status handling
  AMGX_SOLVE_STATUS status;
  /* MPI init (with CUDA GPUs) */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(amgx_mpi_comm, &nranks);
  MPI_Comm_rank(amgx_mpi_comm, &rank);
  // CUDA GPUs
  CUDA_SAFE_CALL(cudaGetDeviceCount(&gpu_count));
  lrank = rank % gpu_count;
  CUDA_SAFE_CALL(cudaSetDevice(lrank));
  printf("Process %d selecting device %d\n", rank, lrank);

  /* check arguments */
  if (argc == 1) {
    printUsageAndExit();
  }

  /* init */
  AMGX_SAFE_CALL(AMGX_initialize());
  /* system */
  AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
  AMGX_SAFE_CALL(AMGX_install_signal_handler());

  /* get api and build info */
  if ((pidx = findParamIndex(argv, argc, "--version")) != -1) {
    AMGX_get_api_version(&major, &minor);
    printf("amgx api version: %d.%d\n", major, minor);
    AMGX_get_build_info_strings(&ver, &date, &time);
    printf("amgx build version: %s\nBuild date and time: %s %s\n", ver, date,
           time);
    AMGX_SAFE_CALL(AMGX_finalize());
    MPI_Finalize();
    exit(0);
  }

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

  sizeof_m_val =
      ((AMGX_GET_MODE_VAL(AMGX_MatPrecision, mode) == AMGX_matDouble))
          ? sizeof(double)
          : sizeof(float);
  sizeof_v_val =
      ((AMGX_GET_MODE_VAL(AMGX_VecPrecision, mode) == AMGX_vecDouble))
          ? sizeof(double)
          : sizeof(float);
  /* get max_it number of outer (non-linear) iteration */
  max_it = 1;

  if ((pidx = findParamIndex(argv, argc, "-it")) != -1) {
    if ((pidy = findParamIndex(argv, argc, "-gpu")) != -1) {
      errAndExit(
          "ERROR: -gpu and -it options are not compatible, you must choose one "
          "or the other option");
    }

    max_it = (int)atol(argv[pidx + 1]);
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
  }

  /* switch on internal error handling (no need to use AMGX_SAFE_CALL after this
   * point) */
  AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));
  /* create resources, matrix, vector and solver */
  AMGX_resources_create(&rsrc, cfg, &amgx_mpi_comm, 1, &lrank);
  AMGX_matrix_create(&A, rsrc, mode);
  AMGX_vector_create(&x, rsrc, mode);
  AMGX_vector_create(&b, rsrc, mode);
  AMGX_solver_create(&solver, rsrc, mode, cfg);

  // read partitioning vector
  if ((pidx = findParamIndex(argv, argc, "-partvec")) != -1) {
    // open the file
    FILE* fin_rowpart = fopen(argv[pidx + 1], "rb");

    if (fin_rowpart == NULL) {
      errAndExit("ERROR: opening the file for the partition vector");
    }

    // find the size of the partition vector
    if (fseek(fin_rowpart, 0L, SEEK_END) != 0) {
      errAndExit("ERROR: reading partition vector");
    }

    partition_vector_size = ftell(fin_rowpart) / sizeof(int);

    if (partition_vector_size == -1L) {
      errAndExit("ERROR: reading partition vector");
    }

    partition_vector = (int*)malloc(partition_vector_size * sizeof(int));
    // reading the partition vector:
    rewind(fin_rowpart);
    int result = fread((void*)partition_vector, sizeof(int),
                       partition_vector_size, fin_rowpart);

    if (result != partition_vector_size) {
      errAndExit("ERROR: reading partition vector");
    }

    printf("Read partition vector, consisting of %d rows\n",
           partition_vector_size);
    fclose(fin_rowpart);
  }

  // read the matrix, [and rhs & solution]
  // WARNING: use 2 rings for classical path
  int nrings;  //=2;
  AMGX_config_get_default_number_of_rings(cfg, &nrings);
  // printf("nrings=%d\n",nrings);

  // Warning: reinterpret_cast is never a good idea here, use with caution!!

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
    CUDA_SAFE_CALL(cudaMemcpy(d_col_indices, h_col_indices,
                              nnz * sizeof(int), cudaMemcpyDefault));
    CUDA_SAFE_CALL(cudaMemcpy(d_row_ptrs, h_row_ptrs,
                              (n + 1) * sizeof(int), cudaMemcpyDefault));
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
    /* pin the memory to improve performance
       WARNING: Even though, internal error handling has been requested,
                AMGX_SAFE_CALL needs to be used on this system call.
                It is an exception to the general rule. */
    // print all input data

    printf("\n");
    printf("Matrix:\n");
    printf("  n=%d\n", n);
    printf("  ncol=%d\n", ncol);
    printf("  nnz=%d\n", nnz);
    printf("  block_size=%d\n", block_size);
    printf("  block_dimx=%d\n", block_dimx);
    printf("  block_dimy=%d\n", block_dimy);
    printf("  sizeof_v_val = %d\n", sizeof_v_val);
    printf("  sizeof_m_val = %d\n", sizeof_m_val);
    printf("  h_x = %p\n", h_x);
    printf("  h_b = %p\n", h_b);
    printf("  h_col_indices = %p\n", h_col_indices);
    printf("  h_row_ptrs = %p\n", h_row_ptrs);
    printf("  h_values = %p\n", h_values);
    printf("  h_diag = %p\n", h_diag);
    printf("\n");

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

  __CHECK__
  /* compute global number of rows */
  int nglobal;
  /* upload the matrix with global indices and compute necessary connectivity
   * information */
  if (partition_vector != NULL) {
    // If no partition vector is given, we assume a partitioning with contiguous
    // blocks (see example above). It is sufficient (and faster/more scalable)
    // to calculate the partition offsets and pass those into the API call
    // instead of creating a full partition vector.
    int64_t* partition_offsets =
        (int64_t*)malloc((nranks + 1) * sizeof(int64_t));
    // gather the number of rows on each rank, and perform an exclusive scan to
    // get the offsets.
    int64_t n64 = n;
    partition_offsets[0] = 0;  // rows of rank 0 always start at index 0
    MPI_Allgather(&n64, 1, MPI_INT64_T, &partition_offsets[1], 1, MPI_INT64_T,
                  amgx_mpi_comm);
    for (int i = 2; i < nranks + 1; ++i) {
      partition_offsets[i] += partition_offsets[i - 1];
    }
    nglobal = partition_offsets[nranks];  // last element always has global
                                          // number of rows

    AMGX_distribution_handle dist;
    AMGX_distribution_create(&dist, cfg);
    if (rank == 0) {
      printf("\n");
      printf("Partitioning:\n");
      printf("  nranks = %d\n", nranks);
      printf("  nglobal = %d\n", nglobal);
      printf("  partition_offsets = [");
      for (int i = 0; i < nranks + 1; ++i) {
        printf("%d ", partition_offsets[i]);
      }
      printf("]\n");
      printf("\n");
      printf("dist:\n");
      for (int i = 0; i < nranks; ++i) {
        printf("%d\n", dist[i]);
      }
      printf("\n");
      __CHECK__
    }
    
    AMGX_distribution_set_partition_data(dist, AMGX_DIST_PARTITION_OFFSETS,
                                         partition_offsets);
    __CHECK__
    AMGX_matrix_upload_distributed(A, nglobal, n, nnz, block_dimx, block_dimy,
                                   row_ptrs, col_indices, values, diag, dist);
    __CHECK__
    AMGX_distribution_destroy(dist);
    free(partition_offsets);
  } else {
    __CHECK__
    // MPI_Allreduce(&n, &nglobal, 1, MPI_INT, MPI_SUM, amgx_mpi_comm);
    nglobal = n * nranks;
    __CHECK__
    // AMGX_matrix_upload_all_global(A, nglobal, n, nnz, block_dimx, block_dimy,
    //                               row_ptrs, col_indices, values, diag, nrings,
    //                               nrings, NULL);
    AMGX_matrix_upload_all(A, n, nnz, block_dimx, block_dimy, row_ptrs,
                           col_indices, values, diag);
    __CHECK__
  }
  __CHECK__

  /* free temporary storage */
  if (partition_vector != NULL) {
    free(partition_vector);
  }

  // AMGX_matrix_upload_all(
  //     A, n, nnz, block_dimx, block_dimy, reinterpret_cast<int*>(mat->rowp),
  //     reinterpret_cast<int*>(mat->cols), mat->Avals.data, mat->diag);
  /* set the connectivity information (for the vector) */
  AMGX_vector_bind(x, A);
  AMGX_vector_bind(b, A);
  /* upload the vector (and the connectivity information) */
  AMGX_vector_upload(x, n, block_dimx, dh_x);
  AMGX_vector_upload(b, n, block_dimx, dh_b);
  // // Set up solver and solve
  // AMGX_solver_setup(solver, A);
  // AMGX_solver_solve(solver, b, x);

  __CHECK__
  /* start outer (non-linear) iterations */
  for (k = 0; k < max_it; k++) {
    /* solver setup */
    // MPI barrier for stability (should be removed in practice to maximize
    // performance)
    // MPI_Barrier(amgx_mpi_comm);
    AMGX_solver_setup(solver, A);
    /* solver solve */
    // MPI barrier for stability (should be removed in practice to maximize
    // performance)
    // MPI_Barrier(amgx_mpi_comm);
    AMGX_solver_solve(solver, b, x);
    /* check the status */
    // MPI_Barrier(amgx_mpi_comm);
    AMGX_solver_get_status(solver, &status);

    /* while not the last iteration */
    if (k + 1 < max_it) {
      /* example of how to change parameters between non-linear iterations */
      // AMGX_config_add_parameters(&cfg, "config_version=2,
      // default:tolerance=1e-12"); AMGX_solver_solve(solver, b, x);

      /* example of how to replace coefficients between non-linear iterations */
      for (i = 0; i < nnz; i++) {
        block_element_random_update(block_dimx, block_dimy, i, values, mode);
      }

      if (diag != NULL) {
        for (i = 0; i < n; i++) {
          block_element_random_update(block_dimx, block_dimy, i, diag, mode);
        }
      }

      // MPI_Barrier(amgx_mpi_comm);
      AMGX_matrix_replace_coefficients(A, n, nnz, values, diag);
      /* upload original vectors (and the connectivity information) */
      AMGX_vector_upload(x, n, block_dimx, dh_x);
      AMGX_vector_upload(b, n, block_dimx, dh_b);
    }
  }

  /* example of how to get (the local part of) the solution */
  // if ((pidx = findParamIndex(argv, argc, "-gpu")) != -1) {
  //     CUDA_SAFE_CALL(cudaMalloc(&d_result, n*block_dimx*sizeof_v_val));
  //     AMGX_vector_download(x, d_result);
  //     CUDA_SAFE_CALL(cudaFree(d_result));
  // }
  // else{
  // void* h_result = malloc(n*block_dimx*sizeof_v_val);
  // AMGX_vector_download(x, h_result);
  // // print h_result
  // printf("\n");
  // printf("h_result:\n");
  // for (int i = 0; i < n; ++i) {
  //   printf("%d ", ((int*)h_result)[i]);
  // }
  // free(h_result);
  AMGX_vector_download(x, h_x);
  // print result h_x
  for (int i = 0; i < block_size; i++) {
    if (((double*)h_x)[i] < 1e-6) {
      ((double*)h_x)[i] = 0.0;
    }
    std::cout << " " << ((double*)h_x)[i];
  }
  std::cout << std::endl;
  // }

  /* example of how to reconstruct the global matrix and write it to a file */
  // AMGX_write_system_distributed(A, b, x, "output_system.mtx", nrings, nranks,
  // partition_sizes, partition_vector_size, partition_vector);

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
    /* unpin the memory
       WARNING: Even though, internal error handling has been requested,
                AMGX_SAFE_CALL needs to be used on this system call.
                It is an exception to the general rule. */
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_x));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_b));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_values));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_row_ptrs));
    AMGX_SAFE_CALL(AMGX_unpin_memory(h_col_indices));

    if (h_diag != NULL) {
      AMGX_SAFE_CALL(AMGX_unpin_memory(h_diag));
    }
  }

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
  MPI_Finalize();
  CUDA_SAFE_CALL(cudaDeviceReset());
  return status;
}

int main(int argc, char** argv) {
  const index_t nx = 32;
  const index_t ny = 32;
  const index_t nz = 32;
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
  // set_dv_3d(nx, ny, nz, *x);
  x->fill(1.0);
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
  // set_force_3d(nx, ny, nz, model, residual);
  // Initialize the solution vector
  auto solution = model->new_solution();
  solution->fill(1.0);
  model->zero_bcs(solution);
  residual->zero();
  BSRMatVecMult(*J, *solution, *residual);
  model->zero_bcs(residual);
  solution->zero();
  // ================================================================

  // // Conjugate gradient solution
  // std::cout << "----------------------------------------------------"
  //           << std::endl;
  // std::cout << "Conjugate gradient methond start" << std::endl;
  // TICK("Conjugate");
  // int num_levels = 3;
  // double omega = 0.6667;
  // double epsilon = 0.0;
  // bool print_info = true;
  // auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);
  // // Compute the solution
  // index_t monitor = 10;
  // index_t max_iters = 80;
  // amg->cg(*residual, *solution, monitor, max_iters);
  // // print solution from solution.begin() to solution.end()
  // // std::cout << "Solution: " << std::endl;
  // // for (index_t i = 0; i < 8; i++) {
  // //   std::cout << " " << solution(i, 0);
  // // }
  // // std::cout << std::endl;

  // TOCK("Conjugate");
  // std::cout << "----------------------------------------------------"
  //           << std::endl;

  std::cout << "----------------------------------------------------"
            << std::endl;
  std::cout << "AMGX solver start" << std::endl;
  // Compute the solution using AMGX
  TICK("AMGX");
  amgx_solve(J, residual, solution, argc, argv);
  TOCK("AMGX");
  std::cout << "----------------------------------------------------"
            << std::endl;

  // ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())>
  // vtk(conn,
  //                                                                        X);
  // vtk.write_mesh();
  // vtk.write_sol("x", *x, 0);
  // vtk.write_sol("ux", *solution, 0);
  // vtk.write_sol("uy", *solution, 1);
  // vtk.write_sol("uz", *solution, 2);
  return 0;
}
