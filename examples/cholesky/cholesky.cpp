#include <iostream>

#include "mpi.h"
#include "sparse/sparse_cholesky.h"

using namespace A2D;

template <typename T>
void build_matrix(int nx, int *_size, int **_colp, int **_rows, T **_kvals,
                  const int *iperm = nullptr) {
  T kmat[][4] = {
      {4.0, 2.0, 2.0, 1.0},
      {2.0, 4.0, 1.0, 2.0},
      {2.0, 1.0, 4.0, 2.0},
      {1.0, 2.0, 2.0, 4.0},
  };

  T ke[64];
  for (int k = 0; k < 64; k++) {
    ke[k] = 0.0;
  }
  for (int ki = 0; ki < 2; ki++) {
    for (int ii = 0; ii < 4; ii++) {
      for (int jj = 0; jj < 4; jj++) {
        ke[8 * (2 * ii + ki) + 2 * jj + ki] = kmat[ii][jj] / 9.0;
      }
    }
  }

  int size = 2 * (nx + 1) * (nx + 1);
  int *colp = new int[size + 1];
  for (int i = 0; i < size + 1; i++) {
    colp[i] = 0;
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nx; j++) {
      int nodes[] = {i + j * (nx + 1), i + 1 + j * (nx + 1),
                     i + (j + 1) * (nx + 1), i + 1 + (j + 1) * (nx + 1)};

      for (int k = 0; k < 2; k++) {
        for (int ii = 0; ii < 4; ii++) {
          int ivar = 2 * nodes[ii] + k;
          if (iperm) {
            ivar = iperm[ivar];
          }
          colp[ivar] += 8;
        }
      }
    }
  }

  int nnz = 0;
  for (int i = 0; i < size; i++) {
    int tmp = colp[i];
    colp[i] = nnz;
    nnz += tmp;
  }
  colp[size] = nnz;

  int *rows = new int[nnz];
  T *kvals = new T[nnz];

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nx; j++) {
      int nodes[] = {i + j * (nx + 1), i + 1 + j * (nx + 1),
                     i + (j + 1) * (nx + 1), i + 1 + (j + 1) * (nx + 1)};

      for (int ki = 0; ki < 2; ki++) {
        for (int kj = 0; kj < 2; kj++) {
          for (int ii = 0; ii < 4; ii++) {
            for (int jj = 0; jj < 4; jj++) {
              int ivar = 2 * nodes[ii] + ki;
              int jvar = 2 * nodes[jj] + kj;
              if (iperm) {
                ivar = iperm[ivar];
                jvar = iperm[jvar];
              }
              rows[colp[ivar]] = jvar;
              kvals[colp[ivar]] = ke[8 * (2 * ii + ki) + (2 * jj + kj)];
              colp[ivar]++;
            }
          }
        }
      }
    }
  }

  for (int i = size - 1; i >= 0; i--) {
    colp[i + 1] = colp[i];
  }
  colp[0] = 0;

  *_size = size;
  *_colp = colp;
  *_rows = rows;
  *_kvals = kvals;
}

int main(int argc, char *argv[]) {
  using T = double;

  int nx = 512;

  int size;
  int *colp;
  int *rows;
  T *kvals;
  build_matrix(nx, &size, &colp, &rows, &kvals, nullptr);

  T *b = new T[size];
  for (int i = 0; i < size; i++) {
    b[i] = 0.0;
  }
  for (int i = 0; i < size; i++) {
    for (int jp = colp[i]; jp < colp[i + 1]; jp++) {
      b[rows[jp]] += kvals[jp];
    }
  }

  std::printf("size = %d\n", size);
  double t0 = MPI_Wtime();
  CholOrderingType order = CholOrderingType::ND;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "ND") == 0) {
      order = CholOrderingType::ND;
    } else if (strcmp(argv[k], "NATURAL") == 0) {
      order = CholOrderingType::NATURAL;
    }
  }
  SparseCholesky<T> *chol = new SparseCholesky<T>(size, colp, rows, order);
  double t1 = MPI_Wtime();
  chol->setValues(size, colp, rows, kvals);

  delete[] colp;
  delete[] rows;
  delete[] kvals;

  double t2 = MPI_Wtime();
  chol->factor();
  double t3 = MPI_Wtime();
  chol->solve(b);
  double t4 = MPI_Wtime();

  std::printf("Setup/order time: %12.5e\n", t1 - t0);
  std::printf("Set values  time: %12.5e\n", t2 - t1);
  std::printf("Factor time:      %12.5e\n", t3 - t2);
  std::printf("Solve time:       %12.5e\n", t4 - t3);

  T err = 0.0;
  for (int i = 0; i < size; i++) {
    err += (1.0 - b[i]) * (1.0 - b[i]);
  }
  std::printf("||x - e||: %25.15e\n", sqrt(err));

  delete chol;

  return 0;
}