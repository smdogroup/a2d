/*
  This file contains implementations of the inline functions.
*/
#ifndef A2D_SPARSE_UTILS_INL_H
#define A2D_SPARSE_UTILS_INL_H

#include "sparse/sparse_utils.h"

namespace A2D {

// Compute y = alpha A * x + beta * y
template <typename T>
void CSRMatVec(double alpha, int nrows, const int *rowp, const int *cols,
               const T *Avals, const T *x, double beta, T *y) {
  if (alpha == 1.0 && beta == 0.0) {
    for (int i = 0; i < nrows; i++) {
      int jp_end = rowp[i + 1];
      T value = 0.0;
      for (int jp = rowp[i]; jp < jp_end; jp++) {
        int j = cols[jp];
        value += Avals[jp] * x[j];
      }
      y[i] = value;
    }
  } else if (alpha == -1.0 && beta == 0.0) {
    for (int i = 0; i < nrows; i++) {
      int jp_end = rowp[i + 1];
      T value = 0.0;
      for (int jp = rowp[i]; jp < jp_end; jp++) {
        int j = cols[jp];
        value += Avals[jp] * x[j];
      }
      y[i] = -value;
    }
  } else {
    if (beta == 0.0) {
      for (int i = 0; i < nrows; i++) {
        int jp_end = rowp[i + 1];
        T value = 0.0;
        for (int jp = rowp[i]; jp < jp_end; jp++) {
          int j = cols[jp];
          value += Avals[jp] * x[j];
        }
        y[i] = alpha * value;
      }
    } else {
      for (int i = 0; i < nrows; i++) {
        int jp_end = rowp[i + 1];
        T value = 0.0;
        for (int jp = rowp[i]; jp < jp_end; jp++) {
          int j = cols[jp];
          value += Avals[jp] * x[j];
        }
        y[i] = beta * y[i] + alpha * value;
      }
    }
  }
}

// Compute y = alpha A * x + beta * y
template <typename T>
void CSCMatVec(double alpha, int nrows, int ncols, const int *colp,
               const int *rows, const T *Avals, const T *x, double beta, T *y) {
  if (beta == 0.0) {
    for (int i = 0; i < nrows; i++) {
      y[i] = 0.0;
    }
  } else {
    for (int i = 0; i < nrows; i++) {
      y[i] = beta * y[i];
    }
  }

  if (alpha == 1.0) {
    for (int i = 0; i < ncols; i++) {
      int jp_end = colp[i + 1];
      T xi = x[i];
      for (int jp = colp[i]; jp < jp_end; jp++) {
        int j = rows[jp];
        y[j] += Avals[jp] * xi;
      }
    }
  } else if (alpha == -1.0) {
    for (int i = 0; i < ncols; i++) {
      int jp_end = colp[i + 1];
      T xi = x[i];
      for (int jp = colp[i]; jp < jp_end; jp++) {
        int j = rows[jp];
        y[j] -= Avals[jp] * xi;
      }
    }
  } else {
    for (int i = 0; i < ncols; i++) {
      int jp_end = colp[i + 1];
      T xi = alpha * x[i];
      for (int jp = colp[i]; jp < jp_end; jp++) {
        int j = rows[jp];
        y[j] += Avals[jp] * xi;
      }
    }
  }
}

template <typename T>
void SparseTranspose(int nrows, int ncols, const int *rowp, const int *cols,
                     const T *Avals, int *colp, int *rows, T *ATvals) {
  for (int j = 0; j < ncols + 1; j++) {
    colp[j] = 0;
  }

  for (int i = 0; i < nrows; i++) {
    int jp_end = rowp[i + 1];
    for (int jp = rowp[i]; jp < jp_end; jp++) {
      int j = cols[jp];
      colp[j + 1]++;
    }
  }

  // Set the colp array to be a pointer into each row
  for (int j = 0; j < ncols; j++) {
    colp[j + 1] += colp[j];
  }

  // Now, add the rows indices
  for (int i = 0; i < nrows; i++) {
    int jp_end = rowp[i + 1];
    for (int jp = rowp[i]; jp < jp_end; jp++) {
      int j = cols[jp];
      rows[colp[j]] = i;
      if (Avals) {
        ATvals[colp[j]] = Avals[jp];
      }
      colp[j]++;
    }
  }

  // Reset the colp array
  for (int j = ncols - 1; j >= 0; j--) {
    colp[j + 1] = colp[j];
  }
  colp[0] = 0;
}

// Compute the number of entries in the matrix product A * A^{T}
template <typename T>
void MatMatTransNumeric(int nrows, int ncols, const int *rowp, const int *cols,
                        const T *Avals, const int *colp, const int *rows,
                        const T *ATvals, const int *Bcolp, int *Brows, T *Bvals,
                        int *flag, T *tmp) {
  for (int i = 0; i < nrows; i++) {
    flag[i] = -1;
  }

  // P_{*j} = A_{*k} * A_{jk}
  for (int j = 0; j < nrows; j++) {
    int nz = 0;

    // Loop over the non-zero columns
    int kp_end = rowp[j + 1];
    for (int kp = rowp[j]; kp < kp_end; kp++) {
      T Ajk = Avals[kp];
      int k = cols[kp];

      // Add the non-zero pattern from column k
      int ip_end = colp[k + 1];
      for (int ip = colp[k]; ip < ip_end; ip++) {
        int i = rows[ip];

        if (flag[i] != j) {
          flag[i] = j;
          tmp[i] = ATvals[ip] * Ajk;
          Brows[Bcolp[j] + nz] = i;
          nz++;
        } else {
          tmp[i] += ATvals[ip] * Ajk;
        }
      }
    }

    // Copy the values from the temporary column
    for (int k = 0; k < nz; k++) {
      Bvals[Bcolp[j] + k] = tmp[Brows[Bcolp[j] + k]];
    }
  }
}

// Compute the matrix C + A * D * A^{T}
template <typename T>
void MatMatTransNumeric(int nrows, int ncols, const T *cvals, const int *rowp,
                        const int *cols, const T *Avals, const T *dvals,
                        const int *colp, const int *rows, const T *ATvals,
                        const int *Bcolp, int *Brows, T *Bvals, int *flag,
                        T *tmp) {
  for (int i = 0; i < nrows; i++) {
    flag[i] = -1;
  }

  // P_{*j} = A_{*k} * A_{jk}
  for (int j = 0; j < nrows; j++) {
    int nz = 0;

    // Loop over the non-zero columns
    int kp_end = rowp[j + 1];
    for (int kp = rowp[j]; kp < kp_end; kp++) {
      int k = cols[kp];
      T dAjk = dvals[k] * Avals[kp];

      // Add the non-zero pattern from column k
      int ip_end = colp[k + 1];
      for (int ip = colp[k]; ip < ip_end; ip++) {
        int i = rows[ip];

        if (flag[i] != j) {
          flag[i] = j;
          if (i == j) {
            tmp[i] = cvals[i] + ATvals[ip] * dAjk;
          } else {
            tmp[i] = ATvals[ip] * dAjk;
          }
          Brows[Bcolp[j] + nz] = i;
          nz++;
        } else {
          tmp[i] += ATvals[ip] * dAjk;
        }
      }
    }

    // Copy the values from the temporary column
    for (int k = 0; k < nz; k++) {
      Bvals[Bcolp[j] + k] = tmp[Brows[Bcolp[j] + k]];
    }
  }
}

// Compute the number of entries in the matrix product A * A^{T}
inline int MatMatTransSymbolic(int nrows, int ncols, const int *rowp,
                               const int *cols, const int *colp,
                               const int *rows, int *Bcolp, int *flag) {
  for (int i = 0; i < nrows; i++) {
    Bcolp[i] = 0;
    flag[i] = -1;
  }

  // P_{*j} = A_{*k} * A_{jk}
  for (int j = 0; j < nrows; j++) {
    int nz = 0;

    // Loop over the non-zero columns
    int kp_end = rowp[j + 1];
    for (int kp = rowp[j]; kp < kp_end; kp++) {
      int k = cols[kp];

      // Add the non-zero pattern from column k
      int ip_end = colp[k + 1];
      for (int ip = colp[k]; ip < ip_end; ip++) {
        int i = rows[ip];

        if (flag[i] != j) {
          flag[i] = j;
          nz++;
        }
      }
    }

    Bcolp[j] = nz;
  }

  int nnz = 0;
  for (int j = 0; j < nrows; j++) {
    int tmp = Bcolp[j];
    Bcolp[j] = nnz;
    nnz += tmp;
  }
  Bcolp[nrows] = nnz;

  return nnz;
}

/*
  Sort an array of length len, then remove duplicate entries and
  entries with values -1.
*/
inline int RemoveDuplicates(int *array, int len, int exclude) {
  std::sort(array, array + len);

  // Remove any negative numbers
  int i = 0;  // location to take entries from
  int j = 0;  // location to place entries

  while (i < len && array[i] < 0) i++;

  if (exclude >= 0) {
    for (; i < len; i++, j++) {
      while ((i < len - 1) && (array[i] == array[i + 1])) i++;

      if (array[i] == exclude) {
        j--;
      } else if (i != j) {
        array[j] = array[i];
      }
    }
  } else {
    for (; i < len; i++, j++) {
      while ((i < len - 1) && (array[i] == array[i + 1])) i++;

      if (i != j) {
        array[j] = array[i];
      }
    }
  }

  return j;  // The new length of the array
}

// Sort and make the data structure unique - remove diagonal
inline void SortAndRemoveDuplicates(int nvars, int *rowp, int *cols,
                                    int remove_diagonal) {
  int begin = rowp[0];
  for (int i = 0; i < nvars; i++) {
    int len = rowp[i + 1] - begin;
    int new_len = -1;
    if (remove_diagonal) {
      new_len = RemoveDuplicates(&cols[begin], len, i);
    } else {
      new_len = RemoveDuplicates(&cols[begin], len);
    }

    if (begin != rowp[i]) {
      for (int k = 0; k < new_len; k++) {
        cols[rowp[i] + k] = cols[begin + k];
      }
    }

    begin = rowp[i + 1];
    rowp[i + 1] = rowp[i] + new_len;
  }
}

// Convert BSRMat to an unblocked, CSR format
template <typename T, index_t M, index_t N>
CSRMat<T> bsr_to_csr(BSRMat<T, M, N> bsr_mat) {
  index_t nbrows = bsr_mat.nbrows;
  index_t nbcols = bsr_mat.nbcols;
  index_t nnz = bsr_mat.nnz;

  index_t nrows = nbrows * M;
  index_t ncols = nbcols * N;
  index_t gnnz = bsr_mat.nnz * M * N;

  CSRMat<T> csr_mat(nrows, ncols, gnnz);

  // Populate rowp
  index_t i, nblocks_this_row;
  for (index_t ib = 0; ib < nbrows; ib++) {
    nblocks_this_row = bsr_mat.rowp(ib + 1) - bsr_mat.rowp(ib);
    for (index_t ii = 0; ii < M; ii++) {
      i = M * ib + ii;
      csr_mat.rowp(i) = (M * bsr_mat.rowp(ib) + ii * nblocks_this_row) * N;
    }
  }
  csr_mat.rowp(nrows) = gnnz;

  // Populate cols and vals
  index_t e = NO_INDEX;
  index_t MN = M * N;
  for (index_t ib = 0; ib < nbrows; ib++) {
    nblocks_this_row = bsr_mat.rowp(ib + 1) - bsr_mat.rowp(ib);
    index_t eb, index;
    for (eb = bsr_mat.rowp(ib), index = 0; eb < bsr_mat.rowp(ib + 1);
         eb++, index++) {
      for (index_t ii = 0; ii < M; ii++) {
        i = M * ib + ii;
        for (index_t jj = 0; jj < M; jj++) {
          e = csr_mat.rowp(i) + index * N + jj;
          csr_mat.vals(e) = bsr_mat.vals(eb, ii, jj);
          csr_mat.cols(e) = N * bsr_mat.cols(eb) + jj;
        }
      }
    }
  }

  return csr_mat;
}

}  // namespace A2D

#endif  // A2D_SPARSE_UTILS_INL_H