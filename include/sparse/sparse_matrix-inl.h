#ifndef A2D_SPARSE_MATRIX_INL_H
#define A2D_SPARSE_MATRIX_INL_H

#include "sparse/sparse_matrix.h"

namespace A2D {

// Find the value index of a block given indices (row, col) of a block
template <typename T, index_t M, index_t N>
index_t BSRMat<T, M, N>::find_value_index(index_t row, index_t col) {
  index_t jp_start = rowp[row];
  index_t jp_end = rowp[row + 1];

  for (index_t jp = jp_start; jp < jp_end; jp++) {
    if (cols[jp] == col) {
      return jp;
    }
  }

  return NO_INDEX;
}

// add values from an element matrix mat of shape (m, n)
template <typename T, index_t M, index_t N>
template <class Mat>
void BSRMat<T, M, N>::add_values(const index_t m, const index_t i[],
                                 const index_t n, const index_t j[], Mat &mat) {
  for (index_t ii = 0; ii < m; ii++) {
    index_t block_row = i[ii] / M;
    index_t local_row = i[ii] % M;

    for (index_t jj = 0; jj < n; jj++) {
      index_t block_col = j[jj] / N;
      index_t local_col = j[jj] % N;

      index_t jp = find_value_index(block_row, block_col);
      if (jp != NO_INDEX) {
        vals(jp, local_row, local_col) += mat(ii, jj);
      }
    }
  }
}

// Zero out rows and set diagonal entry to one for each zeroed row
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::zero_rows(const index_t nbcs, const index_t dof[]) {
  for (index_t ii = 0; ii < nbcs; ii++) {
    index_t block_row = dof[ii] / M;
    index_t local_row = dof[ii] % M;

    for (index_t jp = rowp[block_row]; jp < rowp[block_row + 1]; jp++) {
      for (index_t k = 0; k < N; k++) {
        vals(jp, local_row, k) = 0.0;
      }

      if (cols[jp] == block_row) {
        vals(jp, local_row, local_row) = 1.0;
      }
    }
  }
}

// Convert to a dense matrix
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = M * nbrows;
  index_t n = N * nbcols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t i = 0; i < nbrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];

      for (index_t ii = 0; ii < M; ii++) {
        const index_t irow = M * i + ii;
        for (index_t jj = 0; jj < N; jj++) {
          const index_t jcol = N * j + jj;
          A[n * irow + jcol] = vals(jp, ii, jj);
        }
      }
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

// Export the matrix as mtx format
template <typename T, index_t M, index_t N>
void BSRMat<T, M, N>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  }

  // Write global m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nbrows * M, nbcols * N, nnz * M * N);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t i = 0; i < nbrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];  // (i, j) is the block index pair

      for (index_t ii = 0; ii < M; ii++) {
        const index_t irow = M * i + ii + 1;  // convert to 1-based index
        for (index_t jj = 0; jj < N; jj++) {
          // (irow, jcol) is the entry coo
          const index_t jcol = N * j + jj + 1;  // convert to 1-based index
          T val = vals(jp, ii, jj);
          if constexpr (is_complex<T>::value) {
            std::fprintf(fp, "%d %d %30.20e %30.20e\n", irow, jcol, val.real(),
                         val.imag());
            nnz_mtx++;
          } else if (absfunc(val) >= epsilon) {
            std::fprintf(fp, "%d %d %30.20e\n", irow, jcol, val);
            nnz_mtx++;
          }
        }
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nbrows * M, nbcols * N, nnz_mtx);
  std::fclose(fp);
  return;
}

// Convert to a dense matrix
template <typename T>
void CSRMat<T>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = nrows;
  index_t n = ncols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t i = 0; i < nrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      index_t j = cols[jp];
      A[n * i + j] = vals(jp);
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

// Export the matrix as mtx format
template <typename T>
void CSRMat<T>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }

  // Write m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nrows, ncols, nnz);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t i = 0; i < nrows; i++) {
    for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
      if constexpr (is_complex<T>::value) {
        std::fprintf(fp, "%d %d %30.20e %30.20e\n", i + 1, cols[jp] + 1,
                     vals[jp].real(), vals[jp].imag());
        nnz_mtx++;
      } else if (absfunc(vals[jp]) >= epsilon) {
        std::fprintf(fp, "%d %d %30.20e\n", i + 1, cols[jp] + 1, vals[jp]);
        nnz_mtx++;
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nrows, ncols, nnz_mtx);
  std::fclose(fp);
  return;
}

// Zero out columns and set diagonal entry to one for each zeroed column
template <typename T>
void CSCMat<T>::zero_columns(const index_t nbcs, const index_t dof[]) {
  for (index_t ii = 0; ii < nbcs; ii++) {
    index_t column = dof[ii];

    for (index_t jp = colp[column]; jp < colp[column + 1]; jp++) {
      vals(jp) = 0.0;

      if (rows[jp] == column) {
        vals(jp) = 1.0;
      }
    }
  }
}

// Convert to a dense matrix
template <typename T>
void CSCMat<T>::to_dense(index_t *m_, index_t *n_, T **A_) {
  index_t m = nrows;
  index_t n = ncols;
  index_t size = m * n;

  T *A = new T[size];
  std::fill(A, A + size, T(0.0));

  for (index_t j = 0; j < ncols; j++) {
    for (index_t ip = colp[j]; ip < colp[j + 1]; ip++) {
      index_t i = rows[ip];
      A[n * i + j] = vals(ip);
    }
  }

  *A_ = A;
  *m_ = m;
  *n_ = n;
}

// Export the matrix as mtx format
template <typename T>
void CSCMat<T>::write_mtx(const std::string mtx_name, double epsilon) {
  // Open file and destroy old contents, if any
  std::FILE *fp = std::fopen(mtx_name.c_str(), "w");

  // Write header
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }

  // Write m, n and nnz
  std::fprintf(fp, "%15d%15d%15d\n", nrows, ncols, nnz);

  // Write entries
  index_t nnz_mtx = 0;
  for (index_t j = 0; j < ncols; j++) {
    for (index_t ip = colp[j]; ip < colp[j + 1]; ip++) {
      if constexpr (is_complex<T>::value) {
        std::fprintf(fp, "%d %d %30.20e %30.20e\n", rows[ip] + 1, j + 1,
                     vals[ip].real(), vals[ip].imag());
        nnz_mtx++;
      } else if (absfunc(vals[ip]) >= epsilon) {
        std::fprintf(fp, "%d %d %30.20e\n", rows[ip] + 1, j + 1, vals[ip]);
        nnz_mtx++;
      }
    }
  }
  std::fclose(fp);

  // Modify nnz
  fp = std::fopen(mtx_name.c_str(), "r+");
  if (is_complex<T>::value) {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate complex general\n");
  } else {
    std::fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
  }
  std::fprintf(fp, "%15d%15d%15d", nrows, ncols, nnz_mtx);
  std::fclose(fp);
  return;
}

}  // namespace A2D

#endif  // A2D_SPARSE_MATRIX_INL_H
