#ifndef SPARSE_SYMBOLIC_H
#define SPARSE_SYMBOLIC_H

#include <algorithm>
#include <limits>
#include <set>
#include <vector>

#include "sparse_matrix.h"

namespace A2D {

/*
  Given the CSR data, sort each row
*/
template <typename I>
void SortCSRData(I nrows, std::vector<I>& rowp, std::vector<I>& cols) {
  // Sort the cols array
  typename std::vector<I>::iterator it1;
  typename std::vector<I>::iterator it2 = cols.begin();
  for (I i = 0; i < nrows; i++) {
    it1 = it2;
    it2 += rowp[i + 1] - rowp[i];
    std::sort(it1, it2);
  }
}

/*
  Compute the non-zero pattern of the matrix based on the connectivity pattern
*/
template <typename I, typename T, index_t M, class ConnArray>
BSRMat<I, T, M, M>* BSRMatFromConnectivity(ConnArray& conn) {
  // Set the number of elements
  I nelems = conn.extent(0);

  // Find the number of nodes
  I nnodes = 0;
  for (I i = 0; i < conn.extent(0); i++) {
    for (I j = 0; j < conn.extent(1); j++) {
      if (conn(i, j) > nnodes) {
        nnodes = conn(i, j);
      }
    }
  }
  nnodes++;

  // Insert all the nodes into the node set
  std::set<std::pair<I, I>> node_set;
  for (I i = 0; i < nelems; i++) {
    for (I j1 = 0; j1 < conn.extent(1); j1++) {
      for (I j2 = 0; j2 < conn.extent(1); j2++) {
        node_set.insert(std::pair<I, I>(conn(i, j1), conn(i, j2)));
      }
    }
  }

  // Find the number of nodes referenced by other nodes
  std::vector<I> rowp(nnodes + 1);

  typename std::set<std::pair<I, I>>::iterator it;
  for (it = node_set.begin(); it != node_set.end(); it++) {
    rowp[it->first + 1] += 1;
  }

  // Set the pointer into the rows
  rowp[0] = 0;
  for (I i = 0; i < nnodes; i++) {
    rowp[i + 1] += rowp[i];
  }

  I nnz = rowp[nnodes];
  std::vector<I> cols(nnz);

  for (it = node_set.begin(); it != node_set.end(); it++) {
    cols[rowp[it->first]] = it->second;
    rowp[it->first]++;
  }

  // Reset the pointer into the nodes
  for (I i = nnodes; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Sort the cols array
  SortCSRData(nnodes, rowp, cols);

  BSRMat<I, T, M, M>* A =
      new BSRMat<I, T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return A;
}

template <typename I, typename T, index_t M, class ConnArray>
BSRMat<I, T, M, M>* BSRMatFromConnectivityDeprecated(ConnArray& conn) {
  // Set the number of elements
  I nelems = conn.extent(0);

  // Find the number of nodes
  I nnodes = 0;
  for (I i = 0; i < conn.extent(0); i++) {
    for (I j = 0; j < conn.extent(1); j++) {
      if (conn(i, j) > nnodes) {
        nnodes = conn(i, j);
      }
    }
  }
  nnodes++;

  // Create data to store node -> element connectivity
  std::vector<I> node_to_elem_ptr(nnodes + 1);
  for (I i = 0; i < conn.extent(0); i++) {
    for (I j = 0; j < conn.extent(1); j++) {
      I node = conn(i, j);
      node_to_elem_ptr[node + 1]++;
    }
  }

  for (I i = 0; i < nnodes; i++) {
    node_to_elem_ptr[i + 1] += node_to_elem_ptr[i];
  }

  std::vector<I> node_to_elem(node_to_elem_ptr[nnodes]);
  for (int i = 0; i < conn.extent(0); i++) {
    for (I j = 0; j < conn.extent(1); j++) {
      I node = conn(i, j);
      node_to_elem[node_to_elem_ptr[node]] = i;
      node_to_elem_ptr[node]++;
    }
  }

  // Do an in-place sort of the row and column data
  SortCSRData(nnodes, node_to_elem_ptr, node_to_elem, node_to_elem);

  // Reset the element pointer
  std::vector<I> rowp(nnodes + 1);
  std::vector<I> cols;

  rowp[0] = 0;

  // The set of nodes for each
  std::set<I> node_set;  // The set of nodes
  for (I i = 0; i < nnodes; i++) {
    for (I j = node_to_elem_ptr[i]; j < node_to_elem_ptr[i + 1]; j++) {
      int elem = node_to_elem[j];

      for (I k = 0; k < conn.extent(1); k++) {
        node_set.insert(conn(elem, k));
      }
    }

    // Set the rowp indices for the next set
    rowp[i + 1] = rowp[i] + node_set.size();

    // Push the values
    typename std::set<I>::iterator it;
    for (it = node_set.begin(); it != node_set.end(); it++) {
      cols.push_back(*it);
    }

    node_set.clear();
  }

  I nnz = rowp[nnodes];

  BSRMat<I, T, M, M>* A =
      new BSRMat<I, T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return A;
}

/*
  Symbolic factorization stage
*/
template <typename I, typename T, index_t M>
BSRMat<I, T, M, M>* BSRMatFactorSymbolic(BSRMat<I, T, M, M>& A,
                                         double fill_factor = 5.0) {
  I nrows = A.nbrows;
  I nnz = 0;

  // Column indices associated with the current row
  std::vector<I> rcols(nrows);

  // Row, column and diagonal index data for the new factored matrix
  std::vector<I> rowp(nrows + 1);
  std::vector<I> cols(index_t(fill_factor * A.nnz));
  std::vector<I> diag(nrows);

  rowp[0] = 0;
  for (I i = 0; i < nrows; i++) {
    I nr = 0;  // Number of entries in the current row

    // Add the matrix elements to the current row of the matrix.
    // These new elements are sorted.
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      rcols[nr] = A.cols[jp];
      nr++;
    }

    // Now, perform the symbolic factorization-- this generates new entries
    // Loop over entries in this row, before the diagonal
    I j = 0;
    for (; rcols[j] < i; j++) {
      I p = j + 1;                    // The index into rcols
      I kp_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

      // Start with the first entry after the diagonal in row, cols[j]
      // k is the index into cols for row cols[j]
      for (I kp = diag[rcols[j]] + 1; kp < kp_end; kp++) {
        // Increment p to an entry where we may have cols[k] == rcols[p]
        while (p < nr && rcols[p] < cols[kp]) {
          p++;
        }

        // Add the element into the list of new entries
        if (p >= nr || rcols[p] != cols[kp]) {
          // Insert the new entry into the list and keep the list sorted
          for (I n = nr; n > p; n--) {
            rcols[n] = rcols[n - 1];
          }
          rcols[p] = cols[kp];
          nr++;
        }
      }
    }

    // Make sure that we don't exceed the capacity of the vector
    if (nnz + nr > cols.size()) {
      cols.resize(nnz + nr + A.nnz);
    }

    // Add the new elements
    for (I n = 0; n < nr; n++, nnz++) {
      cols[nnz] = rcols[n];
    }

    rowp[i + 1] = nnz;
    diag[i] = rowp[i] + j;
  }

  BSRMat<I, T, M, M>* bsr =
      new BSRMat<I, T, M, M>(nrows, nrows, nnz, rowp, cols);

  return bsr;
}

/*
  Compute the non-zero pattern for C = A * B
*/
template <typename I, typename T, index_t M, index_t N, index_t P>
BSRMat<I, T, M, P>* BSRMatMatMultSymbolic(BSRMat<I, T, M, N>& A,
                                          BSRMat<I, T, N, P>& B,
                                          double fill_factor = 2.0) {
  I nrows = A.nbrows;
  I ncols = B.nbcols;
  I nnz = 0;

  // Column indices associated with the current row
  const I empty = std::numeric_limits<I>::max();
  const I first_entry = std::numeric_limits<I>::max() - 1;
  std::vector<I> next(ncols, empty);

  // Row, column and diagonal index data for the new factored matrix
  std::vector<I> rowp(nrows + 1);
  std::vector<I> cols(index_t(fill_factor * A.nnz));

  // Compute the non-zero structure of the resulting matrix C = A * B
  // one row at a time
  for (I i = 0; i < A.nbrows; i++) {
    I head = first_entry;
    I num_cols = 0;  // The size of the temporary cols array

    // Add the non-zero pattern to this matrix from each row of B
    // for each column of A.
    I jp_end = A.rowp[i + 1];
    for (I jp = A.rowp[i]; jp < jp_end; jp++) {
      I j = A.cols[jp];

      // Merge the two arrays into cols
      I kp_end = B.rowp[j + 1];
      for (I kp = B.rowp[j]; kp < kp_end; kp++) {
        I k = B.cols[kp];
        if (next[k] == empty) {
          next[k] = head;
          head = k;
          num_cols++;
        }
      }
    }

    if (nnz + num_cols > cols.size()) {
      cols.resize(nnz + num_cols + A.nnz);
    }

    // Reverse through the list
    for (I j = 0; j < num_cols; j++) {
      cols[nnz] = head;
      nnz++;
      I temp = head;
      head = next[head];
      next[temp] = empty;
    }

    rowp[i + 1] = nnz;
  }

  // Sort the CSR data
  SortCSRData(nrows, rowp, cols);

  BSRMat<I, T, M, P>* bsr =
      new BSRMat<I, T, M, P>(nrows, ncols, nnz, rowp, cols);

  return bsr;
}

/*
  Compute the non-zero pattern for C = S + A * B
*/
template <typename I, typename T, index_t M, index_t N, index_t P>
BSRMat<I, T, M, P>* BSRMatMatMultAddSymbolic(BSRMat<I, T, M, P>& S,
                                             BSRMat<I, T, M, N>& A,
                                             BSRMat<I, T, N, P>& B,
                                             double fill_factor = 2.0) {
  I nrows = A.nbrows;
  I ncols = B.nbcols;
  I nnz = 0;

  // Column indices associated with the current row
  const I empty = std::numeric_limits<I>::max();
  const I first_entry = std::numeric_limits<I>::max() - 1;
  std::vector<I> next(ncols, empty);

  // Row, column and diagonal index data for the new factored matrix
  std::vector<I> rowp(nrows + 1);
  std::vector<I> cols(index_t(fill_factor * A.nnz));

  // Compute the non-zero structure of the resulting matrix C = A * B
  // one row at a time
  for (I i = 0; i < A.nbrows; i++) {
    int head = first_entry;
    I num_cols = 0;  // The size of the temporary cols array

    // Merge the two arrays into cols
    I jp_end = S.rowp[i + 1];
    for (I jp = S.rowp[i]; jp < jp_end; jp++) {
      I j = S.cols[jp];
      if (next[j] == empty) {
        next[j] = head;
        head = j;
        num_cols++;
      }
    }

    // Add the non-zero pattern to this matrix from each row of B
    // for each column of A.
    jp_end = A.rowp[i + 1];
    for (I jp = A.rowp[i]; jp < jp_end; jp++) {
      I j = A.cols[jp];

      I kp_end = B.rowp[j + 1];
      for (I kp = B.rowp[j]; kp < kp_end; kp++) {
        I k = B.cols[kp];
        if (next[k] == empty) {
          next[k] = head;
          head = k;
          num_cols++;
        }
      }
    }

    if (nnz + num_cols > cols.size()) {
      cols.resize(nnz + num_cols + A.nnz);
    }

    // Reverse through the list
    for (I j = 0; j < num_cols; j++) {
      cols[nnz] = head;
      nnz++;
      I temp = head;
      head = next[head];
      next[temp] = empty;
    }

    rowp[i + 1] = nnz;
  }

  // Sort the CSR data
  SortCSRData(nrows, rowp, cols);

  BSRMat<I, T, M, P>* bsr =
      new BSRMat<I, T, M, P>(nrows, ncols, nnz, rowp, cols);

  return bsr;
}

/*
  Compute the non-zero pattern of the transpose of the matrix
*/
template <typename I, typename T, index_t M, index_t N>
BSRMat<I, T, N, M>* BSRMatMakeTransposeSymbolic(BSRMat<I, T, M, N>& A) {
  // The number of rows and columns for the transposed matrix
  I nrows = A.nbcols;
  I ncols = A.nbrows;

  std::vector<I> rowp(nrows + 1, 0);

  // Count up the number of references
  for (I i = 0; i < A.nbrows; i++) {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      rowp[j + 1]++;
    }
  }

  for (I i = 0; i < nrows; i++) {
    rowp[i + 1] += rowp[i];
  }

  I nnz = rowp[nrows];
  std::vector<I> cols(nnz);
  for (I i = 0; i < A.nbrows; i++) {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      cols[rowp[j]] = i;
      rowp[j]++;
    }
  }

  // Re-set the rowp array
  for (I i = nrows; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Create the new BSR matrix
  BSRMat<I, T, N, M>* At =
      new BSRMat<I, T, N, M>(nrows, ncols, nnz, rowp, cols);

  return At;
}

/*
  Make a transpose matrix
*/
template <typename I, typename T, index_t M, index_t N>
BSRMat<I, T, N, M>* BSRMatMakeTranspose(BSRMat<I, T, M, N>& A) {
  BSRMat<I, T, N, M>* At = BSRMatMakeTransposeSymbolic(A);

  // Loop over the values in A
  for (I i = 0; i < A.nbrows; i++) {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      auto A0 = MakeSlice(A.Avals, jp);  // Set A0 = A(i, j)

      I* col_ptr = At->find_column_index(j, i);  // Find At(j, i)
      if (col_ptr) {
        I kp = col_ptr - At->cols;
        auto At0 = MakeSlice(At->Avals, kp);

        for (I k1 = 0; k1 < M; k1++) {
          for (I k2 = 0; k2 < N; k2++) {
            At0(k2, k1) = A0(k1, k2);
          }
        }
      } else {
        std::cerr
            << "BSRMatMakeTranspose: Non-zero pattern does not match cannot "
               "copy values"
            << std::endl;
      }
    }
  }

  return At;
}

/*
  Duplicate the non-zero pattern of A without copying the values
*/
template <typename I, typename T, index_t M, index_t N>
BSRMat<I, T, M, N>* BSRMatDuplicate(BSRMat<I, T, M, N>& A) {
  return new BSRMat<I, T, M, N>(A.nbrows, A.nbcols, A.nnz, A.rowp, A.cols);
}

}  // namespace A2D

#endif  // SPARSE_SYMBOLIC_H
