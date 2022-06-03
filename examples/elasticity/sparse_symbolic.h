#ifndef SPARSE_SYMBOLIC_H
#define SPARSE_SYMBOLIC_H

#include <set>
#include <vector>

#include "sparse_matrix.h"

namespace A2D {

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

  // // Do an in-place sort of the row and column data
  // SortCSRData(nnodes, node_to_elem_ptr, node_to_elem, node_to_elem);

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

  BSRMat<I, T, M, M>* bsr =
      new BSRMat<I, T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return bsr;
}

template <typename I, typename T, index_t M, class ConnArray>
BSRMat<I, T, M, M>* BSRMatFromConnectivity2(ConnArray& conn) {
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

  BSRMat<I, T, M, M>* bsr =
      new BSRMat<I, T, M, M>(nnodes, nnodes, nnz, rowp, cols);

  return bsr;
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

//                        I *rowp, I **cols) {
//   // Allocate space for the temporary row info
//   int *rlevs = new int[ncols];
//   int *rcols = new int[ncols];

//   for (I i = 0; i < nrows; i++) {
//     I nr = 0;  // Number of entries in the current row

//     // Add the matrix elements to the current row of the matrix.
//     // These new elements are sorted.
//     int diag_flag = 0;
//     for (I jp = rowp[i]; jp < rowp[i + 1]; jp++) {
//       if (mat->data->cols[j] == i) {
//         diag_flag = 1;
//       }
//       rcols[nr] = mat->data->cols[j];
//       rlevs[nr] = 0;
//       nr++;
//     }

//     // No diagonal element associated with row i, add one!
//     if (!diag_flag) {
//       nr = TacsMergeSortedArrays(nr, rcols, 1, &i);
//     }

//     // Now, perform the symbolic factorization -- this generates new
//     entries I j = 0; for (; rcols[j] < i; j++) {  // For entries in this
//     row, before the diagonal
//       int clev = rlevs[j];       // the level of fill for this entry

//       int p = j + 1;                   // The index into rcols
//       int k_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

//       // Start with the first entry after the diagonal in row, cols[j]
//       // k is the index into cols for row cols[j]
//       for (int k = diag[rcols[j]] + 1; k < k_end; k++) {
//         // Increment p to an entry where we may have cols[k] == rcols[p]
//         while (p < nr && rcols[p] < cols[k]) {
//           p++;
//         }

//         // The element already exists, check if it has a lower level of
//         fill
//         // and update the fill level if necessary
//         if (p < nr && rcols[p] == cols[k]) {
//           if (rlevs[p] > (clev + levs[k] + 1)) {
//             rlevs[p] = clev + levs[k] + 1;
//           }
//         } else if ((clev + levs[k] + 1) <= levFill) {
//           // The element does not exist but should since the level of
//           // fill is low enough. Insert the new entry into the list,
//           // but keep the list sorted
//           for (int n = nr; n > p; n--) {
//             rlevs[n] = rlevs[n - 1];
//             rcols[n] = rcols[n - 1];
//           }

//           rlevs[p] = clev + levs[k] + 1;
//           rcols[p] = cols[k];
//           nr++;
//         }
//       }
//     }

//     // Check if the size will be exceeded by adding the new elements
//     if (size + nr > max_size) {
//       int mat_ext = (int)((fill - 1.0) * mat_size);
//       if (nr > mat_ext) {
//         mat_ext = nr;
//       }
//       max_size = max_size + mat_ext;
//       TacsExtendArray(&cols, size, max_size);
//       TacsExtendArray(&levs, size, max_size);
//     }

//     // Now, put the new entries into the cols/levs arrays
//     for (int k = 0; k < nr; k++) {
//       cols[size] = rcols[k];
//       levs[size] = rlevs[k];
//       size++;
//     }

//     rowp[i + 1] = size;
//     diag[i] = j + rowp[i];
//   }
// }

// template <class I, class T>
// I maximal_independent_set_serial(const I num_rows, const I Ap[],
//                                  const int Ap_size, const I Aj[],
//                                  const int Aj_size, const T active, const T
//                                  C, const T F, T x[], const int x_size) {
//   I N = 0;

//   for (I i = 0; i < num_rows; i++) {
//     if (x[i] != active) continue;

//     x[i] = C;
//     N++;

//     for (I jj = Ap[i]; jj < Ap[i + 1]; jj++) {
//       const I j = Aj[jj];
//       if (x[j] == active) {
//         x[j] = F;
//       }
//     }
//   }

//   return N;
// }

// template <class I, class T, class R>
// I maximal_independent_set_parallel(const I num_rows, const I Ap[],
//                                    const int Ap_size, const I Aj[],
//                                    const int Aj_size, const T active, const
//                                    T C, const T F, T x[], const int x_size,
//                                    const R y[], const int y_size,
//                                    const I max_iters) {
//   I N = 0;
//   I num_iters = 0;

//   bool active_nodes = true;

//   while (active_nodes && (max_iters == -1 || num_iters < max_iters)) {
//     active_nodes = false;

//     num_iters++;

//     for (I i = 0; i < num_rows; i++) {
//       const R yi = y[i];

//       if (x[i] != active) continue;

//       const I row_start = Ap[i];
//       const I row_end = Ap[i + 1];

//       I jj;

//       for (jj = row_start; jj < row_end; jj++) {
//         const I j = Aj[jj];
//         const T xj = x[j];

//         if (xj == C) {
//           x[i] = F;  // neighbor is MIS
//           break;
//         }

//         if (xj == active) {
//           const R yj = y[j];
//           if (yj > yi)
//             break;  // neighbor is larger
//           else if (yj == yi && j > i)
//             break;  // tie breaker goes to neighbor
//         }
//       }

//       if (jj == row_end) {
//         for (jj = row_start; jj < row_end; jj++) {
//           const I j = Aj[jj];
//           if (x[j] == active) x[j] = F;
//         }
//         N++;
//         x[i] = C;
//       } else {
//         active_nodes = true;
//       }
//     }
//   }  // end while

//   return N;
// }

// /*
//   Symbolic factorization stage
// */
// template <class I>
// void BSRFactorSymbolic(int nbrows, int ncols, const I *rowp, const I *cols,
//                        I *rowp, I **cols) {
//   // Allocate space for the temporary row info
//   int *rlevs = new int[ncols];
//   int *rcols = new int[ncols];

//   for (I i = 0; i < nrows; i++) {
//     I nr = 0;  // Number of entries in the current row

//     // Add the matrix elements to the current row of the matrix.
//     // These new elements are sorted.
//     int diag_flag = 0;
//     for (I jp = rowp[i]; jp < rowp[i + 1]; jp++) {
//       if (mat->data->cols[j] == i) {
//         diag_flag = 1;
//       }
//       rcols[nr] = mat->data->cols[j];
//       rlevs[nr] = 0;
//       nr++;
//     }

//     // No diagonal element associated with row i, add one!
//     if (!diag_flag) {
//       nr = TacsMergeSortedArrays(nr, rcols, 1, &i);
//     }

//     // Now, perform the symbolic factorization -- this generates new
//     entries I j = 0; for (; rcols[j] < i; j++) {  // For entries in this
//     row, before the diagonal
//       int clev = rlevs[j];       // the level of fill for this entry

//       int p = j + 1;                   // The index into rcols
//       int k_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

//       // Start with the first entry after the diagonal in row, cols[j]
//       // k is the index into cols for row cols[j]
//       for (int k = diag[rcols[j]] + 1; k < k_end; k++) {
//         // Increment p to an entry where we may have cols[k] == rcols[p]
//         while (p < nr && rcols[p] < cols[k]) {
//           p++;
//         }

//         // The element already exists, check if it has a lower level of
//         fill
//         // and update the fill level if necessary
//         if (p < nr && rcols[p] == cols[k]) {
//           if (rlevs[p] > (clev + levs[k] + 1)) {
//             rlevs[p] = clev + levs[k] + 1;
//           }
//         } else if ((clev + levs[k] + 1) <= levFill) {
//           // The element does not exist but should since the level of
//           // fill is low enough. Insert the new entry into the list,
//           // but keep the list sorted
//           for (int n = nr; n > p; n--) {
//             rlevs[n] = rlevs[n - 1];
//             rcols[n] = rcols[n - 1];
//           }

//           rlevs[p] = clev + levs[k] + 1;
//           rcols[p] = cols[k];
//           nr++;
//         }
//       }
//     }

//     // Check if the size will be exceeded by adding the new elements
//     if (size + nr > max_size) {
//       int mat_ext = (int)((fill - 1.0) * mat_size);
//       if (nr > mat_ext) {
//         mat_ext = nr;
//       }
//       max_size = max_size + mat_ext;
//       TacsExtendArray(&cols, size, max_size);
//       TacsExtendArray(&levs, size, max_size);
//     }

//     // Now, put the new entries into the cols/levs arrays
//     for (int k = 0; k < nr; k++) {
//       cols[size] = rcols[k];
//       levs[size] = rlevs[k];
//       size++;
//     }

//     rowp[i + 1] = size;
//     diag[i] = j + rowp[i];
//   }
// }

}  // namespace A2D

#endif  // SPARSE_SYMBOLIC_H
