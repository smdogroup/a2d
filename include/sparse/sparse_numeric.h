#ifndef A2D_SPARSE_NUMERIC_H
#define A2D_SPARSE_NUMERIC_H

#include <stdexcept>

#include "array.h"
#include "block_numeric.h"
#include "parallel.h"
#include "slice_numeric.h"
#include "sparse_matrix.h"
#include "sparse_symbolic.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/*
  Scatter the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray, class NodeArray, class ElementArray>
void VecElementScatter(ConnArray &conn, NodeArray &X, ElementArray &Xe) {
  // Loop over the elements
  for (index_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (index_t j = 0; j < conn.extent(1); j++) {
      const index_t index = conn(i, j);

      // Loop over the variables
      for (index_t k = 0; k < X.extent(1); k++) {
        Xe(i, j, k) = X(index, k);
      }
    }
  }
}

/*
  Gather the variables stored at the nodes to the a data structure
  that stores the element variables at each element node.
*/
template <class ConnArray, class ElementArray, class NodeArray>
void VecElementGatherAdd(ConnArray &conn, ElementArray &Xe, NodeArray &X) {
  // Loop over the elements
  for (index_t i = 0; i < conn.extent(0); i++) {
    // Loop over each element nodes
    for (index_t j = 0; j < conn.extent(1); j++) {
      const index_t index = conn(i, j);

      // Loop over the variables
      for (index_t k = 0; k < X.extent(1); k++) {
        X(index, k) += Xe(i, j, k);
      }
    }
  }
}

/*
  Add values into the BSR matrix
*/
template <typename T, index_t M, class ConnArray, class JacArray>
void BSRMatAddElementMatrices(ConnArray &conn, JacArray &jac,
                              BSRMat<T, M, M> &A) {
  Timer t("BSRMatAddElementMatrices()");
  for (index_t i = 0; i < conn.extent(0); i++) {
    for (index_t j1 = 0; j1 < conn.extent(1); j1++) {
      index_t row = conn(i, j1);

      for (index_t j2 = 0; j2 < conn.extent(1); j2++) {
        index_t col = conn(i, j2);
        index_t jp = A.find_value_index(row, col);

        if (jp != NO_INDEX) {
          for (index_t k1 = 0; k1 < M; k1++) {
            for (index_t k2 = 0; k2 < M; k2++) {
              A.vals(jp, k1, k2) += jac(i, j1, j2, k1, k2);
            }
          }
        } else {
          std::cerr << "BSRMatAddElementMatrices: Missing entry (" << row
                    << ", " << col << ") " << std::endl;
        }
      }
    }
  }
}

/*
  Compute the matrix-vector product: y = A * x
*/
template <typename T, index_t M, index_t N>
void BSRMatVecMult(BSRMat<T, M, N> &A, MultiArrayNew<T *[N]> &x,
                   MultiArrayNew<T *[M]> &y) {
  parallel_for(
      A.nbrows, KOKKOS_LAMBDA(index_t i)->void {
        for (index_t ii = 0; ii < M; ii++) {
          y(i, ii) = T(0);
        }

        const index_t jp_end = A.rowp[i + 1];
        for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
          index_t j = A.cols[jp];

          blockGemvAddSlice<T, M, N>(A.vals, jp, x, j, y, i);
        }
      });
}

/*
  Compute the matrix-vector product: y += A * x
*/
template <typename T, index_t M, index_t N>
void BSRMatVecMultAdd(BSRMat<T, M, N> &A, MultiArrayNew<T *[N]> &x,
                      MultiArrayNew<T *[M]> &y) {
  parallel_for(
      A.nbrows, KOKKOS_LAMBDA(index_t i)->void {
        const index_t jp_end = A.rowp[i + 1];
        for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
          index_t j = A.cols[jp];

          blockGemvAddSlice<T, M, N>(A.vals, jp, x, j, y, i);
        }
      });
}

/*
  Compute the matrix-vector product: y -= A * x
*/
template <typename T, index_t M, index_t N>
void BSRMatVecMultSub(BSRMat<T, M, N> &A, MultiArrayNew<T *[N]> &x,
                      MultiArrayNew<T *[M]> &y) {
  parallel_for(
      A.nbrows, KOKKOS_LAMBDA(index_t i)->void {
        const index_t jp_end = A.rowp[i + 1];
        for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
          index_t j = A.cols[jp];

          blockGemvSubSlice<T, M, N>(A.vals, jp, x, j, y, i);
        }
      });
}

/*
  Compute the numerical matrix-matrix product

  C = A * B
*/
template <typename T, index_t M, index_t N, index_t P>
void BSRMatMatMult(BSRMat<T, M, N> &A, BSRMat<T, N, P> &B, BSRMat<T, M, P> &C) {
  // C_{ik} = A_{ij} B_{jk}
  // for (index_t i = 0; i < C.nbrows; i++) {
  C.zero();
  parallel_for(
      C.nbrows, KOKKOS_LAMBDA(index_t i)->void {
        for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
          index_t j = A.cols[jp];

          index_t kp = B.rowp[j];
          index_t kp_end = B.rowp[j + 1];

          index_t cp = C.rowp[i];
          index_t cp_end = C.rowp[i + 1];

          for (; kp < kp_end; kp++) {
            while ((cp < cp_end) && (C.cols[cp] < B.cols[kp])) {
              cp++;
            }
            if (cp >= cp_end) {
              break;
            }
            if (B.cols[kp] == C.cols[cp]) {
              blockGemmAddSlice<T, M, N, P>(A.vals, jp, B.vals, kp, C.vals, cp);
            }
          }
        }
      });
}

/*
  Compute the numerical matrix-matrix product

  C += scale * A * B
*/
template <typename T, index_t M, index_t N, index_t P>
void BSRMatMatMultAddScale(T scale, BSRMat<T, M, N> &A, BSRMat<T, N, P> &B,
                           BSRMat<T, M, P> &C) {
  // C_{ik} = A_{ij} B_{jk}
  parallel_for(
      C.nbrows, KOKKOS_LAMBDA(index_t i)->void {
        for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
          index_t j = A.cols[jp];

          index_t kp = B.rowp[j];
          index_t kp_end = B.rowp[j + 1];

          index_t cp = C.rowp[i];
          index_t cp_end = C.rowp[i + 1];

          for (; kp < kp_end; kp++) {
            while ((cp < cp_end) && (C.cols[cp] < B.cols[kp])) {
              cp++;
            }
            if (cp >= cp_end) {
              break;
            }
            if (B.cols[kp] == C.cols[cp]) {
              blockGemmAddScaleSlice<T, M, N, P>(scale, A.vals, jp, B.vals, kp,
                                                 C.vals, cp);
            }
          }
        }
      });
}

/*
  Copy values from the matrix
*/
template <typename T, index_t M, index_t N>
void BSRMatCopy(BSRMat<T, M, N> &src, BSRMat<T, M, N> &dest) {
  // Zero the destination matrix
  dest.zero();

  if (src.nbrows != dest.nbrows || src.nbcols != dest.nbcols) {
    std::cerr << "BSRMatCopy: Matrix dimensions do not match" << std::endl;
    return;
  }

  if (dest.perm.is_allocated() && dest.iperm.is_allocated()) {
    for (index_t i = 0; i < src.nbrows; i++) {
      index_t idest = dest.iperm[i];

      index_t jp = src.rowp[i];
      index_t jp_end = src.rowp[i + 1];

      for (; jp < jp_end; jp++) {
        index_t jdest = dest.iperm[src.cols[jp]];

        index_t kp = dest.find_value_index(idest, jdest);
        if (kp != NO_INDEX) {
          for (index_t k1 = 0; k1 < M; k1++) {
            for (index_t k2 = 0; k2 < N; k2++) {
              dest.vals(kp, k1, k2) = src.vals(jp, k1, k2);
            }
          }
        } else {
          std::cerr << "BSRMatCopy: Non-zero pattern does not match - cannot "
                       "copy values for permuted matrix"
                    << std::endl;
        }
      }
    }
  } else {
    for (index_t i = 0; i < src.nbrows; i++) {
      index_t kp = dest.rowp[i];
      index_t kp_end = dest.rowp[i + 1];

      index_t jp = src.rowp[i];
      index_t jp_end = src.rowp[i + 1];

      for (; (jp < jp_end) && (kp < kp_end); jp++) {
        while (dest.cols[kp] < src.cols[jp] && kp < kp_end) {
          kp++;
        }

        // Copy the matrix if the two entries are equal
        // and the size of both block--matrices is the same
        if (kp < kp_end) {
          if (dest.cols[kp] == src.cols[jp]) {
            for (index_t k1 = 0; k1 < M; k1++) {
              for (index_t k2 = 0; k2 < N; k2++) {
                dest.vals(kp, k1, k2) = src.vals(jp, k1, k2);
              }
            }
          } else {
            std::cerr << "BSRMatCopy: Non-zero pattern does not match - cannot "
                         "copy values"
                      << std::endl;
          }
        }
      }
    }
  }
}

/*
  Apply boundary conditions on a matrix
*/
template <typename T, index_t M, class BCArray>
void BSRMatZeroBCRows(BCArray &bcs, BSRMat<T, M, M> &A) {
  Timer t("BSRMatZeroBCRows()");
  for (index_t i = 0; i < bcs.extent(0); i++) {
    index_t index = bcs(i, 0);
    for (index_t j = 0; j < M; j++) {
      if (bcs(i, 1) & (1U << j)) {
        for (index_t jp = A.rowp[index]; jp < A.rowp[index + 1]; jp++) {
          for (index_t k = 0; k < M; k++) {
            A.vals(jp, j, k) = 0.0;
          }

          // Set the diagonal element to the identity
          if (A.cols[jp] == index) {
            A.vals(jp, j, j) = 1.0;
          }
        }
      }
    }
  }
}

/*
  Apply boundary conditions on a matrix
*/
template <typename T, index_t M, class BCArray>
void VecZeroBCRows(BCArray &bcs, MultiArrayNew<T *[M]> &x) {
  for (index_t i = 0; i < bcs.extent(0); i++) {
    index_t index = bcs(i, 0);
    for (index_t j = 0; j < M; j++) {
      if (bcs(i, 1) & (1U << j)) {
        x(index, j) = 0.0;
      }
    }
  }
}

template <typename T, index_t M, index_t N, class BCArray>
void VecZeroBCRows(BCArray &bcs, MultiArrayNew<T *[M][N]> &x) {
  for (index_t i = 0; i < bcs.extent(0); i++) {
    index_t index = bcs(i, 0);
    for (index_t j = 0; j < M; j++) {
      if (bcs(i, 1) & (1U << j)) {
        for (index_t k = 0; k < N; k++) {
          x(index, j, k) = 0.0;
        }
      }
    }
  }
}

/*
  Factor the matrix in place
*/
template <typename T, index_t M>
void BSRMatFactor(BSRMat<T, M, M> &A) {
  using IdxArray1D_t = MultiArrayNew<index_t *>;

  Vec<index_t, M> ipiv;
  Mat<T, M, M> D;

  // Store the diagonal entries
  IdxArray1D_t diag;
  if (A.diag.is_allocated()) {
    diag = A.diag;
  } else {
    diag = IdxArray1D_t("diag", A.nbrows);
  }

  for (index_t i = 0; i < A.nbrows; i++) {
    // Scan from the first entry in the current row, towards the
    // diagonal
    index_t row_end = A.rowp[i + 1];

    index_t jp = A.rowp[i];
    for (; A.cols[jp] < i; jp++) {
      index_t j = A.cols[jp];

      // D = A[jp] * A[diag[j]]
      blockGemmSlice<T, M, M, M>(A.vals, jp, A.vals, diag[j], D);

      // Scan through the remainder of row i
      index_t kp = jp + 1;

      // The final entry for row j
      index_t pp_end = A.rowp[j + 1];

      // Now, scan through row cj starting at the first entry past the
      // diagonal
      for (index_t pp = diag[j] + 1; (pp < pp_end) && (kp < row_end); pp++) {
        // Determine where the two rows have the same elements
        while (kp < row_end && A.cols[kp] < A.cols[pp]) {
          kp++;
        }

        // A[kp] = A[kp] - D * A[p]
        if (kp < row_end && A.cols[kp] == A.cols[pp]) {
          blockGemmSubSlice<T, M, M, M>(D, A.vals, pp, A.vals, kp);
        }
      }

      // Copy the temporary matrix back
      for (index_t n = 0; n < M; n++) {
        for (index_t m = 0; m < M; m++) {
          A.vals(jp, n, m) = D(n, m);
        }
      }
    }

    if (A.cols[jp] != i) {
      std::cerr << "BSRMatFactor: Failure in factorization of block row " << i
                << " - No diagonal" << std::endl;
      throw std::runtime_error(
          "BSRMatFactor failed");  // TODO: maybe don't do this
    }
    diag[i] = jp;

    // Invert the diagonal matrix component -- Invert( &A[b2*diag[i] )
    int fail = blockInverseSlice<T, M>(A.vals, jp, D, ipiv);

    if (fail) {
      std::cerr << "BSRMatFactor: Failure in factorization of block row " << i
                << " local row " << fail << std::endl;
      throw std::runtime_error(
          "BSRMatFactor failed");  // TODO: maybe don't do this
    } else {
      for (index_t n = 0; n < M; n++) {
        for (index_t m = 0; m < M; m++) {
          A.vals(jp, n, m) = D(n, m);
        }
      }
    }
  }

  // Copy over the diagonal information
  A.diag = diag;
}

/*
  Apply the lower factorization y = L^{-1} y
*/
template <typename T, index_t M>
void BSRMatApplyLower(BSRMat<T, M, M> &A, MultiArrayNew<T *[M]> &y) {
  if (A.perm.is_allocated() && A.iperm.is_allocated()) {
    for (index_t i = 0; i < A.nbrows; i++) {
      index_t end = A.diag[i];
      index_t jp = A.rowp[i];
      for (; jp < end; jp++) {
        index_t j = A.cols[jp];

        blockGemvSubSlice<T, M, M>(A.vals, jp, y, A.perm[j], y, A.perm[i]);
      }
    }
  } else {
    for (index_t i = 0; i < A.nbrows; i++) {
      index_t end = A.diag[i];
      index_t jp = A.rowp[i];
      for (; jp < end; jp++) {
        index_t j = A.cols[jp];

        blockGemvSubSlice<T, M, M>(A.vals, jp, y, j, y, i);
      }
    }
  }
}

/*
  Apply the upper factorization y = U^{-1} y
*/
template <typename T, index_t M>
void BSRMatApplyUpper(BSRMat<T, M, M> &A, MultiArrayNew<T *[M]> &y) {
  Vec<T, M> ty;

  if (A.perm.is_allocated() && A.iperm.is_allocated()) {
    for (index_t i = A.nbrows; i > 0; i--) {
      for (index_t j = 0; j < M; j++) {
        ty(j) = y(A.perm[i - 1], j);
      }

      index_t diag = A.diag[i - 1];
      index_t end = A.rowp[i];
      index_t jp = diag + 1;

      for (; jp < end; jp++) {
        index_t j = A.cols[jp];

        blockGemvSubSlice<T, M, M>(A.vals, jp, y, A.perm[j], ty);
      }

      blockGemvSlice<T, M, M>(A.vals, diag, ty, y, A.perm[i - 1]);
    }
  } else {
    for (index_t i = A.nbrows; i > 0; i--) {
      for (index_t j = 0; j < M; j++) {
        ty(j) = y(i - 1, j);
      }

      index_t diag = A.diag[i - 1];
      index_t end = A.rowp[i];
      index_t jp = diag + 1;

      for (; jp < end; jp++) {
        index_t j = A.cols[jp];

        blockGemvSubSlice<T, M, M>(A.vals, jp, y, j, ty);
      }

      blockGemvSlice<T, M, M>(A.vals, diag, ty, y, i - 1);
    }
  }
}

/*
  Apply the factorization y = U^{-1} L^{-1} x
*/
template <typename T, index_t M>
void BSRMatApplyFactor(BSRMat<T, M, M> &A, MultiArrayNew<T *[M]> &x,
                       MultiArrayNew<T *[M]> &y) {
  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t j = 0; j < M; j++) {
      y(i, j) = x(i, j);
    }
  }

  BSRMatApplyLower<T, M>(A, y);
  BSRMatApplyUpper<T, M>(A, y);
}

/*
  Extract the block-diagonal values and possibly take their inverse
*/
template <typename T, index_t M>
BSRMat<T, M, M> *BSRMatExtractBlockDiagonal(BSRMat<T, M, M> &A,
                                            bool inverse = false) {
  index_t nrows = A.nbrows;
  index_t ncols = A.nbcols;
  index_t nnz = 0;
  std::vector<index_t> rowp(nrows + 1);
  std::vector<index_t> cols(nrows);

  // Create the expected non-zero pattern and create the matrix
  nnz = 0;
  rowp[0] = 0;
  for (index_t i = 0; i < nrows; i++) {
    cols[i] = i;
    nnz++;
    rowp[i + 1] = nnz;
  }

  BSRMat<T, M, M> *D = new BSRMat<T, M, M>(nrows, ncols, nnz, rowp, cols);

  // Now extract the real matrix and set the true non-zero pattern
  Mat<T, M, M> Dinv;
  Vec<index_t, M> ipiv;

  D->nnz = 0;
  D->rowp[0] = 0;
  for (index_t i = 0; i < nrows; i++) {
    index_t jp = A.find_value_index(i, i);
    if (jp != NO_INDEX) {
      // Copy the values
      for (index_t k1 = 0; k1 < M; k1++) {
        for (index_t k2 = 0; k2 < M; k2++) {
          D->vals(D->nnz, k1, k2) = A.vals(jp, k1, k2);
        }
      }

      if (inverse) {
        auto D0 = Kokkos::subview(D->vals, D->nnz, Kokkos::ALL, Kokkos::ALL);
        int fail = blockPseudoInverse(D0, Dinv);
        // int fail = blockInverse<T, M>(D0, Dinv, ipiv);

        if (fail) {
          std::cerr << "BSRMatExtractBlockDiagonal: Failure in factorization "
                       "of block row "
                    << i << " local row " << fail << std::endl;
          throw std::runtime_error(
              "BSRMatExtractBlockDiagonal failed");  // TODO: maybe don't do
                                                     // this
        } else {
          for (index_t k1 = 0; k1 < M; k1++) {
            for (index_t k2 = 0; k2 < M; k2++) {
              D->vals(D->nnz, k1, k2) = Dinv(k1, k2);
            }
          }
        }
      }

      D->cols[D->nnz] = i;
      D->nnz++;
    } else {
      std::cerr << "BSRMatExtractBlockDiagonal: No block diagonal for row " << i
                << std::endl;
    }
    D->rowp[i + 1] = D->nnz;
  }

  return D;
}

/*
  Apply a step of SOR to the system A*x = b for non-zero x.
*/
template <typename T, index_t M>
void BSRApplySOR(BSRMat<T, M, M> &Dinv, BSRMat<T, M, M> &A, T omega,
                 MultiArrayNew<T *[M]> &b, MultiArrayNew<T *[M]> &x) {
  index_t nrows = A.nbrows;

  if (A.perm) {
    for (index_t color = 0, offset = 0; color < A.num_colors; color++) {
      const index_t count = A.color_count[color];

      // for (index_t irow = 0; irow < count; irow++) {
      parallel_for(
          count, KOKKOS_LAMBDA(index_t irow)->void {
            index_t i = A.perm[irow + offset];

            // Copy over the values
            Vec<T, M> t;
            for (index_t m = 0; m < M; m++) {
              t(m) = b(i, m);
            }

            const int jp_end = A.rowp[i + 1];
            for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
              index_t j = A.cols[jp];

              if (i != j) {
                blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
              }
            }

            // x = (1 - omega) * x + omega * D^{-1} * t
            for (index_t m = 0; m < M; m++) {
              x(i, m) = (1.0 - omega) * x(i, m);
            }

            blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
          });

      offset += count;
    }
  } else {
    Vec<T, M> t;

    for (index_t i = 0; i < nrows; i++) {
      // Copy over the values
      for (index_t m = 0; m < M; m++) {
        t(m) = b(i, m);
      }

      const int jp_end = A.rowp[i + 1];
      for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
        index_t j = A.cols[jp];

        if (i != j) {
          blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
        }
      }

      // x = (1 - omega) * x + omega * D^{-1} * t
      for (index_t m = 0; m < M; m++) {
        x(i, m) = (1.0 - omega) * x(i, m);
      }
      blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
    }
  }
}

/*
  Apply a step of Symmetric SOR to the system A*x = b for non-zero x.
*/
template <typename T, index_t M>
void BSRApplySSOR(BSRMat<T, M, M> &Dinv, BSRMat<T, M, M> &A, T omega,
                  MultiArrayNew<T *[M]> &b, MultiArrayNew<T *[M]> &x) {
  index_t nrows = A.nbrows;

  if (A.perm.is_allocated()) {
    for (index_t color = 0, offset = 0; color < A.num_colors; color++) {
      const index_t count = A.color_count[color];

      parallel_for(
          count, KOKKOS_LAMBDA(index_t irow)->void {
            index_t i = A.perm[irow + offset];

            // Copy over the values
            Vec<T, M> t;
            for (index_t m = 0; m < M; m++) {
              t(m) = b(i, m);
            }

            const int jp_end = A.rowp[i + 1];
            for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
              index_t j = A.cols[jp];

              if (i != j) {
                blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
              }
            }

            // x = (1 - omega) * x + omega * D^{-1} * t
            for (index_t m = 0; m < M; m++) {
              x(i, m) = (1.0 - omega) * x(i, m);
            }

            blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
          });

      offset += count;
    }

    index_t offset = A.nbrows - A.color_count[A.num_colors - 1];
    for (index_t color = A.num_colors; color > 0; color--) {
      const index_t count = A.color_count[color - 1];

      parallel_for(
          count, KOKKOS_LAMBDA(index_t irow)->void {
            index_t i = A.perm[irow + offset];

            // Copy over the values
            Vec<T, M> t;
            for (index_t m = 0; m < M; m++) {
              t(m) = b(i, m);
            }

            const int jp_end = A.rowp[i + 1];
            for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
              index_t j = A.cols[jp];

              if (i != j) {
                blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
              }
            }

            // x = (1 - omega) * x + omega * D^{-1} * t
            for (index_t m = 0; m < M; m++) {
              x(i, m) = (1.0 - omega) * x(i, m);
            }
            blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
          });

      if (color >= 2) {
        offset -= A.color_count[color - 2];
      }
    }
  } else {
    Vec<T, M> t;

    for (index_t i = 0; i < nrows; i++) {
      // Copy over the values
      for (index_t m = 0; m < M; m++) {
        t(m) = b(i, m);
      }

      const int jp_end = A.rowp[i + 1];
      for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
        index_t j = A.cols[jp];

        if (i != j) {
          blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
        }
      }

      // x = (1 - omega) * x + omega * D^{-1} * t
      for (index_t m = 0; m < M; m++) {
        x(i, m) = (1.0 - omega) * x(i, m);
      }
      blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
    }

    for (index_t index = nrows; index > 0; index--) {
      index_t i = index - 1;

      // Copy over the values
      for (index_t m = 0; m < M; m++) {
        t(m) = b(i, m);
      }

      const int jp_end = A.rowp[i + 1];
      for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
        index_t j = A.cols[jp];

        if (i != j) {
          blockGemvSubSlice<T, M, M>(A.vals, jp, x, j, t);
        }
      }

      // x = (1 - omega) * x + omega * D^{-1} * t
      for (index_t m = 0; m < M; m++) {
        x(i, m) = (1.0 - omega) * x(i, m);
      }
      blockGemvAddScaleSlice<T, M, M>(omega, Dinv.vals, i, t, x, i);
    }
  }
}

/*
  Estimate the spectral radius using Gerhsgorin's circle theorem
*/
template <typename T, index_t M>
T BSRMatGershgorinSpectralEstimate(BSRMat<T, M, M> &A) {
  // Estimate the spectral radius using Gershgorin and estimate rho
  T rho = 0.0;

  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t k = 0; k < M; k++) {
      T a = 0.0;
      T R = 0.0;

      for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
        for (index_t j = 0; j < M; j++) {
          R += absfunc(A.vals(jp, k, j));
        }

        if (A.cols[jp] == i) {
          a = A.vals(jp, k, k);
        }
      }

      T rho0 = a + (R - absfunc(a));

      if (absfunc(rho0) > absfunc(rho)) {
        rho = rho0;
      }
    }
  }

  return rho;
}

// Declare the LAPACK eigenvalue solver we're about to use
extern "C" {
extern void dgeev_(const char *JOBVL, const char *JOBVR, int *N, double *A,
                   int *LDA, double *WR, double *WI, double *VL, int *LDVL,
                   double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);
}

/*
  Estimate the spectral radius
*/
template <typename T, index_t M>
T BSRMatArnoldiSpectralRadius(BSRMat<T, M, M> &A, index_t size = 15) {
  // Adjust size if matrix is small
  if (A.nbcols * M < size) {
    size = A.nbcols * M;
  }

  // Allocate the Upper Hessenberg matrix of shape (size + 1, size)
  double H[size * (size + 1)];
  std::fill(H, H + size * (size + 1), 0.0);

  // Allocate space for the orthonormal basis vectors
  MultiArrayNew<T *[M]> W[size + 1];

  // Allocate the first vector
  W[0] = MultiArrayNew<T *[M]>("W[0]", A.nbrows);

  // Create an initial random vector
  BLAS::random(W[0]);
  auto norm = RealPart(BLAS::norm(W[0]));
  BLAS::scale(W[0], 1.0 / norm);

  for (index_t i = 0; i < size; i++) {
    // Allocate the next vector
    char label[256];
    snprintf(label, sizeof(label), "W[%d]", i + 1);
    W[i + 1] = MultiArrayNew<T *[M]>(label, A.nbrows);

    // Multiply by the matrix to get the next vector
    BSRMatVecMult(A, W[i], W[i + 1]);

    // Orthogonalize against the existing subspace
    for (index_t j = 0; j <= i; j++) {
      index_t index = i + j * size;  // row-major index for entry H(j, i)
      H[index] = RealPart(BLAS::dot(W[i + 1], W[j]));
      BLAS::axpy(W[i + 1], -H[index], W[j]);
    }

    // Add the term to the matrix
    index_t index = i + 1 + i * size;  // row-major index for entry H(i + 1, i)
    H[index] = RealPart(BLAS::norm(W[i + 1]));
    BLAS::scale(W[i + 1], 1.0 / H[index]);
  }

  // Allocate space for the real/complex eigenvalue
  double eigreal[size];
  double eigimag[size];

  // Compute the eigenspectrum of the reduced Hessenberg matrix
  int hsize = size;
  int lwork = 4 * size;
  double work[lwork];
  int ldv = 1;
  int ldh = size;
  int info = 0;
  dgeev_("N", "N", &hsize, H, &ldh, eigreal, eigimag, NULL, &ldv, NULL, &ldv,
         work, &lwork, &info);

  // Throw runtime error if dgeev failed
  if (info != 0) {
    char msg[256];
    std::snprintf(msg, sizeof(msg), "Eigensolver failed with exit code %d.",
                  info);
    throw std::runtime_error(msg);
  }

  // Find the maximum absolute eigenvalue
  T rho = 0.0;
  for (int i = 0; i < size; i++) {
    double val = sqrt(eigreal[i] * eigreal[i] + eigimag[i] * eigimag[i]);
    if (val > absfunc(rho)) {
      rho = val;
    }
  }

  return rho;
}

}  // namespace A2D

#endif  // A2D_SPARSE_NUMERIC_H
