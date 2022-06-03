#ifndef SPARSE_NUMERIC_H
#define SPARSE_NUMERIC_H

#include "block_numeric.h"
#include "multiarray.h"
#include "sparse_matrix.h"

namespace A2D {

/*
  Add values into the BSR matrix
*/
template <typename I, typename T, index_t M, class ConnArray, class JacArray>
void BSRMatAddElementMatrices(ConnArray &conn, JacArray &jac,
                              BSRMat<I, T, M, M> &A) {
  for (I i = 0; i < conn.extent(0); i++) {
    for (I j1 = 0; j1 < conn.extent(1); j1++) {
      I row = conn(i, j1);

      for (I j2 = 0; j2 < conn.extent(1); j2++) {
        I col = conn(i, j2);
        I *col_ptr = A.find_column_index(row, col);

        if (col_ptr) {
          I jp = col_ptr - A.cols;
          auto Ab = MakeSlice(A.Avals, jp);

          for (I k1 = 0; k1 < M; k1++) {
            for (I k2 = 0; k2 < M; k2++) {
              Ab(k1, k2) += jac(i, j1, j2, k1, k2);
            }
          }
        }
      }
    }
  }
}

/*
  Compute the matrix-vector product: y = A * x
*/
template <typename I, typename T, index_t M, index_t N>
void BSRMatVecMult(BSRMat<I, T, M, N> &A, MultiArray<T, CLayout<N>> &x,
                   MultiArray<T, CLayout<M>> &y) {
  A2D::parallel_for(A.nbrows, [&](A2D::index_t i) -> void {
    auto yb = MakeSlice(y, i);
    // yb.zero();

    I jp = A.rowp[i];
    I jp_end = A.rowp[i + 1];
    for (; jp < jp_end; jp++) {
      I j = A.cols[jp];
      auto xb = MakeSlice(x, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvAdd<T, M, N>(Ab, xb, yb);
    }
  });
}

/*
  Compute the numerical matrix-matrix product

  C = A * B
*/
template <typename I, typename T, index_t M, index_t N, index_t P>
void BSRMatMatMult(BSRMat<I, T, M, N> &A, BSRMat<I, T, N, P> &B,
                   BSRMat<I, T, M, P> &C) {
  // C_{ik} = A_{ij} B_{jk}
  // for (I i = 0; i < C.nbrows; i++) {
  A2D::parallel_for(C.nbrows, [&](A2D::index_t i) -> void {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      auto Ab = MakeSlice(A.Avls, jp);

      I kp = B.rowp[j];
      I kp_end = B.rowp[j + 1];

      I cp = C.rowp[i];
      I cp_end = C.rowp[i + 1];

      for (; kp < kp_end; kp++) {
        while ((cp < cp_end) && (C.cols[cp] < B.cols[kp])) {
          cp++;
        }
        if (cp >= cp_end) {
          break;
        }
        if (B.cols[kp] == C.cols[cp]) {
          auto Bb = MakeSlice(B.Avals, kp);
          auto Cb = MakeSlice(C.Avals, cp);
          blockGemmAdd<T, M, N, P>(Ab, Bb, Cb);
        }
      }
    }
  });
}

/*
  Copy values from the matrix
*/
template <typename I, typename T, index_t M, index_t N>
void BSRMatCopy(BSRMat<I, T, M, N> &src, BSRMat<I, T, M, N> &dest) {
  // Zero the destination matrix
  dest.zero();

  if (src.nbrows != dest.nbrows || src.nbcols != dest.nbcols) {
    std::cerr << "BSRMatCopy: Matrix dimensions do not match" << std::endl;
    return;
  }

  for (I i = 0; i < src.nbrows; i++) {
    I kp = dest.rowp[i];
    I kp_end = dest.rowp[i + 1];

    I jp = src.rowp[i];
    I jp_end = src.rowp[i + 1];

    for (; (jp < jp_end) && (kp < kp_end); jp++) {
      while (dest.cols[kp] < src.cols[jp] && kp < kp_end) {
        kp++;
      }

      // Copy the matrix if the two entries are equal
      // and the size of both block--matrices is the same
      if (kp < kp_end) {
        if (dest.cols[kp] == src.cols[jp]) {
          for (I k1 = 0; k1 < M; k1++) {
            for (I k2 = 0; k2 < N; k2++) {
              dest.Avals(kp, k1, k2) = src.Avals(jp, k1, k2);
            }
          }
        } else {
          std::cerr << "BSRMatCopy: Non-zero pattern does not match cannot "
                       "copy values"
                    << std::endl;
        }
      }
    }
  }
}

/*
  Apply boundary conditions on a matrix
*/
template <typename I, typename T, index_t M, class BCArray>
void BSRMatZeroBCRows(BCArray &bcs, BSRMat<I, T, M, M> &A) {
  for (I i = 0; i < bcs.extent(0); i++) {
    for (I j = 0; j < M; j++) {
      I index = bcs(i, 0);
      if (bcs(i, 1) & (1U << j)) {
        for (I jp = A.rowp[index]; jp < A.rowp[index + 1]; jp++) {
          for (I k = 0; k < M; k++) {
            A.Avals(jp, j, k) = 0.0;
          }

          if (A.cols[jp] == i) {
            A.Avals(jp, j, j) = 1.0;
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
void VecZeroBCRows(BCArray &bcs, MultiArray<T, CLayout<M>> &x) {
  for (index_t i = 0; i < bcs.extent(0); i++) {
    for (index_t j = 0; j < M; j++) {
      index_t index = bcs(i, 0);
      if (bcs(i, 1) & (1U << j)) {
        x(index, j) = 0.0;
      }
    }
  }
}

/*
  Factor the matrix in place
*/
template <typename I, typename T, index_t M>
void BSRMatFactor(BSRMat<I, T, M, M> &A) {
  A2D::Vec<I, M> ipiv;
  A2D::Mat<T, M, M> D;

  // Store the diagonal entries
  I *diag = NULL;
  if (A.diag) {
    diag = A.diag;
  } else {
    diag = new I[A.nbrows];
  }

  for (I i = 0; i < A.nbrows; i++) {
    // Scan from the first entry in the current row, towards the
    // diagonal
    I row_end = A.rowp[i + 1];

    I jp = A.rowp[i];
    for (; A.cols[jp] < i; jp++) {
      I j = A.cols[jp];
      auto Ajp = MakeSlice(A.Avals, jp);
      auto Adiag = MakeSlice(A.Avals, diag[j]);

      // D = A[jp] * A[diag[j]]
      blockGemm<T, M, M, M>(Ajp, Adiag, D);

      // Scan through the remainder of row i
      I kp = jp + 1;

      // The final entry for row j
      I pp_end = A.rowp[j + 1];

      // Now, scan through row cj starting at the first entry past the
      // diagonal
      for (I pp = diag[j] + 1; (pp < pp_end) && (kp < row_end); pp++) {
        // Determine where the two rows have the same elements
        while (kp < row_end && A.cols[kp] < A.cols[pp]) {
          kp++;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (kp < row_end && A.cols[kp] == A.cols[pp]) {
          auto Akp = MakeSlice(A.Avals, kp);
          auto App = MakeSlice(A.Avals, pp);

          blockGemmSub<T, M, M, M>(D, App, Akp);
        }
      }

      // Copy the temporary matrix back
      // auto Ajp = MakeSlice(A.Avals, jp);
      for (I n = 0; n < M; n++) {
        for (I m = 0; m < M; m++) {
          Ajp(n, m) = D(n, m);
        }
      }
    }

    if (A.cols[jp] == i) {
      diag[i] = jp;
      auto Ab = MakeSlice(A.Avals, jp);

      // Invert the diagonal matrix component -- Invert( &A[b2*diag[i] )
      int fail = blockInverse<T, M>(Ab, D, ipiv);

      if (fail) {
        std::cerr << "Failure in factorization of row " << i << " block row "
                  << fail << std::endl;
      } else {
        for (I n = 0; n < M; n++) {
          for (I m = 0; m < M; m++) {
            Ab(n, m) = D(n, m);
          }
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
template <typename I, typename T, index_t M>
void BSRMatApplyLower(BSRMat<I, T, M, M> &A, MultiArray<T, CLayout<M>> &y) {
  for (I i = 0; i < A.nbrows; i++) {
    auto yi = MakeSlice(y, i);

    I end = A.diag[i];
    I jp = A.rowp[i];
    for (; jp < end; jp++) {
      I j = A.cols[jp];
      auto yj = MakeSlice(y, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvSub<T, M, M>(Ab, yj, yi);
    }
  }
}

/*
  Apply the upper factorization y = U^{-1} y
*/
template <typename I, typename T, index_t M>
void BSRMatApplyUpper(BSRMat<I, T, M, M> &A, MultiArray<T, CLayout<M>> &y) {
  A2D::Vec<T, M> ty;

  for (I i = A.nbrows - 1; i >= 0; i--) {
    auto yi = MakeSlice(y, i);
    for (I j = 0; j < M; j++) {
      ty(j) = yi(j);
    }

    I diag = A.diag[i];
    I end = A.rowp[i + 1];
    I jp = diag + 1;

    for (; jp < end; jp++) {
      I j = A.cols[jp];
      auto yj = MakeSlice(y, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvSub<T, M, M>(Ab, yj, ty);
    }

    auto D = MakeSlice(A.Avals, diag);
    blockGemv<T, M, M>(D, ty, yi);
  }
}

/*
  Apply the factorization y = U^{-1} L^{-1} x
*/
template <typename I, typename T, index_t M>
void BSRMatApplyFactor(BSRMat<I, T, M, M> &A, MultiArray<T, CLayout<M>> &x,
                       MultiArray<T, CLayout<M>> &y) {
  for (I i = 0; i < A.nbrows; i++) {
    for (I j = 0; j < M; j++) {
      y(i, j) = x(i, j);
    }
  }

  BSRMatApplyLower<I, T, M>(A, y);
  BSRMatApplyUpper<I, T, M>(A, y);
}

// // /*!
//   Apply a given number of steps of SOR to the system A*x = b.
// */
// void BCSRMatApplySOR(BCSRMatData *Adata, BCSRMatData *Bdata, const
// int start,
//                      const int end, const int var_offset,
//                      const TacsScalar *Adiag, const TacsScalar omega,
//                      const TacsScalar *b, const TacsScalar *xext,
//                      TacsScalar *x) {
//   const int *Arowp = Adata->rowp;
//   const int *Acols = Adata->cols;

//   const int *Browp = NULL;
//   const int *Bcols = NULL;
//   if (Bdata) {
//     Browp = Bdata->rowp;
//     Bcols = Bdata->cols;
//   }

//   int bsize = Adata->bsize;
//   const int b2 = bsize * bsize;

//   TacsScalar *tx = new TacsScalar[bsize];

//   if (start < end) {
//     for (int i = start; i < end; i++) {
//       int bi = bsize * i;

//       // Copy the right-hand-side to the temporary vector
//       // for this row
//       for (int n = 0; n < bsize; n++) {
//         tx[n] = b[bi + n];
//       }

//       // Set the pointer to the beginning of the current
//       // row
//       TacsScalar *a = &Adata->A[b2 * Arowp[i]];

//       // Scan through the row and compute the result:
//       // tx <- b_i - A_{ij}*x_{j} for j != i
//       int end = Arowp[i + 1];
//       for (int k = Arowp[i]; k < end; k++) {
//         int j = Acols[k];
//         int bj = bsize * j;

//         if (i != j) {
//           for (int m = 0; m < bsize; m++) {
//             int bm = bsize * m;
//             for (int n = 0; n < bsize; n++) {
//               tx[m] -= a[bm + n] * x[bj + n];
//             }
//           }
//         }

//         // Increment the block pointer by bsize^2
//         a += b2;
//       }

//       if (Bdata && i >= var_offset) {
//         const int row = i - var_offset;

//         // Set the pointer to the row in B
//         a = &Bdata->A[36 * Browp[row]];
//         end = Browp[row + 1];

//         for (int k = Browp[row]; k < end; k++) {
//           int j = Bcols[k];
//           const TacsScalar *y = &xext[bsize * j];

//           for (int m = 0; m < bsize; m++) {
//             int bm = bsize * m;
//             for (int n = 0; n < bsize; n++) {
//               tx[m] -= a[bm + n] * y[n];
//             }
//           }

//           // Increment the block pointer by bsize^2
//           a += b2;
//         }
//       }

//       // Compute the first term in the update:
//       // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
//       for (int n = 0; n < bsize; n++) {
//         x[bi + n] = (1.0 - omega) * x[bi + n];
//       }

//       // Apply the diagonal inverse and add the result to
//       // the matrix
//       const TacsScalar *adiag = &Adiag[b2 * i];
//       for (int m = 0; m < bsize; m++) {
//         int bm = bsize * m;
//         for (int n = 0; n < bsize; n++) {
//           x[bi + m] += omega * adiag[bm + n] * tx[n];
//         }
//       }
//     }
//   } else {
//     // Go through the matrix with the forward ordering
//     for (int i = start; i < end; i++) {
//       int bi = bsize * i;

//       // Copy the right-hand-side to the temporary vector
//       // for this row
//       for (int n = 0; n < bsize; n++) {
//         tx[n] = b[bi + n];
//       }

//       // Set the pointer to the beginning of the current
//       // row
//       TacsScalar *a = &Adata->A[b2 * Arowp[i]];

//       // Scan through the row and compute the result:
//       // tx <- b_i - A_{ij}*x_{j} for j != i
//       int end = Arowp[i + 1];
//       for (int k = Arowp[i]; k < end; k++) {
//         int j = Acols[k];
//         int bj = bsize * j;

//         if (i != j) {
//           for (int m = 0; m < bsize; m++) {
//             int bm = bsize * m;
//             for (int n = 0; n < bsize; n++) {
//               tx[m] -= a[bm + n] * x[bj + n];
//             }
//           }
//         }

//         // Increment the block pointer by bsize^2
//         a += b2;
//       }

//       if (Bdata && i >= var_offset) {
//         const int row = i - var_offset;

//         // Set the pointer to the row in B
//         a = &Bdata->A[36 * Browp[row]];
//         end = Browp[row + 1];

//         for (int k = Browp[row]; k < end; k++) {
//           int j = Bcols[k];
//           const TacsScalar *y = &xext[bsize * j];

//           for (int m = 0; m < bsize; m++) {
//             int bm = bsize * m;
//             for (int n = 0; n < bsize; n++) {
//               tx[m] -= a[bm + n] * y[n];
//             }
//           }

//           // Increment the block pointer by bsize^2
//           a += b2;
//         }
//       }

//       // Compute the first term in the update:
//       // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
//       for (int n = 0; n < bsize; n++) {
//         x[bi + n] = (1.0 - omega) * x[bi + n];
//       }

//       // Apply the diagonal inverse and add the result to
//       // the matrix
//       const TacsScalar *adiag = &Adiag[b2 * i];
//       for (int m = 0; m < bsize; m++) {
//         int bm = bsize * m;
//         for (int n = 0; n < bsize; n++) {
//           x[bi + m] += omega * adiag[bm + n] * tx[n];
//         }
//       }
//     }
//   }

//   delete[] tx;

}  // namespace A2D

#endif  // SPARSE_NUMERIC_H
