#ifndef SPARSE_NUMERIC_H
#define SPARSE_NUMERIC_H

#include "block_numeric.h"
#include "multiarray.h"
#include "sparse_matrix.h"
#include "sparse_symbolic.h"

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
template <typename I, typename T, index_t M, index_t N>
void BSRMatVecMult(BSRMat<I, T, M, N> &A, MultiArray<T, CLayout<N>> &x,
                   MultiArray<T, CLayout<M>> &y) {
  A2D::parallel_for(A.nbrows, [&](A2D::index_t i) -> void {
    auto yb = MakeSlice(y, i);
    yb.zero();

    const I jp_end = A.rowp[i + 1];
    for (I jp = A.rowp[i]; jp < jp_end; jp++) {
      I j = A.cols[jp];
      auto xb = MakeSlice(x, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvAdd<T, M, N>(Ab, xb, yb);
    }
  });
}

/*
  Compute the matrix-vector product: y += A * x
*/
template <typename I, typename T, index_t M, index_t N>
void BSRMatVecMultAdd(BSRMat<I, T, M, N> &A, MultiArray<T, CLayout<N>> &x,
                      MultiArray<T, CLayout<M>> &y) {
  A2D::parallel_for(A.nbrows, [&](A2D::index_t i) -> void {
    auto yb = MakeSlice(y, i);

    const I jp_end = A.rowp[i + 1];
    for (I jp = A.rowp[i]; jp < jp_end; jp++) {
      I j = A.cols[jp];
      auto xb = MakeSlice(x, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvAdd<T, M, N>(Ab, xb, yb);
    }
  });
}

/*
  Compute the matrix-vector product: y -= A * x
*/
template <typename I, typename T, index_t M, index_t N>
void BSRMatVecMultSub(BSRMat<I, T, M, N> &A, MultiArray<T, CLayout<N>> &x,
                      MultiArray<T, CLayout<M>> &y) {
  A2D::parallel_for(A.nbrows, [&](A2D::index_t i) -> void {
    auto yb = MakeSlice(y, i);

    const I jp_end = A.rowp[i + 1];
    for (I jp = A.rowp[i]; jp < jp_end; jp++) {
      I j = A.cols[jp];
      auto xb = MakeSlice(x, j);
      auto Ab = MakeSlice(A.Avals, jp);

      blockGemvSub<T, M, N>(Ab, xb, yb);
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
  C.zero();
  A2D::parallel_for(C.nbrows, [&](A2D::index_t i) -> void {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      auto Ab = MakeSlice(A.Avals, jp);

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
  Compute the numerical matrix-matrix product

  C += scale * A * B
*/
template <typename I, typename T, index_t M, index_t N, index_t P>
void BSRMatMatMultAddScale(T scale, BSRMat<I, T, M, N> &A,
                           BSRMat<I, T, N, P> &B, BSRMat<I, T, M, P> &C) {
  // C_{ik} = A_{ij} B_{jk}
  A2D::parallel_for(C.nbrows, [&](A2D::index_t i) -> void {
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      I j = A.cols[jp];
      auto Ab = MakeSlice(A.Avals, jp);

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
          blockGemmAddScale<T, M, N, P>(scale, Ab, Bb, Cb);
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

  if (dest.perm && dest.iperm) {
    for (I i = 0; i < src.nbrows; i++) {
      I idest = dest.iperm[i];

      I jp = src.rowp[i];
      I jp_end = src.rowp[i + 1];

      for (; jp < jp_end; jp++) {
        I jdest = dest.iperm[src.cols[jp]];

        I *col_ptr = dest.find_column_index(idest, jdest);
        if (col_ptr) {
          I kp = col_ptr - dest.cols;

          for (I k1 = 0; k1 < M; k1++) {
            for (I k2 = 0; k2 < N; k2++) {
              dest.Avals(kp, k1, k2) = src.Avals(jp, k1, k2);
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
template <typename I, typename T, index_t M, class BCArray>
void BSRMatZeroBCRows(BCArray &bcs, BSRMat<I, T, M, M> &A) {
  for (I i = 0; i < bcs.extent(0); i++) {
    I index = bcs(i, 0);
    for (I j = 0; j < M; j++) {
      if (bcs(i, 1) & (1U << j)) {
        for (I jp = A.rowp[index]; jp < A.rowp[index + 1]; jp++) {
          for (I k = 0; k < M; k++) {
            A.Avals(jp, j, k) = 0.0;
          }

          // Set the diagonal element to the identity
          if (A.cols[jp] == index) {
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
    index_t index = bcs(i, 0);
    for (index_t j = 0; j < M; j++) {
      if (bcs(i, 1) & (1U << j)) {
        x(index, j) = 0.0;
      }
    }
  }
}

template <typename T, index_t M, index_t N, class BCArray>
void VecZeroBCRows(BCArray &bcs, MultiArray<T, CLayout<M, N>> &x) {
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

        // A[kp] = A[kp] - D * A[p]
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

    if (A.cols[jp] != i) {
      std::cerr << "BSRMatFactor: Failure in factorization of block row " << i
                << " - No diagonal" << std::endl;
    }
    diag[i] = jp;
    auto Ab = MakeSlice(A.Avals, jp);

    // Invert the diagonal matrix component -- Invert( &A[b2*diag[i] )
    int fail = blockInverse<T, M>(Ab, D, ipiv);

    if (fail) {
      std::cerr << "BSRMatFactor: Failure in factorization of block row " << i
                << " local row " << fail << std::endl;
    } else {
      for (I n = 0; n < M; n++) {
        for (I m = 0; m < M; m++) {
          Ab(n, m) = D(n, m);
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
  if (A.perm && A.iperm) {
    for (I i = 0; i < A.nbrows; i++) {
      auto yi = MakeSlice(y, A.perm[i]);

      I end = A.diag[i];
      I jp = A.rowp[i];
      for (; jp < end; jp++) {
        I j = A.cols[jp];
        auto yj = MakeSlice(y, A.perm[j]);
        auto Ab = MakeSlice(A.Avals, jp);

        blockGemvSub<T, M, M>(Ab, yj, yi);
      }
    }
  } else {
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
}

/*
  Apply the upper factorization y = U^{-1} y
*/
template <typename I, typename T, index_t M>
void BSRMatApplyUpper(BSRMat<I, T, M, M> &A, MultiArray<T, CLayout<M>> &y) {
  A2D::Vec<T, M> ty;

  if (A.perm && A.iperm) {
    for (I i = A.nbrows; i > 0; i--) {
      auto yi = MakeSlice(y, A.perm[i - 1]);
      for (I j = 0; j < M; j++) {
        ty(j) = yi(j);
      }

      I diag = A.diag[i - 1];
      I end = A.rowp[i];
      I jp = diag + 1;

      for (; jp < end; jp++) {
        I j = A.cols[jp];
        auto yj = MakeSlice(y, A.perm[j]);
        auto Ab = MakeSlice(A.Avals, jp);

        blockGemvSub<T, M, M>(Ab, yj, ty);
      }

      auto D = MakeSlice(A.Avals, diag);
      blockGemv<T, M, M>(D, ty, yi);
    }
  } else {
    for (I i = A.nbrows; i > 0; i--) {
      auto yi = MakeSlice(y, i - 1);
      for (I j = 0; j < M; j++) {
        ty(j) = yi(j);
      }

      I diag = A.diag[i - 1];
      I end = A.rowp[i];
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

/*
  Extract the block-diagonal values and possibly take their inverse
*/
template <typename I, typename T, index_t M>
BSRMat<I, T, M, M> *BSRMatExtractBlockDiagonal(BSRMat<I, T, M, M> &A,
                                               bool inverse = false) {
  I nrows = A.nbrows;
  I ncols = A.nbcols;
  I nnz = 0;
  std::vector<I> rowp(nrows + 1);
  std::vector<I> cols(nrows);

  // Create the expected non-zero pattern and create the matrix
  nnz = 0;
  rowp[0] = 0;
  for (I i = 0; i < nrows; i++) {
    cols[i] = i;
    nnz++;
    rowp[i + 1] = nnz;
  }

  BSRMat<I, T, M, M> *D = new BSRMat<I, T, M, M>(nrows, ncols, nnz, rowp, cols);

  // Now extract the real matrix and set the true non-zero pattern
  A2D::Mat<T, M, M> Dinv;
  A2D::Vec<I, M> ipiv;

  D->nnz = 0;
  D->rowp[0] = 0;
  for (I i = 0; i < nrows; i++) {
    I *col_ptr = A.find_column_index(i, i);
    if (col_ptr) {
      I jp = col_ptr - A.cols;
      auto A0 = MakeSlice(A.Avals, jp);
      auto D0 = MakeSlice(D->Avals, D->nnz);

      // Copy the values
      for (I k1 = 0; k1 < M; k1++) {
        for (I k2 = 0; k2 < M; k2++) {
          D0(k1, k2) = A0(k1, k2);
        }
      }

      if (inverse) {
        int fail = blockInverse<T, M>(D0, Dinv, ipiv);
        if (fail) {
          std::cerr << "BSRMatExtractBlockDiagonal: Failure in factorization "
                       "of block row "
                    << i << " local row " << fail << std::endl;
        } else {
          for (I k1 = 0; k1 < M; k1++) {
            for (I k2 = 0; k2 < M; k2++) {
              D0(k1, k2) = Dinv(k1, k2);
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
template <typename I, typename T, index_t M>
void BSRApplySOR(BSRMat<I, T, M, M> &Dinv, BSRMat<I, T, M, M> &A, T omega,
                 MultiArray<T, CLayout<M>> &b, MultiArray<T, CLayout<M>> &x) {
  I nrows = A.nbrows;

  if (A.perm) {
    for (I color = 0, offset = 0; color < A.num_colors; color++) {
      const index_t count = A.color_count[color];

      // for (I irow = 0; irow < count; irow++) {
      A2D::parallel_for(count, [&](index_t irow) -> void {
        I i = A.perm[irow + offset];

        // Copy over the values
        A2D::Vec<T, M> t;
        for (I m = 0; m < M; m++) {
          t(m) = b(i, m);
        }

        const int jp_end = A.rowp[i + 1];
        for (I jp = A.rowp[i]; jp < jp_end; jp++) {
          I j = A.cols[jp];

          if (i != j) {
            auto xb = MakeSlice(x, j);
            auto Ab = MakeSlice(A.Avals, jp);

            blockGemvSub<T, M, M>(Ab, xb, t);
          }
        }

        // x = (1 - omega) * x + omega * D^{-1} * t
        auto xb = MakeSlice(x, i);
        for (I m = 0; m < M; m++) {
          xb(m) = (1.0 - omega) * xb(m);
        }

        auto D = MakeSlice(Dinv.Avals, i);
        blockGemvAddScale<T, M, M>(omega, D, t, xb);
      });

      offset += count;
    }
  } else {
    A2D::Vec<T, M> t;

    for (I i = 0; i < nrows; i++) {
      // Copy over the values
      for (I m = 0; m < M; m++) {
        t(m) = b(i, m);
      }

      const int jp_end = A.rowp[i + 1];
      for (I jp = A.rowp[i]; jp < jp_end; jp++) {
        I j = A.cols[jp];

        if (i != j) {
          auto xb = MakeSlice(x, j);
          auto Ab = MakeSlice(A.Avals, jp);

          blockGemvSub<T, M, M>(Ab, xb, t);
        }
      }

      // x = (1 - omega) * x + omega * D^{-1} * t
      auto xb = MakeSlice(x, i);
      for (I m = 0; m < M; m++) {
        xb(m) = (1.0 - omega) * xb(m);
      }

      auto D = MakeSlice(Dinv.Avals, i);
      blockGemvAddScale<T, M, M>(omega, D, t, xb);
    }
  }
}

/*
  Estimate the spectral radius using Gerhsgorin's circle theorem
  */
template <typename I, typename T, index_t M>
T BSRMatGershgorinSpectralEstimate(BSRMat<I, T, M, M> &A) {
  // Estimate the spectral radius using Gershgorin and estimate rho
  T rho = 0.0;

  for (I i = 0; i < A.nbrows; i++) {
    for (I k = 0; k < M; k++) {
      T a = 0.0;
      T R = 0.0;

      for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
        for (I j = 0; j < M; j++) {
          R += fabs(A.Avals(jp, k, j));
        }

        if (A.cols[jp] == i) {
          a = A.Avals(jp, k, k);
        }
      }

      T rho0 = a + (R - fabs(a));

      if (fabs(rho0) > fabs(rho)) {
        rho = rho0;
      }
    }
  }

  return rho;
}

}  // namespace A2D

#endif  // SPARSE_NUMERIC_H
