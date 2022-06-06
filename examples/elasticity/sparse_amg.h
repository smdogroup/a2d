#ifndef SPARSE_AMG_H
#define SPARSE_AMG_H

#include "block_numeric.h"
#include "sparse_matrix.h"
#include "sparse_numeric.h"
#include "sparse_symbolic.h"

namespace A2D {

/*
  Compute aggregates for a matrix A stored in CSR format

  - It is assumed that S is symmetric.
  - A may contain diagonal entries (self loops)
  - Unaggregated nodes are marked with a -1
 */
template <class I, class IdxArrayType, class AggArrayType>
I BSRMatStandardAggregation(const I nrows, const IdxArrayType& rowp,
                            const IdxArrayType& cols, AggArrayType& aggr,
                            AggArrayType& cpts) {
  // Bj[n] == -1 means i-th node has not been aggregated
  std::fill(aggr, aggr + nrows, 0);

  I next_aggregate = 1;

  // First pass
  for (I i = 0; i < nrows; i++) {
    if (aggr[i] == 0) {
      const I jp_start = rowp[i];
      const I jp_end = rowp[i + 1];

      // Determine whether all neighbors of this node are free (not already
      // aggregates)
      bool has_aggregated_neighbors = false;
      bool has_neighbors = false;
      for (I jp = jp_start; jp < jp_end; jp++) {
        const I j = cols[jp];
        if (i != j) {
          has_neighbors = true;
          if (aggr[j]) {
            has_aggregated_neighbors = true;
            break;
          }
        }
      }

      if (!has_neighbors) {
        // isolated node, do not aggregate
        aggr[i] = -nrows;
      } else if (!has_aggregated_neighbors) {
        // Make an aggregate out of this node and its neighbors
        aggr[i] = next_aggregate;
        cpts[next_aggregate - 1] = i;
        for (I jp = jp_start; jp < jp_end; jp++) {
          aggr[cols[jp]] = next_aggregate;
        }
        next_aggregate++;
      }
    }
  }

  // Second pass
  // Add unaggregated nodes to any neighboring aggregate
  for (I i = 0; i < nrows; i++) {
    if (aggr[i] == 0) {
      for (I jp = rowp[i]; jp < rowp[i + 1]; jp++) {
        const I j = cols[jp];

        const I xj = aggr[j];
        if (xj > 0) {
          aggr[i] = -xj;
          break;
        }
      }
    }
  }

  next_aggregate--;

  // Third pass
  for (I i = 0; i < nrows; i++) {
    if (aggr[i] != 0) {
      // node i has been aggregated
      if (aggr[i] > 0) {
        aggr[i] = aggr[i] - 1;
      } else if (aggr[i] == -nrows) {
        aggr[i] = -1;
      } else {
        aggr[i] = -aggr[i] - 1;
      }
    } else {
      // node i has not been aggregated
      const I jp_start = rowp[i];
      const I jp_end = rowp[i + 1];

      aggr[i] = next_aggregate;
      cpts[next_aggregate] = i;  // y stores a list of the Cpts

      for (I jp = jp_start; jp < jp_end; jp++) {
        const I j = cols[jp];

        if (aggr[j] == 0) {  // unmarked neighbors
          aggr[j] = next_aggregate;
        }
      }
      next_aggregate++;
    }
  }

  return next_aggregate;
}

/*
  Compute the tentative prolongation operator P.

  This code works by taking the tentative prolongation operators

  At this point, P has the non-zero pattern stored, but does not yet contain
  the numerical values.

  B = P * R
  P^{T} * P = I

  The columns of P are orthonormal.
*/
template <typename I, typename T, index_t M, index_t N>
void BSRMatMakeTentativeProlongation(const I nrows, const I num_aggregates,
                                     const I aggr[],
                                     MultiArray<T, CLayout<M, N>>& B,
                                     MultiArray<T, CLayout<N, N>>& R,
                                     BSRMat<I, T, M, N>** P,
                                     double toler = 1e-10) {
  // Form the non-zero pattern for PT
  std::vector<I> rowp(num_aggregates + 1, 0);

  I nnz = 0;
  for (I i = 0; i < nrows; i++) {
    if (aggr[i] >= 0) {
      rowp[aggr[i] + 1]++;
      nnz++;
    }
  }

  for (I i = 0; i < num_aggregates; i++) {
    rowp[i + 1] += rowp[i];
  }

  std::vector<I> cols(nnz);
  for (I i = 0; i < nrows; i++) {
    if (aggr[i] >= 0) {
      cols[rowp[aggr[i]]] = i;
      rowp[aggr[i]]++;
    }
  }

  for (I i = nrows; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Create the transpose of the matrix
  BSRMat<I, T, N, M> PT(num_aggregates, nrows, nnz, rowp, cols);

  for (I i = 0; i < PT.nbrows; i++) {
    // Copy the block from the near null-space candidates into P
    for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
      I j = PT.cols[jp];
      for (I k1 = 0; k1 < N; k1++) {
        for (I k2 = 0; k2 < M; k2++) {
          PT.Avals(jp, k1, k2) = B(j, k2, k1);
        }
      }
    }

    // Use modified Gram-Schmidt to orthogonalize the rows of the matrix
    for (I k = 0; k < N; k++) {
      // Take the initial norm of the row
      T init_norm = 0.0;
      for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
        for (I m = 0; m < M; m++) {
          init_norm += PT.Avals(jp, k, m) * PT.Avals(jp, k, m);
        }
      }
      init_norm = std::sqrt(init_norm);
      double theta = fabs(init_norm);

      // Take the dot product between rows k and j
      for (I j = 0; j < k; j++) {
        T dot = 0.0;
        for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
          for (I m = 0; m < M; m++) {
            dot += PT.Avals(jp, j, m) * PT.Avals(jp, k, m);
          }
        }

        // Record the result into the R basis
        R(i, j, k) = dot;

        // Subtract the dot product time row j from row k
        for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
          for (I m = 0; m < M; m++) {
            PT.Avals(jp, k, m) -= dot * PT.Avals(jp, j, m);
          }
        }
      }

      // Take the norm of row k
      T norm = 0.0;
      for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
        for (I m = 0; m < M; m++) {
          norm += PT.Avals(jp, k, m) * PT.Avals(jp, k, m);
        }
      }
      norm = std::sqrt(norm);

      // Compute the scalar factor for this row - zero rows
      // that are nearly linearly dependent
      T scale = 0.0;
      if (fabs(norm) > toler * theta) {
        scale = 1.0 / norm;
      }

      // Set the scale value
      R(i, k, k) = scale;
      for (I jp = PT.rowp[i]; jp < PT.rowp[i + 1]; jp++) {
        for (I m = 0; m < M; m++) {
          PT.Avals(jp, k, m) *= scale;
        }
      }
    }
  }

  *P = BSRMatMakeTranspose(PT);
}

/*
  Compute the Jacobi smoothing for the tentative prolongation operator P0

  P = (I - omega/rho(D^{-1} A) * rho(D^{-1} A ) * P0

  Note: This code destroys the entries in A!
*/
template <typename I, typename T, index_t M, index_t N>
BSRMat<I, T, M, N>* BSRJacobiProlongationSmoother(T omega,
                                                  BSRMat<I, T, M, M>& A,
                                                  BSRMat<I, T, M, N>& P0) {
  // Get the diagonal block D^{-1}
  bool inverse = true;
  BSRMat<I, T, M, M>& Dinv = *BSRMatExtractBlockDiagonal(A, inverse);

  // DinvA <- Dinv * A
  BSRMat<I, T, M, M>& DinvA = *BSRMatDuplicate(A);
  for (I i = 0; i < A.nbrows; i++) {
    auto D = MakeSlice(Dinv.Avals, Dinv.rowp[i]);
    for (I jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      auto A0 = MakeSlice(A.Avals, jp);
      auto DinvA0 = MakeSlice(DinvA.Avals, jp);
      blockGemm<T, M, M, M>(D, A0, DinvA0);
    }
  }

  // Estimate the spectral radius using Gerhsgorin
  T rho = BSRMatGershgorinSpectralEstimate(DinvA);

  // Compute the scalar multiple for the matrix-multiplication
  T scale = omega / rho;

  // Compute the non-zero pattern of the smoothed prolongation operator
  BSRMat<I, T, M, N>* P = BSRMatMatMultAddSymbolic(P0, DinvA, P0);

  // Copy values P <- P0
  BSRMatCopy(P0, *P);

  // Compute P = P0 - scale * Dinv * A * P0
  BSRMatMatMultAddScale(-scale, DinvA, P0, *P);

  return P;
}

// template <typename I, typename T, index_t M, index_t N>
// class BSRMatAmgLevelData {
//  public:
//   BSRMatAmgLevelData() : A(NULL), B(NULL), P(NULL), PT(NULL) {}

//   // The near null-space candidates
//   BSRMat<I, T, M, N>* B;

//   // Data for the matrices
//   BSRMat<I, T, M, M>* A;
//   BSRMat<I, T, M, N>* P;
//   BSRMat<I, T, N, M>* PT;

//   // Data for the smoother
//   T omega;
//   BSRMat<I, T, M, M>* Dinv;

//   // Data for the solution
//   MultiArray<T, CLayout<M>>* x;
//   MultiArray<T, CLayout<M>>* b;
//   MultiArray<T, CLayout<M>>* r;
// };

// template <typename I, typename T, index_t M, index_t N>
// void BSRMatSmoothedAmgLevel(T omega, BSRMat<I, T, M, M>& A,
//                             MultiArray<T, CLayout<M, N>>& B,
//                             BSRMat<I, T, M, N>* P, BSRMat<I, T, N, N>* Ar,
//                             BSRMat<I, T, N, N>* Br) {
//   // Compute the strength of connection S - need to fix this
//   // S = BSRMatStrength(A);

//   // Compute the aggregation - based on the strength of connection
//   std::vector<I> aggr(A.nbcols);
//   std::vector<I> cpts(A.nbcols);
//   I num_aggregates =
//       BSRMatStandardAggregation(A.nbrows, A.rowp, A.cols, aggr, cpts);

//   // Based on the aggregates, form a tentative prolongation operator
//   BSRMat<I, T, M, N>* P0;
//   BSRMatMakeTentativeProlongation(A.nbrows, num_aggregates, aggr, B, R, &P0);

//   // Smooth the prolongation operator
//   BSRMat<I, T, M, N>* P = BSRJacobiProlongationSmoother(omega, A, *P0);

//   // Make the transpose operator
//   BSRMat<I, T, N, M>* PT = BSRMatMakeTranspose(P);

//   // AP = A * P
//   BSRMat<I, T, M, N>* AP = BSRMatMatMultSymbolic(A, P);
//   BSRMatMatMult(A, P, AP);

//   // Ar = PT * AP = PT * A * P
//   BSRMat<I, T, N, N>* Ar = BSRMatMatMultSymbolic(PT, AP);
//   BSRMatMatMult(PT, AP, *Ar);
// }

// template <typename I, typename T, index_t M, index_t N>
// class BSRMatSmoothedAmg {
//  public:
//   int nlevels;

//   void applyFactor(MultiArray<T, CLayout<M>>& x, MultiArray<T, CLayout<M>>&
//   y) {

//   }

//   void applyMg(int level) {
//     // Pre-smooth at the current level
//     BSRApplySOR(data[level]->Dinv, data[level]->A, data[level]->b,
//                 data[level]->x);

//     // Compute r = b - A * x

//     // Restrict the residual to the next lowest level

//     // r

//     BSRApplySOR(data[level]->)
//   }

//  private:
// };

}  // namespace A2D

#endif  // SPARSE_AMG_H
