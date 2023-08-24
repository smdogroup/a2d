#ifndef A2D_SPARSE_AMG_H
#define A2D_SPARSE_AMG_H

#include <iomanip>
#include <iostream>
#include <memory>

#include "array.h"
#include "block_numeric.h"
#include "sparse_matrix.h"
#include "sparse_numeric.h"
#include "sparse_symbolic.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/*
  Compute aggregates for a matrix A stored in CSR format

  The output consists of an array with all positive entries. Entries in the
  aggregation array outside the range [0, nrows) indicate an unaggregated
  variable.
*/
template <class IdxArrayType>
index_t BSRMatStandardAggregation(const index_t nrows, const IdxArrayType& rowp,
                                  const IdxArrayType& cols,
                                  std::vector<index_t>& aggr,
                                  std::vector<index_t>& cpts) {
  const index_t not_aggregated = MAX_INDEX;
  const index_t do_not_aggregate = MAX_INDEX - 1;

  // Initially fill the array with the not_aggregated tag
  std::fill(aggr.begin(), aggr.begin() + nrows, not_aggregated);
  index_t num_aggregates = 0;

  // First pass
  for (index_t i = 0; i < nrows; i++) {
    if (aggr[i] == not_aggregated) {
      // Determine whether all neighbors of this node are free (not already
      // aggregates)
      bool has_aggregated_neighbors = false;
      bool has_neighbors = false;

      const index_t jp_end = rowp[i + 1];
      for (index_t jp = rowp[i]; jp < jp_end; jp++) {
        const index_t j = cols[jp];
        if (i != j) {
          has_neighbors = true;
          if (aggr[j] != not_aggregated) {
            has_aggregated_neighbors = true;
            break;
          }
        }
      }

      if (!has_neighbors) {
        // isolated node, do not aggregate
        aggr[i] = do_not_aggregate;  // do not aggregate this node
      } else if (!has_aggregated_neighbors) {
        // Make an aggregate out of this node and its neighbors
        aggr[i] = num_aggregates;
        cpts[num_aggregates] = i;
        for (index_t jp = rowp[i]; jp < jp_end; jp++) {
          aggr[cols[jp]] = num_aggregates;
        }
        num_aggregates++;
      }
    }
  }

  // Second pass
  // Add unaggregated nodes to any neighboring aggregate
  for (index_t i = 0; i < nrows; i++) {
    if (aggr[i] == not_aggregated) {
      for (index_t jp = rowp[i]; jp < rowp[i + 1]; jp++) {
        const index_t j = cols[jp];

        if (aggr[j] != not_aggregated) {
          aggr[i] = aggr[j];
          break;
        }
      }
    }
  }

  // Third pass
  for (index_t i = 0; i < nrows; i++) {
    if (aggr[i] == not_aggregated) {
      // node i has not been aggregated
      aggr[i] = num_aggregates;
      cpts[num_aggregates] = i;

      const index_t jp_end = rowp[i + 1];
      for (index_t jp = rowp[i]; jp < jp_end; jp++) {
        const index_t j = cols[jp];

        // if (aggr[j] == 0) {  // unmarked neighbors
        if (aggr[j] == not_aggregated) {  // unmarked neighbors
          aggr[j] = num_aggregates;
        }
      }
      num_aggregates++;
    } else if (aggr[i] == do_not_aggregate) {
      aggr[i] = nrows + 1;
    }
  }

  return num_aggregates;
}

/*
  Compute the tentative prolongation operator P.

  This code works by taking the tentative prolongation operators

  At this point, P has the non-zero pattern stored, but does not yet contain
  the numerical values.

  B = P * R
  P^{T} * P = index_t

  The columns of P are orthonormal.
*/
template <typename T, index_t M, index_t N>
BSRMat<T, M, N>* BSRMatMakeTentativeProlongation(
    const index_t nrows, const index_t num_aggregates,
    const std::vector<index_t>& aggr, const MultiArrayNew<T* [M][N]>& B,
    MultiArrayNew<T* [N][N]>& R, double toler = 1e-10) {
  // Form the non-zero pattern for PT
  std::vector<index_t> rowp(num_aggregates + 1, 0);

  index_t nnz = 0;
  for (index_t i = 0; i < nrows; i++) {
    if (int(aggr[i]) >= 0 && aggr[i] < num_aggregates) {
      rowp[aggr[i] + 1]++;
      nnz++;
    }
  }

  for (index_t i = 0; i < num_aggregates; i++) {
    rowp[i + 1] += rowp[i];
  }

  std::vector<index_t> cols(nnz);
  for (index_t i = 0; i < nrows; i++) {
    if (int(aggr[i]) >= 0 && aggr[i] < num_aggregates) {
      cols[rowp[aggr[i]]] = i;
      rowp[aggr[i]]++;
    }
  }

  for (index_t i = num_aggregates; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Zero the entries in R
  BLAS::zero(R);

  // Create the transpose of the matrix
  BSRMat<T, N, M> PT(num_aggregates, nrows, nnz, rowp, cols);

  for (index_t i = 0; i < PT.nbrows; i++) {
    // Copy the block from the near null-space candidates into P
    const index_t jp_end = PT.rowp[i + 1];
    for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
      index_t j = PT.cols[jp];
      for (index_t k1 = 0; k1 < N; k1++) {
        for (index_t k2 = 0; k2 < M; k2++) {
          PT.vals(jp, k1, k2) = B(j, k2, k1);
        }
      }
    }

    // Use modified Gram-Schmidt to orthogonalize the rows of the matrix
    for (index_t k = 0; k < N; k++) {
      // Take the initial norm of the row
      T init_norm = 0.0;
      for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
        for (index_t m = 0; m < M; m++) {
          init_norm += PT.vals(jp, k, m) * PT.vals(jp, k, m);
        }
      }
      init_norm = sqrt(init_norm);
      double theta = absfunc(init_norm);

      // Take the dot product between rows k and j
      for (index_t j = 0; j < k; j++) {
        T dot = 0.0;
        for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
          for (index_t m = 0; m < M; m++) {
            dot += PT.vals(jp, j, m) * PT.vals(jp, k, m);
          }
        }

        // Record the result into the R basis
        R(i, j, k) = dot;

        // Subtract the dot product time row j from row k
        for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
          for (index_t m = 0; m < M; m++) {
            PT.vals(jp, k, m) -= dot * PT.vals(jp, j, m);
          }
        }
      }

      // Take the norm of row k
      T norm = 0.0;
      for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
        for (index_t m = 0; m < M; m++) {
          norm += PT.vals(jp, k, m) * PT.vals(jp, k, m);
        }
      }
      norm = sqrt(norm);

      // Compute the scalar factor for this row - zero rows
      // that are nearly linearly dependent
      T scale = 0.0;
      if (absfunc(norm) > toler * theta) {
        scale = 1.0 / norm;
        R(i, k, k) = norm;
      } else {
        std::cerr << "BSRMatMakeTentativeProlongation: Zeroed column"
                  << std::endl;
        throw std::runtime_error(
            "BSRMatMakeTentativeProlongation: Zeroed column");  // TODO: maybe
                                                                // don't do this
      }

      // Set the scale value
      for (index_t jp = PT.rowp[i]; jp < jp_end; jp++) {
        for (index_t m = 0; m < M; m++) {
          PT.vals(jp, k, m) *= scale;
        }
      }
    }
  }

  return BSRMatMakeTranspose(PT);
}

/*
  Compute the Jacobi smoothing for the tentative prolongation operator P0

  P = (index_t - omega/rho(D^{-1} A) * rho(D^{-1} A ) * P0
*/
template <typename T, index_t M, index_t N>
BSRMat<T, M, N>* BSRJacobiProlongationSmoother(T omega, BSRMat<T, M, M>& A,
                                               BSRMat<T, M, M>& Dinv,
                                               BSRMat<T, M, N>& P0, T* rho_) {
  // Compute DinvA <- Dinv * A
  BSRMat<T, M, M>* DinvA = BSRMatDuplicate(A);
  for (index_t i = 0; i < A.nbrows; i++) {
    for (index_t jp = A.rowp[i]; jp < A.rowp[i + 1]; jp++) {
      blockGemmSlice<T, M, M, M>(Dinv.vals, Dinv.rowp[i], A.vals, jp,
                                 DinvA->vals, jp);
    }
  }

  // Estimate the spectral radius using Gerhsgorin
  // T rho = BSRMatGershgorinSpectralEstimate(*DinvA);

  // Spectral estimate using Arnoldi
  T rho = BSRMatArnoldiSpectralRadius(*DinvA, 30u);

  if (rho_) {
    *rho_ = rho;
  }

  // Compute the scalar multiple for the matrix-multiplication
  T scale = omega / rho;

  // Compute the non-zero pattern of the smoothed prolongation operator
  BSRMat<T, M, N>* P = BSRMatMatMultAddSymbolic(P0, *DinvA, P0);

  // Copy values P0 -> P
  BSRMatCopy(P0, *P);

  // Compute (P0 - scale * Dinv * A * P0) -> P
  BSRMatMatMultAddScale(-scale, *DinvA, P0, *P);

  delete DinvA;

  return P;
}

/*
  Compute the strength of connection matrix with the given tolerance
*/
template <typename T, index_t M>
void BSRMatStrengthOfConnection(T epsilon, BSRMat<T, M, M>& A,
                                std::vector<index_t>& rowp,
                                std::vector<index_t>& cols) {
  // Frobenius norm squared for each diagonal entry
  std::vector<T> d(A.nbrows);
  T epsilon4 = epsilon * epsilon * epsilon * epsilon;

  if (A.diag.is_allocated()) {
    for (index_t i = 0; i < A.nbrows; i++) {
      index_t jp = A.diag[i];

      d[i] = 0.0;
      for (index_t ii = 0; ii < M; ii++) {
        for (index_t jj = 0; jj < M; jj++) {
          d[i] += A.vals(jp, ii, jj) * A.vals(jp, ii, jj);
        }
      }
    }
  } else {
    for (index_t i = 0; i < A.nbrows; i++) {
      index_t jp = A.find_value_index(i, i);

      if (jp != NO_INDEX) {
        d[i] = 0.0;
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t jj = 0; jj < M; jj++) {
            d[i] += A.vals(jp, ii, jj) * A.vals(jp, ii, jj);
          }
        }
      }
    }
  }

  rowp[0] = 0;
  for (index_t i = 0, nnz = 0; i < A.nbrows; i++) {
    index_t jp_end = A.rowp[i + 1];
    for (index_t jp = A.rowp[i]; jp < jp_end; jp++) {
      index_t j = A.cols[jp];

      if (i == j) {
        cols[nnz] = j;
        nnz++;
      } else {
        // Compute the Frobenius norm of the entry
        T af = 0.0;
        for (index_t ii = 0; ii < M; ii++) {
          for (index_t jj = 0; jj < M; jj++) {
            af += A.vals(jp, ii, jj) * A.vals(jp, ii, jj);
          }
        }

        if (RealPart(af * af) >= RealPart(epsilon4 * d[i] * d[j])) {
          cols[nnz] = j;
          nnz++;
        }
      }
    }

    rowp[i + 1] = nnz;
  }
}

/*
  Given the matrix A and the near null space basis B compute the block diagonal
  inverse of the matrix A, the prolongation and restriction operators, the
  reduced matrix Ar and the new near null space basis.
*/
template <typename T, index_t M, index_t N>
void BSRMatSmoothedAmgLevel(T omega, T epsilon, BSRMat<T, M, M>& A,
                            MultiArrayNew<T* [M][N]>& B, BSRMat<T, M, M>** Dinv,
                            BSRMat<T, M, N>** P, BSRMat<T, N, M>** PT,
                            BSRMat<T, N, N>** Ar, MultiArrayNew<T* [N][N]>& Br,
                            T* rho_) {
  index_t num_aggregates = 0;
  std::vector<index_t> aggr(A.nbcols);
  std::vector<index_t> cpts(A.nbcols);

  if (absfunc(epsilon) != 0.0) {
    // Compute the strength of connection S
    std::vector<index_t> Srowp(A.nbrows + 1);
    std::vector<index_t> Scols(A.nnz);
    BSRMatStrengthOfConnection(epsilon, A, Srowp, Scols);

    // Compute the aggregation - based on the strength of connection
    num_aggregates =
        BSRMatStandardAggregation(A.nbrows, Srowp, Scols, aggr, cpts);
  } else {
    num_aggregates =
        BSRMatStandardAggregation(A.nbrows, A.rowp, A.cols, aggr, cpts);
  }

  // Based on the aggregates, form a tentative prolongation operator
  MultiArrayNew<T* [N][N]> Br_("Br_", num_aggregates);

  BSRMat<T, M, N>* P0 =
      BSRMatMakeTentativeProlongation(A.nbrows, num_aggregates, aggr, B, Br_);

  // Get the diagonal block D^{-1}
  bool inverse = true;
  BSRMat<T, M, M>* Dinv_ = BSRMatExtractBlockDiagonal(A, inverse);

  // Smooth the prolongation operator
  BSRMat<T, M, N>* P_ =
      BSRJacobiProlongationSmoother(omega, A, *Dinv_, *P0, rho_);

  // Make the transpose operator
  BSRMat<T, N, M>* PT_ = BSRMatMakeTranspose(*P_);

  // AP = A * P
  BSRMat<T, M, N>* AP = BSRMatMatMultSymbolic(A, *P_);
  BSRMatMatMult(A, *P_, *AP);

  // Ar = PT * AP = PT * A * P
  BSRMat<T, N, N>* Ar_ = BSRMatMatMultSymbolic(*PT_, *AP);
  BSRMatMatMult(*PT_, *AP, *Ar_);

  // Copy over the values
  *Dinv = Dinv_;
  *P = P_;
  *PT = PT_;
  *Ar = Ar_;
  Br = Br_;

  delete P0;
  delete AP;
}

template <typename T, index_t M, index_t N>
class BSRMatAmg {
 public:
  BSRMatAmg(int num_levels, T omega, T epsilon,
            std::shared_ptr<BSRMat<T, M, M>> A, MultiArrayNew<T* [M][N]> B,
            bool print_info = false)
      : level(-1),
        A(A),
        B(B),
        P(NULL),
        PT(NULL),
        omega(omega),
        epsilon(epsilon),
        rho(0.0),
        Dinv(NULL),
        Afact(NULL),
        x(NULL),
        b(NULL),
        r(NULL),
        next(NULL) {
    makeAmgLevels(0, num_levels, print_info);
  }
  ~BSRMatAmg() {
    if (P) {
      delete P;
    }
    if (PT) {
      delete PT;
    }
    if (Dinv) {
      delete Dinv;
    }
    if (Afact) {
      delete Afact;
    }
    if (x) {
      delete x;
    }
    if (b) {
      delete b;
    }
    if (r) {
      delete r;
    }
    if (next) {
      delete next;
    }
  }

  /*
    Apply multigrid repeatedly until convergence
  */
  bool mg(MultiArrayNew<T* [M]>& b0, MultiArrayNew<T* [M]>& xk,
          index_t monitor = 0, index_t max_iters = 500, double rtol = 1e-8,
          double atol = 1e-30) {
    Timer timer("BSRMatAmg::mg()");
    // R == the residual
    MultiArrayNew<T* [M]> R("R", b0.layout());

    bool solve_flag = false;
    BLAS::zero(xk);
    BLAS::copy(R, b0);
    T init_norm = BLAS::norm(R);

    if (monitor) {
      std::printf("MG |A * x - b|[  0]: %20.10e\n", fmt(init_norm));
    }

    for (index_t iter = 0; iter < max_iters; iter++) {
      applyMg(b0, xk);

      BLAS::copy(R, b0);
      BSRMatVecMultSub(*A, xk, R);
      T res_norm = BLAS::norm(R);

      if (monitor && (iter + 1) % monitor == 0) {
        std::printf("MG |A * x - b|[%3d]: %20.10e\n", iter + 1, fmt(res_norm));
      }

      if (absfunc(res_norm) < atol ||
          absfunc(res_norm) < rtol * absfunc(init_norm)) {
        if (monitor && !((iter + 1) % monitor == 0)) {
          std::printf("MG |A * x - b|[%3d]: %20.10e\n", iter + 1,
                      fmt(res_norm));
        }
        solve_flag = true;
        break;
      }
    }
    return solve_flag;
  }

  /*
    Apply the preconditioned conjugate gradient method.

    This uses the variant of PCG from the paper "Inexact Preconditioned
    Conjugate Gradient Method with Inner-Outer Iteration" by Golub and Ye.
  */
  bool cg(const std::function<void(MultiArrayNew<T* [M]>&,
                                   MultiArrayNew<T* [M]>&)>& mat_vec,
          MultiArrayNew<T* [M]>& b0, MultiArrayNew<T* [M]>& xk,
          index_t monitor = 0, index_t max_iters = 500, double rtol = 1e-8,
          double atol = 1e-30, index_t iters_per_reset = 100) {
    Timer timer("BSRMatAmg::cg()");
    // R, Z and P and work are temporary vectors
    // R == the residual
    auto b0_layout = b0.layout();
    MultiArrayNew<T* [M]> R("R", b0_layout);
    MultiArrayNew<T* [M]> Z("Z", b0_layout);
    MultiArrayNew<T* [M]> P("P", b0_layout);
    MultiArrayNew<T* [M]> work("work", b0_layout);

    bool solve_flag = false;
    BLAS::zero(xk);
    BLAS::copy(R, b0);  // R = b0
    T init_norm = BLAS::norm(R);

    for (index_t reset = 0, iter = 0; iter < max_iters; reset++) {
      if (reset > 0) {
        BLAS::copy(R, b0);          // R = b0
        mat_vec(xk, work);          // work = A * xk
        BLAS::axpy(R, -1.0, work);  // R = b0 - A * xk
      }

      if (monitor && reset == 0) {
        std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter, fmt(init_norm));
      }

      if (absfunc(init_norm) > atol) {
        // Apply the preconditioner Z = M^{-1} R
        applyFactor(R, Z);

        // Set P = Z
        BLAS::copy(P, Z);

        // Compute rz = (R, Z)
        T rz = BLAS::dot(R, Z);

        for (index_t i = 0; i < iters_per_reset && iter < max_iters;
             i++, iter++) {
          mat_vec(P, work);                   // work = A * P
          T alpha = rz / BLAS::dot(work, P);  // alpha = (R, Z)/(A * P, P)
          BLAS::axpy(xk, alpha, P);           // x = x + alpha * P
          BLAS::axpy(R, -alpha, work);        // R' = R - alpha * A * P

          T res_norm = BLAS::norm(R);

          if (monitor && (iter + 1) % monitor == 0) {
            std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter + 1,
                        fmt(res_norm));
          }

          if (absfunc(res_norm) < atol ||
              absfunc(res_norm) < rtol * absfunc(init_norm)) {
            if (monitor && !((iter + 1) % monitor == 0)) {
              std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter + 1,
                          fmt(res_norm));
            }

            solve_flag = true;
            break;
          }

          applyFactor(R, work);             // work = Z' = M^{-1} * R
          T rz_new = BLAS::dot(R, work);    // rz_new = (R', Z')
          T rz_old = BLAS::dot(R, Z);       // rz_old = (R', Z)
          T beta = (rz_new - rz_old) / rz;  // beta = (R', Z' - Z)/(R, Z)
          BLAS::axpby(P, 1.0, beta, work);  // P' = Z' + beta * P
          BLAS::copy(Z, work);              // Z <- Z'
          rz = rz_new;                      // rz <- (R', Z')
        }
      }

      if (solve_flag) {
        break;
      }
    }
    return solve_flag;
  }

  /*
    Apply one cycle of multigrid with the right-hand-side b and the non-zero
    solution x.
  */
  void applyMg(MultiArrayNew<T* [M]>& b_, MultiArrayNew<T* [M]>& x_) {
    // Set temporary variables
    MultiArrayNew<T* [M]>* bt = b;
    MultiArrayNew<T* [M]>* xt = x;
    b = &b_;
    x = &x_;
    applyMg();
    b = bt;
    x = xt;
  }

  /*
    Apply the multigrid cycle as a preconditioner. This will overwrite
    whatever entries are in x.
  */
  void applyFactor(MultiArrayNew<T* [M]>& b_, MultiArrayNew<T* [M]>& x_) {
    // Set temporary variables
    MultiArrayNew<T* [M]>* bt = b;
    MultiArrayNew<T* [M]>* xt = x;
    b = &b_;
    x = &x_;
    bool zero_solution = true;
    applyMg(zero_solution);
    b = bt;
    x = xt;
  }

  /*
    Update the values of Galerkin projection at each level without
    re-computing the basis
  */
  void update() {
    if (Afact) {
      // Copy values to the matrix
      BSRMatCopy(*A, *Afact);

      // Perform the numerical factorization
      BSRMatFactor(*Afact);
    } else if (next) {
      delete Dinv;
      bool inverse = true;
      Dinv = BSRMatExtractBlockDiagonal(*A, inverse);

      // AP = A * P
      BSRMat<T, M, N>* AP = BSRMatMatMultSymbolic(*A, *P);
      BSRMatMatMult(*A, *P, *AP);

      // next->A = PT * AP = PT * A * P
      BSRMatMatMult(*PT, *AP, *next->A);
      delete AP;

      next->update();
    }
  }

  /*
    Test the accuracy of the Galerkin operator
  */
  void testGalerkin() {
    auto x0 = r->duplicate();
    auto y0 = r->duplicate();
    auto xr = next->r->duplicate();
    auto yr1 = next->r->duplicate();
    auto yr2 = next->r->duplicate();
    BLAS::random(xr);

    // Compute P^{T} * A * P * xr
    BSRMatVecMult(*P, *xr, *x0);
    BSRMatVecMult(*A, *x0, *y0);
    BSRMatVecMult(*PT, *y0, *yr1);

    // Compute Ar * xr
    BSRMatVecMult(*next->A, *xr, *yr2);

    // compute the error
    BLAS::axpy(*yr1, -1.0, *yr2);
    T error = BLAS::norm(*yr1, *yr1);
    T rel_err = BLAS::norm(*yr1, *yr1) / BLAS::norm(*yr2, *yr2);
    std::printf("Galerkin operator check\n");
    std::printf(
        "||Ar * xr - P^{T} * A * P * xr||: %20.10e,  rel. err: %20.10e\n",
        fmt(error), fmt(rel_err));

    delete x0;
    delete y0;
    delete xr;
    delete yr1;
    delete yr2;
  }

 private:
  // Private constructor for initializing the class
  BSRMatAmg(T omega, T epsilon, std::shared_ptr<BSRMat<T, M, M>> A,
            MultiArrayNew<T* [M][N]> B)
      : level(-1),
        A(A),
        B(B),
        P(NULL),
        PT(NULL),
        omega(omega),
        epsilon(epsilon),
        rho(0.0),
        Dinv(NULL),
        Afact(NULL),
        x(NULL),
        b(NULL),
        r(NULL),
        next(NULL) {}

  // Declare a friend since the template parameters may be different at
  // different levels
  template <typename UT, index_t UM, index_t UN>
  friend class BSRMatAmg;

  // Make the different multigrid levels
  void makeAmgLevels(int _level, int num_levels, bool print_info) {
    // Set the multigrid level
    level = _level;

    if (level == num_levels - 1) {
      x = new MultiArrayNew<T* [M]>("x", A->nbrows);
      b = new MultiArrayNew<T* [M]>("b", A->nbrows);

      // Form the sparse factorization - if the matrix is large and sparse, use
      // AMD, otherwise don't bother re-ordering.
      if (A->nbrows >= 20 && A->nnz < 0.25 * A->nbrows * A->nbrows) {
        Afact = BSRMatAMDFactorSymbolic(*A);

      } else {
        Afact = BSRMatFactorSymbolic(*A);
      }

      // Copy values to the matrix
      BSRMatCopy(*A, *Afact);

      // Perform the numerical factorization
      BSRMatFactor(*Afact);

      if (print_info) {
        printf("%10d%15d%15d\n", level, A->nbrows, Afact->nnz);
      }
    } else {
      r = new MultiArrayNew<T* [M]>("r", A->nbrows);
      if (level > 0) {
        r = new MultiArrayNew<T* [M]>("r", A->nbrows);
        x = new MultiArrayNew<T* [M]>("x", A->nbrows);
        b = new MultiArrayNew<T* [M]>("b", A->nbrows);
      }

      // Multicolor order this level
      BSRMatMultiColorOrder(*A);

      // Multi-color the matrix
      BSRMat<T, N, N>* Ar;
      MultiArrayNew<T* [N][N]> Br;

      // Find the new level
      BSRMatSmoothedAmgLevel<T, M, N>(omega, epsilon, *A, B, &Dinv, &P, &PT,
                                      &Ar, Br, &rho);

      // Allocate the next level
      auto Anext = std::shared_ptr<BSRMat<T, N, N>>(Ar);
      auto Bnext = MultiArrayNew<T* [N][N]>(Br);
      next = new BSRMatAmg<T, N, N>(omega, epsilon, Anext, Bnext);

      if (print_info) {
        if (level == 0) {
          printf("%10s%15s%15s%15s%15s\n", "Level", "n(A)", "nnz(A)", "nnz(P)",
                 "rho");
        }
        printf("%10d%15d%15d%15d%15.5f\n", level, A->nbrows, A->nnz, P->nnz,
               fmt(rho));
      }

      next->makeAmgLevels(level + 1, num_levels, print_info);
    }
  }

  void applyMg(bool zero_solution = false) {
    if (Afact) {
      BSRMatApplyFactor(*Afact, *b, *x);
    } else {
      // Pre-smooth with either a zero or non-zero x
      if (zero_solution) {
        BLAS::zero(*x);
      }
      T omega0 = 1.0;
      BSRApplySSOR(*Dinv, *A, omega0, *b, *x);

      // Compute the residuals r = b - A * x
      BLAS::copy(*r, *b);
      BSRMatVecMultSub(*A, *x, *r);

      // Now zero the solution on all subsequent levels
      zero_solution = true;

      // Restrict the residual to the next lowest level
      BSRMatVecMult(*PT, *r, *next->b);

      // Apply multigrid on the next lowest level
      next->applyMg(zero_solution);

      // Interpolate up from the next lowest grid level
      BSRMatVecMultAdd(*P, *next->x, *x);

      // Post-smooth
      BSRApplySSOR(*Dinv, *A, omega0, *b, *x);
    }
  }

  // Integer for the level
  int level;

  // Data for the matrix
  std::shared_ptr<BSRMat<T, M, M>> A;

  // The near null-space candidates
  MultiArrayNew<T* [M][N]> B;

  // Data for the prolongation and restriction
  BSRMat<T, M, N>* P;
  BSRMat<T, N, M>* PT;

  T omega;    // Omega value for constructing the prolongation operator
  T epsilon;  // Strength of connection value
  T rho;      // Estimate of the spectral radius

  BSRMat<T, M, M>* Dinv;  // Block diagonal inverse

  // Data for the full factorization (on the lowest level only)
  BSRMat<T, M, M>* Afact;

  // Data for the solution
  MultiArrayNew<T* [M]>* x;
  MultiArrayNew<T* [M]>* b;
  MultiArrayNew<T* [M]>* r;

  BSRMatAmg<T, N, N>* next;
};

template <typename T, index_t M>
bool conjugate_gradient(
    const std::function<void(MultiArrayNew<T* [M]>&, MultiArrayNew<T* [M]>&)>&
        mat_vec,
    const std::function<void(MultiArrayNew<T* [M]>&, MultiArrayNew<T* [M]>&)>&
        apply_factor,
    MultiArrayNew<T* [M]>& b0, MultiArrayNew<T* [M]>& xk, index_t monitor = 0,
    index_t max_iters = 500, double rtol = 1e-8, double atol = 1e-30,
    index_t iters_per_reset = 100) {
  Timer timer("conjugate_gradient");
  // R, Z and P and work are temporary vectors
  // R == the residual
  auto b0_layout = b0.layout();
  MultiArrayNew<T* [M]> R("R", b0_layout);
  MultiArrayNew<T* [M]> Z("Z", b0_layout);
  MultiArrayNew<T* [M]> P("P", b0_layout);
  MultiArrayNew<T* [M]> work("work", b0_layout);

  bool solve_flag = false;
  BLAS::zero(xk);
  BLAS::copy(R, b0);  // R = b0
  T init_norm = BLAS::norm(R);

  for (index_t reset = 0, iter = 0; iter < max_iters; reset++) {
    if (reset > 0) {
      BLAS::copy(R, b0);          // R = b0
      mat_vec(xk, work);          // work = A * xk
      BLAS::axpy(R, -1.0, work);  // R = b0 - A * xk
    }

    if (monitor && reset == 0) {
      std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter, fmt(init_norm));
    }

    if (absfunc(init_norm) > atol) {
      // Apply the preconditioner Z = M^{-1} R
      apply_factor(R, Z);

      // Set P = Z
      BLAS::copy(P, Z);

      // Compute rz = (R, Z)
      T rz = BLAS::dot(R, Z);

      for (index_t i = 0; i < iters_per_reset && iter < max_iters;
           i++, iter++) {
        mat_vec(P, work);                   // work = A * P
        T alpha = rz / BLAS::dot(work, P);  // alpha = (R, Z)/(A * P, P)
        BLAS::axpy(xk, alpha, P);           // x = x + alpha * P
        BLAS::axpy(R, -alpha, work);        // R' = R - alpha * A * P

        T res_norm = BLAS::norm(R);

        if (monitor && (iter + 1) % monitor == 0) {
          std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter + 1,
                      fmt(res_norm));
        }

        if (absfunc(res_norm) < atol ||
            absfunc(res_norm) < rtol * absfunc(init_norm)) {
          if (monitor && !((iter + 1) % monitor == 0)) {
            std::printf("PCG |A * x - b|[%3d]: %20.10e\n", iter + 1,
                        fmt(res_norm));
          }

          solve_flag = true;
          break;
        }

        apply_factor(R, work);            // work = Z' = M^{-1} * R
        T rz_new = BLAS::dot(R, work);    // rz_new = (R', Z')
        T rz_old = BLAS::dot(R, Z);       // rz_old = (R', Z)
        T beta = (rz_new - rz_old) / rz;  // beta = (R', Z' - Z)/(R, Z)
        BLAS::axpby(P, 1.0, beta, work);  // P' = Z' + beta * P
        BLAS::copy(Z, work);              // Z <- Z'
        rz = rz_new;                      // rz <- (R', Z')
      }
    }

    if (solve_flag) {
      break;
    }
  }
  return solve_flag;
}

template <typename T, index_t M, index_t gmres_size>
bool fgmres(
    const std::function<void(MultiArrayNew<T* [M]>&, MultiArrayNew<T* [M]>&)>&
        mat_vec,
    const std::function<void(MultiArrayNew<T* [M]>&, MultiArrayNew<T* [M]>&)>&
        apply_factor,
    MultiArrayNew<T* [M]>& b0, MultiArrayNew<T* [M]>& x, index_t monitor = 0,
    index_t max_restart = 10, double rtol = 1e-8, double atol = 1e-30) {
  // Initialize the subspaces
  auto b0_layout = b0.layout();
  MultiArrayNew<T* [M]> W[gmres_size + 1];
  MultiArrayNew<T* [M]> Z[gmres_size];
  for (index_t i = 0; i < gmres_size + 1; i++) {
    W[i] = MultiArrayNew<T* [M]>("W", b0_layout);
  }
  for (index_t i = 0; i < gmres_size; i++) {
    Z[i] = MultiArrayNew<T* [M]>("Z", b0_layout);
  }

  // Allocate and initialize data
  index_t Hptr[gmres_size + 1];
  Hptr[0] = 0;
  for (index_t i = 0; i < gmres_size; i++) {
    Hptr[i + 1] = Hptr[i] + i + 2;
  }

  // The Hessenberg matrix
  const index_t hsize = 2 + ((gmres_size + 2) * (gmres_size + 1)) / 2;
  T H[hsize];
  T res[gmres_size + 1];  // The residual

  // The unitary Q matrix in the QR factorixation of H
  T Qsin[gmres_size], Qcos[gmres_size];

  T init_norm = 0.0;
  bool solve_flag = false;

  for (index_t reset = 0, iter = 0; reset < max_restart + 1; reset++) {
    // Compute the residual
    if (reset == 0) {
      BLAS::zero(x);
      BLAS::copy(W[0], b0);  // W[0] = b0

      init_norm = BLAS::norm(W[0]);  // The initial residual
      res[0] = init_norm;
      BLAS::scale(W[0], 1.0 / res[0]);  // W[0] = b/|| b ||
    } else {
      // If the initial guess is non-zero or restarting
      mat_vec(x, W[0]);
      BLAS::axpy(W[0], -1.0, b0);  // W[0] = A*x - b

      res[0] = BLAS::norm(W[0]);
      BLAS::scale(W[0], -1.0 / res[0]);  // W[0] = b/|| b ||
    }

    index_t niters = 0;  // Keep track of the size of the Hessenberg matrix

    if (monitor && reset == 0) {
      std::printf("GMRES |A * x - b|[%3d]: %20.10e\n", iter, fmt(init_norm));
    }

    if (std::fabs(std::real(res[0])) < atol) {
      break;
    }

    for (index_t i = 0; i < gmres_size; i++, iter++) {
      // Apply the preconditioner, Z[i] = M^{-1} W[i]
      apply_factor(W[i], Z[i]);
      mat_vec(Z[i], W[i + 1]);  // W[i+1] = A*Z[i] = A*M^{-1}*W[i]

      // Build the orthogonal basis using MGS
      // orthogonalize(&H[Hptr[i]], W[i + 1], W, i + 1);
      for (index_t j = 0; j < i + 1; j++) {
        H[Hptr[i] + j] = BLAS::dot(W[j], W[i + 1]);
        BLAS::axpy(W[i + 1], -H[Hptr[i] + j], W[j]);
      }

      H[i + 1 + Hptr[i]] = BLAS::norm(W[i + 1]);  // H[i+1,i] = || W[i+1] ||
      BLAS::scale(W[i + 1],
                  1.0 / H[i + 1 + Hptr[i]]);  // W[i+1] = W[i+1]/|| W[i+1] ||

      // Apply the existing part of Q to the new components of
      // the Hessenberg matrix
      T h1, h2;
      for (index_t k = 0; k < i; k++) {
        h1 = H[k + Hptr[i]];
        h2 = H[k + 1 + Hptr[i]];
        H[k + Hptr[i]] = h1 * Qcos[k] + h2 * Qsin[k];
        H[k + 1 + Hptr[i]] = -h1 * Qsin[k] + h2 * Qcos[k];
      }

      // Now, compute the rotation for the new column that was just added
      h1 = H[i + Hptr[i]];
      h2 = H[i + 1 + Hptr[i]];
      T sq = std::sqrt(h1 * h1 + h2 * h2);

      Qcos[i] = h1 / sq;
      Qsin[i] = h2 / sq;
      H[i + Hptr[i]] = h1 * Qcos[i] + h2 * Qsin[i];
      H[i + 1 + Hptr[i]] = -h1 * Qsin[i] + h2 * Qcos[i];

      // Update the residual
      h1 = res[i];
      res[i] = h1 * Qcos[i];
      res[i + 1] = -h1 * Qsin[i];

      if (monitor && (iter + 1) % monitor == 0) {
        std::printf("GMRES |A * x - b|[%3d]: %20.10e\n", iter + 1,
                    fmt(std::fabs(res[i + 1])));
      }

      niters++;

      if (std::fabs(std::real(res[i + 1])) < atol ||
          std::fabs(std::real(res[i + 1])) < rtol * std::real(init_norm)) {
        // Set the solve flag
        solve_flag = true;

        break;
      }
    }

    // Now, compute the solution - the linear combination of the
    // Arnoldi vectors. H is upper triangular

    // Compute the weights
    for (index_t ip = niters; ip > 0; ip--) {
      index_t i = ip - 1;
      for (index_t j = i + 1; j < niters; j++) {  //
        res[i] = res[i] - H[i + Hptr[j]] * res[j];
      }
      res[i] = res[i] / H[i + Hptr[i]];
    }

    // Compute the linear combination
    for (index_t i = 0; i < niters; i++) {
      BLAS::axpy(x, res[i], Z[i]);
    }

    if (solve_flag) {
      break;
    }
  }

  return solve_flag;
}

}  // namespace A2D

#endif  // A2D_SPARSE_AMG_H
