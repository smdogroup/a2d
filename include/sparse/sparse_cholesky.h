#ifndef A2D_SPARSE_CHOLESKY_H
#define A2D_SPARSE_CHOLESKY_H

#include "sparse/sparse_matrix.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dlapack.h"

// Include METIS
extern "C" {
#include "metis.h"
}

namespace A2D {

enum class CholOrderingType { NATURAL, ND };

/*
  Class for the sparse Cholesky factorization.

  This class computes the Cholesky factorization of the matrix A such that

  L * L^{T} = P * A * P^{T}

  where P is the optional permutation matrix.

  This code uses a supernode/supervariable approach in which groups of columns
  with the same nonzero pattern are aggregated into a single block column. This
  enables the use of more level-3 BLAS.

  This is used as one method to solve the sparse systems that arise in the
  interior point method.
*/
template <typename T>
class SparseCholesky {
 public:
  SparseCholesky(int _size, const int *Acolp, const int *Arows,
                 CholOrderingType order = CholOrderingType::ND,
                 const int *_perm = nullptr) {
    construct_cholesky(_size, Acolp, Arows, order, _perm);
  }

  SparseCholesky(CSCMat<T> csc_mat,
                 CholOrderingType order = CholOrderingType::ND,
                 const int *_perm = nullptr, bool set_values = true) {
    construct_cholesky(csc_mat.nrows, (const int *)csc_mat.colp.data(),
                       (const int *)csc_mat.rows.data(), order, _perm);
    setValues(csc_mat);
  }

  ~SparseCholesky();

  // Set values into the Cholesky matrix
  void setValues(const CSCMat<T> &csc_mat) {
    _setValues(csc_mat.ncols, (const int *)csc_mat.colp.data(),
               (const int *)csc_mat.rows.data(), csc_mat.vals.data());
  };
  void setValues(int n, const int Acolp[], const int Arows[], const T Avals[]) {
    _setValues(n, Acolp, Arows, Avals);
  }

  // Factor the matrix
  int factor();

  // Solve the factored system with the specified right-hand-side
  void solve(T *x);

  // Get information about the factorization
  void getInfo(int *_size, int *_num_snodes, int *_nnzL);

 private:
  // Construct the class
  void construct_cholesky(int _size, const int *Acolp, const int *Arows,
                          CholOrderingType order = CholOrderingType::ND,
                          const int *_perm = nullptr);

  // Set values into the Cholesky matrix
  void _setValues(int n, const int Acolp[], const int Arows[], const T Avals[]);

  // Build the elimination tree/forest
  void buildForest(const int Acolp[], const int Arows[], int parent[],
                   int Lnz[]);

  // Initialize the supernodes/supervariables by detecting identical column
  // non-zero patterns
  int initSupernodes(const int parent[], const int Lnz[], int vtn[]);

  // Build the non-zero pattern for the Cholesky factorization
  void buildNonzeroPattern(const int Acolp[], const int Arows[],
                           const int parent[], int Lnz[]);

  // Perform the update to the diagonal matrix
  void updateDiag(const int lsize, const int nlrows, const int lfirst_var,
                  const int *lrows, T *L, const int diag_size, T *diag,
                  T *work);

  // Apply the update to the work column - uses BLAS level 3
  void updateWorkColumn(int lsize, int nl1rows, T *L1, int nl2rows, T *L2,
                        T *temp);

  // Apply the sparse column update
  void updateColumn(const int lwidth, const int nlcols, const int lfirst_var,
                    const int *lrows, int nrows, const int *arows, const T *A,
                    const int *brows, T *B);

  // Perform Cholesky factorization on the diagonal
  int factorDiag(const int diag_size, T *D);

  // Solve L * y = x and output x = y
  void solveDiag(int diag_size, T *L, int nhrs, T *x);

  // Solve L^{T} * y = x and output x = y
  void solveDiagTranspose(int diag_size, T *L, int nhrs, T *x);

  // The following are short cut inline functions.
  // Get the diagonal block index
  inline int get_diag_index(const int i, const int j) {
    if (i >= j) {
      return j + i * (i + 1) / 2;
    } else {
      return i + j * (j + 1) / 2;
    }
  }

  // Given the supernode index, return the pointer to the diagonal matrix
  inline T *get_diag_pointer(const int i) { return &data[data_ptr[i]]; }

  // Given the supernode index, the supernode size and the index into the rows
  // data, return the pointer to the lower factor
  inline T *get_factor_pointer(const int i, const int size, const int index) {
    const int dsize = size * (size + 1) / 2;
    const int offset = index - colp[i];
    return &data[data_ptr[i] + dsize + size * offset];
  }

  // Given the supernode index, the supernode size and the index into the rows
  // data, return the pointer to the lower factor
  inline T *get_factor_pointer(const int i, const int size) {
    const int dsize = size * (size + 1) / 2;
    return &data[data_ptr[i] + dsize];
  }

  // The dimension of the square matrix
  int size;

  // Permutation and inverese permultation for the matrix. Both of these may be
  // nullptr.
  int *perm, *iperm;

  // Temporary vector for solving with the permutation arrays - only allocated
  // if perm and iperm are allocated, otherwise it is nullptr
  T *temp;

  // The row indices for the strict lower-diagonal entries of each super node.
  // This does not contain the row indices for the supernode itself. Only
  // entries below the supernode.
  int *rows;

  // Pointer into the row indices for the strict lower block of the
  // matrix. This does not include the row indices for the supernode.
  int *colp;

  // Number of supernodes
  int num_snodes;

  // Supernode sizes - How many consecutive variables belong to this
  // supernode? sum_{i=1}^{num_snodes} snode_size = size
  int *snode_size;

  // Given the variable index, what is the corresponding supernode?
  int *var_to_snode;

  // Given the supernode, what is the first variable in the node?
  int *snode_to_first_var;

  // Given the supernode index, a pointer into the supernode data
  // This is computed as the following for k = 0 ... num_snodes
  // data_ptr[k] = sum_{i = 0}^{k} snode_size[i] * ((snode_size[i] + 1)/2 +
  // colp[i + 1] - colp[i])
  int *data_ptr;

  // Work_size = max_{i} (snode_size[i] * (colp[i+1] - colp[i]))
  int work_size;

  // The numerical data for all entries size = data_ptr[num_snodes]
  T *data;

  // Optional, if the matrix is given as a block compressed sparse row matrix,
  // we first convert it to a csc format
  CSCMat<T> csc_mat;
};

}  // namespace A2D

#include "sparse/sparse_cholesky-inl.h"

#endif  //  A2D_SPARSE_CHOLESKY_H