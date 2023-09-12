#include <vector>

#include "a2ddefs.h"
#include "sparse/sparse_matrix.h"
#include "sparse/sparse_utils.h"
#include "test_commons.h"

using namespace A2D;

class Environment : public ::testing::Environment {
 public:
  void SetUp() override { Kokkos::initialize(); }
  void TearDown() override { Kokkos::finalize(); }
};

// Create a new environment and initialize kokkos
::testing::Environment *const initialize_kokkos =
    ::testing::AddGlobalTestEnvironment(new Environment);

class BSRMatTest : public ::testing::Test {
 protected:
  static index_t constexpr M = 5;
  static index_t constexpr N = 3;
  using BSRMat_t = BSRMat<double, M, N>;
  using CSRMat_t = CSRMat<double>;
  using CSCMat_t = CSCMat<double>;

  void SetUp() override {
    index_t constexpr nbrows = 10;
    index_t constexpr nbcols = 8;

    srand(0);

    std::vector<index_t> rowp(nbrows + 1);
    for (index_t i = 0; i < nbrows; i++) {
      rowp[i] = rand() % (nbcols + 1);  // [0, nbrows]
    }
    rowp[5] = 0;

    index_t presum = 0, temp = 0;
    for (index_t i = 0; i < nbrows; i++) {
      temp = rowp[i];
      rowp[i] = presum;
      presum += temp;
    }
    index_t nnz = presum;
    rowp[nbrows] = nnz;

    std::vector<index_t> cols(nnz);
    for (index_t n = 0; n < nnz; n++) {
      cols[n] = rand() % nbcols;  // [0, nbcols - 1]
    }
    bsr_mat = BSRMat_t(nbrows, nbcols, nnz, rowp, cols);

    for (index_t n = 0; n < nnz; n++) {
      for (index_t ii = 0; ii < M; ii++) {
        for (index_t jj = 0; jj < N; jj++) {
          bsr_mat.vals(n, ii, jj) = (double)rand() / RAND_MAX;  // [0, 1]
        }
      }
    }
  }

  BSRMat_t bsr_mat;
};

TEST_F(BSRMatTest, BSR_TO_CSR) {
  CSRMat_t csr_mat = bsr_to_csr(bsr_mat);

  index_t m_bsr, n_bsr;
  double *vals_bsr;
  bsr_mat.to_dense(&m_bsr, &n_bsr, &vals_bsr);

  index_t m_csr, n_csr;
  double *vals_csr;
  csr_mat.to_dense(&m_csr, &n_csr, &vals_csr);

  EXPECT_EQ(m_bsr, m_csr);
  EXPECT_EQ(n_bsr, n_csr);
  EXPECT_VEC_EQ(m_bsr * n_bsr, vals_bsr, vals_csr);
}

TEST_F(BSRMatTest, BSR_TO_CSC) {
  CSCMat_t csc_mat = bsr_to_csc(bsr_mat);

  index_t m_bsr, n_bsr;
  double *vals_bsr;
  bsr_mat.to_dense(&m_bsr, &n_bsr, &vals_bsr);

  index_t m_csc, n_csc;
  double *vals_csc;
  csc_mat.to_dense(&m_csc, &n_csc, &vals_csc);

  EXPECT_EQ(m_bsr, m_csc);
  EXPECT_EQ(n_bsr, n_csc);
  EXPECT_VEC_EQ(m_bsr * n_bsr, vals_bsr, vals_csc);
}
