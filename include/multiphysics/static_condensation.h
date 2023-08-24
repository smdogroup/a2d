#ifndef A2D_FE_STATIC_CONDENSATION_H
#define A2D_FE_STATIC_CONDENSATION_H

#include "block_numeric.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "sparse/sparse_amg.h"

namespace A2D {

/**
 * @brief This matrix defines a static condensation matrix that takes the form
 *
 * A = [ B  E ]
 *     [ F  C ]
 *
 * The b-matrix is stored in a global format while the E, F and C matrices are
 * assumed to be element-independent - for instance from an L2 space with only
 * bubble dof.
 *
 * @tparam T
 * @tparam block_size
 * @tparam basis_offset
 * @tparam Basis
 */
template <typename T, index_t block_size, index_t basis_offset, class Basis,
          index_t null_size = 1>
class StaticCondensationMat {
 public:
  // Definition for the index type
  using I = index_t;

  // Define the sizes of the element-wise block matrices
  static const constexpr I bsize =
      Basis::template get_dof_offset<basis_offset>();
  static const constexpr I csize = Basis::ndof - bsize;

  // Definitions for the matrix types
  using BSRMatType = BSRMat<T, block_size, block_size>;
  using BSRMatAmgType = BSRMatAmg<T, block_size, null_size>;
  using BMatType = Mat<T, bsize, bsize>;
  using CMatType = Mat<T, csize, csize>;
  using EMatType = Mat<T, bsize, csize>;
  using FMatType = Mat<T, csize, bsize>;

  StaticCondensationMat(ElementMesh<Basis> &mesh)
      : mesh(mesh),
        bnull("bnull", mesh.get_num_dof() / block_size),
        x("x", mesh.get_num_cumulative_dof(basis_offset - 1)),
        f("f", mesh.get_num_cumulative_dof(basis_offset - 1)),
        B(mesh.get_num_elements()),
        C(mesh.get_num_elements()),
        E(mesh.get_num_elements()),
        F(mesh.get_num_elements()) {
    I nrows;
    std::vector<I> rowp, cols;
    mesh.template create_block_csr<block_size, basis_offset>(nrows, rowp, cols);

    // Create the shared pointer
    mat = std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);

    // Initialize the near null-space to an appropriate vector
    for (I i = 0; i < bnull.extent(0); i++) {
      for (I j = 0; j < bnull.extent(1); j++) {
        bnull(i, j, 0) = 1.0;
      }
    }
  }

  // Required DOF container object (different for each element vector
  // implementation)
  class FEMat {
   public:
    static const index_t size = bsize * bsize;

    FEMat(index_t elem,
          StaticCondensationMat<T, block_size, basis_offset, Basis> &elem_mat)
        : B(elem_mat.get_bmat(elem)),
          C(elem_mat.get_cmat(elem)),
          E(elem_mat.get_emat(elem)),
          F(elem_mat.get_fmat(elem)) {
      B.zero();
      C.zero();
      E.zero();
      F.zero();
    }

    /**
     * @brief Get a reference to the underlying element data
     *
     * @return A reference to the degree of freedom
     */
    T &operator()(const index_t i, const index_t j) {
      if (i < bsize) {
        if (j < bsize) {
          return B(i, j);
        } else {
          return E(i, j - bsize);
        }
      } else {
        if (j < bsize) {
          return F(i - bsize, j);
        } else {
          return C(i - bsize, j - bsize);
        }
      }
    }
    const T &operator()(const index_t i, const index_t j) const {
      if (i < bsize) {
        if (j < bsize) {
          return B(i, j);
        } else {
          return E(i, j - bsize);
        }
      } else {
        if (j < bsize) {
          return F(i - bsize, j);
        } else {
          return C(i - bsize, j - bsize);
        }
      }
    }

   private:
    // Variables for all the basis functions
    BMatType &B;
    CMatType &C;
    EMatType &E;
    FMatType &F;
  };

  /**
   * @brief Get the number of elements
   */
  index_t get_num_elements() const { return mesh.get_num_elements(); }

  /**
   * @brief Add the degree of freedom values to the element vector
   *
   * @param elem the element index
   * @param dof the FEDof object that stores a reference to the degrees of
   * freedom
   *
   * If FEDof contains a pointer to data, this function may do nothing
   */
  void add_element_values(index_t elem, FEMat &elem_mat) {
    I dof[Basis::ndof];
    int sign[Basis::ndof];
    if constexpr (Basis::nbasis > 0) {
      get_dof<0>(elem, dof, sign);
    }

    // Set the correct signs
    for (I i = 0; i < Basis::ndof; i++) {
      for (I j = 0; j < Basis::ndof; j++) {
        elem_mat(i, j) = sign[i] * sign[j] * elem_mat(i, j);
      }
    }
  }

  // Apply boundary conditions
  template <class BcBasisType>
  void zero_bcs(const DirichletBCs<BcBasisType> &bcs) {
    std::vector<I> is_bc(mesh.get_num_dof(), 0);

    const I *bc_dofs = NULL;
    I nbcs = bcs.get_bcs(&bc_dofs);
    for (I i = 0; i < nbcs; i++) {
      is_bc[bc_dofs[i]] = 1;
    }

    // Apply boundary conditions to the element matrices
    for (I elem = 0; elem < mesh.get_num_elements(); elem++) {
      I dof[Basis::ndof];
      int sign[Basis::ndof];
      if constexpr (Basis::nbasis > 0) {
        get_dof<0>(elem, dof, sign);
      }

      for (I i = 0; i < Basis::ndof; i++) {
        if (is_bc[dof[i]]) {
          if (i < bsize) {
            for (I j = 0; j < bsize; j++) {
              B[elem](i, j) = 0.0;
            }
            for (I j = 0; j < csize; j++) {
              E[elem](i, j) = 0.0;
            }
            B[elem](i, i) = 1.0;
          } else {
            for (I j = 0; j < bsize; j++) {
              F[elem](i - bsize, j) = 0.0;
            }
            for (I j = 0; j < csize; j++) {
              C[elem](i - bsize, j) = 0.0;
            }
            C[elem](i - bsize, i - bsize) = 1.0;
          }
        }
      }
    }

    // Apply boundary conditions to the null vector
    for (I i = 0; i < nbcs; i++) {
      I dof = bc_dofs[i];
      for (I j = 0; j < null_size; j++) {
        bnull(dof / block_size, dof % block_size, j) = 0.0;
      }
    }
  }

  void factor() {
    mat->zero();
    for (I elem = 0; elem < mesh.get_num_elements(); elem++) {
      I dof[Basis::ndof];
      int sign[Basis::ndof];
      if constexpr (Basis::nbasis > 0) {
        get_dof<0>(elem, dof, sign);
      }

      // Form B - E * C^{-1} * F before assembling it into the global matrix
      // Compute and store C^{-1} in the C array of matrices
      CMatType copy(C[elem]);
      Vec<I, csize> ipiv;
      blockInverse<T, csize>(copy, C[elem], ipiv);

      // Mulitply temp = C^{-1} * F
      Mat<T, csize, bsize> temp;
      blockGemm<T, csize, csize, bsize>(C[elem], F[elem], temp);

      // Compute B = B - E * temp = B - E * C^{-1} * F
      blockGemmSub<T, bsize, csize, bsize>(E[elem], temp, B[elem]);

      mat->add_values(bsize, dof, bsize, dof, B[elem]);
    }

    // Allocate the solver - we should add some of these as solver options
    I amg_nlevel = 3;
    double omega = 4.0 / 3.0;
    double epsilon = 0.0;
    bool print_info = false;
    amg = std::make_shared<BSRMatAmgType>(amg_nlevel, omega, epsilon, mat,
                                          bnull, print_info);
  }

  /**
   * @brief Apply the factorized matrix
   *
   * [B  E][x] = [b]
   * [F  C][y]   [d]
   *
   * in = [b]   out = [x]
   *      [d]         [y]
   *
   * Form the right-hand side:
   * f = b - E * C^{-1} * d
   *
   * Apply the AMG preconditioner to the system of equations:
   * B x = f
   *
   * Compute the rest of the solution
   * y = C^{-1} * (d - F * x)
   *
   * @param in The input vector
   * @param out The output vector
   */
  void apply_factor(MultiArrayNew<T *[block_size]> &in,
                    MultiArrayNew<T *[block_size]> &out) {
    // Compute f = b - E * C^{-1} * d
    BLAS::zero(f);
    for (I elem = 0; elem < mesh.get_num_elements(); elem++) {
      I dof[Basis::ndof];
      int sign[Basis::ndof];
      if constexpr (Basis::nbasis > 0) {
        get_dof<0>(elem, dof, sign);
      }

      // Extract b and d
      Vec<T, bsize> b;
      for (I i = 0; i < bsize; i++) {
        b(i) = in(dof[i] / block_size, dof[i] % block_size);
      }

      Vec<T, csize> d;
      for (I i = 0; i < csize; i++) {
        d(i) = in(dof[i + bsize] / block_size, dof[i + bsize] % block_size);
      }

      Vec<T, csize> t;  // t = C^{-1} * d
      blockGemv<T, csize, csize>(C[elem], d, t);

      // b = b - E * t = b - E * C^{-1} * d
      blockGemvSub<T, bsize, csize>(E[elem], t, b);

      // Set the values into the right-hand-size
      for (I i = 0; i < bsize; i++) {
        f(dof[i] / block_size, dof[i] % block_size) += b(i);
      }
    }

    // Apply the preconditioner B * x = f
    BLAS::zero(x);
    amg->applyFactor(f, x);

    // Solve for the remaining dof: y = C^{-1} * (d - F * x)
    for (I elem = 0; elem < mesh.get_num_elements(); elem++) {
      I dof[Basis::ndof];
      int sign[Basis::ndof];
      if constexpr (Basis::nbasis > 0) {
        get_dof<0>(elem, dof, sign);
      }

      // Extract the element solution xelem
      Vec<T, bsize> xelem;
      for (I i = 0; i < bsize; i++) {
        xelem(i) = x(dof[i] / block_size, dof[i] % block_size);
      }

      Vec<T, csize> d;
      for (I i = 0; i < csize; i++) {
        d(i) = in(dof[i + bsize] / block_size, dof[i + bsize] % block_size);
      }

      // Compute d = d - F * xe
      blockGemvSub<T, csize, bsize>(F[elem], xelem, d);

      // Compute y = C^{-1} * (d - F * xe)
      Vec<T, csize> yelem;
      blockGemv<T, csize, csize>(C[elem], d, yelem);

      // Set the variables into the solution
      for (I i = 0; i < bsize; i++) {
        out(dof[i] / block_size, dof[i] % block_size) = xelem(i);
      }

      for (I i = 0; i < csize; i++) {
        out(dof[i + bsize] / block_size, dof[i + bsize] % block_size) =
            yelem(i);
      }
    }
  }

  BMatType &get_bmat(I index) { return B[index]; }
  EMatType &get_emat(I index) { return E[index]; }
  FMatType &get_fmat(I index) { return F[index]; }
  CMatType &get_cmat(I index) { return C[index]; }

 private:
  template <index_t basis>
  void get_dof(index_t elem, index_t dof[], int sign[]) {
    for (index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      sign[i + Basis::template get_dof_offset<basis>()] =
          mesh.template get_global_dof_sign<basis>(elem, i);
      dof[i + Basis::template get_dof_offset<basis>()] =
          mesh.template get_global_dof<basis>(elem, i);
    }
    if constexpr (basis + 1 < Basis::nbasis) {
      get_dof<basis + 1>(elem, dof, sign);
    }
  }

  ElementMesh<Basis> &mesh;
  MultiArrayNew<T *[block_size][null_size]> bnull;
  MultiArrayNew<T *[block_size]> x, f;

  std::vector<BMatType> B;
  std::vector<CMatType> C;
  std::vector<EMatType> E;
  std::vector<FMatType> F;

  std::shared_ptr<BSRMatType> mat;
  std::shared_ptr<BSRMatAmgType> amg;
};

}  // namespace A2D

#endif  //  A2D_FE_STATIC_CONDENSATION_H