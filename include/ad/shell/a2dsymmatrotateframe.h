#ifndef A2D_SYM_MAT_ROTATE_FRAME_H
#define A2D_SYM_MAT_ROTATE_FRAME_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dstack.h"
#include "../a2dtest.h"

namespace A2D {

/*
  Define an expression for C = A^T * B * A
*/

template <typename T, int N, bool symA = false, bool symB = false, bool symC = true, bool additive = false>
A2D_FUNCTION void SymMatMatSquareMult(const T A[], const T B[], T C[]) {  
  // C = A * B
  // zero the matrix if not additive
  if constexpr (!additive) {
    if constexpr (symC) {
      std::fill(C, C + N * (N-1)/2, static_cast<T>(0));
    } else {
      std::fill(C, C + N * N, static_cast<T>(0));
    }
  }

  // Precompute the symmetric matrix indices
  std::array<int, N * N> symMatIndices;
  int index = 0;
  for (int i = 0; i < N; ++i) {
      for (int j = i; j < N; ++j, ++index) {
          symMatIndices[i * N + j] = index;
          symMatIndices[j * N + i] = index;
      }
  }


  // Compute Mat A * Mat B => SymMat C (symmetric part only)
  int inner_start;
  if constexpr (symC) {
    for (int irow = 0; irow < N; irow++) {
      for (int icol = irow; icol < N; icol++) { // only populate lower / upper diag once symMat C
        int inner_start;
        if constexpr (symB) { // had extra halfDiag bool here before, only shows up jin reverse mode stuff
          inner_start = icol;
        } else {
          inner_start = 0;
        }
        // int inner_start = 0;
        int ic = symMatIndices[irow * N + icol];
        for (int inner = inner_start; inner < N; inner++) {
          T aVal;
          if constexpr (symA) {
            int ia = symMatIndices[N * irow + inner];
            aVal = A[ia];
          } else {
            aVal = A[N * irow + inner];
          }
          T bVal;
          if constexpr (symB) {
            int ib = symMatIndices[N * inner + icol];
            bVal = B[ib];
          } else {
            bVal = B[N * inner + icol];
          }

          C[ic] += aVal * bVal;
        }
      }
    } 
  } else { // not symC
    for (int irow = 0; irow < N; irow++) {
      for (int icol = 0; icol < N; icol++) { // only populate lower / upper diag once symMat C
        int inner_start;
        if constexpr (symB) { // had extra halfDiag bool here before, only shows up jin reverse mode stuff
          inner_start = icol;
        } else {
          inner_start = 0;
        }
        // int inner_start = 0;
        for (int inner = inner_start; inner < N; inner++) {
          T aVal;
          if constexpr (symA) {
            int ia = symMatIndices[N * irow + inner];
            aVal = A[ia];
          } else {
            aVal = A[N * irow + inner];
          }
          T bVal;
          if constexpr (symB) {
            int ib = symMatIndices[N * inner + icol];
            // printf("ib = %d\n", ib);
            bVal = B[ib];
          } else {
            bVal = B[N * inner + icol];
          }
          
          C[N * irow + icol] += aVal * bVal;
        }
      }
    }
  }
}

template <typename T, int N, bool symA = false, bool symB = true, bool symC = false, bool additive = false>
A2D_FUNCTION void SymMatMatLeftTrSquareMult(const T A[], const T B[], T C[]) {
  // zero the matrix if not additive
  if constexpr (!additive) {
    if constexpr (symC) {
      std::fill(C, C + N * (N-1)/2, static_cast<T>(0));
    } else {
      std::fill(C, C + N * N, static_cast<T>(0));
    }
  }

  // Precompute the symmetric matrix indices
  std::array<int, N * N> symMatIndices;
  int index = 0;
  for (int i = 0; i < N; ++i) {
      for (int j = i; j < N; ++j, ++index) {
          symMatIndices[i * N + j] = index;
          symMatIndices[j * N + i] = index;
      }
  }

  // Compute Mat A^T * SymMat B => C (most of the time)
  // symC is always false so not considering it (just makes it easier to type same # inputs for template)
  for (int inner = 0; inner < N; ++inner) {
    for (int icol = 0; icol < N; ++icol) {
          T bValue;
          if constexpr (symB) {
            // Use precomputed index
            int ib = symMatIndices[inner * N + icol];
            bValue = B[ib]; // Access the value of symmetric matrix A once
          } else {
            bValue = B[N * inner + icol];
          }
          
          // Use a temporary variable to accumulate values for C
          for (int irow = 0; irow < N; ++irow) {
            T aValue;
            if constexpr (symA) {
              int ia = symMatIndices[N * inner + irow];
              aValue = A[ia];
            } else {
              aValue = A[N * inner + irow];
            }
            // printf("W[%d] += S[%d] * T[%d]\n", N * irow + icol, N * inner + irow, ib);
            C[N * irow + icol] += aValue * bValue;
          }
      }
  }
}

template <typename T, int N, bool symA = false, bool symB = false, bool symC = false, bool additive = false>
A2D_FUNCTION void SymMatMatRightTrSquareMult(const T A[], const T B[], T C[]) {
  // zero the matrix if not additive
  if constexpr (!additive) {
    if constexpr (symC) {
      std::fill(C, C + N * (N-1)/2, static_cast<T>(0));
    } else {
      std::fill(C, C + N * N, static_cast<T>(0));
    }
  }

  // Precompute the symmetric matrix indices
  std::array<int, N * N> symMatIndices;
  int index = 0;
  for (int i = 0; i < N; ++i) {
      for (int j = i; j < N; ++j, ++index) {
          symMatIndices[i * N + j] = index;
          symMatIndices[j * N + i] = index;
      }
  }

  // Compute Mat A * SymMat B^T => C
  // symA is not really used here..
  // if C output is symmetric
  if constexpr (symC) {
    for (int inner = 0; inner < N; ++inner) {
      int col_start;
      if constexpr (symB) { // had extra halfDiag bool here before, only shows up jin reverse mode stuff
        col_start = inner;
      } else {
        col_start = 0;
      }
      // int col_start = 0;
      for (int icol = col_start; icol < N; ++icol) {
            // Use precomputed index (because transpose)
            T bValue;
            if constexpr (symB) {
              int ib = symMatIndices[icol * N + inner];
              bValue = B[ib];
            } else {
              bValue = B[N * icol + inner];
            }
            
            // Use a temporary variable to accumulate values for C
            for (int irow = 0; irow < N; ++irow) {
              int ic = symMatIndices[N * irow + icol];
              C[ic] += A[N * irow + inner] * bValue;
            }
        }
    }
  } else  { // C output is not symmetric
    for (int inner = 0; inner < N; ++inner) {
      int col_start;
      if constexpr (symB) { // had extra halfDiag bool here before, only shows up jin reverse mode stuff
        col_start = inner;
      } else {
        col_start = 0;
      }
      // int col_start = 0;
      for (int icol = col_start; icol < N; ++icol) {
            // Use precomputed index (because transpose)
            T bValue;
            if constexpr (symB) {
              int ib = symMatIndices[icol * N + inner];
              bValue = B[ib];
            } else {
              bValue = B[N * icol + inner];
            }

            for (int irow = 0; irow < N; ++irow) {
              C[N * irow + icol] += A[N * irow + inner] * bValue;
            }
        }
    }
  }
}

template <typename T, int N>
A2D_FUNCTION void SymMatRotateFrame(const Mat<T,N,N> &A, const SymMat<T,N> &B, SymMat<T,N> &C) {
  Mat<T, N, N> Ctemp;
  // Ctemp = A^T * B
  SymMatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
  // C = Ctemp * A
  SymMatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
}

template <class Atype, class Btype, class Ctype>
class SymMatRotateFrameExpr {
 public:

  // Extract the numeric type to use
  typedef typename get_object_numeric_type<Ctype>::type T;

  // Extract the dimensions of the matrices
  // if (get_diff_type::)
  // how to get matrix rows for symMat?
  // optional SymMat or Mat here
  // const bool A_issym = get_a2d_object_type<Atype>::value == ADObjType::SYMMAT
  static constexpr int N = get_matrix_rows<Atype>::size;
  static constexpr int M = get_matrix_columns<Atype>::size;
  static constexpr int K = get_symmatrix_size<Btype>::size;
  static constexpr int P = get_symmatrix_size<Ctype>::size;

  // check all square matrices
  static_assert(
    (N == M) && (M == K) && (K == P),
    "all matrices in MatRotateFrameExpr must be same N x N square matrix size."
  );

  // Get the types of the matrices
  static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;
  static constexpr ADiffType adB = get_diff_type<Btype>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<Ctype>::order;

  A2D_FUNCTION SymMatRotateFrameExpr(Atype& A, Btype& B, Ctype& C)
      : A(A), B(B), C(C) {}

  A2D_FUNCTION void eval() {
    Mat<T, N, N> Ctemp;
    // Ctemp = A^T * B
    SymMatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
    // C = Ctemp * A
    SymMatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
  }

  A2D_FUNCTION void bzero() { C.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    // full expression of forward pass:
    //   Cdot = Adot^T * B * A + A^T * B * Adot + A^T * Bdot * A
    // printf("forward\n");

    if constexpr (adA == ADiffType::ACTIVE and adB == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Adot term1
      SymMatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
      // Adot term2
      SymMatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N,true>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));
      // Bdot term
      SymMatMatLeftTrSquareMult<T,N>(get_data(A), GetSeed<seed>::get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N,true>(get_data(Ctemp), get_data(A), get_data(C));

    } else if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Adot term1
      SymMatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
      // Adot term2
      SymMatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N,true>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));

    } else if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Bdot term
      SymMatMatLeftTrSquareMult<T,N>(get_data(A), GetSeed<seed>::get_data(B), get_data(Ctemp));
      SymMatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
    }
  }

  A2D_FUNCTION void reverse() {
    // printf("reverse\n");
    if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expression: Abar += B^T * A * Cbar + B * A * Cbar^T
      // first term B^T * A * Cbar
      SymMatMatLeftTrSquareMult<T,N,true,false,false,false>(get_data(B), get_data(A), get_data(temp));
      SymMatMatSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::b>::get_data(A));

      // second term B * A * Cbar^T added in
      SymMatMatSquareMult<T,N,true,false,false,false>(get_data(B), get_data(A), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::b>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expresion Bbar += A * Cbar * A^T
      SymMatMatSquareMult<T,N,false,true,false,false>(get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,false,true,true>(get_data(temp), get_data(A), GetSeed<ADseed::b>::get_data(B));
    }
  }

  A2D_FUNCTION void hzero() { C.hzero(); }

  A2D_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    // HJP backpropagation based on Aaron's paper and my ppt
    // 
    // Ahat += B^T * A * Chat + B * A * Chat^T +
    //         Bdot^T * A * Cbar + Bdot * A * Cbar^T +
    //         B^T * Adot * Cbar + B * Adot * Cbar^T
    // 
    // Bhat += A * Chat * A^T +
    //         Adot * Cbar * A^T + A * Cbar * Adot^T

    if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // term1 for Ahat : B^T * A * Chat
      SymMatMatLeftTrSquareMult<T,N,true,false,false,false>(get_data(B), get_data(A), get_data(temp));
      SymMatMatSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term2 for Ahat : B * A * Chat^T
      SymMatMatSquareMult<T,N,true,false,false,false>(get_data(B), get_data(A), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term 5 for Ahat : B^T * Adot * Cbar
      SymMatMatLeftTrSquareMult<T,N,true,false,false,false>(get_data(B), GetSeed<ADseed::p>::get_data(A), get_data(temp));
      SymMatMatSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term 6 for Ahat : B * Adot * Cbar^T
      SymMatMatSquareMult<T,N,true,false,false,false>(get_data(B), GetSeed<ADseed::p>::get_data(A), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));
    }

    if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // term 1 for Bhat : A * Chat * A^T
      SymMatMatSquareMult<T,N,false,true,false,false>(get_data(A), GetSeed<ADseed::h>::get_data(C), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,false,true,true>(get_data(temp), get_data(A), GetSeed<ADseed::h>::get_data(B));
    }

    if constexpr (adA == ADiffType::ACTIVE && adB == ADiffType::ACTIVE) {
      // now only remaining terms how up
      Mat<T, N, N> temp;

      // term3 for Ahat : Bdot^T * A * Cbar
      SymMatMatLeftTrSquareMult<T,N,true,false,false,false>(GetSeed<ADseed::p>::get_data(B), get_data(A), get_data(temp));
      SymMatMatSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term4 for Ahat : Bdot * A * Cbar^T
      SymMatMatSquareMult<T,N,true,false,false,false>(GetSeed<ADseed::p>::get_data(B), get_data(A), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,true,false,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term2 for Bhat : Adot * Cbar * A^T
      SymMatMatSquareMult<T,N,false,true,false,false>(GetSeed<ADseed::p>::get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,false,true,true>(get_data(temp), get_data(A), GetSeed<ADseed::h>::get_data(B));

      // term3 for Bhat : A * Cbar * Adot^T
      SymMatMatSquareMult<T,N,false,true,false,false>(get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      SymMatMatRightTrSquareMult<T,N,false,false,true,true>(get_data(temp), GetSeed<ADseed::p>::get_data(A), GetSeed<ADseed::h>::get_data(B));
    }
  }

 private:
  Atype& A;
  Btype& B;
  Ctype& C;
};

// all implementations
template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(ADObj<Atype>& A, ADObj<Btype>& B, ADObj<Ctype>& C) {
  return SymMatRotateFrameExpr<ADObj<Atype>, ADObj<Btype>, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(ADObj<Atype>& A, Btype &B, ADObj<Ctype>& C) {
  return SymMatRotateFrameExpr<ADObj<Atype>, Btype, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(Atype &A, ADObj<Btype>& B, ADObj<Ctype>& C) {
  return SymMatRotateFrameExpr<Atype, ADObj<Btype>, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(A2DObj<Atype>& A, A2DObj<Btype>& B, A2DObj<Ctype>& C) {
  return SymMatRotateFrameExpr<A2DObj<Atype>, A2DObj<Btype>, A2DObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(A2DObj<Atype>& A, Btype &B, A2DObj<Ctype>& C) {
  return SymMatRotateFrameExpr<A2DObj<Atype>, Btype, A2DObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto SymMatRotateFrame(Atype& A, A2DObj<Btype>& B, A2DObj<Ctype>& C) {
  return SymMatRotateFrameExpr<Atype, A2DObj<Btype>, A2DObj<Ctype>>(A, B, C);
}

namespace Test {

template <typename T, int N>
class SymMatRotateFrameTest
    : public A2DTest<T, SymMat<T,N>, Mat<T, N, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>, SymMat<T, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatRotateFrame<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> A;
    SymMat<T,N> B, C;

    x.get_values(A, B);
    SymMatRotateFrame(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) override {
    ADObj<Mat<T, N, N>> A;
    ADObj<SymMat<T,N>> B, C;

    x.get_values(A.value(), B.value());
    auto stack = MakeStack(SymMatRotateFrame(A, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) override {
    A2DObj<Mat<T, N, N>> A;
    A2DObj<SymMat<T,N>> B, C;

    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto stack = MakeStack(SymMatRotateFrame(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

bool SymMatRotateFrameTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  // SymMatRotateFrameTest<Tc, 1> test0;
  // passed = passed && Run(test0, component, write_output);
  SymMatRotateFrameTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatRotateFrameTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  SymMatRotateFrameTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYM_MAT_ROTATE_FRAME_H