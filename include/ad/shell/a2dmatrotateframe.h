#ifndef A2d_MAT_ROTATE_FRAME_H
#define A2d_MAT_ROTATE_FRAME_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dstack.h"
#include "../a2dtest.h"
#include "../core/a2dgemmcore.h"

namespace A2D {

/*
  Define an expression for C = A^T * B * A
*/

template <typename T, int N, bool additive = false>
A2D_FUNCTION void MatMatSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::NORMAL,MatOp::NORMAL, additive>(A,B,C);
}

template <typename T, int N, bool additive = false>
A2D_FUNCTION void MatMatLeftTrSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::TRANSPOSE,MatOp::NORMAL, additive>(A,B,C);
}

template <typename T, int N, bool additive = false>
A2D_FUNCTION void MatMatRightTrSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::NORMAL,MatOp::TRANSPOSE, additive>(A,B,C);
}

template <typename T, int N>
A2D_FUNCTION void MatRotateFrame(const Mat<T,N,N> &A, const Mat<T,N,N> &B, Mat<T,N,N> &C) {
  Mat<T, N, N> Ctemp;
  // Ctemp = A^T * B
  MatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
  // C = Ctemp * A
  MatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
}

template <class Atype, class Btype, class Ctype>
class MatRotateFrameExpr {
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
  static constexpr int K = get_matrix_rows<Btype>::size;
  static constexpr int L = get_matrix_columns<Btype>::size;
  static constexpr int P = get_matrix_rows<Ctype>::size;
  static constexpr int Q = get_matrix_columns<Ctype>::size;

  // check all square matrices
  static_assert(
    (N == M) && (M == K) && (K == L) && (L == P) && (P == Q),
    "all matrices in MatRotateFrameExpr must be same N x N square matrix size."
  );

  // Get the types of the matrices
  static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;
  static constexpr ADiffType adB = get_diff_type<Btype>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<Ctype>::order;

  A2D_FUNCTION MatRotateFrameExpr(Atype& A, Btype& B, Ctype& C)
      : A(A), B(B), C(C) {}

  A2D_FUNCTION void eval() {
    Mat<T, N, N> Ctemp;
    // Ctemp = A^T * B
    MatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
    // C = Ctemp * A
    MatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
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

    if constexpr (adA == ADiffType::ACTIVE and adB == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Adot term1
      MatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
      // Adot term2
      MatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N,true>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));
      // Bdot term
      MatMatLeftTrSquareMult<T,N>(get_data(A), GetSeed<seed>::get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N,true>(get_data(Ctemp), get_data(A), get_data(C));

    } else if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Adot term1
      MatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
      // Adot term2
      MatMatLeftTrSquareMult<T,N>(get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N,true>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));

    } else if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Bdot term
      MatMatLeftTrSquareMult<T,N>(get_data(A), GetSeed<seed>::get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N>(get_data(Ctemp), get_data(A), get_data(C));
    }
  }

  A2D_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expression: Abar += B^T * A * Cbar + B * A * Cbar^T
      // first term B^T * A * Cbar
      MatMatLeftTrSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::b>::get_data(A));

      // second term B * A * Cbar^T added in
      MatMatSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::b>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expresion Bbar += A * Cbar * A^T
      MatMatSquareMult<T,N>(get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), get_data(A), GetSeed<ADseed::b>::get_data(B));
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
      MatMatLeftTrSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term2 for Ahat : B * A * Chat^T
      MatMatSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term 5 for Ahat : B^T * Adot * Cbar
      MatMatLeftTrSquareMult<T,N>(get_data(B), GetSeed<ADseed::p>::get_data(A), get_data(temp));
      MatMatSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term 6 for Ahat : B * Adot * Cbar^T
      MatMatSquareMult<T,N>(get_data(B), GetSeed<ADseed::p>::get_data(A), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));
    }

    if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // term 1 for Bhat : A * Chat * A^T
      MatMatSquareMult<T,N>(get_data(A), GetSeed<ADseed::h>::get_data(C), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), get_data(A), GetSeed<ADseed::h>::get_data(B));
    }

    if constexpr (adA == ADiffType::ACTIVE && adB == ADiffType::ACTIVE) {
      // now only remaining terms how up
      Mat<T, N, N> temp;

      // term3 for Ahat : Bdot^T * A * Cbar
      MatMatLeftTrSquareMult<T,N>(GetSeed<ADseed::p>::get_data(B), get_data(A), get_data(temp));
      MatMatSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term4 for Ahat : Bdot * A * Cbar^T
      MatMatSquareMult<T,N>(GetSeed<ADseed::p>::get_data(B), get_data(A), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::b>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // term2 for Bhat : Adot * Cbar * A^T
      MatMatSquareMult<T,N>(GetSeed<ADseed::p>::get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), get_data(A), GetSeed<ADseed::h>::get_data(B));

      // term3 for Bhat : A * Cbar * Adot^T
      MatMatSquareMult<T,N>(get_data(A), GetSeed<ADseed::b>::get_data(C), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::p>::get_data(A), GetSeed<ADseed::h>::get_data(B));
    }
  }

 private:
  Atype& A;
  Btype& B;
  Ctype& C;
};

// all implementations
template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(ADObj<Atype>& A, ADObj<Btype>& B, ADObj<Ctype>& C) {
  return MatRotateFrameExpr<ADObj<Atype>, ADObj<Btype>, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(ADObj<Atype>& A, Btype &B, ADObj<Ctype>& C) {
  return MatRotateFrameExpr<ADObj<Atype>, Btype, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(Atype &A, ADObj<Btype>& B, ADObj<Ctype>& C) {
  return MatRotateFrameExpr<ADObj<Atype>, Btype, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(A2DObj<Atype>& A, A2DObj<Btype>& B, A2DObj<Ctype>& C) {
  return MatRotateFrameExpr<A2DObj<Atype>, A2DObj<Btype>, A2DObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(A2DObj<Atype>& A, Btype &B, A2DObj<Ctype>& C) {
  return MatRotateFrameExpr<A2DObj<Atype>, Btype, A2DObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatRotateFrame(Atype &A, A2DObj<Btype>& B, A2DObj<Ctype>& C) {
  return MatRotateFrameExpr<A2DObj<Atype>, Btype, A2DObj<Ctype>>(A, B, C);
}

namespace Test {

template <typename T, int N>
class MatRotateFrameTest
    : public A2DTest<T, Mat<T, N, N>, Mat<T, N, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>, Mat<T, N, N>>;
  using Output = VarTuple<T, Mat<T, N, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatRotateFrame<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> A, B, C;

    x.get_values(A, B);
    MatRotateFrame(A, B, C);
    return MakeVarTuple<T>(C);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, N>> A, B, C;

    x.get_values(A.value(), B.value());
    auto stack = MakeStack(MatRotateFrame(A, B, C));
    seed.get_values(C.bvalue());
    stack.reverse();
    g.set_values(A.bvalue(), B.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, N>> A, B, C;

    x.get_values(A.value(), B.value());
    p.get_values(A.pvalue(), B.pvalue());
    auto stack = MakeStack(MatRotateFrame(A, B, C));
    seed.get_values(C.bvalue());
    hval.get_values(C.hvalue());
    stack.hproduct();
    h.set_values(A.hvalue(), B.hvalue());
  }
};

bool MatRotateFrameTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatRotateFrameTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatRotateFrameTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  MatRotateFrameTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2d_MAT_ROTATE_FRAME_H