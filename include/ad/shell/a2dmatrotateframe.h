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
A2D_FUNCTION MatMatSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::NORMAL,MatOp::NORMAL, additive>(A,B,C);
}

template <typename T, int N, bool additive = false>
A2D_FUNCTION MatMatLeftTrSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::TRANSPOSE,MatOp::NORMAL, additive>(A,B,C);
}

template <typename T, int N, bool additive = false>
A2D_FUNCTION MatMatRightTrSquareMult(const T A[], const T B[], T C[]) {
  MatMatMultCore<T,N,N,N,N,N,N,MatOp::NORMAL,MatOp::TRANSPOSE, additive>(A,B,C);
}

template <class Atype, class Btype, class Ctype>
class MatRotateFrameExpr {
 public:

  // Extract the numeric type to use
  typedef typename get_object_numeric_type<Ctype>::type T;

  // Extract the dimensions of the matrices
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
    MatMatLeftTrSquareMult<T,M>(get_data(A), get_data(B), get_data(Ctemp));
    // C = Ctemp * A
    MatMatSquareMult<T,M>(get_data(Ctemp), get_data(A), get_data(C));
  }

  A2D_FUNCTION void bzero() { C.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE and adB == ADiffType:ACTIVE) {
      constexpr bool additive = true;
      Mat<T, N, N> Ctemp;

      // Afwd^T * B * Afwd pass to Cfwd
      MatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));
      // A^T * Bfwd * A added into Cfwd
      MatMatLeftTrSquareMult<T,N>(get_data(A), GetSeed<seed>::get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N,true>(get_data(Ctemp), get_data(A), get_data(C));

    } else if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // Afwd^T * B * Afwd pass to Cfwd
      MatMatLeftTrSquareMult<T,N>(GetSeed<seed>::get_data(A), get_data(B), get_data(Ctemp));
      MatMatSquareMult<T,N>(get_data(Ctemp), GetSeed<seed>::get_data(A), get_data(C));
    } else if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> Ctemp;
      // A^T * Bfwd * A added into Cfwd
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
      MatMatSquareMult<T,N>(get_data(A), GetSeed<ADseed::b>::get_data(C), temp);
      MatMatRightTrSquareMult<T,N,true>(temp, get_data(A), GetSeed<ADseed::b>::get_data(B));
    }
  }

  A2D_FUNCTION void hzero() { C.hzero(); }

  A2D_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    // also need an (adA active and adB active) check?

    // figure out if Hessian second order steps are correct => e.g. in C = AB
    // why is not A2bar += C2bar * B^T * B^T and B2bar += A^2 * C2bar (but only have single products in a2dgemm.h)
    if constexpr (adA == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expression: Abar += B^T * A * Cbar + B * A * Cbar^T
      // first term B^T * A * Cbar
      MatMatLeftTrSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));

      // second term B * A * Cbar^T added in
      MatMatSquareMult<T,N>(get_data(B), get_data(A), get_data(temp));
      MatMatRightTrSquareMult<T,N,true>(get_data(temp), GetSeed<ADseed::h>::get_data(C), GetSeed<ADseed::h>::get_data(A));
    }
    if constexpr (adB == ADiffType::ACTIVE) {
      Mat<T, N, N> temp;
      // full expresion Bbar += A * Cbar * A^T
      MatMatSquareMult<T,N>(get_data(A), GetSeed<ADseed::h>::get_data(C), temp);
      MatMatRightTrSquareMult<T,N,true>(temp, get_data(A), GetSeed<ADseed::h>::get_data(B));
    }
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
A2D_FUNCTION auto MatRotateFrame(A2DObj<Atype>& A, A2DObj<Btype>& B, A2DObj<Ctype>& C) {
  return MatRotateFrameExpr<A2DObj<Atype>, A2DObj<Btype>, A2DObj<Ctype>>(A, B, C);
}

// namespace Test {

// template <MatOp opA, MatOp opB, typename T, int N, int M, int K, int L, int P,
//           int Q>
// class MatMatMultTest
//     : public A2DTest<T, Mat<T, P, Q>, Mat<T, N, M>, Mat<T, K, L>> {
//  public:
//   using Input = VarTuple<T, Mat<T, N, M>, Mat<T, K, L>>;
//   using Output = VarTuple<T, Mat<T, P, Q>>;

//   // Assemble a string to describe the test
//   std::string name() {
//     std::stringstream s;
//     s << "MatMatMult<";
//     if (opA == MatOp::NORMAL) {
//       s << "N,";
//     } else {
//       s << "T,";
//     }
//     if (opB == MatOp::NORMAL) {
//       s << "N,";
//     } else {
//       s << "T,";
//     }
//     s << N << "," << M << "," << K << "," << L << "," << P << "," << Q << ">";

//     return s.str();
//   }

//   // Evaluate the matrix-matrix product
//   Output eval(const Input& x) {
//     Mat<T, N, M> A;
//     Mat<T, K, L> B;
//     Mat<T, P, Q> C;

//     x.get_values(A, B);
//     MatMatMult<opA, opB>(A, B, C);
//     return MakeVarTuple<T>(C);
//   }

//   // Compute the derivative
//   void deriv(const Output& seed, const Input& x, Input& g) {
//     ADObj<Mat<T, N, M>> A;
//     ADObj<Mat<T, K, L>> B;
//     ADObj<Mat<T, P, Q>> C;

//     x.get_values(A.value(), B.value());
//     auto stack = MakeStack(MatMatMult<opA, opB>(A, B, C));
//     seed.get_values(C.bvalue());
//     stack.reverse();
//     g.set_values(A.bvalue(), B.bvalue());
//   }

//   // Compute the second-derivative
//   void hprod(const Output& seed, const Output& hval, const Input& x,
//              const Input& p, Input& h) {
//     A2DObj<Mat<T, N, M>> A;
//     A2DObj<Mat<T, K, L>> B;
//     A2DObj<Mat<T, P, Q>> C;

//     x.get_values(A.value(), B.value());
//     p.get_values(A.pvalue(), B.pvalue());
//     auto stack = MakeStack(MatMatMult<opA, opB>(A, B, C));
//     seed.get_values(C.bvalue());
//     hval.get_values(C.hvalue());
//     stack.hproduct();
//     h.set_values(A.hvalue(), B.hvalue());
//   }
// };

// template <typename T, int N, int M, int K>
// bool MatMatMultTestHelper(bool component = false, bool write_output = true) {
//   const MatOp NORMAL = MatOp::NORMAL;
//   const MatOp TRANSPOSE = MatOp::TRANSPOSE;
//   using Tc = std::complex<T>;

//   bool passed = true;
//   MatMatMultTest<NORMAL, NORMAL, Tc, N, M, M, K, N, K> test1;
//   passed = passed && Run(test1, component, write_output);
//   MatMatMultTest<NORMAL, TRANSPOSE, Tc, N, M, K, M, N, K> test2;
//   passed = passed && Run(test2, component, write_output);
//   MatMatMultTest<TRANSPOSE, NORMAL, Tc, N, M, N, K, M, K> test3;
//   passed = passed && Run(test3, component, write_output);
//   MatMatMultTest<TRANSPOSE, TRANSPOSE, Tc, N, M, K, N, M, K> test4;
//   passed = passed && Run(test4, component, write_output);

//   return passed;
// }

// bool MatMatMultTestAll(bool component = false, bool write_output = true) {
//   bool passed = true;
//   passed =
//       passed && MatMatMultTestHelper<double, 3, 3, 3>(component, write_output);
//   passed =
//       passed && MatMatMultTestHelper<double, 2, 3, 4>(component, write_output);
//   passed =
//       passed && MatMatMultTestHelper<double, 5, 4, 2>(component, write_output);

//   return passed;
// }

// }  // namespace Test

}  // namespace A2D

#endif  // A2d_MAT_ROTATE_FRAME_H