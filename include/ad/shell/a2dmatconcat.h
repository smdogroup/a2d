#ifndef A2D_MAT_CONCAT_H
#define A2D_MAT_CONCAT_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dvec.h"
#include "../a2dstack.h"
#include "../a2dtest.h"

namespace A2D {

template <typename T, int M, int N, int K>
A2D_FUNCTION void MatColConcatCore(const T A[], const T B[], T C[]) {
    // concat A (M x N) with B (M x K) to get (M x (N+K)) matrix C
    for (int i = 0; i < M; i++) {
        // copy A portion
        for (int j = 0; j < N; j++) {
            C[i,j] = A[i,j];
        }
        // copy B portion
        for (int k = 0; k < K; k++) {
            C[i,N+k] = B[i,k];
        }
    }
}

template <typename T, int M, int N, int K>
A2D_FUNCTION void MatColConcatReverseCore(const T Cb, T Ab[], T Bb[]) {
    // concat A (M x N) with B (M x K) to get (M x (N+K)) matrix C
    for (int i = 0; i < M; i++) {
        // backprop to A
        for (int j = 0; j < N; j++) {
            Ab[i,j] += Cb[i,j];
        }
        // backprop to B
        for (int k = 0; k < K; k++) {
            Bb[i,k] += Cb[i,N+k];
        }
    }
}


template <class Atype, class Btype, class Ctype>
class MatColConcatExpr {
    public:
      // extract numeric type to use
      typedef typename get_object_numeric_type<Ctype>::type T;

      // Get the sizes of the matrices, vectors
      static constexpr int M = get_matrix_rows<Atype>::size;
      static constexpr int N = get_matrix_columns<Atype>::size;
      static constexpr int L = get_matrix_rows<Btype>::size;
      static constexpr int K = get_matrix_columns<Btype>::size;
      static constexpr int S = get_matrix_rows<Ctype>::size;
      static constexpr int Q = get_matrix_columns<Ctype>::size;

      // Get the types of the matrices and scalars
      static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;
      static constexpr ADiffType adB = get_diff_type<Btype>::diff_type;
      static constexpr ADiffType adC = get_diff_type<Ctype>::diff_type;

      // assert correct sizes
      static_assert((M == L) && (L == S) && ((N + K) == Q),
                    "MatColConcat requires [A (MxN), B (MxK)] => C (Mx(N+K)) dimensions");

      // get the differentiation order
      static constexpr ADorder order = get_diff_order<Ctype>::order; 

      A2D_FUNCTION MatColConcatExpr(Atype& A, Btype& B, Ctype& C) : A(A), B(B), C(C) {}

      A2D_FUNCTION void eval() {
        MatColConcatCore<T,M,N,K>(get_data(A), get_data(B), get_data(C));
      }

      A2D_FUNCTION void bzero() { C.bzero(); }

      template <ADorder forder>
      A2D_FUNCTION void forward() {
        constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
        
        if constexpr (adA == ADiffType::ACTIVE && adB == ADiffType::ACTIVE) {
            MatColConcatCore<T,M,N,K>(GetSeed<seed>::get_data(A),
                          GetSeed<seed>::get_data(B),
                          GetSeed<seed>::get_data(C));
        } else if constexpr (adA == ADiffType::ACTIVE) {
            Mat<T,M,K> B_void;
            MatColConcatCore<T,M,N,K>(
                GetSeed<seed>::get_data(A),
                get_data(B_void),
                GetSeed<seed>::get_data(C)
            );
        } else if constexpr (adB == ADiffType::ACTIVE) {
            Mat<T,M,N> A_void;
            MatColConcatCore<T,M,N,K>(
                get_data(A_void),
                GetSeed<seed>::get_data(B),
                GetSeed<seed>::get_data(C)
            );
        }
      }

      A2D_FUNCTION void reverse() {
        constexpr ADseed seed = ADseed::b;
        if constexpr (adA == ADiffType::ACTIVE) {
            Mat<T,M,K> B_void;
            MatColConcatReverseCore<T,M,N,K>(
                GetSeed<seed>::get_data(C),
                GetSeed<seed>::get_data(A),
                get_data(B_void)
            );
        }
        if constexpr (adB == ADiffType::ACTIVE) {
            Mat<T,M,K> A_void;
            MatColConcatReverseCore<T,M,N,K>(
                GetSeed<seed>::get_data(C),
                GetSeed<seed>::get_data(A_void),
                get_data(B)
            );
        }
      }

      A2D_FUNCTION void hzero() { C.hzero(); }

      A2D_FUNCTION void hreverse() {
        constexpr ADseed seed = ADseed::h;
        if constexpr (adA == ADiffType::ACTIVE) {
            Mat<T,M,K> B_void;
            MatColConcatReverseCore<T,M,N,K>(
                GetSeed<seed>::get_data(C),
                GetSeed<seed>::get_data(A),
                get_data(B_void)
            );
        }
        if constexpr (adB == ADiffType::ACTIVE) {
            Mat<T,M,K> A_void;
            MatColConcatReverseCore<T,M,N,K>(
                GetSeed<seed>::get_data(C),
                GetSeed<seed>::get_data(A_void),
                get_data(B)
            );
        }
      }

      Atype &A;
      Btype &B;
      Ctype &C;
}; // end of MatColConcatExpr class definition

// Full active variants
template <class Atype, class Btype, class Ctype> // should be const args?
A2D_FUNCTION auto MatColConcat(Atype &A, Btype &B, Ctype &C) {
    return MatColConcatExpr<Atype, Btype, Ctype>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatColConcat(ADObj<Atype> &A, ADObj<Btype> &B, ADObj<Ctype> &C) {
    return MatColConcatExpr<ADObj<Atype>, ADObj<Btype>, ADObj<Ctype>>(A, B, C);
}

template <class Atype, class Btype, class Ctype>
A2D_FUNCTION auto MatColConcat(A2DObj<Atype> &A, A2DObj<Btype> &B, A2DObj<Ctype> &C) {
    return MatColConcatExpr<A2DObj<Atype>, A2DObj<Btype>, A2DObj<Ctype>>(A, B, C);
}


// TODO : add testing here?
// TODO : add MatRowConcat?


} // end of A2D namespace
#endif // A2D_MAT_CONCAT_H