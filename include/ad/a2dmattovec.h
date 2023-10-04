#ifndef A2D_MAT_TO_VEC_H
#define A2D_MAT_TO_VEC_H

#include <type_traits>

#include "a2ddefs.h"
#include "a2dmat.h"

namespace A2D {

template <typename T, typename I, int N, int M>
KOKKOS_FUNCTION void MatColumnToVec(I column, Mat<T, N, M> &A, Vec<T, N> &x) {
  for (int i = 0; i < N; i++) {
    x(i) = A(i, column);
  }
}

template <typename T, typename I, int N>
KOKKOS_FUNCTION void MatColumnToVec(I column, SymMat<T, N> &A, Vec<T, N> &x) {
  for (int i = 0; i < N; i++) {
    x(i) = A(i, column);
  }
}

template <typename I, class Atype, class xtype>
class MatColumnToVecExpr {
 public:
  // Get the dimension of the matrix
  static constexpr int N = conditional_value<
      int, get_a2d_object_type<Atype>::obj_type == ADObjType::SYMMAT,
      get_symmatrix_size<Atype>::size, get_matrix_rows<Atype>::size>::value;
  static_assert(get_vec_size<xtype>::size == N,
                "Matrix and vector dimensions must agree");

  MatColumnToVecExpr(I column, Atype &A, xtype &x)
      : column(column), A(A), x(x) {}

  KOKKOS_FUNCTION void eval() {
    for (int i = 0; i < N; i++) {
      x(i) = A(i, column);
    }
  }

  KOKKOS_FUNCTION void bzero() { x.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    auto xp = GetSeed<seed>::get_obj(x);
    auto Ap = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      xp(i) = Ap(i, column);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    auto xb = GetSeed<seed>::get_obj(x);
    auto Ab = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      Ab(i, column) += xb(i);
    }
  }

  KOKKOS_FUNCTION void hzero() { x.hzero(); }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    auto xh = GetSeed<seed>::get_obj(x);
    auto Ah = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      Ah(i, column) += xh(i);
    }
  }

  const I column;
  Atype &A;
  xtype &x;
};

template <typename I, class Atype, class xtype>
KOKKOS_FUNCTION auto MatColumnToVec(I column, ADObj<Atype> &A,
                                    ADObj<xtype> &x) {
  return MatColumnToVecExpr<I, ADObj<Atype>, ADObj<xtype>>(column, A, x);
}

template <typename I, class Atype, class xtype>
KOKKOS_FUNCTION auto MatColumnToVec(I column, A2DObj<Atype> &A,
                                    A2DObj<xtype> &x) {
  return MatColumnToVecExpr<I, A2DObj<Atype>, A2DObj<xtype>>(column, A, x);
}

template <typename T, typename I, int N, int M>
KOKKOS_FUNCTION void MatRowToVec(I row, Mat<T, N, M> &A, Vec<T, M> &x) {
  for (int i = 0; i < M; i++) {
    x(i) = A(row, i);
  }
}

template <typename T, typename I, int N>
KOKKOS_FUNCTION void MatRowToVec(I row, SymMat<T, N> &A, Vec<T, N> &x) {
  for (int i = 0; i < N; i++) {
    x(i) = A(row, i);
  }
}

template <typename I, class Atype, class xtype>
class MatRowToVecExpr {
 public:
  // Get the dimension of the matrix
  static constexpr int N = conditional_value<
      int, get_a2d_object_type<Atype>::obj_type == ADObjType::SYMMAT,
      get_symmatrix_size<Atype>::size, get_matrix_columns<Atype>::size>::value;
  static_assert(get_vec_size<xtype>::size == N,
                "Matrix and vector dimensions must agree");

  MatRowToVecExpr(I row, Atype &A, xtype &x) : row(row), A(A), x(x) {}

  KOKKOS_FUNCTION void eval() {
    for (int i = 0; i < N; i++) {
      x(i) = A(row, i);
    }
  }

  KOKKOS_FUNCTION void bzero() { x.bzero(); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    auto xp = GetSeed<seed>::get_obj(x);
    auto Ap = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      xp(i) = Ap(row, i);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    auto xb = GetSeed<seed>::get_obj(x);
    auto Ab = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      Ab(row, i) += xb(i);
    }
  }

  KOKKOS_FUNCTION void hzero() { x.hzero(); }

  KOKKOS_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    auto xh = GetSeed<seed>::get_obj(x);
    auto Ah = GetSeed<seed>::get_obj(A);
    for (int i = 0; i < N; i++) {
      Ah(row, i) += xh(i);
    }
  }

  const I row;
  Atype &A;
  xtype &x;
};

template <typename I, class Atype, class xtype>
KOKKOS_FUNCTION auto MatRowToVec(I row, ADObj<Atype> &A, ADObj<xtype> &x) {
  return MatRowToVecExpr<I, ADObj<Atype>, ADObj<xtype>>(row, A, x);
}

template <typename I, class Atype, class xtype>
KOKKOS_FUNCTION auto MatRowToVec(I row, A2DObj<Atype> &A, A2DObj<xtype> &x) {
  return MatRowToVecExpr<I, A2DObj<Atype>, A2DObj<xtype>>(row, A, x);
}

}  // namespace A2D

#endif  // A2D_MAT_TO_VEC_H