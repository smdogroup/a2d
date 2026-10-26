#ifndef A2D_MAT_VEC_MULT_H
#define A2D_MAT_VEC_MULT_H

#include "../a2ddefs.h"
#include "../a2dmat.h"
#include "core/a2dmatveccore.h"
#include "core/a2dsymmatveccore.h"
#include "core/a2dveccore.h"

namespace A2D {

template <typename T, int N, int M>
A2D_FUNCTION void MatVecMult(const Mat<T, N, M>& A, const Vec<T, M>& x,
                             Vec<T, N>& y) {
  MatVecCore<T, N, M>(get_data(A), get_data(x), get_data(y));
}

template <typename T, int N>
A2D_FUNCTION void MatVecMult(const SymMat<T, N>& A, const Vec<T, N>& x,
                             Vec<T, N>& y) {
  SymMatVecCore<T, N>(get_data(A), get_data(x), get_data(y));
}

template <MatOp op, typename T, int N, int M, int K, int P>
A2D_FUNCTION void MatVecMult(const Mat<T, N, M>& A, const Vec<T, K>& x,
                             Vec<T, P>& y) {
  static_assert(((op == MatOp::NORMAL && (M == K && N == P)) ||
                 (op == MatOp::TRANSPOSE && (M == P && N == K))),
                "Matrix and vector dimensions must agree");
  MatVecCore<T, N, M, op>(get_data(A), get_data(x), get_data(y));
}

template <MatOp op, class Atype, class xtype, class ytype>
class MatVecMultExpr {
 public:
  static constexpr MatOp not_op =
      conditional_value<MatOp, op == MatOp::NORMAL, MatOp::TRANSPOSE,
                        MatOp::NORMAL>::value;

  // Extract the numeric type to use
  using T = object_scalar_t<ytype>;

  // Extract the dimensions of the matrices
  static constexpr int N = get_matrix_rows<Atype>::size;
  static constexpr int M = get_matrix_columns<Atype>::size;
  static constexpr int K = get_vec_size<xtype>::size;
  static constexpr int P = get_vec_size<ytype>::size;

  // Get the types of the matrices
  static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;
  static constexpr ADiffType adx = get_diff_type<xtype>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<ytype>::order;

  A2D_FUNCTION MatVecMultExpr(Atype& A, xtype& x, ytype& y) : A(A), x(x), y(y) {
    static_assert(((op == MatOp::NORMAL && (M == K && N == P)) ||
                   (op == MatOp::TRANSPOSE && (M == P && N == K))),
                  "Matrix and vector dimensions must agree");
  }

  A2D_FUNCTION void eval() {
    MatVecCore<T, N, M, op>(get_data(A), get_data(x), get_data(y));
  }

  A2D_FUNCTION void bzero() { y.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adA == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      constexpr bool additive = true;
      MatVecCore<T, N, M, op>(GetSeed<seed>::get_data(A), get_data(x),
                              GetSeed<seed>::get_data(y));
      MatVecCore<T, N, M, op, additive>(get_data(A), GetSeed<seed>::get_data(x),
                                        GetSeed<seed>::get_data(y));

    } else if constexpr (adA == ADiffType::ACTIVE) {
      MatVecCore<T, N, M, op>(GetSeed<seed>::get_data(A), get_data(x),
                              GetSeed<seed>::get_data(y));
    } else if constexpr (adx == ADiffType::ACTIVE) {
      MatVecCore<T, N, M, op>(get_data(A), GetSeed<seed>::get_data(x),
                              GetSeed<seed>::get_data(y));
    }
  }

  A2D_FUNCTION void reverse() {
    constexpr bool additive = true;
    if constexpr (adA == ADiffType::ACTIVE) {
      if constexpr (op == MatOp::NORMAL) {
        VecOuterCore<T, N, M, additive>(GetSeed<ADseed::b>::get_data(y),
                                        get_data(x),
                                        GetSeed<ADseed::b>::get_data(A));
      } else {
        VecOuterCore<T, N, M, additive>(get_data(x),
                                        GetSeed<ADseed::b>::get_data(y),
                                        GetSeed<ADseed::b>::get_data(A));
      }
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      MatVecCore<T, N, M, not_op, additive>(get_data(A),
                                            GetSeed<ADseed::b>::get_data(y),
                                            GetSeed<ADseed::b>::get_data(x));
    }
  }

  A2D_FUNCTION void hzero() { y.hzero(); }

  A2D_FUNCTION void hreverse() {
    constexpr bool additive = true;
    if constexpr (adA == ADiffType::ACTIVE) {
      if constexpr (op == MatOp::NORMAL) {
        VecOuterCore<T, N, M, additive>(GetSeed<ADseed::h>::get_data(y),
                                        get_data(x),
                                        GetSeed<ADseed::h>::get_data(A));
      } else {
        VecOuterCore<T, N, M, additive>(get_data(x),
                                        GetSeed<ADseed::h>::get_data(y),
                                        GetSeed<ADseed::h>::get_data(A));
      }
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      MatVecCore<T, N, M, not_op, additive>(get_data(A),
                                            GetSeed<ADseed::h>::get_data(y),
                                            GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (adA == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      if constexpr (op == MatOp::NORMAL) {
        VecOuterCore<T, N, M, additive>(GetSeed<ADseed::b>::get_data(y),
                                        GetSeed<ADseed::p>::get_data(x),
                                        GetSeed<ADseed::h>::get_data(A));
      } else {
        VecOuterCore<T, N, M, additive>(GetSeed<ADseed::p>::get_data(x),
                                        GetSeed<ADseed::b>::get_data(y),
                                        GetSeed<ADseed::h>::get_data(A));
      }

      MatVecCore<T, N, M, not_op, additive>(GetSeed<ADseed::p>::get_data(A),
                                            GetSeed<ADseed::b>::get_data(y),
                                            GetSeed<ADseed::h>::get_data(x));
    }
  }

 private:
  Atype& A;
  xtype& x;
  ytype& y;
};

template <class Atype, class xtype, class ytype>
class SymMatVecMultExpr {
 public:
  // Extract the numeric type to use
  using T = object_scalar_t<ytype>;

  // Extract the dimensions of the matrices
  static constexpr int N = get_symmatrix_size<Atype>::size;
  static constexpr int K = get_vec_size<xtype>::size;
  static constexpr int P = get_vec_size<ytype>::size;

  static_assert(K == P, "Input vector and output vector must have same size");
  static_assert(N == P, "matrix and vector must have compatible size");

  // Get the types of the matrices
  static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;
  static constexpr ADiffType adx = get_diff_type<xtype>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<ytype>::order;

  A2D_FUNCTION SymMatVecMultExpr(Atype& A, xtype& x, ytype& y)
      : A(A), x(x), y(y) {}

  A2D_FUNCTION void eval() {
    SymMatVecCore<T, N>(get_data(A), get_data(x), get_data(y));
  }

  A2D_FUNCTION void bzero() { y.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adA == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      constexpr bool additive = true;
      SymMatVecCore<T, N>(GetSeed<seed>::get_data(A), get_data(x),
                          GetSeed<seed>::get_data(y));
      SymMatVecCore<T, N, additive>(get_data(A), GetSeed<seed>::get_data(x),
                                    GetSeed<seed>::get_data(y));

    } else if constexpr (adA == ADiffType::ACTIVE) {
      SymMatVecCore<T, N>(GetSeed<seed>::get_data(A), get_data(x),
                          GetSeed<seed>::get_data(y));
    } else if constexpr (adx == ADiffType::ACTIVE) {
      SymMatVecCore<T, N>(get_data(A), GetSeed<seed>::get_data(x),
                          GetSeed<seed>::get_data(y));
    }
  }

  A2D_FUNCTION void reverse() {
    constexpr bool additive = true;
    if constexpr (adA == ADiffType::ACTIVE) {
      DiagonalPreservingVecSymOuterCore<T, N, additive>(
          GetSeed<ADseed::b>::get_data(y), get_data(x),
          GetSeed<ADseed::b>::get_data(A));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      SymMatVecCore<T, N, additive>(get_data(A),
                                    GetSeed<ADseed::b>::get_data(y),
                                    GetSeed<ADseed::b>::get_data(x));
    }
  }

  A2D_FUNCTION void hzero() { y.hzero(); }

  A2D_FUNCTION void hreverse() {
    constexpr bool additive = true;
    if constexpr (adA == ADiffType::ACTIVE) {
      DiagonalPreservingVecSymOuterCore<T, N, additive>(
          GetSeed<ADseed::h>::get_data(y), get_data(x),
          GetSeed<ADseed::h>::get_data(A));
    }
    if constexpr (adx == ADiffType::ACTIVE) {
      SymMatVecCore<T, N, additive>(get_data(A),
                                    GetSeed<ADseed::h>::get_data(y),
                                    GetSeed<ADseed::h>::get_data(x));
    }
    if constexpr (adA == ADiffType::ACTIVE && adx == ADiffType::ACTIVE) {
      DiagonalPreservingVecSymOuterCore<T, N, additive>(
          GetSeed<ADseed::b>::get_data(y), GetSeed<ADseed::p>::get_data(x),
          GetSeed<ADseed::h>::get_data(A));

      SymMatVecCore<T, N, additive>(GetSeed<ADseed::p>::get_data(A),
                                    GetSeed<ADseed::b>::get_data(y),
                                    GetSeed<ADseed::h>::get_data(x));
    }
  }

 private:
  Atype& A;
  xtype& x;
  ytype& y;
};

template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(ADObj<Atype>& A, ADObj<xtype>& x,
                             ADObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<ADObj<Atype>, ADObj<xtype>, ADObj<ytype>>(A, x, y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, ADObj<Atype>, ADObj<xtype>,
                          ADObj<ytype>>(A, x, y);
  }
}
template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(A2DObj<Atype>& A, A2DObj<xtype>& x,
                             A2DObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<A2DObj<Atype>, A2DObj<xtype>, A2DObj<ytype>>(A, x,
                                                                          y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, A2DObj<Atype>, A2DObj<xtype>,
                          A2DObj<ytype>>(A, x, y);
  }
}
template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(ADObj<Atype>& A, const xtype& x, ADObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<ADObj<Atype>, const xtype, ADObj<ytype>>(A, x, y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, ADObj<Atype>, const xtype,
                          ADObj<ytype>>(A, x, y);
  }
}
template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(A2DObj<Atype>& A, const xtype& x,
                             A2DObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<A2DObj<Atype>, const xtype, A2DObj<ytype>>(A, x,
                                                                        y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, A2DObj<Atype>, const xtype,
                          A2DObj<ytype>>(A, x, y);
  }
}
template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(const Atype& A, ADObj<xtype>& x, ADObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<const Atype, ADObj<xtype>, ADObj<ytype>>(A, x, y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, const Atype, ADObj<xtype>,
                          ADObj<ytype>>(A, x, y);
  }
}
template <class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(const Atype& A, A2DObj<xtype>& x,
                             A2DObj<ytype>& y) {
  if constexpr (a2d_obj_type_v<Atype> == ADObjType::SYMMAT) {
    return SymMatVecMultExpr<const Atype, A2DObj<xtype>, A2DObj<ytype>>(A, x,
                                                                        y);
  } else {
    return MatVecMultExpr<MatOp::NORMAL, const Atype, A2DObj<xtype>,
                          A2DObj<ytype>>(A, x, y);
  }
}

template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(ADObj<Atype>& A, ADObj<xtype>& x,
                             ADObj<ytype>& y) {
  return MatVecMultExpr<op, ADObj<Atype>, ADObj<xtype>, ADObj<ytype>>(A, x, y);
}
template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(A2DObj<Atype>& A, A2DObj<xtype>& x,
                             A2DObj<ytype>& y) {
  return MatVecMultExpr<op, A2DObj<Atype>, A2DObj<xtype>, A2DObj<ytype>>(A, x,
                                                                         y);
}
template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(ADObj<Atype>& A, const xtype& x, ADObj<ytype>& y) {
  return MatVecMultExpr<op, ADObj<Atype>, const xtype, ADObj<ytype>>(A, x, y);
}
template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(A2DObj<Atype>& A, const xtype& x,
                             A2DObj<ytype>& y) {
  return MatVecMultExpr<op, A2DObj<Atype>, const xtype, A2DObj<ytype>>(A, x, y);
}

template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(const Atype& A, ADObj<xtype>& x, ADObj<ytype>& y) {
  return MatVecMultExpr<op, const Atype, ADObj<xtype>, ADObj<ytype>>(A, x, y);
}
template <MatOp op, class Atype, class xtype, class ytype>
A2D_FUNCTION auto MatVecMult(const Atype& A, A2DObj<xtype>& x,
                             A2DObj<ytype>& y) {
  return MatVecMultExpr<op, const Atype, A2DObj<xtype>, A2DObj<ytype>>(A, x, y);
}

// now define MatScale
template <typename T, int M, int N>
A2D_FUNCTION void MatScale(const T alpha, const Mat<T, M, N>& x,
                           Mat<T, M, N>& y) {
  MatScaleCore<T, M, N>(alpha, get_data(x), get_data(y));
}

template <class dtype, class Atype, class Btype>
class MatScaleExpr {
 public:
  // Extract the numeric type to use
  using T = object_scalar_t<dtype>;

  // Extract the dimensions of the underlying vectors
  static constexpr int M = get_matrix_rows<Atype>::size;
  static constexpr int N = get_matrix_columns<Atype>::size;
  static constexpr int size = get_num_matrix_entries<Atype>::size;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<Btype>::order;

  // Get the types of the matrices
  static constexpr ADiffType add = get_diff_type<dtype>::diff_type;
  static constexpr ADiffType adA = get_diff_type<Atype>::diff_type;

  // Make sure the matrix dimensions are consistent
  static_assert((a2d_obj_type_v<Atype> == a2d_obj_type_v<Btype>),
                "Matrices are not all of the same type");

  A2D_FUNCTION MatScaleExpr(dtype alpha, Atype& A, Btype& B)
      : alpha(alpha), A(A), B(B) {}

  A2D_FUNCTION void eval() {
    MatScaleCore<T, M, N>(get_data(alpha), get_data(A), get_data(B));
  }

  A2D_FUNCTION void bzero() { B.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (add == ADiffType::ACTIVE && adA == ADiffType::ACTIVE) {
      MatScaleCore<T, M, N>(GetSeed<seed>::get_data(alpha), get_data(A),
                            GetSeed<seed>::get_data(B));
      VecAddCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(A),
                          GetSeed<seed>::get_data(B));
    } else if constexpr (add == ADiffType::ACTIVE) {
      MatScaleCore<T, M, N>(GetSeed<seed>::get_data(alpha), get_data(A),
                            GetSeed<seed>::get_data(B));
    } else if constexpr (adA == ADiffType::ACTIVE) {
      MatScaleCore<T, M, N>(get_data(alpha), GetSeed<seed>::get_data(A),
                            GetSeed<seed>::get_data(B));
    }
  }
  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (add == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(alpha) +=
          VecDotCore<T, size>(GetSeed<seed>::get_data(B), get_data(A));
    }
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(alpha), GetSeed<seed>::get_data(B),
                          GetSeed<seed>::get_data(A));
    }
  }

  A2D_FUNCTION void hzero() { B.hzero(); }

  A2D_FUNCTION void hreverse() {
    if constexpr (add == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(alpha) +=
          VecDotCore<T, size>(GetSeed<ADseed::h>::get_data(B), get_data(A));
    }
    if constexpr (adA == ADiffType::ACTIVE) {
      VecAddCore<T, size>(get_data(alpha), GetSeed<ADseed::h>::get_data(B),
                          GetSeed<ADseed::h>::get_data(B));
    }
    if constexpr (add == ADiffType::ACTIVE && adA == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(alpha) += VecDotCore<T, size>(
          GetSeed<ADseed::b>::get_data(B), GetSeed<ADseed::p>::get_data(A));
      VecAddCore<T, size>(GetSeed<ADseed::p>::get_data(alpha),
                          GetSeed<ADseed::b>::get_data(B),
                          GetSeed<ADseed::h>::get_data(A));
    }
  }

  dtype alpha;
  Atype& A;
  Btype& B;
};

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(ADObj<T>& alpha, ADObj<Atype>& x, ADObj<Btype>& y) {
  return MatScaleExpr<ADObj<T>&, ADObj<Atype>, ADObj<Atype>>(alpha, x, y);
}

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(const T alpha, ADObj<Atype>& x, ADObj<Btype>& y) {
  return MatScaleExpr<const T, ADObj<Atype>, ADObj<Atype>>(alpha, x, y);
}

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(ADObj<T>& alpha, const Atype& x, ADObj<Btype>& y) {
  return MatScaleExpr<ADObj<T>&, const Atype, ADObj<Atype>>(alpha, x, y);
}

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(A2DObj<T>& alpha, A2DObj<Atype>& x,
                           A2DObj<Btype>& y) {
  return MatScaleExpr<A2DObj<T>&, A2DObj<Atype>, A2DObj<Atype>>(alpha, x, y);
}

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(const T alpha, A2DObj<Atype>& x, A2DObj<Btype>& y) {
  return MatScaleExpr<const T, A2DObj<Atype>, A2DObj<Atype>>(alpha, x, y);
}

template <class T, class Atype, class Btype>
A2D_FUNCTION auto MatScale(A2DObj<T>& alpha, const Atype& x, A2DObj<Btype>& y) {
  return MatScaleExpr<A2DObj<T>&, const Atype, A2DObj<Atype>>(alpha, x, y);
}

}  // namespace A2D

#endif  // A2D_MAT_VEC_MULT_H
