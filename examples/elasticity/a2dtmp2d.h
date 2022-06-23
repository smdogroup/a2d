/*
  0: (0, 0)  1: (0, 1)
  2: (1, 0)  3: (1, 1)
*/

#ifndef A2D_TMP_2D_H
#define A2D_TMP_2D_H

#include <stdlib.h>

#include "a2dmatcore2d.h"
#include "a2dobjs.h"
#include "a2dtypes.h"

namespace A2D {

// Mat2x2MatMult
template <typename ScalarType, bool AT = false, bool BT = false>
inline void Mat2x2MatMult(const Mat<ScalarType, 2, 2>& A,
                          const Mat<ScalarType, 2, 2>& B,
                          Mat<ScalarType, 2, 2>& C) {
  if (AT && BT) {
    MatTrans2x2MatTransMultCore(A, B, C);
  } else if (AT) {
    MatTrans2x2MatMultCore(A, B, C);
  } else if (BT) {
    Mat2x2MatTransMultCore(A, B, C);
  } else {
    Mat2x2MatMultCore(A, B, C);
  }
}

template <class AMatType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADMat2x2MatMultExpr
    : public ADExpression<
          ADMat2x2MatMultExpr<AMatType, BMatType, CMatType, AT, BT> > {
 public:
  ADMat2x2MatMultExpr(ADMat<AMatType>& AObj, ADMat<BMatType>& BObj,
                      ADMat<CMatType>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void forward() {
    const AMatType& A = AObj.value();
    const AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    const BMatType& Bb = BObj.bvalue();
    CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ab, B, Cb);
      MatTrans2x2MatTransMultAddCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ab, B, Cb);
      MatTrans2x2MatMultAddCore(A, Bb, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ab, B, Cb);
      Mat2x2MatTransMultAddCore(A, Bb, Cb);
    } else {
      Mat2x2MatMultCore(Ab, B, Cb);
      Mat2x2MatMultAddCore(A, Bb, Cb);
    }
  }

  void reverse() {
    const AMatType& A = AObj.value();
    AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  ADMat<AMatType>& AObj;
  ADMat<BMatType>& BObj;
  ADMat<CMatType>& CObj;
};

template <class AMatType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADMat2x2MatMultExpr<AMatType, BMatType, CMatType, AT, BT> Mat2x2MatMult(
    ADMat<AMatType>& AObj, ADMat<BMatType>& BObj, ADMat<CMatType>& CObj) {
  return ADMat2x2MatMultExpr<AMatType, BMatType, CMatType, AT, BT>(AObj, BObj,
                                                                   CObj);
}

template <int N, class AMatType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DMat2x2MatMultExpr
    : public A2DExpression<
          A2DMat2x2MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT> > {
 public:
  A2DMat2x2MatMultExpr(A2DMat<N, AMatType>& AObj, A2DMat<N, BMatType>& BObj,
                       A2DMat<N, CMatType>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void reverse() {
    const AMatType& A = AObj.value();
    AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  void hforward() {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();

    for (int i = 0; i < N; i++) {
      const AMatType& Ap = AObj.pvalue(i);
      const BMatType& Bp = BObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultCore(Ap, B, Cp);
        MatTrans2x2MatTransMultAddCore(A, Bp, Cp);
      } else if (AT) {
        MatTrans2x2MatMultCore(Ap, B, Cp);
        MatTrans2x2MatMultAddCore(A, Bp, Cp);
      } else if (BT) {
        Mat2x2MatTransMultCore(Ap, B, Cp);
        Mat2x2MatTransMultAddCore(A, Bp, Cp);
      } else {
        Mat2x2MatMultCore(Ap, B, Cp);
        Mat2x2MatMultAddCore(A, Bp, Cp);
      }
    }
  }

  void hreverse() {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();

    for (int i = 0; i < N; i++) {
      AMatType& Ah = AObj.hvalue(i);
      BMatType& Bh = BObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);
      const CMatType& Cb = CObj.bvalue();
      const AMatType& Ap = AObj.pvalue(i);
      const BMatType& Bp = BObj.pvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultAddCore(B, Ch, Ah);
        MatTrans2x2MatTransMultAddCore(Ch, A, Bh);

        for (int ii = 0; ii < 2; ii++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int kk = 0; kk < 2; kk++) {
              Ah(jj, ii) += Cb(ii, kk) * Bp(kk, jj);
              Bh(jj, ii) += Cb(kk, jj) * Ap(ii, kk);
            }
          }
        }
      } else if (AT) {
        Mat2x2MatTransMultAddCore(B, Ch, Ah);
        Mat2x2MatMultAddCore(A, Ch, Bh);

        for (int ii = 0; ii < 2; ii++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int kk = 0; kk < 2; kk++) {
              Ah(jj, ii) += Cb(ii, kk) * Bp(jj, kk);
              Bh(ii, jj) += Cb(kk, jj) * Ap(ii, kk);
            }
          }
        }
      } else if (BT) {
        Mat2x2MatMultAddCore(Ch, B, Ah);
        MatTrans2x2MatMultAddCore(Ch, A, Bh);

        for (int ii = 0; ii < 2; ii++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int kk = 0; kk < 2; kk++) {
              Ah(ii, jj) += Cb(ii, kk) * Bp(kk, jj);
              Bh(jj, ii) += Cb(kk, jj) * Ap(kk, ii);
            }
          }
        }
      } else {
        Mat2x2MatTransMultAddCore(Ch, B, Ah);
        MatTrans2x2MatMultAddCore(A, Ch, Bh);

        for (int ii = 0; ii < 2; ii++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int kk = 0; kk < 2; kk++) {
              Ah(ii, jj) += Cb(ii, kk) * Bp(jj, kk);
              Bh(ii, jj) += Cb(kk, jj) * Ap(kk, ii);
            }
          }
        }
      }
    }
  }

  A2DMat<N, AMatType>& AObj;
  A2DMat<N, BMatType>& BObj;
  A2DMat<N, CMatType>& CObj;
};

template <int N, class AMatType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DMat2x2MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT>
Mat2x2MatMult(A2DMat<N, AMatType>& AObj, A2DMat<N, BMatType>& BObj,
              A2DMat<N, CMatType>& CObj) {
  return A2DMat2x2MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT>(
      AObj, BObj, CObj);
}

template <typename ScalarType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADpMat2x2MatMultExpr
    : public ADExpression<
          ADpMat2x2MatMultExpr<ScalarType, BMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  ADpMat2x2MatMultExpr(Mat2x2& A, ADMat<BMatType>& BObj, ADMat<CMatType>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void forward() {
    const BMatType& B = BObj.value();
    const BMatType& Bb = BObj.bvalue();
    CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, Bb, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, Bb, Cb);
    } else {
      Mat2x2MatMultCore(A, Bb, Cb);
    }
  }

  void reverse() {
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  Mat2x2& A;
  ADMat<BMatType>& BObj;
  ADMat<CMatType>& CObj;
};

template <class ScalarType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADpMat2x2MatMultExpr<ScalarType, BMatType, CMatType, AT, BT>
Mat2x2MatMult(Mat<ScalarType, 2, 2>& A, ADMat<BMatType>& BObj,
              ADMat<CMatType>& CObj) {
  return ADpMat2x2MatMultExpr<ScalarType, BMatType, CMatType, AT, BT>(A, BObj,
                                                                      CObj);
}

template <typename ScalarType, class AMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADMat2x2pMatMultExpr
    : public ADExpression<
          ADMat2x2pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  ADMat2x2pMatMultExpr(ADMat<AMatType>& AObj, Mat2x2& B, ADMat<CMatType>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const AMatType& A = AObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void forward() {
    const AMatType& A = AObj.value();
    const AMatType& Ab = AObj.bvalue();
    CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(Ab, B, Cb);
    } else if (AT) {
      MatTrans2x2MatMultCore(Ab, B, Cb);
    } else if (BT) {
      Mat2x2MatTransMultCore(Ab, B, Cb);
    } else {
      Mat2x2MatMultCore(Ab, B, Cb);
    }
  }

  void reverse() {
    AMatType& Ab = AObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
    }
  }

  ADMat<AMatType>& AObj;
  Mat2x2& B;
  ADMat<CMatType>& CObj;
};

template <typename ScalarType, class AMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADMat2x2pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT>
Mat2x2MatMult(ADMat<AMatType>& AObj, Mat<ScalarType, 2, 2>& B,
              ADMat<CMatType>& CObj) {
  return ADMat2x2pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT>(AObj, B,
                                                                      CObj);
}

template <int N, typename ScalarType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DpMat2x2MatMultExpr
    : public A2DExpression<
          A2DpMat2x2MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  A2DpMat2x2MatMultExpr(Mat2x2& A, A2DMat<N, BMatType>& BObj,
                        A2DMat<N, CMatType>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void reverse() {
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat2x2MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans2x2MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans2x2MatMultAddCore(A, Cb, Bb);
    }
  }

  void hforward() {
    for (int i = 0; i < N; i++) {
      const BMatType& Bp = BObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultCore(A, Bp, Cp);
      } else if (AT) {
        MatTrans2x2MatMultCore(A, Bp, Cp);
      } else if (BT) {
        Mat2x2MatTransMultCore(A, Bp, Cp);
      } else {
        Mat2x2MatMultCore(A, Bp, Cp);
      }
    }
  }

  void hreverse() {
    for (int i = 0; i < N; i++) {
      BMatType& Bh = BObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultAddCore(Ch, A, Bh);
      } else if (AT) {
        Mat2x2MatMultAddCore(A, Ch, Bh);
      } else if (BT) {
        MatTrans2x2MatMultAddCore(Ch, A, Bh);
      } else {
        MatTrans2x2MatMultAddCore(A, Ch, Bh);
      }
    }
  }

  Mat2x2& A;
  A2DMat<N, BMatType>& BObj;
  A2DMat<N, CMatType>& CObj;
};

template <int N, class ScalarType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DpMat2x2MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT>
Mat2x2MatMult(Mat<ScalarType, 2, 2>& A, A2DMat<N, BMatType>& BObj,
              A2DMat<N, CMatType>& CObj) {
  return A2DpMat2x2MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT>(
      A, BObj, CObj);
}

template <int N, typename ScalarType, class AMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DMat2x2pMatMultExpr
    : public ADExpression<
          A2DMat2x2pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 2, 2> Mat2x2;

  A2DMat2x2pMatMultExpr(A2DMat<N, AMatType>& AObj, Mat2x2& B,
                        A2DMat<N, CMatType>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const AMatType& A = AObj.value();
    CMatType& C = CObj.value();

    if (AT && BT) {
      MatTrans2x2MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans2x2MatMultCore(A, B, C);
    } else if (BT) {
      Mat2x2MatTransMultCore(A, B, C);
    } else {
      Mat2x2MatMultCore(A, B, C);
    }
  }

  void reverse() {
    AMatType& Ab = AObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat2x2MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat2x2MatMultAddCore(Cb, B, Ab);
    } else {
      Mat2x2MatTransMultAddCore(Cb, B, Ab);
    }
  }

  void hforward() {
    for (int i = 0; i < N; i++) {
      const AMatType& Ap = AObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultCore(Ap, B, Cp);
      } else if (AT) {
        MatTrans2x2MatMultCore(Ap, B, Cp);
      } else if (BT) {
        Mat2x2MatTransMultCore(Ap, B, Cp);
      } else {
        Mat2x2MatMultCore(Ap, B, Cp);
      }
    }
  }

  void hreverse() {
    for (int i = 0; i < N; i++) {
      AMatType& Ah = AObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);

      if (AT && BT) {
        MatTrans2x2MatTransMultAddCore(B, Ch, Ah);
      } else if (AT) {
        Mat2x2MatTransMultAddCore(B, Ch, Ah);
      } else if (BT) {
        Mat2x2MatMultAddCore(Ch, B, Ah);
      } else {
        Mat2x2MatTransMultAddCore(Ch, B, Ah);
      }
    }
  }

  A2DMat<N, AMatType>& AObj;
  Mat2x2& B;
  A2DMat<N, CMatType>& CObj;
};

template <int N, typename ScalarType, class AMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DMat2x2pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT>
Mat2x2MatMult(A2DMat<N, AMatType>& AObj, Mat<ScalarType, 2, 2>& B,
              A2DMat<N, CMatType>& CObj) {
  return A2DMat2x2pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT>(
      AObj, B, CObj);
}

// Mat2x2Det
template <typename ScalarType>
inline void Mat2x2Det(const Mat<ScalarType, 2, 2>& A, ScalarType& det) {
  det = A(1, 1) * A(0, 0) - A(0, 1) * A(1, 0);
}

template <class MatType, class ScalarType>
class ADMat2x2DetExpr
    : public ADExpression<ADMat2x2DetExpr<MatType, ScalarType> > {
 public:
  ADMat2x2DetExpr(ADMat<MatType>& AObj, ADScalar<ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const MatType& A = AObj.value();

    detObj.value = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
  }

  void forward() {
    const MatType& A = AObj.value();
    const MatType& Ad = AObj.bvalue();

    detObj.bvalue = (Ad(0, 0) * A(1, 1) + Ad(1, 1) * A(0, 0) -
                     Ad(0, 1) * A(1, 0) - Ad(1, 0) * A(0, 1));
  }

  void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    MatType& Ab = AObj.bvalue();

    Ab(0, 0) += A(1, 1) * bdet;
    Ab(1, 1) += A(0, 0) * bdet;
    Ab(0, 1) -= A(1, 0) * bdet;
    Ab(1, 0) -= A(0, 1) * bdet;
  }

  ADMat<MatType>& AObj;
  ADScalar<ScalarType>& detObj;
};

template <class MatType, typename ScalarType>
inline ADMat2x2DetExpr<MatType, ScalarType> Mat2x2Det(
    ADMat<MatType>& A, ADScalar<ScalarType>& det) {
  return ADMat2x2DetExpr<MatType, ScalarType>(A, det);
}

template <int N, class MatType, class ScalarType>
class A2DMat2x2DetExpr
    : public ADExpression<A2DMat2x2DetExpr<N, MatType, ScalarType> > {
 public:
  A2DMat2x2DetExpr(A2DMat<N, MatType>& AObj, A2DScalar<N, ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const MatType& A = AObj.value();

    detObj.value = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));
  }

  void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    MatType& Ab = AObj.bvalue();

    Ab(0, 0) += A(1, 1) * bdet;
    Ab(1, 1) += A(0, 0) * bdet;
    Ab(0, 1) -= A(1, 0) * bdet;
    Ab(1, 0) -= A(0, 1) * bdet;
  }

  void hforward() {
    const MatType& A = AObj.value();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
      detObj.pvalue[i] = (Ap(0, 0) * A(1, 1) + Ap(1, 1) * A(0, 0) -
                          Ap(0, 1) * A(1, 0) - Ap(1, 0) * A(0, 1));
    }
  }

  void hreverse() {
    const ScalarType bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    const MatType& Ab = AObj.bvalue();

    for (int i = 0; i < N; i++) {
      const ScalarType hdet = detObj.hvalue[i];
      const MatType& Ap = AObj.pvalue(i);
      MatType& Ah = AObj.hvalue(i);

      Ah(0, 0) += Ap(1, 1) * bdet;
      Ah(0, 1) -= Ap(1, 0) * bdet;
      Ah(1, 0) -= Ap(0, 1) * bdet;
      Ah(1, 1) += Ap(0, 0) * bdet;

      Ah(0, 0) += A(1, 1) * hdet;
      Ah(1, 1) += A(0, 0) * hdet;
      Ah(0, 1) -= A(1, 0) * hdet;
      Ah(1, 0) -= A(0, 1) * hdet;
    }
  }

  A2DMat<N, MatType>& AObj;
  A2DScalar<N, ScalarType>& detObj;
};

template <int N, class MatType, typename ScalarType>
inline A2DMat2x2DetExpr<N, MatType, ScalarType> Mat2x2Det(
    A2DMat<N, MatType>& A, A2DScalar<N, ScalarType>& det) {
  return A2DMat2x2DetExpr<N, MatType, ScalarType>(A, det);
}

// Mat2x2Inverse
template <typename ScalarType>
inline void Mat2x2Inverse(const Mat<ScalarType, 2, 2>& A,
                          Mat<ScalarType, 2, 2>& Ainv) {
  ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
  ScalarType detinv = 1.0 / det;

  Ainv(0, 0) = A(1, 1) * detinv;
  Ainv(0, 1) = -A(0, 1) * detinv;
  Ainv(1, 0) = -A(1, 0) * detinv;
  Ainv(1, 1) = A(0, 0) * detinv;
}

template <class MatType, class InvMatType>
class ADMat2x2InverseExpr
    : public ADExpression<ADMat2x2InverseExpr<MatType, InvMatType> > {
 public:
  typedef typename MatType::type ScalarType;

  ADMat2x2InverseExpr(ADMat<MatType>& AObj, ADMat<InvMatType>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const MatType& A = AObj.value();
    InvMatType& Ainv = AinvObj.value();

    ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = A(1, 1) * detinv;
    Ainv(0, 1) = -A(0, 1) * detinv;
    Ainv(1, 0) = -A(1, 0) * detinv;
    Ainv(1, 1) = A(0, 0) * detinv;
  }

  void forward() {
    const InvMatType& Ainv = AinvObj.value();
    const MatType& Ad = AObj.bvalue();
    InvMatType& Ainvd = AinvObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    Mat2x2MatMultCore(Ainv, Ad, tmp);
    Mat2x2MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvd);
  }

  void reverse() {
    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvb = AinvObj.bvalue();
    MatType& Ab = AObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    MatTrans2x2MatMultCore(Ainv, Ainvb, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ab);
  }

  ADMat<MatType>& AObj;
  ADMat<InvMatType>& AinvObj;
};

template <class MatType, class InvMatType>
inline ADMat2x2InverseExpr<MatType, InvMatType> Mat2x2Inverse(
    ADMat<MatType>& AObj, ADMat<InvMatType>& AinvObj) {
  return ADMat2x2InverseExpr<MatType, InvMatType>(AObj, AinvObj);
}

template <int N, class MatType, class InvMatType>
class A2DMat2x2InverseExpr
    : public A2DExpression<A2DMat2x2InverseExpr<N, MatType, InvMatType> > {
 public:
  typedef typename MatType::type ScalarType;

  A2DMat2x2InverseExpr(A2DMat<N, MatType>& AObj, A2DMat<N, InvMatType>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const MatType& A = AObj.value();
    InvMatType& Ainv = AinvObj.value();

    ScalarType det = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = A(1, 1) * detinv;
    Ainv(0, 1) = -A(0, 1) * detinv;
    Ainv(1, 0) = -A(1, 0) * detinv;
    Ainv(1, 1) = A(0, 0) * detinv;
  }

  void reverse() {
    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvd = AinvObj.bvalue();
    MatType& Ad = AObj.bvalue();

    Mat<ScalarType, 2, 2> tmp;
    MatTrans2x2MatMultCore(Ainv, Ainvd, tmp);
    Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ad);
  }

  void hforward() {
    Mat<ScalarType, 2, 2> tmp;
    const InvMatType& Ainv = AinvObj.value();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
      InvMatType& Ainvp = AinvObj.pvalue(i);

      Mat2x2MatMultCore(Ainv, Ap, tmp);
      Mat2x2MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvp);
    }
  }

  // hA = A^{-T} * Ap^{T} * A^{-T} * Ainvb * A^{-T} +
  //      A^{-T} * Ainvb * A^{-T} * Ap^{T} * A^{-T} =
  //    = - (A^{-T} * Ap^{T} * Ab + Ab * Ap^{T} * A^{-T})

  void hreverse() {
    // Temporary matrix
    Mat<ScalarType, 2, 2> tmp, tmp2;

    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvb = AinvObj.bvalue();
    const MatType& Ab = AObj.bvalue();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
      const InvMatType& Ainvh = AinvObj.hvalue(i);
      MatType& Ah = AObj.hvalue(i);

      // Ainv^{T} * Ap^{T} * Ab
      MatTrans2x2MatTransMultCore(Ainv, Ap, tmp);
      Mat2x2MatMultAddScaleCore(ScalarType(-1.0), tmp, Ab, Ah);

      // Ab * Ap^{T} * A^{-T}
      Mat2x2MatTransMultCore(Ab, Ap, tmp);
      Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);

      MatTrans2x2MatMultCore(Ainv, Ainvh, tmp);
      Mat2x2MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);
    }
  }

  A2DMat<N, MatType>& AObj;
  A2DMat<N, InvMatType>& AinvObj;
};

template <int N, class MatType, class InvMatType>
inline A2DMat2x2InverseExpr<N, MatType, InvMatType> Mat2x2Inverse(
    A2DMat<N, MatType>& AObj, A2DMat<N, InvMatType>& AinvObj) {
  return A2DMat2x2InverseExpr<N, MatType, InvMatType>(AObj, AinvObj);
}

// Symm2x2SymmMultTrace
template <typename ScalarType>
inline void Symm2x2SymmMultTrace(const SymmMat<ScalarType, 2>& S,
                                 const SymmMat<ScalarType, 2>& E,
                                 ScalarType& trace) {
  trace = S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
}

template <class SMatType, class EMatType, class ScalarType>
class ADSymm2x2SymmMultTraceExpr
    : public ADExpression<
          ADSymm2x2SymmMultTraceExpr<SMatType, EMatType, ScalarType> > {
 public:
  ADSymm2x2SymmMultTraceExpr(ADMat<SMatType>& SObj, ADMat<EMatType>& EObj,
                             ADScalar<ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SMatType& S = SObj.value();
    const EMatType& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
  }

  void forward() {
    const EMatType& E = EObj.value();
    const EMatType& Ed = EObj.bvalue();
    const SMatType& S = SObj.value();
    const SMatType& Sd = SObj.bvalue();

    output.bvalue = S(0, 0) * Ed(0, 0) + S(1, 1) * Ed(1, 1) +
                    2.0 * (S(0, 1) * Ed(0, 1)) + Sd(0, 0) * E(0, 0) +
                    Sd(1, 1) * E(1, 1) + 2.0 * (Sd(0, 1) * E(0, 1));
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const SMatType& S = SObj.value();
    SMatType& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
  }

  ADMat<SMatType>& SObj;
  ADMat<EMatType>& EObj;
  ADScalar<ScalarType>& output;
};

template <class SMatType, class EMatType, class ScalarType>
inline ADSymm2x2SymmMultTraceExpr<SMatType, EMatType, ScalarType>
Symm2x2SymmMultTrace(ADMat<SMatType>& S, ADMat<EMatType>& E,
                     ADScalar<ScalarType>& trace) {
  return ADSymm2x2SymmMultTraceExpr<SMatType, EMatType, ScalarType>(S, E,
                                                                    trace);
}

template <int N, class SMatType, class EMatType, class ScalarType>
class A2DSymm2x2SymmMultTraceExpr
    : public A2DExpression<
          A2DSymm2x2SymmMultTraceExpr<N, SMatType, EMatType, ScalarType> > {
 public:
  A2DSymm2x2SymmMultTraceExpr(A2DMat<N, SMatType>& SObj,
                              A2DMat<N, EMatType>& EObj,
                              A2DScalar<N, ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SMatType& S = SObj.value();
    const EMatType& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + 2.0 * (S(0, 1) * E(0, 1));
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const SMatType& S = SObj.value();
    SMatType& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  void hforward() {
    const EMatType& E = EObj.value();
    const SMatType& S = SObj.value();

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      const SMatType& Sp = SObj.pvalue(i);

      output.pvalue[i] = S(0, 0) * Ep(0, 0) + S(1, 1) * Ep(1, 1) +
                         2.0 * S(0, 1) * Ep(0, 1) + Sp(0, 0) * E(0, 0) +
                         Sp(1, 1) * E(1, 1) + 2.0 * Sp(0, 1) * E(0, 1);
    }
  }

  void hreverse() {
    const EMatType& E = EObj.value();
    const SMatType& S = SObj.value();

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      const SMatType& Sp = SObj.pvalue(i);
      EMatType& Eh = EObj.hvalue(i);
      SMatType& Sh = SObj.hvalue(i);

      Eh(0, 0) += output.bvalue * Sp(0, 0);
      Eh(1, 1) += output.bvalue * Sp(1, 1);
      Eh(0, 1) += 2.0 * output.bvalue * Sp(0, 1);

      Sh(0, 0) += output.bvalue * Ep(0, 0);
      Sh(1, 1) += output.bvalue * Ep(1, 1);
      Sh(0, 1) += 2.0 * output.bvalue * Ep(0, 1);

      Eh(0, 0) += output.hvalue[i] * S(0, 0);
      Eh(1, 1) += output.hvalue[i] * S(1, 1);
      Eh(0, 1) += 2.0 * output.hvalue[i] * S(0, 1);

      Sh(0, 0) += output.hvalue[i] * E(0, 0);
      Sh(1, 1) += output.hvalue[i] * E(1, 1);
      Sh(0, 1) += 2.0 * output.hvalue[i] * E(0, 1);
    }
  }

  A2DMat<N, SMatType>& SObj;
  A2DMat<N, EMatType>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class SMatType, class EMatType, class ScalarType>
inline A2DSymm2x2SymmMultTraceExpr<N, SMatType, EMatType, ScalarType>
Symm2x2SymmMultTrace(A2DMat<N, SMatType>& S, A2DMat<N, EMatType>& E,
                     A2DScalar<N, ScalarType>& trace) {
  return A2DSymm2x2SymmMultTraceExpr<N, SMatType, EMatType, ScalarType>(S, E,
                                                                        trace);
}

template <class ScalarType>
inline void Symm2x2IsotropicConstitutive(const ScalarType& mu,
                                         const ScalarType& lambda,
                                         const SymmMat<ScalarType, 2>& E,
                                         SymmMat<ScalarType, 2>& S) {
  ScalarType tr = lambda * (E(0, 0) + E(1, 1));
  ScalarType mu2 = 2.0 * mu;
  S(0, 0) = mu2 * E(0, 0) + tr;
  S(0, 1) = mu2 * E(0, 1);
  S(1, 1) = mu2 * E(1, 1) + tr;
}

template <class ScalarType, class EMatType, class SMatType>
class ADSymm2x2IsotropicConstitutiveExpr
    : public ADExpression<
          ADSymm2x2IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType> > {
 public:
  ADSymm2x2IsotropicConstitutiveExpr(const ScalarType& mu,
                                     const ScalarType& lambda,
                                     ADMat<EMatType>& EObj,
                                     ADMat<SMatType>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const EMatType& E = EObj.value();
    SMatType& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(1, 1) = mu2 * E(1, 1) + tr;
  }

  void forward() {
    const EMatType& Ed = EObj.bvalue();
    SMatType& Sd = SObj.bvalue();

    ScalarType tr = lambda * (Ed(0, 0) + Ed(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Sd(0, 0) = mu2 * Ed(0, 0) + tr;
    Sd(0, 1) = mu2 * Ed(0, 1);
    Sd(1, 1) = mu2 * Ed(1, 1) + tr;
  }

  void reverse() {
    const SMatType& Sb = SObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<EMatType>& EObj;
  ADMat<SMatType>& SObj;
};

template <class ScalarType, class EMatType, class SMatType>
inline ADSymm2x2IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType>
Symm2x2IsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                             ADMat<EMatType>& E, ADMat<SMatType>& S) {
  return ADSymm2x2IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType>(
      mu, lambda, E, S);
}

template <int N, class ScalarType, class EMatType, class SMatType>
class A2DSymm2x2IsotropicConstitutiveExpr
    : public A2DExpression<A2DSymm2x2IsotropicConstitutiveExpr<
          N, ScalarType, EMatType, SMatType> > {
 public:
  A2DSymm2x2IsotropicConstitutiveExpr(const ScalarType& mu,
                                      const ScalarType& lambda,
                                      A2DMat<N, EMatType>& EObj,
                                      A2DMat<N, SMatType>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const EMatType& E = EObj.value();
    SMatType& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(1, 1) = mu2 * E(1, 1) + tr;
  }

  void reverse() {
    const SMatType& Sb = SObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
  }

  void hforward() {
    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      SMatType& Sp = SObj.pvalue(i);

      ScalarType tr = lambda * (Ep(0, 0) + Ep(1, 1));
      ScalarType mu2 = 2.0 * mu;
      Sp(0, 0) = mu2 * Ep(0, 0) + tr;
      Sp(0, 1) = mu2 * Ep(0, 1);
      Sp(1, 1) = mu2 * Ep(1, 1) + tr;
    }
  }

  void hreverse() {
    for (int i = 0; i < N; i++) {
      const SMatType& Sh = SObj.hvalue(i);
      EMatType& Eh = EObj.hvalue(i);

      ScalarType tr = lambda * (Sh(0, 0) + Sh(1, 1));
      ScalarType mu2 = 2.0 * mu;
      Eh(0, 0) += mu2 * Sh(0, 0) + tr;
      Eh(0, 1) += mu2 * Sh(0, 1);
      Eh(1, 1) += mu2 * Sh(1, 1) + tr;
    }
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<N, EMatType>& EObj;
  A2DMat<N, SMatType>& SObj;
};

template <int N, class ScalarType, class EMatType, class SMatType>
inline A2DSymm2x2IsotropicConstitutiveExpr<N, ScalarType, EMatType, SMatType>
Symm2x2IsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                             A2DMat<N, EMatType>& E, A2DMat<N, SMatType>& S) {
  return A2DSymm2x2IsotropicConstitutiveExpr<N, ScalarType, EMatType, SMatType>(
      mu, lambda, E, S);
}

template <class ScalarType>
inline void Symm2x2IsotropicEnergy(const ScalarType& mu,
                                   const ScalarType& lambda,
                                   const SymmMat<ScalarType, 2>& E,
                                   ScalarType& output) {
  ScalarType tr = E(0, 0) + E(1, 1);
  ScalarType trE =
      E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

  output = mu * trE + 0.5 * lambda * tr * tr;
}

template <class ScalarType, class EMatType>
class ADSymm2x2IsotropicEnergyExpr
    : public ADExpression<ADSymm2x2IsotropicEnergyExpr<ScalarType, EMatType> > {
 public:
  ADSymm2x2IsotropicEnergyExpr(const ScalarType& mu, const ScalarType& lambda,
                               ADMat<EMatType>& EObj,
                               ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  void forward() {
    const EMatType& E = EObj.value();
    const EMatType& Ed = EObj.bvalue();
    ScalarType tr = (E(0, 0) + E(1, 1));
    ScalarType trd = (Ed(0, 0) + Ed(1, 1));
    ScalarType trEd = 2.0 * (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) +
                             2.0 * E(0, 1) * Ed(0, 1));

    output.bvalue = mu * trEd + lambda * tr * trd;
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = (E(0, 0) + E(1, 1));
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<EMatType>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType, class EMatType>
inline ADSymm2x2IsotropicEnergyExpr<ScalarType, EMatType>
Symm2x2IsotropicEnergy(const ScalarType& mu, const ScalarType& lambda,
                       ADMat<EMatType>& E, ADScalar<ScalarType>& output) {
  return ADSymm2x2IsotropicEnergyExpr<ScalarType, EMatType>(mu, lambda, E,
                                                            output);
}

template <class ScalarType, class EMatType>
class ADSymm2x2ADIsotropicEnergyExpr
    : public ADExpression<
          ADSymm2x2ADIsotropicEnergyExpr<ScalarType, EMatType> > {
 public:
  ADSymm2x2ADIsotropicEnergyExpr(ADScalar<ScalarType>& mu,
                                 ADScalar<ScalarType>& lambda,
                                 ADMat<EMatType>& EObj,
                                 ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  void forward() {
    const EMatType& E = EObj.value();
    const EMatType& Ed = EObj.bvalue();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trd = Ed(0, 0) + Ed(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);
    ScalarType trEd = 2.0 * (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) +
                             2.0 * E(0, 1) * Ed(0, 1));

    output.bvalue = mu.value * trEd + lambda.value * tr * trd +
                    mu.bvalue * trE + 0.5 * lambda.bvalue * tr * tr;
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  ADScalar<ScalarType>& mu;
  ADScalar<ScalarType>& lambda;
  ADMat<EMatType>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType, class EMatType>
inline ADSymm2x2ADIsotropicEnergyExpr<ScalarType, EMatType>
Symm2x2IsotropicEnergy(ADScalar<ScalarType>& mu, ADScalar<ScalarType>& lambda,
                       ADMat<EMatType>& E, ADScalar<ScalarType>& output) {
  return ADSymm2x2ADIsotropicEnergyExpr<ScalarType, EMatType>(mu, lambda, E,
                                                              output);
}

template <int N, class ScalarType, class EMatType>
class A2DSymm2x2IsotropicEnergyExpr
    : public A2DExpression<
          A2DSymm2x2IsotropicEnergyExpr<N, ScalarType, EMatType> > {
 public:
  A2DSymm2x2IsotropicEnergyExpr(const ScalarType& mu, const ScalarType& lambda,
                                A2DMat<N, EMatType>& EObj,
                                A2DScalar<N, ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = E(0, 0) + E(1, 1);
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
  }

  void hforward() {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      ScalarType trd = Ep(0, 0) + Ep(1, 1);
      ScalarType trEd = 2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) +
                               2.0 * E(0, 1) * Ep(0, 1));

      output.pvalue[i] = mu * trEd + lambda * tr * trd;
    }
  }

  void hreverse() {
    const EMatType& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu;
    ScalarType tr = E(0, 0) + E(1, 1);

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      EMatType& Eh = EObj.hvalue(i);

      // by * (d^2y/dx^2 * px)
      ScalarType trp = Ep(0, 0) + Ep(1, 1);
      Eh(0, 0) += (mu2 * Ep(0, 0) + lambda * trp) * output.bvalue;
      Eh(1, 1) += (mu2 * Ep(1, 1) + lambda * trp) * output.bvalue;

      Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;

      // hy * (dy/dx)
      Eh(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.hvalue[i];
      Eh(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.hvalue[i];

      Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue[i];
    }
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<N, EMatType>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType, class EMatType>
inline A2DSymm2x2IsotropicEnergyExpr<N, ScalarType, EMatType>
Symm2x2IsotropicEnergy(const ScalarType& mu, const ScalarType& lambda,
                       A2DMat<N, EMatType>& E,
                       A2DScalar<N, ScalarType>& output) {
  return A2DSymm2x2IsotropicEnergyExpr<N, ScalarType, EMatType>(mu, lambda, E,
                                                                output);
}

template <int N, class ScalarType, class EMatType>
class A2DSymm2x2A2DIsotropicEnergyExpr
    : public A2DExpression<
          A2DSymm2x2A2DIsotropicEnergyExpr<N, ScalarType, EMatType> > {
 public:
  A2DSymm2x2A2DIsotropicEnergyExpr(A2DScalar<N, ScalarType>& mu,
                                   A2DScalar<N, ScalarType>& lambda,
                                   A2DMat<N, EMatType>& EObj,
                                   A2DScalar<N, ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  void hforward() {
    const EMatType& E = EObj.value();
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      ScalarType trd = Ep(0, 0) + Ep(1, 1);
      ScalarType trEd = 2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) +
                               +2.0 * E(0, 1) * Ep(0, 1));

      output.pvalue[i] = mu.value * trEd + lambda.value * tr * trd +
                         mu.pvalue[i] * trE + 0.5 * lambda.pvalue[i] * tr * tr;
    }
  }

  void hreverse() {
    const EMatType& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu.value;
    ScalarType tr = E(0, 0) + E(1, 1);
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + 2.0 * E(0, 1) * E(0, 1);

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      EMatType& Eh = EObj.hvalue(i);

      // by * (d^2y/dx^2 * px)
      ScalarType trp = Ep(0, 0) + Ep(1, 1);
      Eh(0, 0) += (mu2 * Ep(0, 0) + lambda.value * trp) * output.bvalue;
      Eh(1, 1) += (mu2 * Ep(1, 1) + lambda.value * trp) * output.bvalue;

      Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;

      // hy * (dy/dx)
      Eh(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.hvalue[i];
      Eh(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.hvalue[i];

      Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue[i];

      mu.hvalue[i] += trE * output.hvalue[i];
      lambda.hvalue[i] += 0.5 * tr * tr * output.hvalue[i];

      // account for Hessian blocks w.r.t. E and mu or lambda
      Eh(0, 0) +=
          (2 * E(0, 0) * mu.pvalue[i] + tr * lambda.pvalue[i]) * output.bvalue;
      Eh(1, 1) +=
          (2 * E(1, 1) * mu.pvalue[i] + tr * lambda.pvalue[i]) * output.bvalue;
      Eh(0, 1) += 4 * E(0, 1) * mu.pvalue[i] * output.bvalue;

      mu.hvalue[i] += (2 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1)) +
                       4 * E(0, 1) * Ep(0, 1)) *
                      output.bvalue;
      lambda.hvalue[i] += tr * trp * output.bvalue;
    }
  }

  A2DScalar<N, ScalarType>& mu;
  A2DScalar<N, ScalarType>& lambda;
  A2DMat<N, EMatType>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType, class EMatType>
inline A2DSymm2x2A2DIsotropicEnergyExpr<N, ScalarType, EMatType>
Symm2x2IsotropicEnergy(A2DScalar<N, ScalarType>& mu,
                       A2DScalar<N, ScalarType>& lambda, A2DMat<N, EMatType>& E,
                       A2DScalar<N, ScalarType>& output) {
  return A2DSymm2x2A2DIsotropicEnergyExpr<N, ScalarType, EMatType>(mu, lambda,
                                                                   E, output);
}

template <class ScalarType>
inline void Mat2x2GreenStrain(const Mat<ScalarType, 2, 2>& Ux,
                              SymmMat<ScalarType, 2>& E) {
  E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
  E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

  E(0, 1) =
      0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
}

template <class UxMatType, class EMatType>
class ADMat2x2GreenStrainExpr
    : public ADExpression<ADMat2x2GreenStrainExpr<UxMatType, EMatType> > {
 public:
  ADMat2x2GreenStrainExpr(ADMat<UxMatType>& UxObj, ADMat<EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

    E(0, 1) =
        0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
  }

  void forward() {
    const UxMatType& Ux = UxObj.value();
    const UxMatType& Uxd = UxObj.bvalue();
    EMatType& Ed = EObj.bvalue();

    Ed(0, 0) = Uxd(0, 0) + Ux(0, 0) * Uxd(0, 0) + Ux(1, 0) * Uxd(1, 0);
    Ed(1, 1) = Uxd(1, 1) + Ux(0, 1) * Uxd(0, 1) + Ux(1, 1) * Uxd(1, 1);

    Ed(0, 1) = 0.5 * (Uxd(0, 1) + Uxd(1, 0) + Ux(0, 0) * Uxd(0, 1) +
                      Ux(1, 0) * Uxd(1, 1) + Uxd(0, 0) * Ux(0, 1) +
                      Uxd(1, 0) * Ux(1, 1));
  }

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1);
  }

  ADMat<UxMatType>& UxObj;
  ADMat<EMatType>& EObj;
};

template <class UxMatType, class EMatType>
inline ADMat2x2GreenStrainExpr<UxMatType, EMatType> Mat2x2GreenStrain(
    ADMat<UxMatType>& Ux, ADMat<EMatType>& E) {
  return ADMat2x2GreenStrainExpr<UxMatType, EMatType>(Ux, E);
}

template <int N, class UxMatType, class EMatType>
class A2DMat2x2GreenStrainExpr
    : public A2DExpression<A2DMat2x2GreenStrainExpr<N, UxMatType, EMatType> > {
 public:
  A2DMat2x2GreenStrainExpr(A2DMat<N, UxMatType>& UxObj,
                           A2DMat<N, EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1));

    E(0, 1) =
        0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) + Ux(1, 0) * Ux(1, 1));
  }

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1);
  }

  void hforward() {
    const UxMatType& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      EMatType& Ep = EObj.pvalue(i);
      Ep(0, 0) = Uxp(0, 0) + Ux(0, 0) * Uxp(0, 0) + Ux(1, 0) * Uxp(1, 0);
      Ep(1, 1) = Uxp(1, 1) + Ux(0, 1) * Uxp(0, 1) + Ux(1, 1) * Uxp(1, 1);

      Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0) + Ux(0, 0) * Uxp(0, 1) +
                        Ux(1, 0) * Uxp(1, 1) + Uxp(0, 0) * Ux(0, 1) +
                        Uxp(1, 0) * Ux(1, 1));
    }
  }

  void hreverse() {
    const UxMatType& Eb = EObj.bvalue();
    const UxMatType& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      const EMatType& Eh = EObj.hvalue(i);
      UxMatType& Uxh = UxObj.hvalue(i);

      Uxh(0, 0) += Uxp(0, 0) * Eb(0, 0) + 0.5 * Uxp(0, 1) * Eb(0, 1);
      Uxh(0, 1) += 0.5 * Uxp(0, 0) * Eb(0, 1) + Uxp(0, 1) * Eb(1, 1);

      Uxh(1, 0) += Uxp(1, 0) * Eb(0, 0) + 0.5 * Uxp(1, 1) * Eb(0, 1);
      Uxh(1, 1) += 0.5 * Uxp(1, 0) * Eb(0, 1) + Uxp(1, 1) * Eb(1, 1);

      Uxh(0, 0) += (Ux(0, 0) + 1.0) * Eh(0, 0) + 0.5 * Ux(0, 1) * Eh(0, 1);
      Uxh(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eh(0, 1) + Ux(0, 1) * Eh(1, 1);

      Uxh(1, 0) += Ux(1, 0) * Eh(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eh(0, 1);
      Uxh(1, 1) += 0.5 * Ux(1, 0) * Eh(0, 1) + (Ux(1, 1) + 1.0) * Eh(1, 1);
    }
  }

  A2DMat<N, UxMatType>& UxObj;
  A2DMat<N, EMatType>& EObj;
};

template <int N, class UxMatType, class EMatType>
inline A2DMat2x2GreenStrainExpr<N, UxMatType, EMatType> Mat2x2GreenStrain(
    A2DMat<N, UxMatType>& Ux, A2DMat<N, EMatType>& E) {
  return A2DMat2x2GreenStrainExpr<N, UxMatType, EMatType>(Ux, E);
}

template <class ScalarType>
inline void Mat2x2LinearGreenStrain(const Mat<ScalarType, 2, 2>& Ux,
                                    SymmMat<ScalarType, 2>& E) {
  E(0, 0) = Ux(0, 0);
  E(1, 1) = Ux(1, 1);

  E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
}

template <class UxMatType, class EMatType>
class ADMat2x2LinearGreenStrainExpr
    : public ADExpression<ADMat2x2LinearGreenStrainExpr<UxMatType, EMatType> > {
 public:
  ADMat2x2LinearGreenStrainExpr(ADMat<UxMatType>& UxObj, ADMat<EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  }

  void forward() {
    const UxMatType& Uxb = UxObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    Eb(0, 0) = Uxb(0, 0);
    Eb(1, 1) = Uxb(1, 1);

    Eb(0, 1) = 0.5 * (Uxb(0, 1) + Uxb(1, 0));
  }

  void reverse() {
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
  }

  ADMat<UxMatType>& UxObj;
  ADMat<EMatType>& EObj;
};

template <class UxMatType, class EMatType>
inline ADMat2x2LinearGreenStrainExpr<UxMatType, EMatType>
Mat2x2LinearGreenStrain(ADMat<UxMatType>& Ux, ADMat<EMatType>& E) {
  return ADMat2x2LinearGreenStrainExpr<UxMatType, EMatType>(Ux, E);
}

template <int N, class UxMatType, class EMatType>
class A2DMat2x2LinearGreenStrainExpr
    : public A2DExpression<
          A2DMat2x2LinearGreenStrainExpr<N, UxMatType, EMatType> > {
 public:
  A2DMat2x2LinearGreenStrainExpr(A2DMat<N, UxMatType>& UxObj,
                                 A2DMat<N, EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  }

  void reverse() {
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

    // Uxb = Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
  }

  void hforward() {
    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      EMatType& Ep = EObj.pvalue(i);

      Ep(0, 0) = Uxp(0, 0);
      Ep(1, 1) = Uxp(1, 1);

      Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0));
    }
  }

  void hreverse() {
    for (int i = 0; i < N; i++) {
      const EMatType& Eh = EObj.hvalue(i);
      UxMatType& Uxh = UxObj.hvalue(i);

      Uxh(0, 0) += Eh(0, 0);
      Uxh(0, 1) += 0.5 * Eh(0, 1);

      Uxh(1, 0) += 0.5 * Eh(0, 1);
      Uxh(1, 1) += Eh(1, 1);
    }
  }

  A2DMat<N, UxMatType>& UxObj;
  A2DMat<N, EMatType>& EObj;
};

template <int N, class UxMatType, class EMatType>
inline A2DMat2x2LinearGreenStrainExpr<N, UxMatType, EMatType>
Mat2x2LinearGreenStrain(A2DMat<N, UxMatType>& Ux, A2DMat<N, EMatType>& E) {
  return A2DMat2x2LinearGreenStrainExpr<N, UxMatType, EMatType>(Ux, E);
}

}  // namespace A2D

#endif  // A2D_TMP_2D_H
