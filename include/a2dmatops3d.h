#ifndef A2D_TMP_3D_H
#define A2D_TMP_3D_H

#include <stdlib.h>

#include "a2dmatcore3d.h"
#include "a2dobjs.h"
#include "a2dtypes.h"

namespace A2D {

// MatMatMult
template <typename ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION void MatMatMult(const Mat<ScalarType, 3, 3>& A,
                                    const Mat<ScalarType, 3, 3>& B,
                                    Mat<ScalarType, 3, 3>& C) {
  if (AT && BT) {
    MatTrans3x3MatTransMultCore(A, B, C);
  } else if (AT) {
    MatTrans3x3MatMultCore(A, B, C);
  } else if (BT) {
    Mat3x3MatTransMultCore(A, B, C);
  } else {
    Mat3x3MatMultCore(A, B, C);
  }
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADMat3x3MatMultExpr
    : public ADExpression<ADMat3x3MatMultExpr<ScalarType, AT, BT>> {
 public:
  A2D_INLINE_FUNCTION ADMat3x3MatMultExpr(ADMat<Mat<ScalarType, 3, 3>>& AObj,
                                          ADMat<Mat<ScalarType, 3, 3>>& BObj,
                                          ADMat<Mat<ScalarType, 3, 3>>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    const Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(Ab, B, Cb);
      MatTrans3x3MatTransMultAddCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans3x3MatMultCore(Ab, B, Cb);
      MatTrans3x3MatMultAddCore(A, Bb, Cb);
    } else if (BT) {
      Mat3x3MatTransMultCore(Ab, B, Cb);
      Mat3x3MatTransMultAddCore(A, Bb, Cb);
    } else {
      Mat3x3MatMultCore(Ab, B, Cb);
      Mat3x3MatMultAddCore(A, Bb, Cb);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(B, Cb, Ab);
      MatTrans3x3MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat3x3MatTransMultAddCore(B, Cb, Ab);
      Mat3x3MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat3x3MatMultAddCore(Cb, B, Ab);
      MatTrans3x3MatMultAddCore(Cb, A, Bb);
    } else {
      Mat3x3MatTransMultAddCore(Cb, B, Ab);
      MatTrans3x3MatMultAddCore(A, Cb, Bb);
    }
  }

  ADMat<Mat<ScalarType, 3, 3>>& AObj;
  ADMat<Mat<ScalarType, 3, 3>>& BObj;
  ADMat<Mat<ScalarType, 3, 3>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION ADMat3x3MatMultExpr<ScalarType, AT, BT> MatMatMult(
    ADMat<Mat<ScalarType, 3, 3>>& AObj, ADMat<Mat<ScalarType, 3, 3>>& BObj,
    ADMat<Mat<ScalarType, 3, 3>>& CObj) {
  return ADMat3x3MatMultExpr<ScalarType, AT, BT>(AObj, BObj, CObj);
}

template <int N, typename ScalarType, bool AT = false, bool BT = false>
class A2DMat3x3MatMultExpr
    : public A2DExpression<A2DMat3x3MatMultExpr<N, ScalarType, AT, BT>> {
 public:
  A2D_INLINE_FUNCTION A2DMat3x3MatMultExpr(
      A2DMat<N, Mat<ScalarType, 3, 3>>& AObj,
      A2DMat<N, Mat<ScalarType, 3, 3>>& BObj,
      A2DMat<N, Mat<ScalarType, 3, 3>>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(B, Cb, Ab);
      MatTrans3x3MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat3x3MatTransMultAddCore(B, Cb, Ab);
      Mat3x3MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      Mat3x3MatMultAddCore(Cb, B, Ab);
      MatTrans3x3MatMultAddCore(Cb, A, Bb);
    } else {
      Mat3x3MatTransMultAddCore(Cb, B, Ab);
      MatTrans3x3MatMultAddCore(A, Cb, Bb);
    }
  }

  A2D_INLINE_FUNCTION void hforward() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& B = BObj.value();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      const Mat<ScalarType, 3, 3>& Bp = BObj.pvalue(i);
      Mat<ScalarType, 3, 3>& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultCore(Ap, B, Cp);
        MatTrans3x3MatTransMultAddCore(A, Bp, Cp);
      } else if (AT) {
        MatTrans3x3MatMultCore(Ap, B, Cp);
        MatTrans3x3MatMultAddCore(A, Bp, Cp);
      } else if (BT) {
        Mat3x3MatTransMultCore(Ap, B, Cp);
        Mat3x3MatTransMultAddCore(A, Bp, Cp);
      } else {
        Mat3x3MatMultCore(Ap, B, Cp);
        Mat3x3MatMultAddCore(A, Bp, Cp);
      }
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& B = BObj.value();

    for (int i = 0; i < N; i++) {
      Mat<ScalarType, 3, 3>& Ah = AObj.hvalue(i);
      Mat<ScalarType, 3, 3>& Bh = BObj.hvalue(i);
      const Mat<ScalarType, 3, 3>& Ch = CObj.hvalue(i);
      const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      const Mat<ScalarType, 3, 3>& Bp = BObj.pvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultAddCore(B, Ch, Ah);
        MatTrans3x3MatTransMultAddCore(Ch, A, Bh);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++) {
              Ah(jj, ii) += Cb(ii, kk) * Bp(kk, jj);
              Bh(jj, ii) += Cb(kk, jj) * Ap(ii, kk);
            }
          }
        }
      } else if (AT) {
        Mat3x3MatTransMultAddCore(B, Ch, Ah);
        Mat3x3MatMultAddCore(A, Ch, Bh);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++) {
              Ah(jj, ii) += Cb(ii, kk) * Bp(jj, kk);
              Bh(ii, jj) += Cb(kk, jj) * Ap(ii, kk);
            }
          }
        }
      } else if (BT) {
        Mat3x3MatMultAddCore(Ch, B, Ah);
        MatTrans3x3MatMultAddCore(Ch, A, Bh);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++) {
              Ah(ii, jj) += Cb(ii, kk) * Bp(kk, jj);
              Bh(jj, ii) += Cb(kk, jj) * Ap(kk, ii);
            }
          }
        }
      } else {
        Mat3x3MatTransMultAddCore(Ch, B, Ah);
        MatTrans3x3MatMultAddCore(A, Ch, Bh);

        for (int ii = 0; ii < 3; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++) {
              Ah(ii, jj) += Cb(ii, kk) * Bp(jj, kk);
              Bh(ii, jj) += Cb(kk, jj) * Ap(kk, ii);
            }
          }
        }
      }
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& AObj;
  A2DMat<N, Mat<ScalarType, 3, 3>>& BObj;
  A2DMat<N, Mat<ScalarType, 3, 3>>& CObj;
};

template <int N, typename ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION A2DMat3x3MatMultExpr<N, ScalarType, AT, BT> MatMatMult(
    A2DMat<N, Mat<ScalarType, 3, 3>>& AObj,
    A2DMat<N, Mat<ScalarType, 3, 3>>& BObj,
    A2DMat<N, Mat<ScalarType, 3, 3>>& CObj) {
  return A2DMat3x3MatMultExpr<N, ScalarType, AT, BT>(AObj, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADpMat3x3MatMultExpr
    : public ADExpression<ADpMat3x3MatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2D_INLINE_FUNCTION ADpMat3x3MatMultExpr(Mat3x3& A,
                                           ADMat<Mat<ScalarType, 3, 3>>& BObj,
                                           ADMat<Mat<ScalarType, 3, 3>>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, Bb, Cb);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, Bb, Cb);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, Bb, Cb);
    } else {
      Mat3x3MatMultCore(A, Bb, Cb);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat3x3MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans3x3MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans3x3MatMultAddCore(A, Cb, Bb);
    }
  }

  Mat3x3& A;
  ADMat<Mat<ScalarType, 3, 3>>& BObj;
  ADMat<Mat<ScalarType, 3, 3>>& CObj;
};

template <class ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION ADpMat3x3MatMultExpr<ScalarType, AT, BT> MatMatMult(
    Mat<ScalarType, 3, 3>& A, ADMat<Mat<ScalarType, 3, 3>>& BObj,
    ADMat<Mat<ScalarType, 3, 3>>& CObj) {
  return ADpMat3x3MatMultExpr<ScalarType, AT, BT>(A, BObj, CObj);
}

template <typename ScalarType, bool AT = false, bool BT = false>
class ADMat3x3pMatMultExpr
    : public ADExpression<ADMat3x3pMatMultExpr<ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2D_INLINE_FUNCTION ADMat3x3pMatMultExpr(ADMat<Mat<ScalarType, 3, 3>>& AObj,
                                           Mat3x3& B,
                                           ADMat<Mat<ScalarType, 3, 3>>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(Ab, B, Cb);
    } else if (AT) {
      MatTrans3x3MatMultCore(Ab, B, Cb);
    } else if (BT) {
      Mat3x3MatTransMultCore(Ab, B, Cb);
    } else {
      Mat3x3MatMultCore(Ab, B, Cb);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat3x3MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat3x3MatMultAddCore(Cb, B, Ab);
    } else {
      Mat3x3MatTransMultAddCore(Cb, B, Ab);
    }
  }

  ADMat<Mat<ScalarType, 3, 3>>& AObj;
  Mat3x3& B;
  ADMat<Mat<ScalarType, 3, 3>>& CObj;
};

template <typename ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION ADMat3x3pMatMultExpr<ScalarType, AT, BT> MatMatMult(
    ADMat<Mat<ScalarType, 3, 3>>& AObj, Mat<ScalarType, 3, 3>& B,
    ADMat<Mat<ScalarType, 3, 3>>& CObj) {
  return ADMat3x3pMatMultExpr<ScalarType, AT, BT>(AObj, B, CObj);
}

template <int N, typename ScalarType, bool AT = false, bool BT = false>
class A2DpMat3x3MatMultExpr
    : public A2DExpression<A2DpMat3x3MatMultExpr<N, ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2D_INLINE_FUNCTION A2DpMat3x3MatMultExpr(
      Mat3x3& A, A2DMat<N, Mat<ScalarType, 3, 3>>& BObj,
      A2DMat<N, Mat<ScalarType, 3, 3>>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& B = BObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    Mat<ScalarType, 3, 3>& Bb = BObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(Cb, A, Bb);
    } else if (AT) {
      Mat3x3MatMultAddCore(A, Cb, Bb);
    } else if (BT) {
      MatTrans3x3MatMultAddCore(Cb, A, Bb);
    } else {
      MatTrans3x3MatMultAddCore(A, Cb, Bb);
    }
  }

  A2D_INLINE_FUNCTION void hforward() {
    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Bp = BObj.pvalue(i);
      Mat<ScalarType, 3, 3>& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultCore(A, Bp, Cp);
      } else if (AT) {
        MatTrans3x3MatMultCore(A, Bp, Cp);
      } else if (BT) {
        Mat3x3MatTransMultCore(A, Bp, Cp);
      } else {
        Mat3x3MatMultCore(A, Bp, Cp);
      }
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    for (int i = 0; i < N; i++) {
      Mat<ScalarType, 3, 3>& Bh = BObj.hvalue(i);
      const Mat<ScalarType, 3, 3>& Ch = CObj.hvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultAddCore(Ch, A, Bh);
      } else if (AT) {
        Mat3x3MatMultAddCore(A, Ch, Bh);
      } else if (BT) {
        MatTrans3x3MatMultAddCore(Ch, A, Bh);
      } else {
        MatTrans3x3MatMultAddCore(A, Ch, Bh);
      }
    }
  }

  Mat3x3& A;
  A2DMat<N, Mat<ScalarType, 3, 3>>& BObj;
  A2DMat<N, Mat<ScalarType, 3, 3>>& CObj;
};

template <int N, class ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION A2DpMat3x3MatMultExpr<N, ScalarType, AT, BT> MatMatMult(
    Mat<ScalarType, 3, 3>& A, A2DMat<N, Mat<ScalarType, 3, 3>>& BObj,
    A2DMat<N, Mat<ScalarType, 3, 3>>& CObj) {
  return A2DpMat3x3MatMultExpr<N, ScalarType, AT, BT>(A, BObj, CObj);
}

template <int N, typename ScalarType, bool AT = false, bool BT = false>
class A2DMat3x3pMatMultExpr
    : public ADExpression<A2DMat3x3pMatMultExpr<N, ScalarType, AT, BT>> {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2D_INLINE_FUNCTION A2DMat3x3pMatMultExpr(
      A2DMat<N, Mat<ScalarType, 3, 3>>& AObj, Mat3x3& B,
      A2DMat<N, Mat<ScalarType, 3, 3>>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& C = CObj.value();

    if (AT && BT) {
      MatTrans3x3MatTransMultCore(A, B, C);
    } else if (AT) {
      MatTrans3x3MatMultCore(A, B, C);
    } else if (BT) {
      Mat3x3MatTransMultCore(A, B, C);
    } else {
      Mat3x3MatMultCore(A, B, C);
    }
  }

  A2D_INLINE_FUNCTION void reverse() {
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();
    const Mat<ScalarType, 3, 3>& Cb = CObj.bvalue();

    if (AT && BT) {
      MatTrans3x3MatTransMultAddCore(B, Cb, Ab);
    } else if (AT) {
      Mat3x3MatTransMultAddCore(B, Cb, Ab);
    } else if (BT) {
      Mat3x3MatMultAddCore(Cb, B, Ab);
    } else {
      Mat3x3MatTransMultAddCore(Cb, B, Ab);
    }
  }

  A2D_INLINE_FUNCTION void hforward() {
    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      Mat<ScalarType, 3, 3>& Cp = CObj.pvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultCore(Ap, B, Cp);
      } else if (AT) {
        MatTrans3x3MatMultCore(Ap, B, Cp);
      } else if (BT) {
        Mat3x3MatTransMultCore(Ap, B, Cp);
      } else {
        Mat3x3MatMultCore(Ap, B, Cp);
      }
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    for (int i = 0; i < N; i++) {
      Mat<ScalarType, 3, 3>& Ah = AObj.hvalue(i);
      const Mat<ScalarType, 3, 3>& Ch = CObj.hvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultAddCore(B, Ch, Ah);
      } else if (AT) {
        Mat3x3MatTransMultAddCore(B, Ch, Ah);
      } else if (BT) {
        Mat3x3MatMultAddCore(Ch, B, Ah);
      } else {
        Mat3x3MatTransMultAddCore(Ch, B, Ah);
      }
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& AObj;
  Mat3x3& B;
  A2DMat<N, Mat<ScalarType, 3, 3>>& CObj;
};

template <int N, typename ScalarType, bool AT = false, bool BT = false>
A2D_INLINE_FUNCTION A2DMat3x3pMatMultExpr<N, ScalarType, AT, BT> MatMatMult(
    A2DMat<N, Mat<ScalarType, 3, 3>>& AObj, Mat<ScalarType, 3, 3>& B,
    A2DMat<N, Mat<ScalarType, 3, 3>>& CObj) {
  return A2DMat3x3pMatMultExpr<N, ScalarType, AT, BT>(AObj, B, CObj);
}

// Mat3x3Det
template <typename ScalarType>
A2D_INLINE_FUNCTION void MatDet(const Mat<ScalarType, 3, 3>& A,
                                ScalarType& det) {
  det = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
         A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
         A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
}

template <class ScalarType>
class ADMat3x3DetExpr : public ADExpression<ADMat3x3DetExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADMat3x3DetExpr(ADMat<Mat<ScalarType, 3, 3>>& AObj,
                                      ADScalar<ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();

    detObj.value = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                    A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                    A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    const Mat<ScalarType, 3, 3>& Ad = AObj.bvalue();

    detObj.bvalue = (Ad(0, 0) * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) +
                     Ad(0, 1) * (A(2, 0) * A(1, 2) - A(2, 2) * A(1, 0)) +
                     Ad(0, 2) * (A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) +
                     Ad(1, 0) * (A(2, 1) * A(0, 2) - A(2, 2) * A(0, 1)) +
                     Ad(1, 1) * (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) +
                     Ad(1, 2) * (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)) +
                     Ad(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) +
                     Ad(2, 1) * (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) +
                     Ad(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();

    Ab(0, 0) += (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) * bdet;
    Ab(0, 1) += (A(2, 0) * A(1, 2) - A(2, 2) * A(1, 0)) * bdet;
    Ab(0, 2) += (A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) * bdet;
    Ab(1, 0) += (A(2, 1) * A(0, 2) - A(2, 2) * A(0, 1)) * bdet;
    Ab(1, 1) += (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) * bdet;
    Ab(1, 2) += (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)) * bdet;
    Ab(2, 0) += (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * bdet;
    Ab(2, 1) += (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * bdet;
    Ab(2, 2) += (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * bdet;
  }

  ADMat<Mat<ScalarType, 3, 3>>& AObj;
  ADScalar<ScalarType>& detObj;
};

template <typename ScalarType>
A2D_INLINE_FUNCTION ADMat3x3DetExpr<ScalarType> MatDet(
    ADMat<Mat<ScalarType, 3, 3>>& A, ADScalar<ScalarType>& det) {
  return ADMat3x3DetExpr<ScalarType>(A, det);
}

template <int N, class ScalarType>
class A2DMat3x3DetExpr : public ADExpression<A2DMat3x3DetExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DMat3x3DetExpr(A2DMat<N, Mat<ScalarType, 3, 3>>& AObj,
                                       A2DScalar<N, ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();

    detObj.value = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                    A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                    A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();

    Ab(0, 0) += (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) * bdet;
    Ab(0, 1) += (A(2, 0) * A(1, 2) - A(2, 2) * A(1, 0)) * bdet;
    Ab(0, 2) += (A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) * bdet;
    Ab(1, 0) += (A(2, 1) * A(0, 2) - A(2, 2) * A(0, 1)) * bdet;
    Ab(1, 1) += (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) * bdet;
    Ab(1, 2) += (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)) * bdet;
    Ab(2, 0) += (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * bdet;
    Ab(2, 1) += (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * bdet;
    Ab(2, 2) += (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * bdet;
  }

  A2D_INLINE_FUNCTION void hforward() {
    const Mat<ScalarType, 3, 3>& A = AObj.value();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      detObj.pvalue[i] = (Ap(0, 0) * (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) +
                          Ap(0, 1) * (A(2, 0) * A(1, 2) - A(2, 2) * A(1, 0)) +
                          Ap(0, 2) * (A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) +
                          Ap(1, 0) * (A(2, 1) * A(0, 2) - A(2, 2) * A(0, 1)) +
                          Ap(1, 1) * (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) +
                          Ap(1, 2) * (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)) +
                          Ap(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) +
                          Ap(2, 1) * (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) +
                          Ap(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)));
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const ScalarType bdet = detObj.bvalue;
    const Mat<ScalarType, 3, 3>& A = AObj.value();

    for (int i = 0; i < N; i++) {
      const ScalarType hdet = detObj.hvalue[i];
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      Mat<ScalarType, 3, 3>& Ah = AObj.hvalue(i);

      Ah(0, 0) += (A(2, 2) * Ap(1, 1) - A(2, 1) * Ap(1, 2) +
                   Ap(2, 2) * A(1, 1) - Ap(2, 1) * A(1, 2)) *
                  bdet;
      Ah(0, 1) += (A(2, 0) * Ap(1, 2) - A(2, 2) * Ap(1, 0) +
                   Ap(2, 0) * A(1, 2) - Ap(2, 2) * A(1, 0)) *
                  bdet;
      Ah(0, 2) += (A(2, 1) * Ap(1, 0) - A(2, 0) * Ap(1, 1) +
                   Ap(2, 1) * A(1, 0) - Ap(2, 0) * A(1, 1)) *
                  bdet;
      Ah(1, 0) += (A(2, 1) * Ap(0, 2) - A(2, 2) * Ap(0, 1) +
                   Ap(2, 1) * A(0, 2) - Ap(2, 2) * A(0, 1)) *
                  bdet;
      Ah(1, 1) += (A(2, 2) * Ap(0, 0) - A(2, 0) * Ap(0, 2) +
                   Ap(2, 2) * A(0, 0) - Ap(2, 0) * A(0, 2)) *
                  bdet;
      Ah(1, 2) += (A(2, 0) * Ap(0, 1) - A(2, 1) * Ap(0, 0) +
                   Ap(2, 0) * A(0, 1) - Ap(2, 1) * A(0, 0)) *
                  bdet;
      Ah(2, 0) += (A(0, 1) * Ap(1, 2) - A(0, 2) * Ap(1, 1) +
                   Ap(0, 1) * A(1, 2) - Ap(0, 2) * A(1, 1)) *
                  bdet;
      Ah(2, 1) += (A(1, 0) * Ap(0, 2) - A(0, 0) * Ap(1, 2) +
                   Ap(1, 0) * A(0, 2) - Ap(0, 0) * A(1, 2)) *
                  bdet;
      Ah(2, 2) += (A(0, 0) * Ap(1, 1) - A(1, 0) * Ap(0, 1) +
                   Ap(0, 0) * A(1, 1) - Ap(1, 0) * A(0, 1)) *
                  bdet;

      Ah(0, 0) += (A(2, 2) * A(1, 1) - A(2, 1) * A(1, 2)) * hdet;
      Ah(0, 1) += (A(2, 0) * A(1, 2) - A(2, 2) * A(1, 0)) * hdet;
      Ah(0, 2) += (A(2, 1) * A(1, 0) - A(2, 0) * A(1, 1)) * hdet;
      Ah(1, 0) += (A(2, 1) * A(0, 2) - A(2, 2) * A(0, 1)) * hdet;
      Ah(1, 1) += (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)) * hdet;
      Ah(1, 2) += (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)) * hdet;
      Ah(2, 0) += (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * hdet;
      Ah(2, 1) += (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * hdet;
      Ah(2, 2) += (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * hdet;
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& AObj;
  A2DScalar<N, ScalarType>& detObj;
};

template <int N, typename ScalarType>
A2D_INLINE_FUNCTION A2DMat3x3DetExpr<N, ScalarType> MatDet(
    A2DMat<N, Mat<ScalarType, 3, 3>>& A, A2DScalar<N, ScalarType>& det) {
  return A2DMat3x3DetExpr<N, ScalarType>(A, det);
}

// Mat3x3Inverse
template <typename ScalarType>
A2D_INLINE_FUNCTION void MatInverse(const Mat<ScalarType, 3, 3>& A,
                                    Mat<ScalarType, 3, 3>& Ainv) {
  ScalarType det = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                    A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                    A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
  ScalarType detinv = 1.0 / det;

  Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * detinv;
  Ainv(0, 1) = -(A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) * detinv;
  Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * detinv;

  Ainv(1, 0) = -(A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) * detinv;
  Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * detinv;
  Ainv(1, 2) = -(A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0)) * detinv;

  Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) * detinv;
  Ainv(2, 1) = -(A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0)) * detinv;
  Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) * detinv;
}

template <typename ScalarType>
class ADMat3x3InverseExpr
    : public ADExpression<ADMat3x3InverseExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADMat3x3InverseExpr(ADMat<Mat<ScalarType, 3, 3>>& AObj,
                                          ADMat<Mat<ScalarType, 3, 3>>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();

    ScalarType det = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                      A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                      A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * detinv;
    Ainv(0, 1) = -(A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) * detinv;
    Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * detinv;

    Ainv(1, 0) = -(A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) * detinv;
    Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * detinv;
    Ainv(1, 2) = -(A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0)) * detinv;

    Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) * detinv;
    Ainv(2, 1) = -(A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0)) * detinv;
    Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) * detinv;
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();
    const Mat<ScalarType, 3, 3>& Ad = AObj.bvalue();
    Mat<ScalarType, 3, 3>& Ainvd = AinvObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    Mat3x3MatMultCore(Ainv, Ad, tmp);
    Mat3x3MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvd);
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();
    const Mat<ScalarType, 3, 3>& Ainvb = AinvObj.bvalue();
    Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    MatTrans3x3MatMultCore(Ainv, Ainvb, tmp);
    Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ab);
  }

  ADMat<Mat<ScalarType, 3, 3>>& AObj;
  ADMat<Mat<ScalarType, 3, 3>>& AinvObj;
};

template <typename ScalarType>
A2D_INLINE_FUNCTION ADMat3x3InverseExpr<ScalarType> MatInverse(
    ADMat<Mat<ScalarType, 3, 3>>& AObj, ADMat<Mat<ScalarType, 3, 3>>& AinvObj) {
  return ADMat3x3InverseExpr<ScalarType>(AObj, AinvObj);
}

template <int N, typename ScalarType>
class A2DMat3x3InverseExpr
    : public A2DExpression<A2DMat3x3InverseExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DMat3x3InverseExpr(
      A2DMat<N, Mat<ScalarType, 3, 3>>& AObj,
      A2DMat<N, Mat<ScalarType, 3, 3>>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const Mat<ScalarType, 3, 3>& A = AObj.value();
    Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();

    ScalarType det = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                      A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                      A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
    ScalarType detinv = 1.0 / det;

    Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) * detinv;
    Ainv(0, 1) = -(A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) * detinv;
    Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * detinv;

    Ainv(1, 0) = -(A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) * detinv;
    Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * detinv;
    Ainv(1, 2) = -(A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0)) * detinv;

    Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) * detinv;
    Ainv(2, 1) = -(A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0)) * detinv;
    Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) * detinv;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();
    const Mat<ScalarType, 3, 3>& Ainvd = AinvObj.bvalue();
    Mat<ScalarType, 3, 3>& Ad = AObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    MatTrans3x3MatMultCore(Ainv, Ainvd, tmp);
    Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ad);
  }

  A2D_INLINE_FUNCTION void hforward() {
    Mat<ScalarType, 3, 3> tmp;
    const Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      Mat<ScalarType, 3, 3>& Ainvp = AinvObj.pvalue(i);

      Mat3x3MatMultCore(Ainv, Ap, tmp);
      Mat3x3MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvp);
    }
  }

  // hA = A^{-T} * Ap^{T} * A^{-T} * Ainvb * A^{-T} +
  //      A^{-T} * Ainvb * A^{-T} * Ap^{T} * A^{-T} =
  //    = - (A^{-T} * Ap^{T} * Ab + Ab * Ap^{T} * A^{-T})

  A2D_INLINE_FUNCTION void hreverse() {
    // Temporary matrix
    Mat<ScalarType, 3, 3> tmp, tmp2;

    const Mat<ScalarType, 3, 3>& Ainv = AinvObj.value();
    const Mat<ScalarType, 3, 3>& Ab = AObj.bvalue();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Ap = AObj.pvalue(i);
      const Mat<ScalarType, 3, 3>& Ainvh = AinvObj.hvalue(i);
      Mat<ScalarType, 3, 3>& Ah = AObj.hvalue(i);

      // Ainv^{T} * Ap^{T} * Ab
      MatTrans3x3MatTransMultCore(Ainv, Ap, tmp);
      Mat3x3MatMultAddScaleCore(ScalarType(-1.0), tmp, Ab, Ah);

      // Ab * Ap^{T} * A^{-T}
      Mat3x3MatTransMultCore(Ab, Ap, tmp);
      Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);

      MatTrans3x3MatMultCore(Ainv, Ainvh, tmp);
      Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ah);
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& AObj;
  A2DMat<N, Mat<ScalarType, 3, 3>>& AinvObj;
};

template <int N, typename ScalarType>
A2D_INLINE_FUNCTION A2DMat3x3InverseExpr<N, ScalarType> MatInverse(
    A2DMat<N, Mat<ScalarType, 3, 3>>& AObj,
    A2DMat<N, Mat<ScalarType, 3, 3>>& AinvObj) {
  return A2DMat3x3InverseExpr<N, ScalarType>(AObj, AinvObj);
}

// SymmTrace
template <typename ScalarType>
A2D_INLINE_FUNCTION void SymmTrace(const SymmMat<ScalarType, 3>& S,
                                   ScalarType& trace) {
  trace = S(0, 0) + S(1, 1) + S(2, 2);
}

template <class ScalarType>
class ADSymm3x3TraceExpr : public ADExpression<ADSymm3x3TraceExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADSymm3x3TraceExpr(ADMat<SymmMat<ScalarType, 3>>& SObj,
                                         ADScalar<ScalarType>& output)
      : SObj(SObj), output(output) {
    const SymmMat<ScalarType, 3>& S = SObj.value();
    output.value = S(0, 0) + S(1, 1) + S(2, 2);
  }

  A2D_INLINE_FUNCTION void forward() {
    const SymmMat<ScalarType, 3>& Sd = SObj.bvalue();
    output.bvalue = Sd(0, 0) + Sd(1, 1) + Sd(2, 2);
  }

  A2D_INLINE_FUNCTION void reverse() {
    SymmMat<ScalarType, 3>& Sb = SObj.bvalue();

    Sb(0, 0) += output.bvalue;
    Sb(1, 1) += output.bvalue;
    Sb(2, 2) += output.bvalue;
  }

  ADMat<SymmMat<ScalarType, 3>>& SObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
A2D_INLINE_FUNCTION ADSymm3x3TraceExpr<ScalarType> SymmTrace(
    ADMat<SymmMat<ScalarType, 3>>& S, ADScalar<ScalarType>& trace) {
  return ADSymm3x3TraceExpr<ScalarType>(S, trace);
}

template <int N, class ScalarType>
class A2DSymm3x3TraceExpr
    : public A2DExpression<A2DSymm3x3TraceExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DSymm3x3TraceExpr(
      A2DMat<N, SymmMat<ScalarType, 3>>& SObj, A2DScalar<N, ScalarType>& output)
      : SObj(SObj), output(output) {
    const SymmMat<ScalarType, 3>& S = SObj.value();
    output.value = S(0, 0) + S(1, 1) + S(2, 2);
  }

  A2D_INLINE_FUNCTION void reverse() {
    SymmMat<ScalarType, 3>& Sb = SObj.bvalue();

    Sb(0, 0) += output.bvalue;
    Sb(1, 1) += output.bvalue;
    Sb(2, 2) += output.bvalue;
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  A2D_INLINE_FUNCTION void hforward() {
    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Sp = SObj.pvalue(i);
      output.pvalue[i] = Sp(0, 0) + Sp(1, 1) + Sp(2, 2);
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    for (int i = 0; i < N; i++) {
      SymmMat<ScalarType, 3>& Sh = SObj.hvalue(i);

      Sh(0, 0) += output.hvalue[i];
      Sh(1, 1) += output.hvalue[i];
      Sh(2, 2) += output.hvalue[i];
    }
  }

  A2DMat<N, SymmMat<ScalarType, 3>>& SObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType>
A2D_INLINE_FUNCTION A2DSymm3x3TraceExpr<N, ScalarType> SymmTrace(
    A2DMat<N, SymmMat<ScalarType, 3>>& S, A2DScalar<N, ScalarType>& trace) {
  return A2DSymm3x3TraceExpr<N, ScalarType>(S, trace);
}

// Symm3x3SymmMultTrace
template <typename ScalarType>
A2D_INLINE_FUNCTION void SymmSymmMultTrace(const SymmMat<ScalarType, 3>& S,
                                           const SymmMat<ScalarType, 3>& E,
                                           ScalarType& trace) {
  trace = (S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
           2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2)));
}

template <class ScalarType>
class ADSymm3x3SymmMultTraceExpr
    : public ADExpression<ADSymm3x3SymmMultTraceExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADSymm3x3SymmMultTraceExpr(
      ADMat<SymmMat<ScalarType, 3>>& SObj, ADMat<SymmMat<ScalarType, 3>>& EObj,
      ADScalar<ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& S = SObj.value();
    const SymmMat<ScalarType, 3>& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
        2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2));
  }

  A2D_INLINE_FUNCTION void forward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const SymmMat<ScalarType, 3>& Ed = EObj.bvalue();
    const SymmMat<ScalarType, 3>& S = SObj.value();
    const SymmMat<ScalarType, 3>& Sd = SObj.bvalue();

    output.bvalue =
        S(0, 0) * Ed(0, 0) + S(1, 1) * Ed(1, 1) + S(2, 2) * Ed(2, 2) +
        2.0 * (S(0, 1) * Ed(0, 1) + S(0, 2) * Ed(0, 2) + S(1, 2) * Ed(1, 2)) +
        Sd(0, 0) * E(0, 0) + Sd(1, 1) * E(1, 1) + Sd(2, 2) * E(2, 2) +
        2.0 * (Sd(0, 1) * E(0, 1) + Sd(0, 2) * E(0, 2) + Sd(1, 2) * E(1, 2));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const SymmMat<ScalarType, 3>& S = SObj.value();
    SymmMat<ScalarType, 3>& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(2, 2) += output.bvalue * S(2, 2);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);
    Eb(0, 2) += 2.0 * output.bvalue * S(0, 2);
    Eb(1, 2) += 2.0 * output.bvalue * S(1, 2);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(2, 2) += output.bvalue * E(2, 2);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
    Sb(0, 2) += 2.0 * output.bvalue * E(0, 2);
    Sb(1, 2) += 2.0 * output.bvalue * E(1, 2);
  }

  ADMat<SymmMat<ScalarType, 3>>& SObj;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
A2D_INLINE_FUNCTION ADSymm3x3SymmMultTraceExpr<ScalarType> SymmSymmMultTrace(
    ADMat<SymmMat<ScalarType, 3>>& S, ADMat<SymmMat<ScalarType, 3>>& E,
    ADScalar<ScalarType>& trace) {
  return ADSymm3x3SymmMultTraceExpr<ScalarType>(S, E, trace);
}

template <int N, class ScalarType>
class A2DSymm3x3SymmMultTraceExpr
    : public A2DExpression<A2DSymm3x3SymmMultTraceExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DSymm3x3SymmMultTraceExpr(
      A2DMat<N, SymmMat<ScalarType, 3>>& SObj,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj, A2DScalar<N, ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& S = SObj.value();
    const SymmMat<ScalarType, 3>& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
        2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const SymmMat<ScalarType, 3>& S = SObj.value();
    SymmMat<ScalarType, 3>& Sb = SObj.bvalue();

    Eb(0, 0) += output.bvalue * S(0, 0);
    Eb(1, 1) += output.bvalue * S(1, 1);
    Eb(2, 2) += output.bvalue * S(2, 2);
    Eb(0, 1) += 2.0 * output.bvalue * S(0, 1);
    Eb(0, 2) += 2.0 * output.bvalue * S(0, 2);
    Eb(1, 2) += 2.0 * output.bvalue * S(1, 2);

    Sb(0, 0) += output.bvalue * E(0, 0);
    Sb(1, 1) += output.bvalue * E(1, 1);
    Sb(2, 2) += output.bvalue * E(2, 2);
    Sb(0, 1) += 2.0 * output.bvalue * E(0, 1);
    Sb(0, 2) += 2.0 * output.bvalue * E(0, 2);
    Sb(1, 2) += 2.0 * output.bvalue * E(1, 2);
  }

  // Compute E.pvalue() = J * Ux.pvalue()
  A2D_INLINE_FUNCTION void hforward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const SymmMat<ScalarType, 3>& S = SObj.value();

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      const SymmMat<ScalarType, 3>& Sp = SObj.pvalue(i);

      output.pvalue[i] =
          S(0, 0) * Ep(0, 0) + S(1, 1) * Ep(1, 1) + S(2, 2) * Ep(2, 2) +
          2.0 * (S(0, 1) * Ep(0, 1) + S(0, 2) * Ep(0, 2) + S(1, 2) * Ep(1, 2)) +
          Sp(0, 0) * E(0, 0) + Sp(1, 1) * E(1, 1) + Sp(2, 2) * E(2, 2) +
          2.0 * (Sp(0, 1) * E(0, 1) + Sp(0, 2) * E(0, 2) + Sp(1, 2) * E(1, 2));
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const SymmMat<ScalarType, 3>& S = SObj.value();

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      const SymmMat<ScalarType, 3>& Sp = SObj.pvalue(i);
      SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);
      SymmMat<ScalarType, 3>& Sh = SObj.hvalue(i);

      Eh(0, 0) += output.bvalue * Sp(0, 0);
      Eh(1, 1) += output.bvalue * Sp(1, 1);
      Eh(2, 2) += output.bvalue * Sp(2, 2);
      Eh(0, 1) += 2.0 * output.bvalue * Sp(0, 1);
      Eh(0, 2) += 2.0 * output.bvalue * Sp(0, 2);
      Eh(1, 2) += 2.0 * output.bvalue * Sp(1, 2);

      Sh(0, 0) += output.bvalue * Ep(0, 0);
      Sh(1, 1) += output.bvalue * Ep(1, 1);
      Sh(2, 2) += output.bvalue * Ep(2, 2);
      Sh(0, 1) += 2.0 * output.bvalue * Ep(0, 1);
      Sh(0, 2) += 2.0 * output.bvalue * Ep(0, 2);
      Sh(1, 2) += 2.0 * output.bvalue * Ep(1, 2);

      Eh(0, 0) += output.hvalue[i] * S(0, 0);
      Eh(1, 1) += output.hvalue[i] * S(1, 1);
      Eh(2, 2) += output.hvalue[i] * S(2, 2);
      Eh(0, 1) += 2.0 * output.hvalue[i] * S(0, 1);
      Eh(0, 2) += 2.0 * output.hvalue[i] * S(0, 2);
      Eh(1, 2) += 2.0 * output.hvalue[i] * S(1, 2);

      Sh(0, 0) += output.hvalue[i] * E(0, 0);
      Sh(1, 1) += output.hvalue[i] * E(1, 1);
      Sh(2, 2) += output.hvalue[i] * E(2, 2);
      Sh(0, 1) += 2.0 * output.hvalue[i] * E(0, 1);
      Sh(0, 2) += 2.0 * output.hvalue[i] * E(0, 2);
      Sh(1, 2) += 2.0 * output.hvalue[i] * E(1, 2);
    }
  }

  A2DMat<N, SymmMat<ScalarType, 3>>& SObj;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType>
A2D_INLINE_FUNCTION A2DSymm3x3SymmMultTraceExpr<N, ScalarType>
SymmSymmMultTrace(A2DMat<N, SymmMat<ScalarType, 3>>& S,
                  A2DMat<N, SymmMat<ScalarType, 3>>& E,
                  A2DScalar<N, ScalarType>& trace) {
  return A2DSymm3x3SymmMultTraceExpr<N, ScalarType>(S, E, trace);
}

template <class ScalarType>
A2D_INLINE_FUNCTION void SymmIsotropicConstitutive(
    const ScalarType& mu, const ScalarType& lambda,
    const SymmMat<ScalarType, 3>& E, SymmMat<ScalarType, 3>& S) {
  ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
  ScalarType mu2 = 2.0 * mu;
  S(0, 0) = mu2 * E(0, 0) + tr;
  S(0, 1) = mu2 * E(0, 1);
  S(0, 2) = mu2 * E(0, 2);
  S(1, 1) = mu2 * E(1, 1) + tr;
  S(1, 2) = mu2 * E(1, 2);
  S(2, 2) = mu2 * E(2, 2) + tr;
}

template <class ScalarType>
class ADSymm3x3IsotropicConstitutiveExpr
    : public ADExpression<ADSymm3x3IsotropicConstitutiveExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADSymm3x3IsotropicConstitutiveExpr(
      const ScalarType& mu, const ScalarType& lambda,
      ADMat<SymmMat<ScalarType, 3>>& EObj, ADMat<SymmMat<ScalarType, 3>>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(0, 2) = mu2 * E(0, 2);
    S(1, 1) = mu2 * E(1, 1) + tr;
    S(1, 2) = mu2 * E(1, 2);
    S(2, 2) = mu2 * E(2, 2) + tr;
  }

  A2D_INLINE_FUNCTION void forward() {
    const SymmMat<ScalarType, 3>& Ed = EObj.bvalue();
    SymmMat<ScalarType, 3>& Sd = SObj.bvalue();

    ScalarType tr = lambda * (Ed(0, 0) + Ed(1, 1) + Ed(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Sd(0, 0) = mu2 * Ed(0, 0) + tr;
    Sd(0, 1) = mu2 * Ed(0, 1);
    Sd(0, 2) = mu2 * Ed(0, 2);
    Sd(1, 1) = mu2 * Ed(1, 1) + tr;
    Sd(1, 2) = mu2 * Ed(1, 2);
    Sd(2, 2) = mu2 * Ed(2, 2) + tr;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& Sb = SObj.bvalue();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1) + Sb(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(0, 2) += mu2 * Sb(0, 2);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
    Eb(1, 2) += mu2 * Sb(1, 2);
    Eb(2, 2) += mu2 * Sb(2, 2) + tr;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
  ADMat<SymmMat<ScalarType, 3>>& SObj;
};

template <class ScalarType>
A2D_INLINE_FUNCTION ADSymm3x3IsotropicConstitutiveExpr<ScalarType>
SymmIsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                          ADMat<SymmMat<ScalarType, 3>>& E,
                          ADMat<SymmMat<ScalarType, 3>>& S) {
  return ADSymm3x3IsotropicConstitutiveExpr<ScalarType>(mu, lambda, E, S);
}

template <int N, class ScalarType>
class A2DSymm3x3IsotropicConstitutiveExpr
    : public A2DExpression<A2DSymm3x3IsotropicConstitutiveExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DSymm3x3IsotropicConstitutiveExpr(
      const ScalarType& mu, const ScalarType& lambda,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj,
      A2DMat<N, SymmMat<ScalarType, 3>>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(0, 2) = mu2 * E(0, 2);
    S(1, 1) = mu2 * E(1, 1) + tr;
    S(1, 2) = mu2 * E(1, 2);
    S(2, 2) = mu2 * E(2, 2) + tr;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& Sb = SObj.bvalue();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1) + Sb(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(0, 2) += mu2 * Sb(0, 2);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
    Eb(1, 2) += mu2 * Sb(1, 2);
    Eb(2, 2) += mu2 * Sb(2, 2) + tr;
  }

  A2D_INLINE_FUNCTION void hforward() {
    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      SymmMat<ScalarType, 3>& Sp = SObj.pvalue(i);

      ScalarType tr = lambda * (Ep(0, 0) + Ep(1, 1) + Ep(2, 2));
      ScalarType mu2 = 2.0 * mu;
      Sp(0, 0) = mu2 * Ep(0, 0) + tr;
      Sp(0, 1) = mu2 * Ep(0, 1);
      Sp(0, 2) = mu2 * Ep(0, 2);
      Sp(1, 1) = mu2 * Ep(1, 1) + tr;
      Sp(1, 2) = mu2 * Ep(1, 2);
      Sp(2, 2) = mu2 * Ep(2, 2) + tr;
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Sh = SObj.hvalue(i);
      SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);

      ScalarType tr = lambda * (Sh(0, 0) + Sh(1, 1) + Sh(2, 2));
      ScalarType mu2 = 2.0 * mu;
      Eh(0, 0) += mu2 * Sh(0, 0) + tr;
      Eh(0, 1) += mu2 * Sh(0, 1);
      Eh(0, 2) += mu2 * Sh(0, 2);
      Eh(1, 1) += mu2 * Sh(1, 1) + tr;
      Eh(1, 2) += mu2 * Sh(1, 2);
      Eh(2, 2) += mu2 * Sh(2, 2) + tr;
    }
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
  A2DMat<N, SymmMat<ScalarType, 3>>& SObj;
};

template <int N, class ScalarType>
A2D_INLINE_FUNCTION A2DSymm3x3IsotropicConstitutiveExpr<N, ScalarType>
SymmIsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                          A2DMat<N, SymmMat<ScalarType, 3>>& E,
                          A2DMat<N, SymmMat<ScalarType, 3>>& S) {
  return A2DSymm3x3IsotropicConstitutiveExpr<N, ScalarType>(mu, lambda, E, S);
}

template <class ScalarType>
A2D_INLINE_FUNCTION void SymmIsotropicEnergy(const ScalarType& mu,
                                             const ScalarType& lambda,
                                             const SymmMat<ScalarType, 3>& E,
                                             ScalarType& output) {
  ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
  ScalarType trE =
      E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
      2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

  output = mu * trE + 0.5 * lambda * tr * tr;
}

template <class ScalarType>
class ADSymm3x3IsotropicEnergyExpr
    : public ADExpression<ADSymm3x3IsotropicEnergyExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADSymm3x3IsotropicEnergyExpr(
      const ScalarType& mu, const ScalarType& lambda,
      ADMat<SymmMat<ScalarType, 3>>& EObj, ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  A2D_INLINE_FUNCTION void forward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const SymmMat<ScalarType, 3>& Ed = EObj.bvalue();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trd = (Ed(0, 0) + Ed(1, 1) + Ed(2, 2));
    ScalarType trEd =
        2.0 *
        (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) + E(2, 2) * Ed(2, 2) +
         2.0 * (E(0, 1) * Ed(0, 1) + E(0, 2) * Ed(0, 2) + E(1, 2) * Ed(1, 2)));

    output.bvalue = mu * trEd + lambda * tr * trd;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;
    Eb(2, 2) += (mu2 * E(2, 2) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
    Eb(0, 2) += 2.0 * mu2 * E(0, 2) * output.bvalue;
    Eb(1, 2) += 2.0 * mu2 * E(1, 2) * output.bvalue;
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
A2D_INLINE_FUNCTION ADSymm3x3IsotropicEnergyExpr<ScalarType>
SymmIsotropicEnergy(const ScalarType& mu, const ScalarType& lambda,
                    ADMat<SymmMat<ScalarType, 3>>& E,
                    ADScalar<ScalarType>& output) {
  return ADSymm3x3IsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
class ADSymm3x3ADIsotropicEnergyExpr
    : public ADExpression<ADSymm3x3ADIsotropicEnergyExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADSymm3x3ADIsotropicEnergyExpr(
      ADScalar<ScalarType>& mu, ADScalar<ScalarType>& lambda,
      ADMat<SymmMat<ScalarType, 3>>& EObj, ADScalar<ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  A2D_INLINE_FUNCTION void forward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const SymmMat<ScalarType, 3>& Ed = EObj.bvalue();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trd = (Ed(0, 0) + Ed(1, 1) + Ed(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));
    ScalarType trEd =
        2.0 *
        (E(0, 0) * Ed(0, 0) + E(1, 1) * Ed(1, 1) + E(2, 2) * Ed(2, 2) +
         2.0 * (E(0, 1) * Ed(0, 1) + E(0, 2) * Ed(0, 2) + E(1, 2) * Ed(1, 2)));

    output.bvalue = mu.value * trEd + lambda.value * tr * trd +
                    mu.bvalue * trE + 0.5 * lambda.bvalue * tr * tr;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;
    Eb(2, 2) += (mu2 * E(2, 2) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
    Eb(0, 2) += 2.0 * mu2 * E(0, 2) * output.bvalue;
    Eb(1, 2) += 2.0 * mu2 * E(1, 2) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  ADScalar<ScalarType>& mu;
  ADScalar<ScalarType>& lambda;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
  ADScalar<ScalarType>& output;
};

template <class ScalarType>
A2D_INLINE_FUNCTION ADSymm3x3ADIsotropicEnergyExpr<ScalarType>
SymmIsotropicEnergy(ADScalar<ScalarType>& mu, ADScalar<ScalarType>& lambda,
                    ADMat<SymmMat<ScalarType, 3>>& E,
                    ADScalar<ScalarType>& output) {
  return ADSymm3x3ADIsotropicEnergyExpr<ScalarType>(mu, lambda, E, output);
}

template <int N, class ScalarType>
class A2DSymm3x3IsotropicEnergyExpr
    : public A2DExpression<A2DSymm3x3IsotropicEnergyExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DSymm3x3IsotropicEnergyExpr(
      const ScalarType& mu, const ScalarType& lambda,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj, A2DScalar<N, ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    output.value = mu * trE + 0.5 * lambda * tr * tr;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu;

    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    Eb(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.bvalue;
    Eb(2, 2) += (mu2 * E(2, 2) + lambda * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
    Eb(0, 2) += 2.0 * mu2 * E(0, 2) * output.bvalue;
    Eb(1, 2) += 2.0 * mu2 * E(1, 2) * output.bvalue;
  }

  A2D_INLINE_FUNCTION void hforward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      ScalarType trd = (Ep(0, 0) + Ep(1, 1) + Ep(2, 2));
      ScalarType trEd =
          2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) + E(2, 2) * Ep(2, 2) +
                 2.0 * (E(0, 1) * Ep(0, 1) + E(0, 2) * Ep(0, 2) +
                        E(1, 2) * Ep(1, 2)));

      output.pvalue[i] = mu * trEd + lambda * tr * trd;
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu;
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);

      // by * (d^2y/dx^2 * px)
      ScalarType trp = (Ep(0, 0) + Ep(1, 1) + Ep(2, 2));
      Eh(0, 0) += (mu2 * Ep(0, 0) + lambda * trp) * output.bvalue;
      Eh(1, 1) += (mu2 * Ep(1, 1) + lambda * trp) * output.bvalue;
      Eh(2, 2) += (mu2 * Ep(2, 2) + lambda * trp) * output.bvalue;

      Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;
      Eh(0, 2) += 2.0 * mu2 * Ep(0, 2) * output.bvalue;
      Eh(1, 2) += 2.0 * mu2 * Ep(1, 2) * output.bvalue;

      // hy * (dy/dx)
      Eh(0, 0) += (mu2 * E(0, 0) + lambda * tr) * output.hvalue[i];
      Eh(1, 1) += (mu2 * E(1, 1) + lambda * tr) * output.hvalue[i];
      Eh(2, 2) += (mu2 * E(2, 2) + lambda * tr) * output.hvalue[i];

      Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue[i];
      Eh(0, 2) += 2.0 * mu2 * E(0, 2) * output.hvalue[i];
      Eh(1, 2) += 2.0 * mu2 * E(1, 2) * output.hvalue[i];
    }
  }

  const ScalarType& mu;
  const ScalarType& lambda;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType>
A2D_INLINE_FUNCTION A2DSymm3x3IsotropicEnergyExpr<N, ScalarType>
SymmIsotropicEnergy(const ScalarType& mu, const ScalarType& lambda,
                    A2DMat<N, SymmMat<ScalarType, 3>>& E,
                    A2DScalar<N, ScalarType>& output) {
  return A2DSymm3x3IsotropicEnergyExpr<N, ScalarType>(mu, lambda, E, output);
}

template <int N, class ScalarType>
class A2DSymm3x3A2DIsotropicEnergyExpr
    : public A2DExpression<A2DSymm3x3A2DIsotropicEnergyExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DSymm3x3A2DIsotropicEnergyExpr(
      A2DScalar<N, ScalarType>& mu, A2DScalar<N, ScalarType>& lambda,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj, A2DScalar<N, ScalarType>& output)
      : mu(mu), lambda(lambda), EObj(EObj), output(output) {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    output.value = mu.value * trE + 0.5 * lambda.value * tr * tr;
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    const ScalarType mu2 = 2.0 * mu.value;

    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    Eb(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.bvalue;
    Eb(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.bvalue;
    Eb(2, 2) += (mu2 * E(2, 2) + lambda.value * tr) * output.bvalue;

    Eb(0, 1) += 2.0 * mu2 * E(0, 1) * output.bvalue;
    Eb(0, 2) += 2.0 * mu2 * E(0, 2) * output.bvalue;
    Eb(1, 2) += 2.0 * mu2 * E(1, 2) * output.bvalue;

    mu.bvalue = trE * output.bvalue;
    lambda.bvalue = 0.5 * tr * tr * output.bvalue;
  }

  A2D_INLINE_FUNCTION void hforward() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      ScalarType trd = (Ep(0, 0) + Ep(1, 1) + Ep(2, 2));
      ScalarType trEd =
          2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) + E(2, 2) * Ep(2, 2) +
                 2.0 * (E(0, 1) * Ep(0, 1) + E(0, 2) * Ep(0, 2) +
                        E(1, 2) * Ep(1, 2)));

      output.pvalue[i] = mu.value * trEd + lambda.value * tr * trd +
                         mu.pvalue[i] * trE + 0.5 * lambda.pvalue[i] * tr * tr;
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const SymmMat<ScalarType, 3>& E = EObj.value();
    const ScalarType mu2 = 2.0 * mu.value;
    ScalarType tr = (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType trE =
        E(0, 0) * E(0, 0) + E(1, 1) * E(1, 1) + E(2, 2) * E(2, 2) +
        2.0 * (E(0, 1) * E(0, 1) + E(0, 2) * E(0, 2) + E(1, 2) * E(1, 2));

    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);

      // by * (d^2y/dx^2 * px)
      ScalarType trp = (Ep(0, 0) + Ep(1, 1) + Ep(2, 2));
      ScalarType trEp =
          2.0 * (E(0, 0) * Ep(0, 0) + E(1, 1) * Ep(1, 1) + E(2, 2) * Ep(2, 2) +
                 2.0 * (E(0, 1) * Ep(0, 1) + E(0, 2) * Ep(0, 2) +
                        E(1, 2) * Ep(1, 2)));

      Eh(0, 0) += (mu2 * Ep(0, 0) + lambda.value * trp) * output.bvalue;
      Eh(1, 1) += (mu2 * Ep(1, 1) + lambda.value * trp) * output.bvalue;
      Eh(2, 2) += (mu2 * Ep(2, 2) + lambda.value * trp) * output.bvalue;

      Eh(0, 1) += 2.0 * mu2 * Ep(0, 1) * output.bvalue;
      Eh(0, 2) += 2.0 * mu2 * Ep(0, 2) * output.bvalue;
      Eh(1, 2) += 2.0 * mu2 * Ep(1, 2) * output.bvalue;

      mu.hvalue[i] += trEp * output.bvalue;
      lambda.hvalue[i] += tr * trp * output.bvalue;

      // hy * (dy/dx)
      Eh(0, 0) += (mu2 * E(0, 0) + lambda.value * tr) * output.hvalue[i];
      Eh(1, 1) += (mu2 * E(1, 1) + lambda.value * tr) * output.hvalue[i];
      Eh(2, 2) += (mu2 * E(2, 2) + lambda.value * tr) * output.hvalue[i];

      Eh(0, 1) += 2.0 * mu2 * E(0, 1) * output.hvalue[i];
      Eh(0, 2) += 2.0 * mu2 * E(0, 2) * output.hvalue[i];
      Eh(1, 2) += 2.0 * mu2 * E(1, 2) * output.hvalue[i];

      mu.hvalue[i] += trE * output.hvalue[i];
      lambda.hvalue[i] += 0.5 * tr * tr * output.hvalue[i];

      // account for Hessian blocks w.r.t. E and mu or lambda
      Eh(0, 0) += (2.0 * E(0, 0) * mu.pvalue[i] + tr * lambda.pvalue[i]) *
                  output.bvalue;
      Eh(1, 1) += (2.0 * E(1, 1) * mu.pvalue[i] + tr * lambda.pvalue[i]) *
                  output.bvalue;
      Eh(2, 2) += (2.0 * E(2, 2) * mu.pvalue[i] + tr * lambda.pvalue[i]) *
                  output.bvalue;
      Eh(0, 1) += 4.0 * E(0, 1) * mu.pvalue[i] * output.bvalue;
      Eh(0, 2) += 4.0 * E(0, 2) * mu.pvalue[i] * output.bvalue;
      Eh(1, 2) += 4.0 * E(1, 2) * mu.pvalue[i] * output.bvalue;
    }
  }

  A2DScalar<N, ScalarType>& mu;
  A2DScalar<N, ScalarType>& lambda;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class ScalarType>
A2D_INLINE_FUNCTION A2DSymm3x3A2DIsotropicEnergyExpr<N, ScalarType>
SymmIsotropicEnergy(A2DScalar<N, ScalarType>& mu,
                    A2DScalar<N, ScalarType>& lambda,
                    A2DMat<N, SymmMat<ScalarType, 3>>& E,
                    A2DScalar<N, ScalarType>& output) {
  return A2DSymm3x3A2DIsotropicEnergyExpr<N, ScalarType>(mu, lambda, E, output);
}

template <class ScalarType>
A2D_INLINE_FUNCTION void MatGreenStrain(const Mat<ScalarType, 3, 3>& Ux,
                                        SymmMat<ScalarType, 3>& E) {
  E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0) +
                              Ux(2, 0) * Ux(2, 0));
  E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1) +
                              Ux(2, 1) * Ux(2, 1));
  E(2, 2) = Ux(2, 2) + 0.5 * (Ux(0, 2) * Ux(0, 2) + Ux(1, 2) * Ux(1, 2) +
                              Ux(2, 2) * Ux(2, 2));

  E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) +
                   Ux(1, 0) * Ux(1, 1) + Ux(2, 0) * Ux(2, 1));
  E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0) + Ux(0, 0) * Ux(0, 2) +
                   Ux(1, 0) * Ux(1, 2) + Ux(2, 0) * Ux(2, 2));
  E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1) + Ux(0, 1) * Ux(0, 2) +
                   Ux(1, 1) * Ux(1, 2) + Ux(2, 1) * Ux(2, 2));
}

template <typename ScalarType>
class ADMat3x3GreenStrainExpr
    : public ADExpression<ADMat3x3GreenStrainExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADMat3x3GreenStrainExpr(
      ADMat<Mat<ScalarType, 3, 3>>& UxObj, ADMat<SymmMat<ScalarType, 3>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    SymmMat<ScalarType, 3>& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0) +
                                Ux(2, 0) * Ux(2, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1) +
                                Ux(2, 1) * Ux(2, 1));
    E(2, 2) = Ux(2, 2) + 0.5 * (Ux(0, 2) * Ux(0, 2) + Ux(1, 2) * Ux(1, 2) +
                                Ux(2, 2) * Ux(2, 2));

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) +
                     Ux(1, 0) * Ux(1, 1) + Ux(2, 0) * Ux(2, 1));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0) + Ux(0, 0) * Ux(0, 2) +
                     Ux(1, 0) * Ux(1, 2) + Ux(2, 0) * Ux(2, 2));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1) + Ux(0, 1) * Ux(0, 2) +
                     Ux(1, 1) * Ux(1, 2) + Ux(2, 1) * Ux(2, 2));
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    const Mat<ScalarType, 3, 3>& Uxd = UxObj.bvalue();
    SymmMat<ScalarType, 3>& Ed = EObj.bvalue();

    Ed(0, 0) = Uxd(0, 0) + Ux(0, 0) * Uxd(0, 0) + Ux(1, 0) * Uxd(1, 0) +
               Ux(2, 0) * Uxd(2, 0);
    Ed(1, 1) = Uxd(1, 1) + Ux(0, 1) * Uxd(0, 1) + Ux(1, 1) * Uxd(1, 1) +
               Ux(2, 1) * Uxd(2, 1);
    Ed(2, 2) = Uxd(2, 2) + Ux(0, 2) * Uxd(0, 2) + Ux(1, 2) * Uxd(1, 2) +
               Ux(2, 2) * Uxd(2, 2);

    Ed(0, 1) = 0.5 * (Uxd(0, 1) + Uxd(1, 0) + Ux(0, 0) * Uxd(0, 1) +
                      Ux(1, 0) * Uxd(1, 1) + Ux(2, 0) * Uxd(2, 1) +
                      Uxd(0, 0) * Ux(0, 1) + Uxd(1, 0) * Ux(1, 1) +
                      Uxd(2, 0) * Ux(2, 1));
    Ed(0, 2) = 0.5 * (Uxd(0, 2) + Uxd(2, 0) + Ux(0, 0) * Uxd(0, 2) +
                      Ux(1, 0) * Uxd(1, 2) + Ux(2, 0) * Uxd(2, 2) +
                      Uxd(0, 0) * Ux(0, 2) + Uxd(1, 0) * Ux(1, 2) +
                      Uxd(2, 0) * Ux(2, 2));
    Ed(1, 2) = 0.5 * (Uxd(1, 2) + Uxd(2, 1) + Ux(0, 1) * Uxd(0, 2) +
                      Ux(1, 1) * Uxd(1, 2) + Ux(2, 1) * Uxd(2, 2) +
                      Uxd(0, 1) * Ux(0, 2) + Uxd(1, 1) * Ux(1, 2) +
                      Uxd(2, 1) * Ux(2, 2));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    const SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    Mat<ScalarType, 3, 3>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) +
                 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1) +
                 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +
                 Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) +
                 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1) +
                 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +
                 Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) += Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) + Ux(2, 1) * Eb(1, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +
                 (Ux(2, 2) + 1.0) * Eb(2, 2);
  }

  ADMat<Mat<ScalarType, 3, 3>>& UxObj;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
};

template <typename ScalarType>
A2D_INLINE_FUNCTION ADMat3x3GreenStrainExpr<ScalarType> MatGreenStrain(
    ADMat<Mat<ScalarType, 3, 3>>& Ux, ADMat<SymmMat<ScalarType, 3>>& E) {
  return ADMat3x3GreenStrainExpr<ScalarType>(Ux, E);
}

template <int N, typename ScalarType>
class A2DMat3x3GreenStrainExpr
    : public A2DExpression<A2DMat3x3GreenStrainExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DMat3x3GreenStrainExpr(
      A2DMat<N, Mat<ScalarType, 3, 3>>& UxObj,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    SymmMat<ScalarType, 3>& E = EObj.value();
    E(0, 0) = Ux(0, 0) + 0.5 * (Ux(0, 0) * Ux(0, 0) + Ux(1, 0) * Ux(1, 0) +
                                Ux(2, 0) * Ux(2, 0));
    E(1, 1) = Ux(1, 1) + 0.5 * (Ux(0, 1) * Ux(0, 1) + Ux(1, 1) * Ux(1, 1) +
                                Ux(2, 1) * Ux(2, 1));
    E(2, 2) = Ux(2, 2) + 0.5 * (Ux(0, 2) * Ux(0, 2) + Ux(1, 2) * Ux(1, 2) +
                                Ux(2, 2) * Ux(2, 2));

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0) + Ux(0, 0) * Ux(0, 1) +
                     Ux(1, 0) * Ux(1, 1) + Ux(2, 0) * Ux(2, 1));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0) + Ux(0, 0) * Ux(0, 2) +
                     Ux(1, 0) * Ux(1, 2) + Ux(2, 0) * Ux(2, 2));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1) + Ux(0, 1) * Ux(0, 2) +
                     Ux(1, 1) * Ux(1, 2) + Ux(2, 1) * Ux(2, 2));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    const SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    Mat<ScalarType, 3, 3>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += (Ux(0, 0) + 1.0) * Eb(0, 0) + 0.5 * Ux(0, 1) * Eb(0, 1) +
                 0.5 * Ux(0, 2) * Eb(0, 2);
    Uxb(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 1) + Ux(0, 1) * Eb(1, 1) +
                 0.5 * Ux(0, 2) * Eb(1, 2);
    Uxb(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eb(0, 2) + 0.5 * Ux(0, 1) * Eb(1, 2) +
                 Ux(0, 2) * Eb(2, 2);

    Uxb(1, 0) += Ux(1, 0) * Eb(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eb(0, 1) +
                 0.5 * Ux(1, 2) * Eb(0, 2);
    Uxb(1, 1) += 0.5 * Ux(1, 0) * Eb(0, 1) + (Ux(1, 1) + 1.0) * Eb(1, 1) +
                 0.5 * Ux(1, 2) * Eb(1, 2);
    Uxb(1, 2) += 0.5 * Ux(1, 0) * Eb(0, 2) + 0.5 * (Ux(1, 1) + 1.0) * Eb(1, 2) +
                 Ux(1, 2) * Eb(2, 2);

    Uxb(2, 0) += Ux(2, 0) * Eb(0, 0) + 0.5 * Ux(2, 1) * Eb(0, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Ux(2, 0) * Eb(0, 1) + Ux(2, 1) * Eb(1, 1) +
                 0.5 * (Ux(2, 2) + 1.0) * Eb(1, 2);
    Uxb(2, 2) += 0.5 * Ux(2, 0) * Eb(0, 2) + 0.5 * Ux(2, 1) * Eb(1, 2) +
                 (Ux(2, 2) + 1.0) * Eb(2, 2);
  }

  A2D_INLINE_FUNCTION void hforward() {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Uxp = UxObj.pvalue(i);
      SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);
      Ep(0, 0) = Uxp(0, 0) + Ux(0, 0) * Uxp(0, 0) + Ux(1, 0) * Uxp(1, 0) +
                 Ux(2, 0) * Uxp(2, 0);
      Ep(1, 1) = Uxp(1, 1) + Ux(0, 1) * Uxp(0, 1) + Ux(1, 1) * Uxp(1, 1) +
                 Ux(2, 1) * Uxp(2, 1);
      Ep(2, 2) = Uxp(2, 2) + Ux(0, 2) * Uxp(0, 2) + Ux(1, 2) * Uxp(1, 2) +
                 Ux(2, 2) * Uxp(2, 2);

      Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0) + Ux(0, 0) * Uxp(0, 1) +
                        Ux(1, 0) * Uxp(1, 1) + Ux(2, 0) * Uxp(2, 1) +
                        Uxp(0, 0) * Ux(0, 1) + Uxp(1, 0) * Ux(1, 1) +
                        Uxp(2, 0) * Ux(2, 1));
      Ep(0, 2) = 0.5 * (Uxp(0, 2) + Uxp(2, 0) + Ux(0, 0) * Uxp(0, 2) +
                        Ux(1, 0) * Uxp(1, 2) + Ux(2, 0) * Uxp(2, 2) +
                        Uxp(0, 0) * Ux(0, 2) + Uxp(1, 0) * Ux(1, 2) +
                        Uxp(2, 0) * Ux(2, 2));
      Ep(1, 2) = 0.5 * (Uxp(1, 2) + Uxp(2, 1) + Ux(0, 1) * Uxp(0, 2) +
                        Ux(1, 1) * Uxp(1, 2) + Ux(2, 1) * Uxp(2, 2) +
                        Uxp(0, 1) * Ux(0, 2) + Uxp(1, 1) * Ux(1, 2) +
                        Uxp(2, 1) * Ux(2, 2));
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    const Mat<ScalarType, 3, 3>& Eb = EObj.bvalue();
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Uxp = UxObj.pvalue(i);
      const SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);
      Mat<ScalarType, 3, 3>& Uxh = UxObj.hvalue(i);

      Uxh(0, 0) += Uxp(0, 0) * Eb(0, 0) + 0.5 * Uxp(0, 1) * Eb(0, 1) +
                   0.5 * Uxp(0, 2) * Eb(0, 2);
      Uxh(0, 1) += 0.5 * Uxp(0, 0) * Eb(0, 1) + Uxp(0, 1) * Eb(1, 1) +
                   0.5 * Uxp(0, 2) * Eb(1, 2);
      Uxh(0, 2) += 0.5 * Uxp(0, 0) * Eb(0, 2) + 0.5 * Uxp(0, 1) * Eb(1, 2) +
                   Uxp(0, 2) * Eb(2, 2);

      Uxh(1, 0) += Uxp(1, 0) * Eb(0, 0) + 0.5 * Uxp(1, 1) * Eb(0, 1) +
                   0.5 * Uxp(1, 2) * Eb(0, 2);
      Uxh(1, 1) += 0.5 * Uxp(1, 0) * Eb(0, 1) + Uxp(1, 1) * Eb(1, 1) +
                   0.5 * Uxp(1, 2) * Eb(1, 2);
      Uxh(1, 2) += 0.5 * Uxp(1, 0) * Eb(0, 2) + 0.5 * Uxp(1, 1) * Eb(1, 2) +
                   Uxp(1, 2) * Eb(2, 2);

      Uxh(2, 0) += Uxp(2, 0) * Eb(0, 0) + 0.5 * Uxp(2, 1) * Eb(0, 1) +
                   0.5 * Uxp(2, 2) * Eb(0, 2);
      Uxh(2, 1) += 0.5 * Uxp(2, 0) * Eb(0, 1) + Uxp(2, 1) * Eb(1, 1) +
                   0.5 * Uxp(2, 2) * Eb(1, 2);
      Uxh(2, 2) += 0.5 * Uxp(2, 0) * Eb(0, 2) + 0.5 * Uxp(2, 1) * Eb(1, 2) +
                   Uxp(2, 2) * Eb(2, 2);

      Uxh(0, 0) += (Ux(0, 0) + 1.0) * Eh(0, 0) + 0.5 * Ux(0, 1) * Eh(0, 1) +
                   0.5 * Ux(0, 2) * Eh(0, 2);
      Uxh(0, 1) += 0.5 * (Ux(0, 0) + 1.0) * Eh(0, 1) + Ux(0, 1) * Eh(1, 1) +
                   0.5 * Ux(0, 2) * Eh(1, 2);
      Uxh(0, 2) += 0.5 * (Ux(0, 0) + 1.0) * Eh(0, 2) +
                   0.5 * Ux(0, 1) * Eh(1, 2) + Ux(0, 2) * Eh(2, 2);

      Uxh(1, 0) += Ux(1, 0) * Eh(0, 0) + 0.5 * (Ux(1, 1) + 1.0) * Eh(0, 1) +
                   0.5 * Ux(1, 2) * Eh(0, 2);
      Uxh(1, 1) += 0.5 * Ux(1, 0) * Eh(0, 1) + (Ux(1, 1) + 1.0) * Eh(1, 1) +
                   0.5 * Ux(1, 2) * Eh(1, 2);
      Uxh(1, 2) += 0.5 * Ux(1, 0) * Eh(0, 2) +
                   0.5 * (Ux(1, 1) + 1.0) * Eh(1, 2) + Ux(1, 2) * Eh(2, 2);

      Uxh(2, 0) += Ux(2, 0) * Eh(0, 0) + 0.5 * Ux(2, 1) * Eh(0, 1) +
                   0.5 * (Ux(2, 2) + 1.0) * Eh(0, 2);
      Uxh(2, 1) += 0.5 * Ux(2, 0) * Eh(0, 1) + Ux(2, 1) * Eh(1, 1) +
                   0.5 * (Ux(2, 2) + 1.0) * Eh(1, 2);
      Uxh(2, 2) += 0.5 * Ux(2, 0) * Eh(0, 2) + 0.5 * Ux(2, 1) * Eh(1, 2) +
                   (Ux(2, 2) + 1.0) * Eh(2, 2);
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& UxObj;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
};

template <int N, typename ScalarType>
A2D_INLINE_FUNCTION A2DMat3x3GreenStrainExpr<N, ScalarType> MatGreenStrain(
    A2DMat<N, Mat<ScalarType, 3, 3>>& Ux,
    A2DMat<N, SymmMat<ScalarType, 3>>& E) {
  return A2DMat3x3GreenStrainExpr<N, ScalarType>(Ux, E);
}

template <class ScalarType>
A2D_INLINE_FUNCTION void MatLinearGreenStrain(const Mat<ScalarType, 3, 3>& Ux,
                                              SymmMat<ScalarType, 3>& E) {
  E(0, 0) = Ux(0, 0);
  E(1, 1) = Ux(1, 1);
  E(2, 2) = Ux(2, 2);

  E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
  E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
}

template <typename ScalarType>
class ADMat3x3LinearGreenStrainExpr
    : public ADExpression<ADMat3x3LinearGreenStrainExpr<ScalarType>> {
 public:
  A2D_INLINE_FUNCTION ADMat3x3LinearGreenStrainExpr(
      ADMat<Mat<ScalarType, 3, 3>>& UxObj, ADMat<SymmMat<ScalarType, 3>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    SymmMat<ScalarType, 3>& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);
    E(2, 2) = Ux(2, 2);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
  }

  A2D_INLINE_FUNCTION void forward() {
    const Mat<ScalarType, 3, 3>& Uxb = UxObj.bvalue();
    SymmMat<ScalarType, 3>& Eb = EObj.bvalue();

    Eb(0, 0) = Uxb(0, 0);
    Eb(1, 1) = Uxb(1, 1);
    Eb(2, 2) = Uxb(2, 2);

    Eb(0, 1) = 0.5 * (Uxb(0, 1) + Uxb(1, 0));
    Eb(0, 2) = 0.5 * (Uxb(0, 2) + Uxb(2, 0));
    Eb(1, 2) = 0.5 * (Uxb(1, 2) + Uxb(2, 1));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    Mat<ScalarType, 3, 3>& Uxb = UxObj.bvalue();

    // Uxb = (I + Ux) * Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);
    Uxb(0, 2) += 0.5 * Eb(0, 2);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
    Uxb(1, 2) += 0.5 * Eb(1, 2);

    Uxb(2, 0) += 0.5 * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Eb(1, 2);
    Uxb(2, 2) += Eb(2, 2);
  }

  ADMat<Mat<ScalarType, 3, 3>>& UxObj;
  ADMat<SymmMat<ScalarType, 3>>& EObj;
};

template <typename ScalarType>
A2D_INLINE_FUNCTION ADMat3x3LinearGreenStrainExpr<ScalarType>
MatLinearGreenStrain(ADMat<Mat<ScalarType, 3, 3>>& Ux,
                     ADMat<SymmMat<ScalarType, 3>>& E) {
  return ADMat3x3LinearGreenStrainExpr<ScalarType>(Ux, E);
}

template <int N, typename ScalarType>
class A2DMat3x3LinearGreenStrainExpr
    : public A2DExpression<A2DMat3x3LinearGreenStrainExpr<N, ScalarType>> {
 public:
  A2D_INLINE_FUNCTION A2DMat3x3LinearGreenStrainExpr(
      A2DMat<N, Mat<ScalarType, 3, 3>>& UxObj,
      A2DMat<N, SymmMat<ScalarType, 3>>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const Mat<ScalarType, 3, 3>& Ux = UxObj.value();
    SymmMat<ScalarType, 3>& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);
    E(2, 2) = Ux(2, 2);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
  }

  A2D_INLINE_FUNCTION void reverse() {
    const SymmMat<ScalarType, 3>& Eb = EObj.bvalue();
    Mat<ScalarType, 3, 3>& Uxb = UxObj.bvalue();

    // Uxb = Eb
    Uxb(0, 0) += Eb(0, 0);
    Uxb(0, 1) += 0.5 * Eb(0, 1);
    Uxb(0, 2) += 0.5 * Eb(0, 2);

    Uxb(1, 0) += 0.5 * Eb(0, 1);
    Uxb(1, 1) += Eb(1, 1);
    Uxb(1, 2) += 0.5 * Eb(1, 2);

    Uxb(2, 0) += 0.5 * Eb(0, 2);
    Uxb(2, 1) += 0.5 * Eb(1, 2);
    Uxb(2, 2) += Eb(2, 2);
  }

  A2D_INLINE_FUNCTION void hforward() {
    for (int i = 0; i < N; i++) {
      const Mat<ScalarType, 3, 3>& Uxp = UxObj.pvalue(i);
      SymmMat<ScalarType, 3>& Ep = EObj.pvalue(i);

      Ep(0, 0) = Uxp(0, 0);
      Ep(1, 1) = Uxp(1, 1);
      Ep(2, 2) = Uxp(2, 2);

      Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0));
      Ep(0, 2) = 0.5 * (Uxp(0, 2) + Uxp(2, 0));
      Ep(1, 2) = 0.5 * (Uxp(1, 2) + Uxp(2, 1));
    }
  }

  A2D_INLINE_FUNCTION void hreverse() {
    for (int i = 0; i < N; i++) {
      const SymmMat<ScalarType, 3>& Eh = EObj.hvalue(i);
      Mat<ScalarType, 3, 3>& Uxh = UxObj.hvalue(i);

      Uxh(0, 0) += Eh(0, 0);
      Uxh(0, 1) += 0.5 * Eh(0, 1);
      Uxh(0, 2) += 0.5 * Eh(0, 2);

      Uxh(1, 0) += 0.5 * Eh(0, 1);
      Uxh(1, 1) += Eh(1, 1);
      Uxh(1, 2) += 0.5 * Eh(1, 2);

      Uxh(2, 0) += 0.5 * Eh(0, 2);
      Uxh(2, 1) += 0.5 * Eh(1, 2);
      Uxh(2, 2) += Eh(2, 2);
    }
  }

  A2DMat<N, Mat<ScalarType, 3, 3>>& UxObj;
  A2DMat<N, SymmMat<ScalarType, 3>>& EObj;
};

template <int N, typename ScalarType>
A2D_INLINE_FUNCTION A2DMat3x3LinearGreenStrainExpr<N, ScalarType>
MatLinearGreenStrain(A2DMat<N, Mat<ScalarType, 3, 3>>& Ux,
                     A2DMat<N, SymmMat<ScalarType, 3>>& E) {
  return A2DMat3x3LinearGreenStrainExpr<N, ScalarType>(Ux, E);
}

}  // namespace A2D

#endif  // A2D_TMP_3D_H
