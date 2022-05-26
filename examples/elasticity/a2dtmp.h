/*
  0: (0, 0)  1: (0, 1)  2: (0, 2)
  3: (1, 0)  4: (1, 1)  5: (1, 2)
  6: (2, 0)  7: (2, 1)  8: (2, 2)
*/

#ifndef A2D_TMP_H
#define A2D_TMP_H

#include <stdlib.h>

#include "a2dmatcore.h"
#include "a2dobjs.h"
#include "a2dtypes.h"

namespace A2D {

// Mat3x3MatMult
template <typename ScalarType, bool AT = false, bool BT = false>
inline void Mat3x3MatMult(const Mat<ScalarType, 3, 3>& A,
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

template <class AMatType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADMat3x3MatMultExpr
    : public ADExpression<
          ADMat3x3MatMultExpr<AMatType, BMatType, CMatType, AT, BT> > {
 public:
  ADMat3x3MatMultExpr(ADMat<AMatType>& AObj, ADMat<BMatType>& BObj,
                      ADMat<CMatType>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

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

  void forward() {
    const AMatType& A = AObj.value();
    const AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    const BMatType& Bb = BObj.bvalue();
    CMatType& Cb = CObj.bvalue();

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

  void reverse() {
    const AMatType& A = AObj.value();
    AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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

  ADMat<AMatType>& AObj;
  ADMat<BMatType>& BObj;
  ADMat<CMatType>& CObj;
};

template <class AMatType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADMat3x3MatMultExpr<AMatType, BMatType, CMatType, AT, BT> Mat3x3MatMult(
    ADMat<AMatType>& AObj, ADMat<BMatType>& BObj, ADMat<CMatType>& CObj) {
  return ADMat3x3MatMultExpr<AMatType, BMatType, CMatType, AT, BT>(AObj, BObj,
                                                                   CObj);
}

template <int N, class AMatType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DMat3x3MatMultExpr
    : public A2DExpression<
          A2DMat3x3MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT> > {
 public:
  A2DMat3x3MatMultExpr(A2DMat<N, AMatType>& AObj, A2DMat<N, BMatType>& BObj,
                       A2DMat<N, CMatType>& CObj)
      : AObj(AObj), BObj(BObj), CObj(CObj) {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

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

  void reverse() {
    const AMatType& A = AObj.value();
    AMatType& Ab = AObj.bvalue();
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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

  void hforward() {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();

    for (int i = 0; i < N; i++) {
      const AMatType& Ap = AObj.pvalue(i);
      const BMatType& Bp = BObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

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

  void hreverse() {
    const AMatType& A = AObj.value();
    const BMatType& B = BObj.value();

    for (int i = 0; i < N; i++) {
      AMatType& Ah = AObj.hvalue(i);
      BMatType& Bh = BObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);

      if (AT && BT) {
        MatTrans3x3MatTransMultAddCore(B, Ch, Ah);
        MatTrans3x3MatTransMultAddCore(Ch, A, Bh);
      } else if (AT) {
        Mat3x3MatTransMultAddCore(B, Ch, Ah);
        Mat3x3MatMultAddCore(A, Ch, Bh);
      } else if (BT) {
        Mat3x3MatMultAddCore(Ch, B, Ah);
        MatTrans3x3MatMultAddCore(Ch, A, Bh);
      } else {
        Mat3x3MatTransMultAddCore(Ch, B, Ah);
        MatTrans3x3MatMultAddCore(A, Ch, Bh);
      }
    }
  }

  A2DMat<N, AMatType>& AObj;
  A2DMat<N, BMatType>& BObj;
  A2DMat<N, CMatType>& CObj;
};

template <int N, class AMatType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DMat3x3MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT>
Mat3x3MatMult(A2DMat<N, AMatType>& AObj, A2DMat<N, BMatType>& BObj,
              A2DMat<N, CMatType>& CObj) {
  return A2DMat3x3MatMultExpr<N, AMatType, BMatType, CMatType, AT, BT>(
      AObj, BObj, CObj);
}

template <typename ScalarType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADpMat3x3MatMultExpr
    : public ADExpression<
          ADpMat3x3MatMultExpr<ScalarType, BMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  ADpMat3x3MatMultExpr(Mat3x3& A, ADMat<BMatType>& BObj, ADMat<CMatType>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

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

  void forward() {
    const BMatType& B = BObj.value();
    const BMatType& Bb = BObj.bvalue();
    CMatType& Cb = CObj.bvalue();

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

  void reverse() {
    const BMatType& B = BObj.value();
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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
  ADMat<BMatType>& BObj;
  ADMat<CMatType>& CObj;
};

template <class ScalarType, class BMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADpMat3x3MatMultExpr<ScalarType, BMatType, CMatType, AT, BT>
Mat3x3MatMult(Mat<ScalarType, 3, 3>& A, ADMat<BMatType>& BObj,
              ADMat<CMatType>& CObj) {
  return ADpMat3x3MatMultExpr<ScalarType, BMatType, CMatType, AT, BT>(A, BObj,
                                                                      CObj);
}

template <typename ScalarType, class AMatType, class CMatType, bool AT = false,
          bool BT = false>
class ADMat3x3pMatMultExpr
    : public ADExpression<
          ADMat3x3pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  ADMat3x3pMatMultExpr(ADMat<AMatType>& AObj, Mat3x3& B, ADMat<CMatType>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const AMatType& A = AObj.value();
    CMatType& C = CObj.value();

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

  void forward() {
    const AMatType& A = AObj.value();
    const AMatType& Ab = AObj.bvalue();
    CMatType& Cb = CObj.bvalue();

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

  void reverse() {
    const AMatType& A = AObj.value();
    AMatType& Ab = AObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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

  ADMat<AMatType>& AObj;
  Mat3x3& B;
  ADMat<CMatType>& CObj;
};

template <typename ScalarType, class AMatType, class CMatType, bool AT = false,
          bool BT = false>
inline ADMat3x3pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT>
Mat3x3MatMult(ADMat<AMatType>& AObj, Mat<ScalarType, 3, 3>& B,
              ADMat<CMatType>& CObj) {
  return ADMat3x3pMatMultExpr<ScalarType, AMatType, CMatType, AT, BT>(AObj, B,
                                                                      CObj);
}

template <int N, typename ScalarType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DpMat3x3MatMultExpr
    : public A2DExpression<
          A2DpMat3x3MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2DpMat3x3MatMultExpr(Mat3x3& A, A2DMat<N, BMatType>& BObj,
                        A2DMat<N, CMatType>& CObj)
      : A(A), BObj(BObj), CObj(CObj) {
    const BMatType& B = BObj.value();
    CMatType& C = CObj.value();

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

  void reverse() {
    BMatType& Bb = BObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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

  void hforward() {
    for (int i = 0; i < N; i++) {
      const BMatType& Bp = BObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

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

  void hreverse() {
    for (int i = 0; i < N; i++) {
      BMatType& Bh = BObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);

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
  A2DMat<N, BMatType>& BObj;
  A2DMat<N, CMatType>& CObj;
};

template <int N, class ScalarType, class BMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DpMat3x3MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT>
Mat3x3MatMult(Mat<ScalarType, 3, 3>& A, A2DMat<N, BMatType>& BObj,
              A2DMat<N, CMatType>& CObj) {
  return A2DpMat3x3MatMultExpr<N, ScalarType, BMatType, CMatType, AT, BT>(
      A, BObj, CObj);
}

template <int N, typename ScalarType, class AMatType, class CMatType,
          bool AT = false, bool BT = false>
class A2DMat3x3pMatMultExpr
    : public ADExpression<
          A2DMat3x3pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT> > {
 public:
  typedef Mat<ScalarType, 3, 3> Mat3x3;

  A2DMat3x3pMatMultExpr(A2DMat<N, AMatType>& AObj, Mat3x3& B,
                        A2DMat<N, CMatType>& CObj)
      : AObj(AObj), B(B), CObj(CObj) {
    const AMatType& A = AObj.value();
    CMatType& C = CObj.value();

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

  void reverse() {
    AMatType& Ab = AObj.bvalue();
    const CMatType& Cb = CObj.bvalue();

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

  void hforward() {
    for (int i = 0; i < N; i++) {
      const AMatType& Ap = AObj.pvalue(i);
      CMatType& Cp = CObj.pvalue(i);

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

  void hreverse() {
    for (int i = 0; i < N; i++) {
      AMatType& Ah = AObj.hvalue(i);
      const CMatType& Ch = CObj.hvalue(i);

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

  A2DMat<N, AMatType>& AObj;
  Mat3x3& B;
  A2DMat<N, CMatType>& CObj;
};

template <int N, typename ScalarType, class AMatType, class CMatType,
          bool AT = false, bool BT = false>
inline A2DMat3x3pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT>
Mat3x3MatMult(A2DMat<N, AMatType>& AObj, Mat<ScalarType, 3, 3>& B,
              A2DMat<N, CMatType>& CObj) {
  return A2DMat3x3pMatMultExpr<N, ScalarType, AMatType, CMatType, AT, BT>(
      AObj, B, CObj);
}

// Mat3x3Det
template <typename ScalarType>
inline void Mat3x3Det(const Mat<ScalarType, 3, 3>& A, ScalarType& det) {
  det = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
         A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
         A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
}

template <class MatType, class ScalarType>
class ADMat3x3DetExpr
    : public ADExpression<ADMat3x3DetExpr<MatType, ScalarType> > {
 public:
  ADMat3x3DetExpr(ADMat<MatType>& AObj, ADScalar<ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const MatType& A = AObj.value();

    detObj.value = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                    A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                    A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
  }

  void forward() {
    const MatType& A = AObj.value();
    const MatType& Ad = AObj.bvalue();

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

  void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    MatType& Ab = AObj.bvalue();

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

  ADMat<MatType>& AObj;
  ADScalar<ScalarType>& detObj;
};

template <class MatType, typename ScalarType>
inline ADMat3x3DetExpr<MatType, ScalarType> Mat3x3Det(
    ADMat<MatType>& A, ADScalar<ScalarType>& det) {
  return ADMat3x3DetExpr<MatType, ScalarType>(A, det);
}

template <int N, class MatType, class ScalarType>
class A2DMat3x3DetExpr
    : public ADExpression<A2DMat3x3DetExpr<N, MatType, ScalarType> > {
 public:
  A2DMat3x3DetExpr(A2DMat<N, MatType>& AObj, A2DScalar<N, ScalarType>& detObj)
      : AObj(AObj), detObj(detObj) {
    const MatType& A = AObj.value();

    detObj.value = (A(2, 2) * (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) -
                    A(2, 1) * (A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) +
                    A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)));
  }

  void reverse() {
    const ScalarType& bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    MatType& Ab = AObj.bvalue();

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

  void hforward() {
    const MatType& A = AObj.value();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
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

  void hreverse() {
    const ScalarType bdet = detObj.bvalue;
    const MatType& A = AObj.value();
    const MatType& Ab = AObj.bvalue();

    for (int i = 0; i < N; i++) {
      const ScalarType hdet = detObj.hvalue[i];
      const MatType& Ap = AObj.pvalue(i);
      MatType& Ah = AObj.hvalue(i);

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

  A2DMat<N, MatType>& AObj;
  A2DScalar<N, ScalarType>& detObj;
};

template <int N, class MatType, typename ScalarType>
inline A2DMat3x3DetExpr<N, MatType, ScalarType> Mat3x3Det(
    A2DMat<N, MatType>& A, A2DScalar<N, ScalarType>& det) {
  return A2DMat3x3DetExpr<N, MatType, ScalarType>(A, det);
}

// Mat3x3Inverse
template <typename ScalarType>
inline void Mat3x3Inverse(const Mat<ScalarType, 3, 3>& A,
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

template <class MatType, class InvMatType>
class ADMat3x3InverseExpr
    : public ADExpression<ADMat3x3InverseExpr<MatType, InvMatType> > {
 public:
  typedef typename MatType::type ScalarType;

  ADMat3x3InverseExpr(ADMat<MatType>& AObj, ADMat<InvMatType>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const MatType& A = AObj.value();
    InvMatType& Ainv = AinvObj.value();

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

  void forward() {
    const InvMatType& Ainv = AinvObj.value();
    const MatType& Ad = AObj.bvalue();
    InvMatType& Ainvd = AinvObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    Mat3x3MatMultCore(Ainv, Ad, tmp);
    Mat3x3MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvd);
  }

  void reverse() {
    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvb = AinvObj.bvalue();
    MatType& Ab = AObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    MatTrans3x3MatMultCore(Ainv, Ainvb, tmp);
    Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ab);
  }

  ADMat<MatType>& AObj;
  ADMat<InvMatType>& AinvObj;
};

template <class MatType, class InvMatType>
inline ADMat3x3InverseExpr<MatType, InvMatType> Mat3x3Inverse(
    ADMat<MatType>& AObj, ADMat<InvMatType>& AinvObj) {
  return ADMat3x3InverseExpr<MatType, InvMatType>(AObj, AinvObj);
}

template <int N, class MatType, class InvMatType>
class A2DMat3x3InverseExpr
    : public A2DExpression<A2DMat3x3InverseExpr<N, MatType, InvMatType> > {
 public:
  typedef typename MatType::type ScalarType;

  A2DMat3x3InverseExpr(A2DMat<N, MatType>& AObj, A2DMat<N, InvMatType>& AinvObj)
      : AObj(AObj), AinvObj(AinvObj) {
    const MatType& A = AObj.value();
    InvMatType& Ainv = AinvObj.value();

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

  void reverse() {
    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvd = AinvObj.bvalue();
    MatType& Ad = AObj.bvalue();

    Mat<ScalarType, 3, 3> tmp;
    MatTrans3x3MatMultCore(Ainv, Ainvd, tmp);
    Mat3x3MatTransMultAddScaleCore(ScalarType(-1.0), tmp, Ainv, Ad);
  }

  void hforward() {
    Mat<ScalarType, 3, 3> tmp;
    const InvMatType& Ainv = AinvObj.value();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
      InvMatType& Ainvp = AinvObj.pvalue(i);

      Mat3x3MatMultCore(Ainv, Ap, tmp);
      Mat3x3MatMultScaleCore(ScalarType(-1.0), tmp, Ainv, Ainvp);
    }
  }

  // hA = A^{-T} * Ap^{T} * A^{-T} * Ainvb * A^{-T} +
  //      A^{-T} * Ainvb * A^{-T} * Ap^{T} * A^{-T} =
  //    = - (A^{-T} * Ap^{T} * Ab + Ab * Ap^{T} * A^{-T})

  void hreverse() {
    // Temporary matrix
    Mat<ScalarType, 3, 3> tmp, tmp2;

    const InvMatType& Ainv = AinvObj.value();
    const InvMatType& Ainvb = AinvObj.bvalue();
    const MatType& Ab = AObj.bvalue();

    for (int i = 0; i < N; i++) {
      const MatType& Ap = AObj.pvalue(i);
      const InvMatType& Ainvh = AinvObj.hvalue(i);
      MatType& Ah = AObj.hvalue(i);

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

  A2DMat<N, MatType>& AObj;
  A2DMat<N, InvMatType>& AinvObj;
};

template <int N, class MatType, class InvMatType>
inline A2DMat3x3InverseExpr<N, MatType, InvMatType> Mat3x3Inverse(
    A2DMat<N, MatType>& AObj, A2DMat<N, InvMatType>& AinvObj) {
  return A2DMat3x3InverseExpr<N, MatType, InvMatType>(AObj, AinvObj);
}

// Symm3x3SymmMultTrace
template <typename ScalarType>
inline void Symm3x3SymmMultTrace(const SymmMat<ScalarType, 3>& S,
                                 const SymmMat<ScalarType, 3>& E,
                                 ScalarType& trace) {
  trace = (S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
           2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2)));
}

template <class SMatType, class EMatType, class ScalarType>
class ADSymm3x3SymmMultTraceExpr
    : public ADExpression<
          ADSymm3x3SymmMultTraceExpr<SMatType, EMatType, ScalarType> > {
 public:
  ADSymm3x3SymmMultTraceExpr(ADMat<SMatType>& SObj, ADMat<EMatType>& EObj,
                             ADScalar<ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SMatType& S = SObj.value();
    const EMatType& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
        2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2));
  }

  void forward() {
    const EMatType& E = EObj.value();
    const EMatType& Ed = EObj.bvalue();
    const SMatType& S = SObj.value();
    const SMatType& Sd = SObj.bvalue();

    output.pvalue =
        S(0, 0) * Ed(0, 0) + S(1, 1) * Ed(1, 1) + S(2, 2) * Ed(2, 2) +
        2.0 * (S(0, 1) * Ed(0, 1) + S(0, 2) * Ed(0, 2) + S(1, 2) * Ed(1, 2)) +
        Sd(0, 0) * E(0, 0) + Sd(1, 1) * E(1, 1) + Sd(2, 2) * E(2, 2) +
        2.0 * (Sd(0, 1) * E(0, 1) + Sd(0, 2) * E(0, 2) + Sd(1, 2) * E(1, 2));
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const SMatType& S = SObj.value();
    SMatType& Sb = SObj.bvalue();

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

  ADMat<SMatType>& SObj;
  ADMat<EMatType>& EObj;
  ADScalar<ScalarType>& output;
};

template <class SMatType, class EMatType, class ScalarType>
inline ADSymm3x3SymmMultTraceExpr<SMatType, EMatType, ScalarType>
Symm3x3SymmMultTrace(ADMat<SMatType>& S, ADMat<EMatType>& E,
                     ADScalar<ScalarType>& trace) {
  return ADSymm3x3SymmMultTraceExpr<SMatType, EMatType, ScalarType>(S, E,
                                                                    trace);
}

template <int N, class SMatType, class EMatType, class ScalarType>
class A2DSymm3x3SymmMultTraceExpr
    : public A2DExpression<
          A2DSymm3x3SymmMultTraceExpr<N, SMatType, EMatType, ScalarType> > {
 public:
  A2DSymm3x3SymmMultTraceExpr(A2DMat<N, SMatType>& SObj,
                              A2DMat<N, EMatType>& EObj,
                              A2DScalar<N, ScalarType>& output)
      : SObj(SObj), EObj(EObj), output(output) {
    const SMatType& S = SObj.value();
    const EMatType& E = EObj.value();

    output.value =
        S(0, 0) * E(0, 0) + S(1, 1) * E(1, 1) + S(2, 2) * E(2, 2) +
        2.0 * (S(0, 1) * E(0, 1) + S(0, 2) * E(0, 2) + S(1, 2) * E(1, 2));
  }

  void reverse() {
    const EMatType& E = EObj.value();
    EMatType& Eb = EObj.bvalue();
    const SMatType& S = SObj.value();
    SMatType& Sb = SObj.bvalue();

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
  void hforward() {
    const EMatType& E = EObj.value();
    const SMatType& S = SObj.value();

    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      const SMatType& Sp = SObj.pvalue(i);

      output.pvalue[i] =
          S(0, 0) * Ep(0, 0) + S(1, 1) * Ep(1, 1) + S(2, 2) * Ep(2, 2) +
          2.0 * (S(0, 1) * Ep(0, 1) + S(0, 2) * Ep(0, 2) + S(1, 2) * Ep(1, 2)) +
          Sp(0, 0) * E(0, 0) + Sp(1, 1) * E(1, 1) + Sp(2, 2) * E(2, 2) +
          2.0 * (Sp(0, 1) * E(0, 1) + Sp(0, 2) * E(0, 2) + Sp(1, 2) * E(1, 2));
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

  A2DMat<N, SMatType>& SObj;
  A2DMat<N, EMatType>& EObj;
  A2DScalar<N, ScalarType>& output;
};

template <int N, class SMatType, class EMatType, class ScalarType>
inline A2DSymm3x3SymmMultTraceExpr<N, SMatType, EMatType, ScalarType>
Symm3x3SymmMultTrace(A2DMat<N, SMatType>& S, A2DMat<N, EMatType>& E,
                     A2DScalar<N, ScalarType>& trace) {
  return A2DSymm3x3SymmMultTraceExpr<N, SMatType, EMatType, ScalarType>(S, E,
                                                                        trace);
}

template <class ScalarType>
inline void Symm3x3IsotropicConstitutive(const ScalarType& mu,
                                         const ScalarType& lambda,
                                         const SymmMat<ScalarType, 3>& E,
                                         SymmMat<ScalarType, 3>& S) {
  ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
  ScalarType mu2 = 2.0 * mu;
  S(0, 0) = mu2 * E(0, 0) + tr;
  S(0, 1) = mu2 * E(0, 1);
  S(0, 2) = mu2 * E(0, 2);
  S(1, 1) = mu2 * E(1, 1) + tr;
  S(1, 2) = mu2 * E(1, 2);
  S(2, 2) = mu2 * E(2, 2) + tr;
}

template <class ScalarType, class EMatType, class SMatType>
class ADSymm3x3IsotropicConstitutiveExpr
    : public ADExpression<
          ADSymm3x3IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType> > {
 public:
  ADSymm3x3IsotropicConstitutiveExpr(const ScalarType& mu,
                                     const ScalarType& lambda,
                                     ADMat<EMatType>& EObj,
                                     ADMat<SMatType>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const EMatType& E = EObj.value();
    SMatType& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(0, 2) = mu2 * E(0, 2);
    S(1, 1) = mu2 * E(1, 1) + tr;
    S(1, 2) = mu2 * E(1, 2);
    S(2, 2) = mu2 * E(2, 2) + tr;
  }

  void forward() {
    const EMatType& Ed = EObj.bvalue();
    SMatType& Sd = SObj.pvalue();

    ScalarType tr = lambda * (Ed(0, 0) + Ed(1, 1) + Ed(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Sd(0, 0) = mu2 * Ed(0, 0) + tr;
    Sd(0, 1) = mu2 * Ed(0, 1);
    Sd(0, 2) = mu2 * Ed(0, 2);
    Sd(1, 1) = mu2 * Ed(1, 1) + tr;
    Sd(1, 2) = mu2 * Ed(1, 2);
    Sd(2, 2) = mu2 * Ed(2, 2) + tr;
  }

  void reverse() {
    const SMatType& Sb = SObj.bvalue();
    EMatType& Eb = EObj.bvalue();

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
  ADMat<EMatType>& EObj;
  ADMat<SMatType>& SObj;
};

template <class ScalarType, class EMatType, class SMatType>
inline ADSymm3x3IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType>
Symm3x3IsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                             ADMat<EMatType>& E, ADMat<SMatType>& S) {
  return ADSymm3x3IsotropicConstitutiveExpr<ScalarType, EMatType, SMatType>(
      mu, lambda, E, S);
}

template <int N, class ScalarType, class EMatType, class SMatType>
class A2DSymm3x3IsotropicConstitutiveExpr
    : public A2DExpression<A2DSymm3x3IsotropicConstitutiveExpr<
          N, ScalarType, EMatType, SMatType> > {
 public:
  A2DSymm3x3IsotropicConstitutiveExpr(const ScalarType& mu,
                                      const ScalarType& lambda,
                                      A2DMat<N, EMatType>& EObj,
                                      A2DMat<N, SMatType>& SObj)
      : mu(mu), lambda(lambda), EObj(EObj), SObj(SObj) {
    const EMatType& E = EObj.value();
    SMatType& S = SObj.value();
    ScalarType tr = lambda * (E(0, 0) + E(1, 1) + E(2, 2));
    ScalarType mu2 = 2.0 * mu;
    S(0, 0) = mu2 * E(0, 0) + tr;
    S(0, 1) = mu2 * E(0, 1);
    S(0, 2) = mu2 * E(0, 2);
    S(1, 1) = mu2 * E(1, 1) + tr;
    S(1, 2) = mu2 * E(1, 2);
    S(2, 2) = mu2 * E(2, 2) + tr;
  }

  void reverse() {
    const SMatType& Sb = SObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    ScalarType tr = lambda * (Sb(0, 0) + Sb(1, 1) + Sb(2, 2));
    ScalarType mu2 = 2.0 * mu;
    Eb(0, 0) += mu2 * Sb(0, 0) + tr;
    Eb(0, 1) += mu2 * Sb(0, 1);
    Eb(0, 2) += mu2 * Sb(0, 2);
    Eb(1, 1) += mu2 * Sb(1, 1) + tr;
    Eb(1, 2) += mu2 * Sb(1, 2);
    Eb(2, 2) += mu2 * Sb(2, 2) + tr;
  }

  void hforward() {
    for (int i = 0; i < N; i++) {
      const EMatType& Ep = EObj.pvalue(i);
      SMatType& Sp = SObj.pvalue(i);

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

  void hreverse() {
    for (int i = 0; i < N; i++) {
      const SMatType& Sh = SObj.hvalue(i);
      EMatType& Eh = EObj.hvalue(i);

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
  A2DMat<N, EMatType>& EObj;
  A2DMat<N, SMatType>& SObj;
};

template <int N, class ScalarType, class EMatType, class SMatType>
inline A2DSymm3x3IsotropicConstitutiveExpr<N, ScalarType, EMatType, SMatType>
Symm3x3IsotropicConstitutive(const ScalarType& mu, const ScalarType& lambda,
                             A2DMat<N, EMatType>& E, A2DMat<N, SMatType>& S) {
  return A2DSymm3x3IsotropicConstitutiveExpr<N, ScalarType, EMatType, SMatType>(
      mu, lambda, E, S);
}

template <class ScalarType>
inline void Mat3x3GreenStrain(const Mat<ScalarType, 3, 3>& Ux,
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

template <class UxMatType, class EMatType>
class ADMat3x3GreenStrainExpr
    : public ADExpression<ADMat3x3GreenStrainExpr<UxMatType, EMatType> > {
 public:
  ADMat3x3GreenStrainExpr(ADMat<UxMatType>& UxObj, ADMat<EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
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

  void forward() {
    const UxMatType& Ux = UxObj.value();
    const UxMatType& Uxd = UxObj.bvalue();
    EMatType& Ed = EObj.bvalue();

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

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

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

  ADMat<UxMatType>& UxObj;
  ADMat<EMatType>& EObj;
};

template <class UxMatType, class EMatType>
inline ADMat3x3GreenStrainExpr<UxMatType, EMatType> Mat3x3GreenStrain(
    ADMat<UxMatType>& Ux, ADMat<EMatType>& E) {
  return ADMat3x3GreenStrainExpr<UxMatType, EMatType>(Ux, E);
}

template <int N, class UxMatType, class EMatType>
class A2DMat3x3GreenStrainExpr
    : public A2DExpression<A2DMat3x3GreenStrainExpr<N, UxMatType, EMatType> > {
 public:
  A2DMat3x3GreenStrainExpr(A2DMat<N, UxMatType>& UxObj,
                           A2DMat<N, EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
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

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

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

  void hforward() {
    const UxMatType& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      EMatType& Ep = EObj.pvalue(i);
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

  void hreverse() {
    const UxMatType& Eb = EObj.bvalue();
    const UxMatType& Ux = UxObj.value();

    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      const EMatType& Eh = EObj.hvalue(i);
      UxMatType& Uxh = UxObj.hvalue(i);

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

  A2DMat<N, UxMatType>& UxObj;
  A2DMat<N, EMatType>& EObj;
};

template <int N, class UxMatType, class EMatType>
inline A2DMat3x3GreenStrainExpr<N, UxMatType, EMatType> Mat3x3GreenStrain(
    A2DMat<N, UxMatType>& Ux, A2DMat<N, EMatType>& E) {
  return A2DMat3x3GreenStrainExpr<N, UxMatType, EMatType>(Ux, E);
}

template <class ScalarType>
inline void Mat3x3LinearGreenStrain(const Mat<ScalarType, 3, 3>& Ux,
                                    SymmMat<ScalarType, 3>& E) {
  E(0, 0) = Ux(0, 0);
  E(1, 1) = Ux(1, 1);
  E(2, 2) = Ux(2, 2);

  E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
  E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
  E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
}

template <class UxMatType, class EMatType>
class ADMat3x3LinearGreenStrainExpr
    : public ADExpression<ADMat3x3LinearGreenStrainExpr<UxMatType, EMatType> > {
 public:
  ADMat3x3LinearGreenStrainExpr(ADMat<UxMatType>& UxObj, ADMat<EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);
    E(2, 2) = Ux(2, 2);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
  }

  void forward() {
    const UxMatType& Ux = UxObj.value();
    const UxMatType& Uxb = UxObj.bvalue();
    EMatType& Eb = EObj.bvalue();

    Eb(0, 0) = Uxb(0, 0);
    Eb(1, 1) = Uxb(1, 1);
    Eb(2, 2) = Uxb(2, 2);

    Eb(0, 1) = 0.5 * (Uxb(0, 1) + Uxb(1, 0));
    Eb(0, 2) = 0.5 * (Uxb(0, 2) + Uxb(2, 0));
    Eb(1, 2) = 0.5 * (Uxb(1, 2) + Uxb(2, 1));
  }

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

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

  ADMat<UxMatType>& UxObj;
  ADMat<EMatType>& EObj;
};

template <class UxMatType, class EMatType>
inline ADMat3x3LinearGreenStrainExpr<UxMatType, EMatType>
Mat3x3LinearGreenStrain(ADMat<UxMatType>& Ux, ADMat<EMatType>& E) {
  return ADMat3x3LinearGreenStrainExpr<UxMatType, EMatType>(Ux, E);
}

template <int N, class UxMatType, class EMatType>
class A2DMat3x3LinearGreenStrainExpr
    : public A2DExpression<
          A2DMat3x3LinearGreenStrainExpr<N, UxMatType, EMatType> > {
 public:
  A2DMat3x3LinearGreenStrainExpr(A2DMat<N, UxMatType>& UxObj,
                                 A2DMat<N, EMatType>& EObj)
      : UxObj(UxObj), EObj(EObj) {
    const UxMatType& Ux = UxObj.value();
    EMatType& E = EObj.value();
    E(0, 0) = Ux(0, 0);
    E(1, 1) = Ux(1, 1);
    E(2, 2) = Ux(2, 2);

    E(0, 1) = 0.5 * (Ux(0, 1) + Ux(1, 0));
    E(0, 2) = 0.5 * (Ux(0, 2) + Ux(2, 0));
    E(1, 2) = 0.5 * (Ux(1, 2) + Ux(2, 1));
  }

  void reverse() {
    const UxMatType& Ux = UxObj.value();
    const EMatType& Eb = EObj.bvalue();
    UxMatType& Uxb = UxObj.bvalue();

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

  void hforward() {
    for (int i = 0; i < N; i++) {
      const UxMatType& Uxp = UxObj.pvalue(i);
      EMatType& Ep = EObj.pvalue(i);

      Ep(0, 0) = Uxp(0, 0);
      Ep(1, 1) = Uxp(1, 1);
      Ep(2, 2) = Uxp(2, 2);

      Ep(0, 1) = 0.5 * (Uxp(0, 1) + Uxp(1, 0));
      Ep(0, 2) = 0.5 * (Uxp(0, 2) + Uxp(2, 0));
      Ep(1, 2) = 0.5 * (Uxp(1, 2) + Uxp(2, 1));
    }
  }

  void hreverse() {
    for (int i = 0; i < N; i++) {
      const EMatType& Eh = EObj.hvalue(i);
      UxMatType& Uxh = UxObj.hvalue(i);

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

  A2DMat<N, UxMatType>& UxObj;
  A2DMat<N, EMatType>& EObj;
};

template <int N, class UxMatType, class EMatType>
inline A2DMat3x3LinearGreenStrainExpr<N, UxMatType, EMatType>
Mat3x3LinearGreenStrain(A2DMat<N, UxMatType>& Ux, A2DMat<N, EMatType>& E) {
  return A2DMat3x3LinearGreenStrainExpr<N, UxMatType, EMatType>(Ux, E);
}

}  // namespace A2D

#endif  // A2D_TMP_H
