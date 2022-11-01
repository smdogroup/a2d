#ifndef A2D_MAT_VEC_OPS_3D_H
#define A2D_MAT_VEC_OPS_3D_H

#include "a2dmatveccore3d.h"
#include "a2dobjs.h"
#include "a2dtypes.h"

namespace A2D {

// TODO: remove disassemblies that store to zero vectors and instead replace them with implementations of calls to NULL
//  in place of that vector
// TODO: remove assemblies that use zero vectors and instead replace them with implementations of calls to NULL in
//  place of that vector

/*
  Mat3x2ToVec3
  This represents the disassembly of a 3x2 matrix into two 3 vectors
  u(i) = A(i, 0)
  v(i) = A(i, 1)
*/
template <typename T>
inline void Mat3x2ToVec3(const Mat<T, 3, 2>& A,
                         Vec<T, 3>& u,
                         Vec<T, 3>& v) {
  Mat3x2ToVec3Core(A, u, v);
}

template <typename T>
class ADMat3x2ToVec3Expr {
 public:
  ADMat3x2ToVec3Expr(ADMat<Mat<T, 3, 2>>& AObj,
                     ADVec<Vec<T, 3>>& uObj,
                     ADVec<Vec<T, 3>>& vObj)
      : AObj(AObj), uObj(uObj), vObj(vObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Vec<T, 3>& u = uObj.value();
    Vec<T, 3>& v = vObj.value();
    Mat3x2ToVec3Core(A, u, v);
  }

  void forward() {
//            input:
    const Mat<T, 3, 2>& Ab = AObj.bvalue();
//            output:
    Vec<T, 3>& ub = uObj.bvalue();
    Vec<T, 3>& vb = vObj.bvalue();
//            operations:
    Mat3x2ToVec3Core(Ab, ub, vb);
  };

  void reverse() {
//            input:
    const Vec<T, 3>& ub = uObj.bvalue();
    const Vec<T, 3>& vb = vObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
//            operations:
    Mat3x2FromTwoVec3Core(ub, vb, Ab);
  };

  ADMat<Mat<T, 3, 2>>& AObj;
  ADVec<Vec<T, 3>>& uObj;
  ADVec<Vec<T, 3>>& vObj;
};

template <typename T>
inline ADMat3x2ToVec3Expr<T> Mat3x2ToVec3(ADMat<Mat<T, 3, 2>>& A,
                                          ADVec<Vec<T, 3>>& u,
                                          ADVec<Vec<T, 3>>& v) {
  return ADMat3x2ToVec3Expr<T>(A, u, v);
}

template <int N, typename T>
class A2DMat3x2ToVec3Expr {
 public:
  A2DMat3x2ToVec3Expr(A2DMat<N, Mat<T, 3, 2>>& AObj,
                      A2DVec<N, Vec<T, 3>>& uObj,
                      A2DVec<N, Vec<T, 3>>& vObj)
      : AObj(AObj), uObj(uObj), vObj(vObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Vec<T, 3>& u = uObj.value();
    Vec<T, 3>& v = vObj.value();
    Mat3x2ToVec3Core(A, u, v);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& ub = uObj.bvalue();
    const Vec<T, 3>& vb = vObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
//            operations:
    Mat3x2FromTwoVec3Core(ub, vb, Ab);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 2>& Ap = AObj.pvalue(i);
//                output:
      Vec<T, 3>& up = uObj.pvalue(i);
      Vec<T, 3>& vp = vObj.pvalue(i);
//                operations:
      Mat3x2ToVec3Core(Ap, up, vp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3>& uh = uObj.hvalue(i);
      const Vec<T, 3>& vh = vObj.hvalue(i);
//                output:
      Mat<T, 3, 2>& Ah = AObj.hvalue(i);
//                operations:
      Mat3x2FromTwoVec3Core(uh, vh, Ah);
    }
  };

  A2DMat<N, Mat<T, 3, 2>>& AObj;
  A2DVec<N, Vec<T, 3>>& uObj;
  A2DVec<N, Vec<T, 3>>& vObj;
};

template <int N, typename T>
inline A2DMat3x2ToVec3Expr<N, T> Mat3x2ToVec3(A2DMat<N, Mat<T, 3, 2>>& A,
                                              A2DVec<N, Vec<T, 3>>& u,
                                              A2DVec<N, Vec<T, 3>>& v) {
  return A2DMat3x2ToVec3Expr<N, T>(A, u, v);
}

/*
  Mat3x3FromThreeVec3
  This represents the assembly of three 3 vectors into a 3x3 matrix
  C(i, 0) = x(i)
  C(i, 1) = y(i)
  C(i, 2) = z(i)
*/
template <typename T>
inline void Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                const Vec<T, 3>& y,
                                const Vec<T, 3>& z,
                                Mat<T, 3, 3>& C) {
  Mat3x3FromThreeVec3Core(x, y, z, C);
}

template <typename T>
class ADVecVecVecToMat3x3Expr {
 public:
  ADVecVecVecToMat3x3Expr(ADVec<Vec<T, 3>>& xObj,
                          const Vec<T, 3>& y,
                          const Vec<T, 3>& z,
                          ADMat<Mat<T, 3, 3>>& CObj)
      : xObj(xObj), y(y), z(z), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3> zeros{};
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(xb, zeros, zeros, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, zeros, zeros);
  };

  ADVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  const Vec<T, 3>& z;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADVecVecVecToMat3x3Expr<T> Mat3x3FromThreeVec3(ADVec<Vec<T, 3>>& x,
                                                      const Vec<T, 3>& y,
                                                      const Vec<T, 3>& z,
                                                      ADMat<Mat<T, 3, 3>>& C) {
  return ADVecVecVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class VecADVecVecToMat3x3Expr {
 public:
  VecADVecVecToMat3x3Expr(const Vec<T, 3>& x,
                          ADVec<Vec<T, 3>>& yObj,
                          const Vec<T, 3>& z,
                          ADMat<Mat<T, 3, 3>>& CObj)
      : x(x), yObj(yObj), z(z), CObj(CObj) {
    const Vec<T, 3>& y = yObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3> zeros{};
    const Vec<T, 3>& yb = yObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(zeros, yb, zeros, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& yb = yObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, yb, zeros);
  };

  const Vec<T, 3>& x;
  ADVec<Vec<T, 3>>& yObj;
  const Vec<T, 3>& z;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline VecADVecVecToMat3x3Expr<T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                      ADVec<Vec<T, 3>>& y,
                                                      const Vec<T, 3>& z,
                                                      ADMat<Mat<T, 3, 3>>& C) {
  return VecADVecVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class VecVecADVecToMat3x3Expr {
 public:
  VecVecADVecToMat3x3Expr(const Vec<T, 3>& x,
                          const Vec<T, 3>& y,
                          ADVec<Vec<T, 3>>& zObj,
                          ADMat<Mat<T, 3, 3>>& CObj)
      : x(x), y(y), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3> zeros{};
    const Vec<T, 3>& zb = zObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(zeros, zeros, zb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, zeros, zb);
  };

  const Vec<T, 3>& x;
  const Vec<T, 3>& y;
  ADVec<Vec<T, 3>>& zObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline VecVecADVecToMat3x3Expr<T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                      const Vec<T, 3>& y,
                                                      ADVec<Vec<T, 3>>& z,
                                                      ADMat<Mat<T, 3, 3>>& C) {
  return VecVecADVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class VecADVecADVecToMat3x3Expr {
 public:
  VecADVecADVecToMat3x3Expr(const Vec<T, 3>& x,
                            ADVec<Vec<T, 3>>& yObj,
                            ADVec<Vec<T, 3>>& zObj,
                            ADMat<Mat<T, 3, 3>>& CObj)
      : x(x), yObj(yObj), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3> zeros{};
    const Vec<T, 3>& yb = yObj.bvalue();
    const Vec<T, 3>& zb = zObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(zeros, yb, zb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, yb, zb);
  };

  const Vec<T, 3>& x;
  ADVec<Vec<T, 3>>& yObj;
  ADVec<Vec<T, 3>>& zObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline VecADVecADVecToMat3x3Expr<T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                        ADVec<Vec<T, 3>>& y,
                                                        ADVec<Vec<T, 3>>& z,
                                                        ADMat<Mat<T, 3, 3>>& C) {
  return VecADVecADVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class ADVecVecADVecToMat3x3Expr {
 public:
  ADVecVecADVecToMat3x3Expr(ADVec<Vec<T, 3>>& xObj,
                            const Vec<T, 3>& y,
                            ADVec<Vec<T, 3>>& zObj,
                            ADMat<Mat<T, 3, 3>>& CObj)
      : xObj(xObj), y(y), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3> zeros{};
    const Vec<T, 3>& zb = zObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(xb, zeros, zb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3> zeros{};
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, zeros, zb);
  };

  ADVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  ADVec<Vec<T, 3>>& zObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADVecVecADVecToMat3x3Expr<T> Mat3x3FromThreeVec3(ADVec<Vec<T, 3>>& x,
                                                        const Vec<T, 3>& y,
                                                        ADVec<Vec<T, 3>>& z,
                                                        ADMat<Mat<T, 3, 3>>& C) {
  return ADVecVecADVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class ADVecADVecVecToMat3x3Expr {
 public:
  ADVecADVecVecToMat3x3Expr(ADVec<Vec<T, 3>>& xObj,
                            ADVec<Vec<T, 3>>& yObj,
                            const Vec<T, 3>& z,
                            ADMat<Mat<T, 3, 3>>& CObj)
      : xObj(xObj), yObj(yObj), z(z), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3>& yb = yObj.bvalue();
    const Vec<T, 3> zeros{};
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(xb, yb, zeros, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, yb, zeros);
  };

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& yObj;
  const Vec<T, 3>& z;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADVecADVecVecToMat3x3Expr<T> Mat3x3FromThreeVec3(ADVec<Vec<T, 3>>& x,
                                                        ADVec<Vec<T, 3>>& y,
                                                        const Vec<T, 3>& z,
                                                        ADMat<Mat<T, 3, 3>>& C) {
  return ADVecADVecVecToMat3x3Expr<T>(x, y, z, C);
}

template <typename T>
class ADVecADVecADVecToMat3x3Expr {
 public:
  ADVecADVecADVecToMat3x3Expr(ADVec<Vec<T, 3>>& xObj,
                              ADVec<Vec<T, 3>>& yObj,
                              ADVec<Vec<T, 3>>& zObj,
                              ADMat<Mat<T, 3, 3>>& CObj)
      : xObj(xObj), yObj(yObj), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3>& yb = yObj.bvalue();
    const Vec<T, 3>& zb = zObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromThreeVec3Core(xb, yb, zb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, yb, zb);
  };

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& yObj;
  ADVec<Vec<T, 3>>& zObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADVecADVecADVecToMat3x3Expr<T> Mat3x3FromThreeVec3(ADVec<Vec<T, 3>>& x,
                                                          ADVec<Vec<T, 3>>& y,
                                                          ADVec<Vec<T, 3>>& z,
                                                          ADMat<Mat<T, 3, 3>>& C) {
  return ADVecADVecADVecToMat3x3Expr<T>(x, y, z, C);
}

template <int N, typename T>
class A2DVecVecVecToMat3x3Expr {
 public:
  A2DVecVecVecToMat3x3Expr(A2DVec<N, Vec<T, 3>>& xObj,
                           const Vec<T, 3>& y,
                           const Vec<T, 3>& z,
                           A2DMat<N, Mat<T, 3, 3>>& CObj)
      : xObj(xObj), y(y), z(z), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, zeros, zeros);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3> zeros{};
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(xp, zeros, zeros, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3> zeros{};
//                operations:
      Mat3x3ToThreeVec3Core(Ch, xh, zeros, zeros);
    }
  };

  A2DVec<N, Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  const Vec<T, 3>& z;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DVecVecVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(A2DVec<N, Vec<T, 3>>& x,
                                                          const Vec<T, 3>& y,
                                                          const Vec<T, 3>& z,
                                                          A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DVecVecVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class VecA2DVecVecToMat3x3Expr {
 public:
  VecA2DVecVecToMat3x3Expr(const Vec<T, 3>& x,
                           A2DVec<N, Vec<T, 3>>& yObj,
                           const Vec<T, 3>& z,
                           A2DMat<N, Mat<T, 3, 3>>& CObj)
      : x(x), yObj(yObj), z(z), CObj(CObj) {
    const Vec<T, 3>& y = yObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& yb = yObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, yb, zeros);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3> zeros{};
      const Vec<T, 3>& yp = yObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(zeros, yp, zeros, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3> zeros{};
      Vec<T, 3>& yh = yObj.hvalue(i);
//                operations:
      Mat3x3ToThreeVec3Core(Ch, zeros, yh, zeros);
    }
  };

  const Vec<T, 3>& x;
  A2DVec<N, Vec<T, 3>>& yObj;
  const Vec<T, 3>& z;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline VecA2DVecVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                          A2DVec<N, Vec<T, 3>>& y,
                                                          const Vec<T, 3>& z,
                                                          A2DMat<N, Mat<T, 3, 3>>& C) {
  return VecA2DVecVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class VecVecA2DVecToMat3x3Expr {
 public:
  VecVecA2DVecToMat3x3Expr(const Vec<T, 3>& x,
                           const Vec<T, 3>& y,
                           A2DVec<N, Vec<T, 3>>& zObj,
                           A2DMat<N, Mat<T, 3, 3>>& CObj)
      : x(x), y(y), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, zeros, zb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3> zeros{};
      const Vec<T, 3>& zp = zObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(zeros, zeros, zp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3> zeros{};
      Vec<T, 3>& zh = zObj.hvalue(i);
//                operations:
      Mat3x3ToThreeVec3Core(Ch, zeros, zeros, zh);
    }
  };

  const Vec<T, 3>& x;
  const Vec<T, 3>& y;
  A2DVec<N, Vec<T, 3>>& zObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline VecVecA2DVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                          const Vec<T, 3>& y,
                                                          A2DVec<N, Vec<T, 3>>& z,
                                                          A2DMat<N, Mat<T, 3, 3>>& C) {
  return VecVecA2DVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class VecA2DVecA2DVecToMat3x3Expr {
 public:
  VecA2DVecA2DVecToMat3x3Expr(const Vec<T, 3>& x,
                              A2DVec<N, Vec<T, 3>>& yObj,
                              A2DVec<N, Vec<T, 3>>& zObj,
                              A2DMat<N, Mat<T, 3, 3>>& CObj)
      : x(x), yObj(yObj), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3> zeros{};
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, zeros, yb, zb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3> zeros{};
      const Vec<T, 3>& yp = yObj.pvalue(i);
      const Vec<T, 3>& zp = zObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(zeros, yp, zp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3> zeros{};
      Vec<T, 3>& yh = yObj.hvalue(i);
      Vec<T, 3>& zh = zObj.hvalue(i);
//                operations:
      Mat3x3ToThreeVec3Core(Ch, zeros, yh, zh);
    }
  };

  const Vec<T, 3>& x;
  A2DVec<N, Vec<T, 3>>& yObj;
  A2DVec<N, Vec<T, 3>>& zObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline VecA2DVecA2DVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(const Vec<T, 3>& x,
                                                             A2DVec<N, Vec<T, 3>>& y,
                                                             A2DVec<N, Vec<T, 3>>& z,
                                                             A2DMat<N, Mat<T, 3, 3>>& C) {
  return VecA2DVecA2DVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class A2DVecVecA2DVecToMat3x3Expr {
 public:
  A2DVecVecA2DVecToMat3x3Expr(A2DVec<N, Vec<T, 3>>& xObj,
                              const Vec<T, 3>& y,
                              A2DVec<N, Vec<T, 3>>& zObj,
                              A2DMat<N, Mat<T, 3, 3>>& CObj)
      : xObj(xObj), y(y), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3> zeros{};
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, zeros, zb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3> zeros{};
      const Vec<T, 3>& zp = zObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(xp, zeros, zp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3> zeros{};
      Vec<T, 3>& zh = zObj.hvalue(i);
//                operations:
      Mat3x3ToThreeVec3Core(Ch, xh, zeros, zh);
    }
  };

  A2DVec<N, Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  A2DVec<N, Vec<T, 3>>& zObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DVecVecA2DVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(A2DVec<N, Vec<T, 3>>& x,
                                                             const Vec<T, 3>& y,
                                                             A2DVec<N, Vec<T, 3>>& z,
                                                             A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DVecVecA2DVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class A2DVecA2DVecVecToMat3x3Expr {
 public:
  A2DVecA2DVecVecToMat3x3Expr(A2DVec<N, Vec<T, 3>>& xObj,
                              A2DVec<N, Vec<T, 3>>& yObj,
                              const Vec<T, 3>& z,
                              A2DMat<N, Mat<T, 3, 3>>& CObj)
      : xObj(xObj), yObj(yObj), z(z), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, yb, zeros);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      const Vec<T, 3> zeros{};
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(xp, yp, zeros, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      Vec<T, 3> zeros{};
//                operations:
      Mat3x3ToThreeVec3Core(Ch, xh, yh, zeros);
    }
  };

  A2DVec<N, Vec<T, 3>>& xObj;
  A2DVec<N, Vec<T, 3>>& yObj;
  const Vec<T, 3>& z;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DVecA2DVecVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(A2DVec<N, Vec<T, 3>>& x,
                                                             A2DVec<N, Vec<T, 3>>& y,
                                                             const Vec<T, 3>& z,
                                                             A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DVecA2DVecVecToMat3x3Expr<N, T>(x, y, z, C);
}

template <int N, typename T>
class A2DVecA2DVecA2DVecToMat3x3Expr {
 public:
  A2DVecA2DVecA2DVecToMat3x3Expr(A2DVec<N, Vec<T, 3>>& xObj,
                                 A2DVec<N, Vec<T, 3>>& yObj,
                                 A2DVec<N, Vec<T, 3>>& zObj,
                                 A2DMat<N, Mat<T, 3, 3>>& CObj)
      : xObj(xObj), yObj(yObj), zObj(zObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& z = zObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromThreeVec3Core(x, y, z, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    Vec<T, 3>& zb = zObj.bvalue();
//            operations:
    Mat3x3ToThreeVec3Core(Cb, xb, yb, zb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      const Vec<T, 3>& zp = zObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromThreeVec3Core(xp, yp, zp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      Vec<T, 3>& zh = zObj.hvalue(i);
//                operations:
      Mat3x3ToThreeVec3Core(Ch, xh, yh, zh);
    }
  };

  A2DVec<N, Vec<T, 3>>& xObj;
  A2DVec<N, Vec<T, 3>>& yObj;
  A2DVec<N, Vec<T, 3>>& zObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DVecA2DVecA2DVecToMat3x3Expr<N, T> Mat3x3FromThreeVec3(A2DVec<N, Vec<T, 3>>& x,
                                                                A2DVec<N, Vec<T, 3>>& y,
                                                                A2DVec<N, Vec<T, 3>>& z,
                                                                A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DVecA2DVecA2DVecToMat3x3Expr<N, T>(x, y, z, C);
}

/*
  Mat3x3FromMat3x2AndVec3
  This represents the assembly of a 3x2 matrix and a 3 vector into a 3x3 matrix
  C(i, 0) = A(i, 0)
  C(i, 1) = A(i, 1)
  C(i, 2) = x(i)
*/
template <typename T>
inline void Mat3x3FromMat3x2AndVec3(const Mat<T, 3, 2>& A,
                                    const Vec<T, 3>& x,
                                    Mat<T, 3, 3>& C) {
  Mat3x3FromMat3x2AndVec3Core(A, x, C);
}

template <typename T>
class ADMatAndVecToMat3x3Expr {
 public:
  ADMatAndVecToMat3x3Expr(ADMat<Mat<T, 3, 2>>& AObj,
                          const Vec<T, 3>& x,
                          ADMat<Mat<T, 3, 3>>& CObj)
      : AObj(AObj), x(x), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void forward() {
//            input:
    const Mat<T, 3, 2>& Ab = AObj.bvalue();
    const Vec<T, 3> zeros{};
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromMat3x2AndVec3Core(Ab, zeros, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, Ab, zeros);
  };

  ADMat<Mat<T, 3, 2>>& AObj;
  const Vec<T, 3>& x;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADMatAndVecToMat3x3Expr<T> Mat3x3FromMat3x2AndVec3(ADMat<Mat<T, 3, 2>>& A,
                                                          const Vec<T, 3>& x,
                                                          ADMat<Mat<T, 3, 3>>& C) {
  return ADMatAndVecToMat3x3Expr<T>(A, x, C);
}

template <typename T>
class MatAndADVecToMat3x3Expr {
 public:
  MatAndADVecToMat3x3Expr(const Mat<T, 3, 2>& A,
                          ADVec<Vec<T, 3>>& xObj,
                          ADMat<Mat<T, 3, 3>>& CObj)
      : A(A), xObj(xObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void forward() {
//            input:
    const Mat<T, 3, 2> zeros{};
    const Vec<T, 3>& xb = xObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromMat3x2AndVec3Core(zeros, xb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2> zeros{};
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, zeros, xb);
  };

  const Mat<T, 3, 2>& A;
  ADVec<Vec<T, 3>>& xObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline MatAndADVecToMat3x3Expr<T> Mat3x3FromMat3x2AndVec3(const Mat<T, 3, 2>& A,
                                                          ADVec<Vec<T, 3>>& x,
                                                          ADMat<Mat<T, 3, 3>>& C) {
  return MatAndADVecToMat3x3Expr<T>(A, x, C);
}

template <typename T>
class ADMatAndADVecToMat3x3Expr {
 public:
  ADMatAndADVecToMat3x3Expr(ADMat<Mat<T, 3, 2>>& AObj,
                            ADVec<Vec<T, 3>>& xObj,
                            ADMat<Mat<T, 3, 3>>& CObj)
      : AObj(AObj), xObj(xObj), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void forward() {
//            input:
    const Mat<T, 3, 2>& Ab = AObj.bvalue();
    const Vec<T, 3>& xb = xObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromMat3x2AndVec3Core(Ab, xb, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, Ab, xb);
  };

  ADMat<Mat<T, 3, 2>>& AObj;
  ADVec<Vec<T, 3>>& xObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADMatAndADVecToMat3x3Expr<T> Mat3x3FromMat3x2AndVec3(ADMat<Mat<T, 3, 2>>& A,
                                                            ADVec<Vec<T, 3>>& x,
                                                            ADMat<Mat<T, 3, 3>>& C) {
  return ADMatAndADVecToMat3x3Expr<T>(A, x, C);
}

template <int N, typename T>
class A2DMatAndVecToMat3x3Expr {
 public:
  A2DMatAndVecToMat3x3Expr(A2DMat<N, Mat<T, 3, 2>>& AObj,
                           const Vec<T, 3>& x,
                           A2DMat<N, Mat<T, 3, 3>>& CObj)
      : AObj(AObj), x(x), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
    Vec<T, 3> zeros{};
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, Ab, zeros);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 2>& Ap = AObj.pvalue(i);
      const Vec<T, 3> zeros{};
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromMat3x2AndVec3Core(Ap, zeros, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Mat<T, 3, 2>& Ah = AObj.hvalue(i);
      Vec<T, 3> zeros{};
//                operations:
      Mat3x3ToMat3x2AndVec3Core(Ch, Ah, zeros);
    }
  };

  A2DMat<N, Mat<T, 3, 2>>& AObj;
  const Vec<T, 3>& x;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DMatAndVecToMat3x3Expr<N, T> Mat3x3FromMat3x2AndVec3(A2DMat<N, Mat<T, 3, 2>>& A,
                                                              const Vec<T, 3>& x,
                                                              A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DMatAndVecToMat3x3Expr<N, T>(A, x, C);
}

template <int N, typename T>
class MatAndA2DVecToMat3x3Expr {
 public:
  MatAndA2DVecToMat3x3Expr(const Mat<T, 3, 2>& A,
                           A2DVec<N, Vec<T, 3>>& xObj,
                           A2DMat<N, Mat<T, 3, 3>>& CObj)
      : A(A), xObj(xObj), CObj(CObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2> zeros{};
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, zeros, xb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 2> zeros{};
      const Vec<T, 3>& xp = xObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromMat3x2AndVec3Core(zeros, xp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Mat<T, 3, 2> zeros{};
      Vec<T, 3>& xh = xObj.hvalue(i);
//                operations:
      Mat3x3ToMat3x2AndVec3Core(Ch, zeros, xh);
    }
  };

  const Mat<T, 3, 2>& A;
  A2DVec<N, Vec<T, 3>>& xObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline MatAndA2DVecToMat3x3Expr<N, T> Mat3x3FromMat3x2AndVec3(const Mat<T, 3, 2>& A,
                                                              A2DVec<N, Vec<T, 3>>& x,
                                                              A2DMat<N, Mat<T, 3, 3>>& C) {
  return MatAndA2DVecToMat3x3Expr<N, T>(A, x, C);
}

template <int N, typename T>
class A2DMatAndA2DVecToMat3x3Expr {
 public:
  A2DMatAndA2DVecToMat3x3Expr(A2DMat<N, Mat<T, 3, 2>>& AObj,
                              A2DVec<N, Vec<T, 3>>& xObj,
                              A2DMat<N, Mat<T, 3, 3>>& CObj)
      : AObj(AObj), xObj(xObj), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2AndVec3Core(A, x, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Mat3x3ToMat3x2AndVec3Core(Cb, Ab, xb);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 2>& Ap = AObj.pvalue(i);
      const Vec<T, 3>& xp = xObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromMat3x2AndVec3Core(Ap, xp, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Mat<T, 3, 2>& Ah = AObj.hvalue(i);
      Vec<T, 3>& xh = xObj.hvalue(i);
//                operations:
      Mat3x3ToMat3x2AndVec3Core(Ch, Ah, xh);
    }
  };

  A2DMat<N, Mat<T, 3, 2>>& AObj;
  A2DVec<N, Vec<T, 3>>& xObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DMatAndA2DVecToMat3x3Expr<N, T> Mat3x3FromMat3x2AndVec3(A2DMat<N, Mat<T, 3, 2>>& A,
                                                                 A2DVec<N, Vec<T, 3>>& x,
                                                                 A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DMatAndA2DVecToMat3x3Expr<N, T>(A, x, C);
}

/*
  Mat3x3FromMat3x2
  This represents the casting of a 3x2 matrix to a 3x3 matrix
  C(i, 0) = A(i, 0)
  C(i, 1) = A(i, 1)
  C(i, 2) = 0
*/
template <typename T>
inline void Mat3x3FromMat3x2(const Mat<T, 3, 2>& A,
                             Mat<T, 3, 3>& C) {
  Mat3x3FromMat3x2Core(A, C);
}

template <typename T>
class ADMat3x2ToMat3x3Expr {
 public:
  ADMat3x2ToMat3x3Expr(ADMat<Mat<T, 3, 2>>& AObj,
                       ADMat<Mat<T, 3, 3>>& CObj)
      : AObj(AObj), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2Core(A, C);
  }

  void forward() {
//            input:
    const Mat<T, 3, 2>& Ab = AObj.bvalue();
//            output:
    Mat<T, 3, 3>& Cb = CObj.bvalue();
//            operations:
    Mat3x3FromMat3x2Core(Ab, Cb);
  };

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
//            operations:
    Mat3x3ToMat3x2Core(Cb, Ab);
  };

  ADMat<Mat<T, 3, 2>>& AObj;
  ADMat<Mat<T, 3, 3>>& CObj;
};

template <typename T>
inline ADMat3x2ToMat3x3Expr<T> Mat3x3FromMat3x2(ADMat<Mat<T, 3, 2>>& A,
                                                ADMat<Mat<T, 3, 3>>& C) {
  return ADMat3x2ToMat3x3Expr<T>(A, C);
}

template <int N, typename T>
class A2DMat3x2ToMat3x3Expr {
 public:
  A2DMat3x2ToMat3x3Expr(A2DMat<N, Mat<T, 3, 2>>& AObj,
                        A2DMat<N, Mat<T, 3, 3>>& CObj)
      : AObj(AObj), CObj(CObj) {
    const Mat<T, 3, 2>& A = AObj.value();
    Mat<T, 3, 3>& C = CObj.value();
    Mat3x3FromMat3x2Core(A, C);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Cb = CObj.bvalue();
//            output:
    Mat<T, 3, 2>& Ab = AObj.bvalue();
//            operations:
    Mat3x3ToMat3x2Core(Cb, Ab);
  };

  void hforward() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 2>& Ap = AObj.pvalue(i);
//                output:
      Mat<T, 3, 3>& Cp = CObj.pvalue(i);
//                operations:
      Mat3x3FromMat3x2Core(Ap, Cp);
    }
  };

  void hreverse() {
//            main loop:
    for (int i = 0; i < N; ++i) {
//                input:
      const Mat<T, 3, 3>& Ch = CObj.hvalue(i);
//                output:
      Mat<T, 3, 2>& Ah = AObj.hvalue(i);
//                operations:
      Mat3x3ToMat3x2Core(Ch, Ah);
    }
  };

  A2DMat<N, Mat<T, 3, 2>>& AObj;
  A2DMat<N, Mat<T, 3, 3>>& CObj;
};

template <int N, typename T>
inline A2DMat3x2ToMat3x3Expr<N, T> Mat3x3FromMat3x2(A2DMat<N, Mat<T, 3, 2>>& A,
                                                    A2DMat<N, Mat<T, 3, 3>>& C) {
  return A2DMat3x2ToMat3x3Expr<N, T>(A, C);
}

}  // namespace A2D

#endif //A2D_MAT_VEC_OPS_3D_H
