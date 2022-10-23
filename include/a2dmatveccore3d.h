#ifndef A2D_MAT_VEC_CORE_3D_H
#define A2D_MAT_VEC_CORE_3D_H

#include "a2dobjs.h"

namespace A2D {

template <typename T>
inline void Mat3x2ToVec3Core(const Mat<T, 3, 2>& A,
                             Vec<T, 3>& u,
                             Vec<T, 3>& v) {
  u(0) = A(0, 0);
  u(1) = A(1, 0);
  u(2) = A(2, 0);
  v(0) = A(0, 1);
  v(1) = A(1, 1);
  v(2) = A(2, 1);
}

template <typename T>
inline void Mat3x2FromTwoVec3Core(const Vec<T, 3>& x,
                                  const Vec<T, 3>& y,
                                  Mat<T, 3, 2>& C) {
  C(0, 0) = x(0);
  C(1, 0) = x(1);
  C(2, 0) = x(2);
  C(0, 1) = y(0);
  C(1, 1) = y(1);
  C(2, 1) = y(2);
}

template <typename T>
inline void Mat3x3FromThreeVec3Core(const Vec<T, 3>& x,
                                    const Vec<T, 3>& y,
                                    const Vec<T, 3>& z,
                                    Mat<T, 3, 3>& C) {
  C(0, 0) = x(0);
  C(1, 0) = x(1);
  C(2, 0) = x(2);
  C(0, 1) = y(0);
  C(1, 1) = y(1);
  C(2, 1) = y(2);
  C(0, 2) = z(0);
  C(1, 2) = z(1);
  C(2, 2) = z(2);
}

template <typename T>
inline void Mat3x3ToThreeVec3Core(const Mat<T, 3, 3>& A,
                                  Vec<T, 3>& u,
                                  Vec<T, 3>& v,
                                  Vec<T, 3>& w) {
  u(0) = A(0, 0);
  u(1) = A(1, 0);
  u(2) = A(2, 0);
  v(0) = A(0, 1);
  v(1) = A(1, 1);
  v(2) = A(2, 1);
  w(0) = A(0, 2);
  w(1) = A(1, 2);
  w(2) = A(2, 2);
}

template <typename T>
inline void Mat3x3FromMat3x2AndVec3Core(const Mat<T, 3, 2>& A,
                                        const Vec<T, 3>& x,
                                        Mat<T, 3, 3>& C) {
  C(0, 0) = A(0, 0);
  C(1, 0) = A(1, 0);
  C(2, 0) = A(2, 0);
  C(0, 1) = A(0, 1);
  C(1, 1) = A(1, 1);
  C(2, 1) = A(2, 1);
  C(0, 2) = x(0);
  C(1, 2) = x(1);
  C(2, 2) = x(2);
}

template <typename T>
inline void Mat3x3ToMat3x2AndVec3Core(const Mat<T, 3, 3>& A,
                                      Mat<T, 3, 2>& C,
                                      Vec<T, 3>& u) {
  C(0, 0) = A(0, 0);
  C(1, 0) = A(1, 0);
  C(2, 0) = A(2, 0);
  C(0, 1) = A(0, 1);
  C(1, 1) = A(1, 1);
  C(2, 1) = A(2, 1);
  u(0) = A(0, 2);
  u(1) = A(1, 2);
  u(2) = A(2, 2);
}

template <typename T>
inline void Mat3x3FromMat3x2Core(const Mat<T, 3, 2>& A,
                                 Mat<T, 3, 3>& C) {
  C(0, 0) = A(0, 0);
  C(1, 0) = A(1, 0);
  C(2, 0) = A(2, 0);
  C(0, 1) = A(0, 1);
  C(1, 1) = A(1, 1);
  C(2, 1) = A(2, 1);
  C(0, 2) = 0;
  C(1, 2) = 0;
  C(2, 2) = 0;
}

template <typename T>
inline void Mat3x3ToMat3x2Core(const Mat<T, 3, 3>& A,
                               Mat<T, 3, 2>& C) {
  C(0, 0) = A(0, 0);
  C(1, 0) = A(1, 0);
  C(2, 0) = A(2, 0);
  C(0, 1) = A(0, 1);
  C(1, 1) = A(1, 1);
  C(2, 1) = A(2, 1);
}

}  // namespace A2D

#endif //A2D_MAT_VEC_CORE_3D_H
