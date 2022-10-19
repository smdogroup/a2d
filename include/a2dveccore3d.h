#ifndef A2D_VEC_CORE_H
#define A2D_VEC_CORE_H

namespace A2D {

/**
 * @brief vector addition:  v = x + y
 * @param x    :  3-Vector to add
 * @param y    :  3-Vector to add
 * @param v    :  3-Vector destination
 */
template <typename T, class VecType>
inline void Vec3AddCore(const VecType& x,
                        const VecType& y,
                        VecType& v) {
  v(0) = x(0) + y(0);
  v(1) = x(1) + y(1);
  v(2) = x(2) + y(2);
}

/**
 * @brief vector addition in place:  v += x
 * @param x    :  3-Vector to add
 * @param y    :  3-Vector to add
 * @param v    :  3-Vector destination
 */
template <class VecType>
inline void Vec3AddInPlaceCore(const VecType& x,
                               VecType& v) {
  v(0) += x(0);
  v(1) += x(1);
  v(2) += x(2);
}

/**
 * @brief scalar-vector multiply and add:  v = a * x + b * y
 * @param alpha:  Scalar to multiply x
 * @param x    :  3-Vector
 * @param beta :  Scalar to multiply y
 * @param y    :  3-Vector
 * @param v    :  3-Vector to store result
 */
template <typename T, class VecType>
inline void Vec3AXPBYCore(const T alpha,
                          const VecType& x,
                          const T beta,
                          const VecType& y,
                          VecType& v) {
  v(0) = (alpha * x(0)) + (beta * y(0));
  v(1) = (alpha * x(1)) + (beta * y(1));
  v(2) = (alpha * x(2)) + (beta * y(2));
}

/**
 * @brief vector add then multiply (in place):  v = (v + x) * a
 * @param alpha:  Scalar to multiply result
 * @param x    :  3-Vector
 * @param v    :  3-Vector to add to then multiply
 */
template <typename T, class VecType>
inline void Vec3AddThenScaleCore(const T alpha,
                                 const VecType& x,
                                 VecType& v) {
  v(0) = alpha * (x(0) + v(0));
  v(1) = alpha * (x(1) + v(1));
  v(2) = alpha * (x(2) + v(2));
}

/**
 * @brief vector scale then add (in place):  v = v + (a * x)
 * @param alpha:  Scalar to multiply x vector
 * @param x    :  3-Vector to scale
 * @param v    :  3-Vector to add scaled vector
 */
template <typename T, class VecType>
inline void Vec3ScaleAndAddCore(const T alpha,
                                const VecType& x,
                                VecType& v) {
  v(0) += alpha * x(0);
  v(1) += alpha * x(1);
  v(2) += alpha * x(2);
}

template <typename T, class VecType>
inline T Vec3DotCore(const VecType& x, const VecType& y) {
  return (x(0) * y(0) + x(1) * y(1) + x(2) * y(2));
}

template <typename T, class VecType>
inline T Vec3NormCore(const VecType& x) {
  return sqrt(Vec3DotCore<T>(x, x));
}

template <typename T, class VecType>
inline void Vec3ScaleCore(const T alpha, const VecType& x, VecType& v) {
  v(0) = alpha * x(0);
  v(1) = alpha * x(1);
  v(2) = alpha * x(2);
}

template <typename T, class VecType>
inline void Vec3AXPYCore(const T alpha, const VecType& x, const VecType& y,
                         VecType& v) {
  v(0) = alpha * x(0) + y(0);
  v(1) = alpha * x(1) + y(1);
  v(2) = alpha * x(2) + y(2);
}

template <class VecType>
inline void Vec3CrossProductCore(const VecType& x, const VecType& y,
                                 VecType& v) {
  v(0) = x(1) * y(2) - x(2) * y(1);
  v(1) = x(2) * y(0) - x(0) * y(2);
  v(2) = x(0) * y(1) - x(1) * y(0);
}

template <class VecType>
inline void Vec3CrossProductAddCore(const VecType& x, const VecType& y,
                                    VecType& v) {
  v(0) += x(1) * y(2) - x(2) * y(1);
  v(1) += x(2) * y(0) - x(0) * y(2);
  v(2) += x(0) * y(1) - x(1) * y(0);
}

template <class VecType, class MatType>
inline void Vec3OuterProductCore(const VecType& x, const VecType& y,
                                 MatType& A) {
  A(0, 0) = x(0) * y(0);
  A(0, 1) = x(0) * y(1);
  A(0, 2) = x(0) * y(2);
  A(1, 0) = x(1) * y(0);
  A(1, 1) = x(1) * y(1);
  A(1, 2) = x(1) * y(2);
  A(2, 0) = x(2) * y(0);
  A(2, 1) = x(2) * y(1);
  A(2, 2) = x(2) * y(2);
}

template <class VecType, class MatType>
inline void Vec3OuterProductAddCore(const VecType& x, const VecType& y,
                                    MatType& A) {
  A(0, 0) += x(0) * y(0);
  A(0, 1) += x(0) * y(1);
  A(0, 2) += x(0) * y(2);
  A(1, 0) += x(1) * y(0);
  A(1, 1) += x(1) * y(1);
  A(1, 2) += x(1) * y(2);
  A(2, 0) += x(2) * y(0);
  A(2, 1) += x(2) * y(1);
  A(2, 2) += x(2) * y(2);
}

template <typename T, class VecType, class MatType>
inline void Vec3OuterProductScaleCore(const T scale, const VecType& x,
                                      const VecType& y, MatType& A) {
  A(0, 0) = scale * x(0) * y(0);
  A(0, 1) = scale * x(0) * y(1);
  A(0, 2) = scale * x(0) * y(2);
  A(1, 0) = scale * x(1) * y(0);
  A(1, 1) = scale * x(1) * y(1);
  A(1, 2) = scale * x(1) * y(2);
  A(2, 0) = scale * x(2) * y(0);
  A(2, 1) = scale * x(2) * y(1);
  A(2, 2) = scale * x(2) * y(2);
}

template <typename T, class VecType, class MatType>
inline void Vec3OuterProductAddScaleCore(const T scale, const VecType& x,
                                         const VecType& y, MatType& A) {
  A(0, 0) += scale * x(0) * y(0);
  A(0, 1) += scale * x(0) * y(1);
  A(0, 2) += scale * x(0) * y(2);
  A(1, 0) += scale * x(1) * y(0);
  A(1, 1) += scale * x(1) * y(1);
  A(1, 2) += scale * x(1) * y(2);
  A(2, 0) += scale * x(2) * y(0);
  A(2, 1) += scale * x(2) * y(1);
  A(2, 2) += scale * x(2) * y(2);
}

template <class VecType, class MatType>
inline void Mat3x3VecMultCore(const MatType& A, const VecType& x, VecType& y) {
  y(0) = A(0, 0) * x(0) + A(0, 1) * x(1) + A(0, 2) * x(2);
  y(1) = A(1, 0) * x(0) + A(1, 1) * x(1) + A(1, 2) * x(2);
  y(2) = A(2, 0) * x(0) + A(2, 1) * x(1) + A(2, 2) * x(2);
}

template <typename T, class VecType, class MatType>
inline void Mat3x3VecMultScaleCore(const T scale, const MatType& A,
                                   const VecType& x, VecType& y) {
  y(0) = scale * (A(0, 0) * x(0) + A(0, 1) * x(1) + A(0, 2) * x(2));
  y(1) = scale * (A(1, 0) * x(0) + A(1, 1) * x(1) + A(1, 2) * x(2));
  y(2) = scale * (A(2, 0) * x(0) + A(2, 1) * x(1) + A(2, 2) * x(2));
}

template <typename T, class VecType, class MatType>
inline void Mat3x3VecMultAddCore(const MatType& A, const VecType& x,
                                 VecType& y) {
  y(0) += A(0, 0) * x(0) + A(0, 1) * x(1) + A(0, 2) * x(2);
  y(1) += A(1, 0) * x(0) + A(1, 1) * x(1) + A(1, 2) * x(2);
  y(2) += A(2, 0) * x(0) + A(2, 1) * x(1) + A(2, 2) * x(2);
}

template <typename T, class VecType, class MatType>
inline void Mat3x3VecMultAddScaleCore(const T scale, const MatType& A,
                                      const VecType& x, VecType& y) {
  y(0) += scale * (A(0, 0) * x(0) + A(0, 1) * x(1) + A(0, 2) * x(2));
  y(1) += scale * (A(1, 0) * x(0) + A(1, 1) * x(1) + A(1, 2) * x(2));
  y(2) += scale * (A(2, 0) * x(0) + A(2, 1) * x(1) + A(2, 2) * x(2));
}

template <class VecType, class MatType>
inline void MatTrans3x3VecMultCore(const MatType& A, const VecType& x,
                                   VecType& y) {
  y(0) = A(0, 0) * x(0) + A(1, 0) * x(1) + A(2, 0) * x(2);
  y(1) = A(0, 1) * x(0) + A(1, 1) * x(1) + A(2, 1) * x(2);
  y(2) = A(0, 2) * x(0) + A(1, 2) * x(1) + A(2, 2) * x(2);
}

template <typename T, class VecType, class MatType>
inline void MatTrans3x3VecMultScaleCore(const T scale, const MatType& A,
                                        const VecType& x, VecType& y) {
  y(0) = scale * (A(0, 0) * x(0) + A(1, 0) * x(1) + A(2, 0) * x(2));
  y(1) = scale * (A(0, 1) * x(0) + A(1, 1) * x(1) + A(2, 1) * x(2));
  y(2) = scale * (A(0, 2) * x(0) + A(1, 2) * x(1) + A(2, 2) * x(2));
}

template <typename T, class VecType, class MatType>
inline void MatTrans3x3VecMultAddCore(const MatType& A, const VecType& x,
                                      VecType& y) {
  y(0) += A(0, 0) * x(0) + A(1, 0) * x(1) + A(2, 0) * x(2);
  y(1) += A(0, 1) * x(0) + A(1, 1) * x(1) + A(2, 1) * x(2);
  y(2) += A(0, 2) * x(0) + A(1, 2) * x(1) + A(2, 2) * x(2);
}

template <typename T, class VecType, class MatType>
inline void MatTrans3x3VecMultAddScaleCore(const T scale, const MatType& A,
                                           const VecType& x, VecType& y) {
  y(0) += scale * (A(0, 0) * x(0) + A(1, 0) * x(1) + A(2, 0) * x(2));
  y(1) += scale * (A(0, 1) * x(0) + A(1, 1) * x(1) + A(2, 1) * x(2));
  y(2) += scale * (A(0, 2) * x(0) + A(1, 2) * x(1) + A(2, 2) * x(2));
}

template <typename T, class VecType, class MatType>
inline T Mat3x3InnerProductCore(const MatType& A, const VecType& x,
                                const VecType& y) {
  return x(0) * (A(0, 0) * y(0) + A(0, 1) * y(1) + A(0, 2) * y(2)) +
      x(1) * (A(1, 0) * y(0) + A(1, 1) * y(1) + A(1, 2) * y(2)) +
      x(2) * (A(2, 0) * y(0) + A(2, 1) * y(1) + A(2, 2) * y(2));
}

}  // namespace A2D

#endif  // A2D_VEC_CORE_H
