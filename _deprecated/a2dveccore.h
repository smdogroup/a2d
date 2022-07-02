#ifndef A2D_VEC_CORE_H
#define A2D_VEC_CORE_H

namespace A2D {

  template<typename T>
  inline T Vec3DotCore( const T x[],
                        const T y[] ){
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
  }

  template<typename T>
  inline T Vec3NormCore( const T x[] ){
    return sqrt(Vec3DotCore(x, x));
  }

  template<typename T>
  inline void Vec3ScaleCore( const T alpha, const T x[], T v[] ){
    v[0] = alpha * x[0];
    v[1] = alpha * x[1];
    v[2] = alpha * x[2];
  }

  template<typename T>
  inline void Vec3AXPYCore( const T alpha, const T x[],
                            const T y[], T v[] ){
    v[0] = alpha * x[0] + y[0];
    v[1] = alpha * x[1] + y[1];
    v[2] = alpha * x[2] + y[2];
  }

  template<typename T>
  inline void Vec3CrossProductCore( const T x[], const T y[],
                                    T v[] ){

    v[0] = x[1] * y[2] - x[2] * y[1];
    v[1] = x[2] * y[0] - x[0] * y[2];
    v[2] = x[0] * y[1] - x[1] * y[0];
  }

  template<typename T>
  inline void Vec3CrossProductAddCore( const T x[], const T y[],
                                       T v[] ){

    v[0] += x[1] * y[2] - x[2] * y[1];
    v[1] += x[2] * y[0] - x[0] * y[2];
    v[2] += x[0] * y[1] - x[1] * y[0];
  }

  template<typename T>
  inline void Vec3OuterProductCore( const T x[], const T y[],
                                    T A[] ){
    A[0] = x[0] * y[0];
    A[1] = x[0] * y[1];
    A[2] = x[0] * y[2];
    A[3] = x[1] * y[0];
    A[4] = x[1] * y[1];
    A[5] = x[1] * y[2];
    A[6] = x[2] * y[0];
    A[7] = x[2] * y[1];
    A[8] = x[2] * y[2];
  }

  template<typename T>
  inline void Vec3OuterProductAddCore( const T x[], const T y[],
                                       T A[] ){
    A[0] += x[0] * y[0];
    A[1] += x[0] * y[1];
    A[2] += x[0] * y[2];
    A[3] += x[1] * y[0];
    A[4] += x[1] * y[1];
    A[5] += x[1] * y[2];
    A[6] += x[2] * y[0];
    A[7] += x[2] * y[1];
    A[8] += x[2] * y[2];
  }

  template<typename T>
  inline void Vec3OuterProductScaleCore( const T scale,
                                         const T x[], const T y[],
                                         T A[] ){
    A[0] = scale * x[0] * y[0];
    A[1] = scale * x[0] * y[1];
    A[2] = scale * x[0] * y[2];
    A[3] = scale * x[1] * y[0];
    A[4] = scale * x[1] * y[1];
    A[5] = scale * x[1] * y[2];
    A[6] = scale * x[2] * y[0];
    A[7] = scale * x[2] * y[1];
    A[8] = scale * x[2] * y[2];
  }

  template<typename T>
  inline void Vec3OuterProductAddScaleCore( const T scale,
                                            const T x[], const T y[],
                                            T A[] ){
    A[0] += scale * x[0] * y[0];
    A[1] += scale * x[0] * y[1];
    A[2] += scale * x[0] * y[2];
    A[3] += scale * x[1] * y[0];
    A[4] += scale * x[1] * y[1];
    A[5] += scale * x[1] * y[2];
    A[6] += scale * x[2] * y[0];
    A[7] += scale * x[2] * y[1];
    A[8] += scale * x[2] * y[2];
  }

  template<typename T>
  inline void Mat3x3VecMultCore( const T A[],
                                 const T x[],
                                 T y[] ){
    y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
    y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
    y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
  }

  template<typename T>
  inline void Mat3x3VecMultScaleCore( const T scale,
                                      const T A[],
                                      const T x[],
                                      T y[] ){
    y[0] = scale * (A[0] * x[0] + A[1] * x[1] + A[2] * x[2]);
    y[1] = scale * (A[3] * x[0] + A[4] * x[1] + A[5] * x[2]);
    y[2] = scale * (A[6] * x[0] + A[7] * x[1] + A[8] * x[2]);
  }

  template<typename T>
  inline void Mat3x3VecMultAddCore( const T A[],
                                    const T x[],
                                    T y[] ){
    y[0] += A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
    y[1] += A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
    y[2] += A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
  }

  template<typename T>
  inline void Mat3x3VecMultAddScaleCore( const T scale,
                                         const T A[],
                                         const T x[],
                                         T y[] ){
    y[0] += scale * (A[0] * x[0] + A[1] * x[1] + A[2] * x[2]);
    y[1] += scale * (A[3] * x[0] + A[4] * x[1] + A[5] * x[2]);
    y[2] += scale * (A[6] * x[0] + A[7] * x[1] + A[8] * x[2]);
  }

  template<typename T>
  inline void MatTrans3x3VecMultCore( const T A[],
                                      const T x[],
                                      T y[] ){
    y[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
    y[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
    y[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
  }

  template<typename T>
  inline void MatTrans3x3VecMultScaleCore( const T scale,
                                           const T A[],
                                           const T x[],
                                           T y[] ){
    y[0] = scale * (A[0] * x[0] + A[3] * x[1] + A[6] * x[2]);
    y[1] = scale * (A[1] * x[0] + A[4] * x[1] + A[7] * x[2]);
    y[2] = scale * (A[2] * x[0] + A[5] * x[1] + A[8] * x[2]);
  }

  template<typename T>
  inline void MatTrans3x3VecMultAddCore( const T A[],
                                         const T x[],
                                         T y[] ){
    y[0] += A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
    y[1] += A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
    y[2] += A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
  }

  template<typename T>
  inline void MatTrans3x3VecMultAddScaleCore( const T scale,
                                              const T A[],
                                              const T x[],
                                              T y[] ){
    y[0] += scale * (A[0] * x[0] + A[3] * x[1] + A[6] * x[2]);
    y[1] += scale * (A[1] * x[0] + A[4] * x[1] + A[7] * x[2]);
    y[2] += scale * (A[2] * x[0] + A[5] * x[1] + A[8] * x[2]);
  }

  template<typename T>
  inline T Mat3x3InnerProductCore( const T A[],
                                   const T x[],
                                   const T y[] ){
    return
      x[0] * (A[0] * y[0] + A[1] * y[1] + A[2] * y[2]) +
      x[1] * (A[3] * y[0] + A[4] * y[1] + A[5] * y[2]) +
      x[2] * (A[6] * y[0] + A[7] * y[1] + A[8] * y[2]);
  }

} // namespace AD

#endif // A2D_VEC_CORE_H
