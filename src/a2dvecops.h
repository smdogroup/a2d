#ifndef A2D_VEC_OPS_H
#define A2D_VEC_OPS_H

#include "a2dobjs.h"
#include "a2dveccore.h"

namespace A2D {

  /*
    Vec3Norm
  */
  template<typename T>
  class Vec3Norm {
  public:
    Vec3Norm( const Vec3<T>& x, Scalar<T>& alpha ){
      alpha.value = Vec3NormCore(x.x);
    }
  };

  template<typename T>
  class ADVec3Norm {
  public:
    ADVec3Norm( ADVec3<T>& x, ADScalar<T>& alpha ) : x(x), alpha(alpha) {
      alpha.value = Vec3NormCore(x.x);
    }
    void forward(){
      if (alpha.value != 0.0){
        alpha.valued = Vec3DotCore(x.x, x.xd)/alpha.value;
      }
    }
    void reverse(){
      if (alpha.value != 0.0){
        T ainv = alpha.valued/alpha.value;
        x.xd[0] += ainv * x.x[0];
        x.xd[1] += ainv * x.x[1];
        x.xd[2] += ainv * x.x[2];
      }
    }

    ADVec3<T>& x;
    ADScalar<T>& alpha;
  };

  /*
    Vec3Scale
  */
  template<typename T>
  class Vec3Scale {
  public:
    Vec3Scale( const Scalar<T>& alpha, Vec3<T>& x, Vec3<T>& v ){
      Vec3ScaleCore(alpha.value, x.x, v.x);
    }
  };

  template<typename T>
  class ADVec3Scale {
  public:
    ADVec3Scale( ADScalar<T>& alpha, ADVec3<T>& x, ADVec3<T>& v ) : alpha(alpha), x(x), v(v) {
      Vec3ScaleCore(alpha.value, x.x, v.x);
    }
    void forward(){
      v.xd[0] = alpha.valued * x.x[0] + alpha.value * x.xd[0];
      v.xd[1] = alpha.valued * x.x[1] + alpha.value * x.xd[1];
      v.xd[2] = alpha.valued * x.x[2] + alpha.value * x.xd[2];
    }
    void reverse(){
      alpha.valued += Vec3DotCore(v.xd, x.x);
      x.xd[0] += alpha.value * v.xd[0];
      x.xd[1] += alpha.value * v.xd[1];
      x.xd[2] += alpha.value * v.xd[2];
    }

    ADScalar<T>& alpha;
    ADVec3<T>& x;
    ADVec3<T>& v;
  };

  /*
    Vec3Axpy
  */
  template<typename T>
  class Vec3Axpy {
  public:
    Vec3Axpy( const Scalar<T>& alpha, const Vec3<T>& x, const Vec3<T>& y, Vec3<T>& v ){
      Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
    }
    Vec3Axpy( const T scale, const Scalar<T>& alpha, const Vec3<T>& x, const Vec3<T>& y, Vec3<T>& v ){
      Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
    }
  };

  template<typename T>
  class Vec3VecADScalarAxpy {
  public:
    Vec3VecADScalarAxpy( ADScalar<T>& alpha, const Vec3<T>& x, const Vec3<T>& y, ADVec3<T>& v ) : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
    }
    Vec3VecADScalarAxpy( const T scale, ADScalar<T>& alpha, const Vec3<T>& x, const Vec3<T>& y, ADVec3<T>& v ) : scale(scale), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
    }
    void forward(){
      v.xd[0] = scale * alpha.valued * x.x[0];
      v.xd[1] = scale * alpha.valued * x.x[1];
      v.xd[2] = scale * alpha.valued * x.x[2];
    }
    void reverse(){
      alpha.valued += scale * Vec3DotCore(x.x, v.xd);
    }

    const T scale;
    ADScalar<T>& alpha;
    const Vec3<T>& x;
    const Vec3<T>& y;
    ADVec3<T>& v;
  };

  template<typename T>
  class ADVec3VecADScalarAxpy {
  public:
    ADVec3VecADScalarAxpy( ADScalar<T>& alpha, ADVec3<T>& x, const Vec3<T>& y, ADVec3<T>& v ) : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
    }
    ADVec3VecADScalarAxpy( const T scale, ADScalar<T>& alpha, ADVec3<T>& x, const Vec3<T>& y, ADVec3<T>& v ) : scale(scale), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
    }
    void forward(){
      v.xd[0] = scale * (alpha.valued * x.x[0] + alpha.value * x.xd[0]);
      v.xd[1] = scale * (alpha.valued * x.x[1] + alpha.value * x.xd[1]);
      v.xd[2] = scale * (alpha.valued * x.x[2] + alpha.value * x.xd[2]);
    }
    void reverse(){
      alpha.valued += scale * Vec3DotCore(x.x, v.xd);
      x.xd[0] += scale * alpha.value * v.xd[0];
      x.xd[1] += scale * alpha.value * v.xd[1];
      x.xd[2] += scale * alpha.value * v.xd[2];
    }

    const T scale;
    ADScalar<T>& alpha;
    ADVec3<T>& x;
    const Vec3<T>& y;
    ADVec3<T>& v;
  };

  template<typename T>
  class ADVec3ADVecScalarAxpy {
  public:
    ADVec3ADVecScalarAxpy( const Scalar<T>& alpha, ADVec3<T>& x, ADVec3<T>& y, ADVec3<T>& v ) : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
    }
    ADVec3ADVecScalarAxpy( const T scale, const Scalar<T>& alpha, ADVec3<T>& x, ADVec3<T>& y, ADVec3<T>& v ) : scale(scale), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
    }
    void forward(){
      v.xd[0] = scale * (alpha.value * x.xd[0]) + y.xd[0];
      v.xd[1] = scale * (alpha.value * x.xd[1]) + y.xd[1];
      v.xd[2] = scale * (alpha.value * x.xd[2]) + y.xd[2];
    }
    void reverse(){
      x.xd[0] += scale * alpha.value * v.xd[0];
      x.xd[1] += scale * alpha.value * v.xd[1];
      x.xd[2] += scale * alpha.value * v.xd[2];
      y.xd[0] += v.xd[0];
      y.xd[1] += v.xd[1];
      y.xd[2] += v.xd[2];
    }

    const T scale;
    const Scalar<T>& alpha;
    ADVec3<T>& x;
    ADVec3<T>& y;
    ADVec3<T>& v;
  };

  template<typename T>
  class ADVec3Axpy {
  public:
    ADVec3Axpy( ADScalar<T>& alpha, ADVec3<T>& x, ADVec3<T>& y, ADVec3<T>& v ) : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
    }
    ADVec3Axpy( const T scale, ADScalar<T>& alpha, ADVec3<T>& x, ADVec3<T>& y, ADVec3<T>& v ) : scale(scale), alpha(alpha), x(x), y(y), v(v) {
      Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
    }
    void forward(){
      v.xd[0] = scale * (alpha.valued * x.x[0] + alpha.value * x.xd[0]) + y.xd[0];
      v.xd[1] = scale * (alpha.valued * x.x[1] + alpha.value * x.xd[1]) + y.xd[1];
      v.xd[2] = scale * (alpha.valued * x.x[2] + alpha.value * x.xd[2]) + y.xd[2];
    }
    void reverse(){
      alpha.valued += scale * Vec3DotCore(x.x, v.xd);
      x.xd[0] += scale * alpha.value * v.xd[0];
      x.xd[1] += scale * alpha.value * v.xd[1];
      x.xd[2] += scale * alpha.value * v.xd[2];
      y.xd[0] += v.xd[0];
      y.xd[1] += v.xd[1];
      y.xd[2] += v.xd[2];
    }

    const T scale;
    ADScalar<T>& alpha;
    ADVec3<T>& x;
    ADVec3<T>& y;
    ADVec3<T>& v;
  };

  /*
    Vec3Dot
  */
  template<typename T>
  class Vec3Dot {
  public:
    Vec3Dot( const Vec3<T>& x, const Vec3<T>& y, Scalar<T>& alpha ){
      alpha.value = Vec3DotCore(x.x, y.x);
    }
    Vec3Dot( const T scale, const Vec3<T>& x, const Vec3<T>& y, Scalar<T>& alpha ){
      alpha.value = scale * Vec3DotCore(x.x, y.x);
    }
  };

  template<typename T>
  class Vec3ADVecDot {
  public:
    Vec3ADVecDot( const Vec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : scale(1.0), x(x), y(y), alpha(alpha) {
      alpha.value = Vec3DotCore(x.x, y.x);
    }
    Vec3ADVecDot( const T scale, const Vec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : scale(scale), x(x), y(y), alpha(alpha) {
      alpha.value = scale * Vec3DotCore(x.x, y.x);
    }
    void forward(){
      alpha.valued = scale * Vec3DotCore(x.x, y.xd);
    }
    void reverse(){
      T s = scale * alpha.valued;
      y.xd[0] += s * x.x[0];
      y.xd[1] += s * x.x[1];
      y.xd[2] += s * x.x[2];
    }

    const T scale;
    const Vec3<T>& x;
    ADVec3<T>& y;
    ADScalar<T>& alpha;
  };

  template<typename T>
  class ADVec3Dot {
  public:
    ADVec3Dot( ADVec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : scale(1.0), x(x), y(y), alpha(alpha) {
      alpha.value = Vec3DotCore(x.x, y.x);
    }
    ADVec3Dot( const T scale, ADVec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : scale(scale), x(x), y(y), alpha(alpha) {
      alpha.value = scale * Vec3DotCore(x.x, y.x);
    }
    void forward(){
      alpha.valued = scale * (Vec3DotCore(x.x, y.xd) + Vec3DotCore(x.xd, y.x));
    }
    void reverse(){
      T s = scale * alpha.valued;
      x.xd[0] += s * y.x[0];
      x.xd[1] += s * y.x[1];
      x.xd[2] += s * y.x[2];

      y.xd[0] += s * x.x[0];
      y.xd[1] += s * x.x[1];
      y.xd[2] += s * x.x[2];
    }

    const T scale;
    ADVec3<T>& x;
    ADVec3<T>& y;
    ADScalar<T>& alpha;
  };

  /*
    Vec3CrossProduct
  */
  template<typename T>
  class Vec3CrossProduct {
  public:
    Vec3CrossProduct( const Vec3<T>& x, const Vec3<T>& y, Vec3<T>& v ){
      Vec3CrossProductCore(x.x, y.x, v.x);
    }
  };

  template<typename T>
  class ADVec3CrossProduct {
  public:
    ADVec3CrossProduct( ADVec3<T>& x, ADVec3<T>& y, ADVec3<T>& v ) : x(x), y(y), v(v) {
      Vec3CrossProductCore(x.x, y.x, v.x);
    }
    void forward(){
      Vec3CrossProductCore(x.xd, y.x, v.xd);
      Vec3CrossProductAddCore(x.x, y.xd, v.xd);
    }
    void reverse(){
      Vec3CrossProductAddCore(y.x, v.xd, x.xd);
      Vec3CrossProductAddCore(v.xd, x.x, y.xd);
    }

    ADVec3<T>& x;
    ADVec3<T>& y;
    ADVec3<T>& v;
  };

  /*
    Vec3 normalize
  */
  template<typename T>
  class Vec3Normalize {
  public:
    Vec3Normalize( const Vec3<T>& x, Vec3<T>& y ){
      T alpha = Vec3DotCore(x.x, x.x);
      if (alpha != 0.0){
        T inv = 1.0/sqrt(alpha);
        y.x[0] = inv * x.x[0];
        y.x[1] = inv * x.x[1];
        y.x[2] = inv * x.x[2];
      }
    }
  };

  template<typename T>
  class ADVec3Normalize {
  public:
    ADVec3Normalize( ADVec3<T>& x, ADVec3<T>& y ) : x(x), y(y) {
      alpha = Vec3DotCore(x.x, x.x);
      if (alpha != 0.0){
        inv = 1.0/sqrt(alpha);
        y.x[0] = inv * x.x[0];
        y.x[1] = inv * x.x[1];
        y.x[2] = inv * x.x[2];
      }
      else {
        inv = 0.0;
      }
    }
    void forward(){
      T beta = Vec3DotCore(x.x, x.xd);
      T scale = inv * inv * inv;
      y.xd[0] = (alpha * x.xd[0] - beta * x.x[0]) * scale;
      y.xd[1] = (alpha * x.xd[1] - beta * x.x[1]) * scale;
      y.xd[2] = (alpha * x.xd[2] - beta * x.x[2]) * scale;
    }
    void reverse(){
      T beta = Vec3DotCore(x.x, y.xd);
      T scale = inv * inv * inv;
      x.xd[0] += (alpha * y.xd[0] - beta * x.x[0]) * scale;
      x.xd[1] += (alpha * y.xd[1] - beta * x.x[1]) * scale;
      x.xd[2] += (alpha * y.xd[2] - beta * x.x[2]) * scale;
    }

    T alpha, inv;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class Mat3x2ToVec3 {
  public:
    Mat3x2ToVec3( const Mat3x2<T>& A, Vec3<T>& x, Vec3<T>& y ){
      x.x[0] = A.A[0];
      x.x[1] = A.A[2];
      x.x[2] = A.A[4];

      y.x[0] = A.A[1];
      y.x[1] = A.A[3];
      y.x[2] = A.A[5];
    }
  };

  template<typename T>
  class ADMat3x2ToADVec3 {
  public:
    ADMat3x2ToADVec3( ADMat3x2<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      x.x[0] = A.A[0];
      x.x[1] = A.A[2];
      x.x[2] = A.A[4];

      y.x[0] = A.A[1];
      y.x[1] = A.A[3];
      y.x[2] = A.A[5];
    }
    void forward(){
      x.xd[0] = A.Ad[0];
      x.xd[1] = A.Ad[2];
      x.xd[2] = A.Ad[4];

      y.xd[0] = A.Ad[1];
      y.xd[1] = A.Ad[3];
      y.xd[2] = A.Ad[5];
    }
    void reverse(){
      A.Ad[0] += x.xd[0];
      A.Ad[2] += x.xd[1];
      A.Ad[4] += x.xd[2];

      A.Ad[1] += y.xd[0];
      A.Ad[3] += y.xd[1];
      A.Ad[5] += y.xd[2];
    }

    ADMat3x2<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  /*
    Matrix-vector product
  */
  template<typename T>
  class Mat3x3VecMult {
  public:
    Mat3x3VecMult( const Mat3x3<T>& A, const Vec3<T>& x, Vec3<T>& y ){
      Mat3x3VecMultCore(A.A, x.x, y.x);
    }
  };

  template<typename T>
  class ADMat3x3VecMult {
  public:
    ADMat3x3VecMult( ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      Mat3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultCore(A.Ad, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddCore(y.xd, x.x, A.Ad);
    }

    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class Mat3x3ADVecMult {
  public:
    Mat3x3ADVecMult( const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      Mat3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultAddCore(A.A, x.xd, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddCore(A.A, y.xd, x.xd);
    }

    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMat3x3ADVecMult {
  public:
    ADMat3x3ADVecMult( ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      Mat3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultCore(A.Ad, x.x, y.xd);
      Mat3x3VecMultAddCore(A.A, x.xd, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddCore(A.A, y.xd, x.xd);
      Vec3OuterProductAddCore(y.xd, x.x, A.Ad);
    }

    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  /*
    Transpose matrix-vector product
  */
  template<typename T>
  class MatTrans3x3VecMult {
  public:
    MatTrans3x3VecMult( const Mat3x3<T>& A, const Vec3<T>& x, Vec3<T>& y ){
      MatTrans3x3VecMultCore(A.A, x.x, y.x);
    }
  };

  template<typename T>
  class ADMatTrans3x3VecMult {
  public:
    ADMatTrans3x3VecMult( ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      MatTrans3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultCore(A.Ad, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddCore(x.x, y.xd, A.Ad);
    }

    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class MatTrans3x3ADVecMult {
  public:
    MatTrans3x3ADVecMult( const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      MatTrans3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultAddCore(A.A, x.xd, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddCore(A.A, y.xd, x.xd);
    }

    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMatTrans3x3ADVecMult {
  public:
    ADMatTrans3x3ADVecMult( ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : A(A), x(x), y(y) {
      MatTrans3x3VecMultCore(A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultCore(A.Ad, x.x, y.xd);
      MatTrans3x3VecMultAddCore(A.A, x.xd, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddCore(A.A, y.xd, x.xd);
      Vec3OuterProductAddCore(x.x, y.xd, A.Ad);
    }

    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  /*
    Matrix-vector product
  */
  template<typename T>
  class Mat3x3VecMultScale {
  public:
    Mat3x3VecMultScale( const Scalar<T>& scale, const Mat3x3<T>& A, const Vec3<T>& x, Vec3<T>& y ){
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
  };

  template<typename T>
  class ADMat3x3VecMultScale {
  public:
    ADMat3x3VecMultScale( const Scalar<T>& scale, ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
    }

    const Scalar<T>& scale;
    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class Mat3x3ADVecMultScale {
  public:
    Mat3x3ADVecMultScale( const Scalar<T>& scale, const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    }

    const Scalar<T>& scale;
    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMat3x3ADVecMultScale {
  public:
    ADMat3x3ADVecMultScale( const Scalar<T>&scale, ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
    }

    const Scalar<T>& scale;
    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class Mat3x3VecMultADScale {
  public:
    Mat3x3VecMultADScale( ADScalar<T>& scale, const Mat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
    }

    ADScalar<T>& scale;
    const Mat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMat3x3VecMultADScale {
  public:
    ADMat3x3VecMultADScale( ADScalar<T>& scale, ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
      scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
    }

    ADScalar<T>& scale;
    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class Mat3x3ADVecMultADScale {
  public:
    Mat3x3ADVecMultADScale( ADScalar<T>& scale, const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.value, A.A, x.xd, y.xd);
      Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
    }

    ADScalar<T>& scale;
    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMat3x3ADVecMultADScale {
  public:
    ADMat3x3ADVecMultADScale( ADScalar<T>&scale, ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
      Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
      scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
    }

    ADScalar<T>& scale;
    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  /*
    Transpose matrix-vector product
  */
  template<typename T>
  class MatTrans3x3VecMultScale {
  public:
    MatTrans3x3VecMultScale( const Scalar<T>& scale, const Mat3x3<T>& A, const Vec3<T>& x, Vec3<T>& y ){
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
  };

  template<typename T>
  class ADMatTrans3x3VecMultScale {
  public:
    ADMatTrans3x3VecMultScale( const Scalar<T>& scale, ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
    }

    const Scalar<T>& scale;
    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class MatTrans3x3ADVecMultScale {
  public:
    MatTrans3x3ADVecMultScale( const Scalar<T>& scale, const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    }

    const Scalar<T>& scale;
    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMatTrans3x3ADVecMultScale {
  public:
    ADMatTrans3x3ADVecMultScale( const Scalar<T>& scale, ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
    }

    const Scalar<T>& scale;
    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class MatTrans3x3VecMultADScale {
  public:
    MatTrans3x3VecMultADScale( ADScalar<T>& scale, const Mat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }

    ADScalar<T>& scale;
    const Mat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMatTrans3x3VecMultADScale {
  public:
    ADMatTrans3x3VecMultADScale( ADScalar<T>& scale, ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
      scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }

    ADScalar<T>& scale;
    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class MatTrans3x3ADVecMultADScale {
  public:
    MatTrans3x3ADVecMultADScale( ADScalar<T>& scale, const Mat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.xd, y.xd);
      MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }

    ADScalar<T>& scale;
    const Mat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  template<typename T>
  class ADMatTrans3x3ADVecMultADScale {
  public:
    ADMatTrans3x3ADVecMultADScale( ADScalar<T>& scale, ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y ) : scale(scale), A(A), x(x), y(y) {
      MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
    }
    void forward(){
      MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
      MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
      MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
      Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
      scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }

    ADScalar<T>& scale;
    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
  };

  /*
    Inner products alpha = x^{T} * A * y
  */
  template<typename T>
  class Mat3x3VecVecInnerProduct {
  public:
    Mat3x3VecVecInnerProduct( const Mat3x3<T>& A, const Vec3<T>& x, const Vec3<T>& y, Scalar<T>& alpha ){
      alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
    }
  };

  template<typename T>
  class ADMat3x3VecVecInnerProduct {
  public:
    ADMat3x3VecVecInnerProduct( ADMat3x3<T>& A, const Vec3<T>& x, const Vec3<T>& y, ADScalar<T>& alpha ) : A(A), x(x), y(y), alpha(alpha) {
      alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
    }
    void forward(){
      alpha.valued = Mat3x3InnerProductCore(A.Ad, x.x, y.x);
    }
    void reverse(){
      Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
    }

    ADMat3x3<T>& A;
    const Vec3<T>& x;
    const Vec3<T>& y;
    ADScalar<T>& alpha;
  };

  template<typename T>
  class ADMat3x3VecADVecInnerProduct {
  public:
    ADMat3x3VecADVecInnerProduct( ADMat3x3<T>& A, const Vec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : A(A), x(x), y(y), alpha(alpha) {
      alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
    }
    void forward(){
      alpha.valued =
        Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
        Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }
    void reverse(){
      MatTrans3x3VecMultAddScaleCore(alpha.valued, A.A, x.x, y.xd);
      Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
    }

    ADMat3x3<T>& A;
    const Vec3<T>& x;
    ADVec3<T>& y;
    ADScalar<T>& alpha;
  };

  template<typename T>
  class ADMat3x3ADVecVecInnerProduct {
  public:
    ADMat3x3ADVecVecInnerProduct( ADMat3x3<T>& A, ADVec3<T>& x, const Vec3<T>& y, ADScalar<T>& alpha ) : A(A), x(x), y(y), alpha(alpha) {
      alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
    }
    void forward(){
      alpha.valued =
        Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
        Mat3x3InnerProductCore(A.A, x.xd, y.x);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(alpha.valued, A.A, y.x, x.xd);
      Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
    }

    ADMat3x3<T>& A;
    ADVec3<T>& x;
    const Vec3<T>& y;
    ADScalar<T>& alpha;
  };

  template<typename T>
  class ADMat3x3ADVecADVecInnerProduct {
  public:
    ADMat3x3ADVecADVecInnerProduct( ADMat3x3<T>& A, ADVec3<T>& x, ADVec3<T>& y, ADScalar<T>& alpha ) : A(A), x(x), y(y), alpha(alpha) {
      alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
    }
    void forward(){
      alpha.valued =
        Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
        Mat3x3InnerProductCore(A.A, x.xd, y.x) +
        Mat3x3InnerProductCore(A.A, x.x, y.xd);
    }
    void reverse(){
      Mat3x3VecMultAddScaleCore(alpha.valued, A.A, y.x, x.xd);
      MatTrans3x3VecMultAddScaleCore(alpha.valued, A.A, x.x, y.xd);
      Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
    }

    ADMat3x3<T>& A;
    ADVec3<T>& x;
    ADVec3<T>& y;
    ADScalar<T>& alpha;
  };

} // namespace AD

#endif // A2D_VEC_OPS_H
