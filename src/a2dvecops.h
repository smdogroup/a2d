#ifndef A2D_VEC_OPS_H
#define A2D_VEC_OPS_H

#include "a2dobjs.h"
#include "a2dveccore.h"

/*
  Vec3Norm
*/
class Vec3Norm {
public:
  Vec3Norm( Vec3& x, Scalar& alpha ){
    alpha.value = Vec3NormCore(x.x);
  }
};

class ADVec3Norm {
public:
  ADVec3Norm( ADVec3& x, ADScalar& alpha ) : x(x), alpha(alpha) {
    alpha.value = Vec3NormCore(x.x);
  }
  void forward(){
    if (alpha.value != 0.0){
      alpha.valued = Vec3DotCore(x.x, x.xd)/alpha.value;
    }
  }
  void reverse(){
    if (alpha.value != 0.0){
      TacsScalar ainv = alpha.valued/alpha.value;
      x.xd[0] += ainv * x.x[0];
      x.xd[1] += ainv * x.x[1];
      x.xd[2] += ainv * x.x[2];
    }
  }

  ADVec3& x;
  ADScalar& alpha;
};

/*
  Vec3Scale
*/
class Vec3Scale {
public:
  Vec3Scale( Scalar& alpha, Vec3& x, Vec3& v ){
    Vec3ScaleCore(alpha.value, x.x, v.x);
  }
};

class ADVec3Scale {
public:
  ADVec3Scale( ADScalar& alpha, ADVec3& x, ADVec3& v ) : alpha(alpha), x(x), v(v) {
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

  ADScalar& alpha;
  ADVec3& x;
  ADVec3& v;
};

/*
  Vec3Axpy
*/
class Vec3Axpy {
public:
  Vec3Axpy( Scalar& alpha, Vec3& x, Vec3& y, Vec3& v ){
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
};

class ADVec3Axpy {
public:
  ADVec3Axpy( ADScalar& alpha, ADVec3& x, ADVec3& y, ADVec3& v ) : alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  void forward(){
    v.xd[0] = alpha.valued * x.x[0] + alpha.value * x.xd[0] + y.xd[0];
    v.xd[1] = alpha.valued * x.x[1] + alpha.value * x.xd[1] + y.xd[1];
    v.xd[2] = alpha.valued * x.x[2] + alpha.value * x.xd[2] + y.xd[2];
  }
  void reverse(){
    alpha.valued = Vec3DotCore(x.x, v.xd);
    x.xd[0] += alpha.value * v.xd[0];
    x.xd[1] += alpha.value * v.xd[1];
    x.xd[2] += alpha.value * v.xd[2];
    y.xd[0] += v.xd[0];
    y.xd[1] += v.xd[1];
    y.xd[2] += v.xd[2];
  }

  ADScalar& alpha;
  ADVec3& x;
  ADVec3& y;
  ADVec3& v;
};

/*
  Vec3Dot
*/
class Vec3Dot {
public:
  Vec3Dot( Vec3& x, Vec3& y, Scalar& alpha ){
    alpha.value = Vec3DotCore(x.x, y.x);
  }
};

class ADVec3Dot {
public:
  ADVec3Dot( ADVec3& x, ADVec3& y, ADScalar& alpha ) : x(x), y(y), alpha(alpha) {
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  void forward(){
    alpha.valued = Vec3DotCore(x.x, y.xd) + Vec3DotCore(x.xd, y.x);
  }
  void reverse(){
    x.xd[0] += alpha.valued * y.x[0];
    x.xd[1] += alpha.valued * y.x[1];
    x.xd[2] += alpha.valued * y.x[2];

    y.xd[0] += alpha.valued * x.x[0];
    y.xd[1] += alpha.valued * x.x[1];
    y.xd[2] += alpha.valued * x.x[2];
  }

  ADVec3& x;
  ADVec3& y;
  ADScalar& alpha;
};

/*
  Vec3CrossProduct
*/
class Vec3CrossProduct {
public:
  Vec3CrossProduct( Vec3& x, Vec3& y, Vec3& v ){
    Vec3CrossProductCore(x.x, y.x, v.x);
  }
};

class ADVec3CrossProduct {
public:
  ADVec3CrossProduct( ADVec3& x, ADVec3& y, ADVec3& v ) : x(x), y(y), v(v) {
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

  ADVec3& x;
  ADVec3& y;
  ADVec3& v;
};

// /*
//   Vec3SymmOuterProduct
// */
// class Vec3SymmOuterProduct {
// public:
//   Vec3SymmOuterProduct( Scalar& alpha, Vec3& x, Symm3x3& S ){
//     Vec3SymmOuterProductCore(alpha.value, x.x, S.A);
//   }
// };

// class ADVec3SymmOuterProduct {
// public:
//   ADVec3SymmOuterProduct( ADScalar& alpha, ADVec3& x, ADSymm3x3& S ){
//     Vec3SymmOuterProductCore(alpha.value, x.x, S.A);
//   }
//   void computeDeriv(){

//   }

//   ADScalar& alpha;
//   ADVec3& x;
//   ADSymm3x3& S;
// };

// /*
//   Vec3OuterProduct
// */
// class Vec3OuterProduct {
// public:
//   Vec3OuterProduct( ADScalar& alpha, Vec3& x, Vec3& y, Mat3x3& A ){
//     Vec3OuterPproductCore(alpha.value, x.x, y.x, A.A);
//   }
// };

// class ADVec3OuterProduct {
// public:
//   ADVec3OuterProduct( ADScalar& alpha, Vec3& x, Vec3& y, Mat3x3& A ) : alpha(alpha), x(x), y(y), A(A) {
//     Vec3OuterPproductCore(alpha.value, x.x, y.x, A.A);
//   }
//   void computeDeriv(){

//   }
//   ADScalar& alpha;
//   Vec3& x;
//   Vec3& y;
//   Mat3x3& A;
// };

#endif // A2D_VEC_OPS_H