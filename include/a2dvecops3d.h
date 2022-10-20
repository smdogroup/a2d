#ifndef A2D_VEC_OPS_H
#define A2D_VEC_OPS_H

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2dveccore3d.h"

namespace A2D {

/*
  Vec3Norm
  This represents the L2 norm (or 2 norm) of a 3 vector
  normObj = ||xObj||
*/
template <typename T>
inline void Vec3Norm(const Vec<T, 3>& x,
                     T& norm) {
  norm = Vec3NormCore<T>(x);
}

template <typename T>
class ADVec3NormExpr {
 public:
  ADVec3NormExpr(ADVec<Vec<T, 3>>& xObj,
                 ADScalar<T>& normObj)
      : xObj(xObj), normObj(normObj) {
    const Vec<T, 3>& x = xObj.value();
    normObj.value = Vec3NormCore<T>(x);
  }

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
    const T& n = normObj.value;
// //            output:
//            T& nb = normObj.bvalue;
// //            operations:
//            nb = Vec3DotCore<T, Vec<T, 3>>(x, xb) / n;
    normObj.bvalue = Vec3DotCore<T, Vec<T, 3>>(x, xb) / n;
  }

  void reverse() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const T& n = normObj.value;
    const T& nb = normObj.bvalue;
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Vec3ScaleCore(nb / n, x, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADScalar<T>& normObj;
};

template <typename T>
inline ADVec3NormExpr<T> Vec3Norm(ADVec<Vec<T, 3>>& x,
                                  ADScalar<T>& norm) {
  return ADVec3NormExpr<T>(x, norm);
}

/*
  Vec3Scale
  This represents scalar multiplication of a 3 vector
  yObj = aObj * xObj
*/
template <typename T>
inline void Vec3Scale(const Vec<T, 3>& x,
                      const T& a,
                      Vec<T, 3>& v) {
  Vec3ScaleCore(a, x, v);
}

template <typename T>
class ADVec3ScaleExpr {
 public:
  ADVec3ScaleExpr(ADVec<Vec<T, 3>>& xObj,
                  const T& a,
                  ADVec<Vec<T, 3>>& vObj)
      : xObj(xObj), a(a), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  };

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
//            output:
    Vec<T, 3>& vb = vObj.bvalue();
//            operations:
    Vec3ScaleCore(a, xb, vb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Vec3ScaleCore(a, vb, xb);
  };

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& vObj;
  const T& a;
};

template <typename T>
inline ADVec3ScaleExpr<T> Vec3Scale(ADVec<Vec<T, 3>>& x,
                                    const T& a,
                                    ADVec<Vec<T, 3>>& v) {
  return ADVec3ScaleExpr<T>(x, a, v);
}

template <typename T>
class Vec3ADScaleExpr {
 public:
  Vec3ADScaleExpr(const Vec<T, 3>& x,
                  ADScalar<T>& aObj,
                  ADVec<Vec<T, 3>>& vObj)
      : x(x), aObj(aObj), vObj(vObj) {
    const T& a = aObj.value;
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  }

  void forward() {
//            input:
    const T& ab = aObj.bvalue;
//            output:
    Vec<T, 3>& vb = vObj.bvalue();
//            operations:
    Vec3ScaleCore(ab, x, vb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& vb = vObj.bvalue();
//            output:
    T& ab = aObj.bvalue;
//            operations:
    ab = Vec3DotCore<T, Vec<T, 3>>(vb, x);
  };

  const Vec<T, 3>& x;
  ADVec<Vec<T, 3>>& vObj;
  ADScalar<T>& aObj;
};

template <typename T>
inline Vec3ADScaleExpr<T> Vec3Scale(const Vec<T, 3>& x,
                                    ADScalar<T>& a,
                                    ADVec<Vec<T, 3>>& v) {
  return Vec3ADScaleExpr<T>(x, a, v);
}

template <typename T>
class ADVec3ADScaleExpr {
 public:
  ADVec3ADScaleExpr(ADVec<Vec<T, 3>>& xObj,
                    ADScalar<T>& aObj,
                    ADVec<Vec<T, 3>>& vObj)
      : xObj(xObj), aObj(aObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  }

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
    const T& a = aObj.value;
    const T& ab = aObj.bvalue;
//            output:
    Vec<T, 3>& vb = vObj.bvalue();
//            operations:
    Vec3AXPBYCore(ab, x, a, xb, vb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    T& ab = aObj.bvalue;
//            operations:
    ab = Vec3DotCore<T, Vec<T, 3>>(vb, x);
    Vec3ScaleCore(a, vb, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& vObj;
  ADScalar<T>& aObj;
};

template <typename T>
inline ADVec3ADScaleExpr<T> Vec3Scale(ADVec<Vec<T, 3>>& x,
                                      ADScalar<T>& a,
                                      ADVec<Vec<T, 3>>& v) {
  return ADVec3ADScaleExpr<T>(x, a, v);
}

/*
  Vec3Axpy
*/
template <typename T>
class Vec3VecADScalarAxpyExpr {
 public:
  Vec3VecADScalarAxpyExpr(ADScalar<T>& aobj, const Vec<T, 3>& x,
                          const Vec<T, 3>& y, ADVec<Vec<T, 3>>& vobj)
      : scale(1.0), aobj(aobj), x(x), y(y), vobj(vobj) {
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(aobj.value, x, y, v);
  }
  Vec3VecADScalarAxpyExpr(const T& scale, ADScalar<T>& aobj, const Vec<T, 3>& x,
                          const Vec<T, 3>& y, ADVec<Vec<T, 3>>& vobj)
      : scale(scale), aobj(aobj), x(x), y(y), vobj(vobj) {
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(scale * aobj.value, x, y, v);
  }
  void forward() {
    Vec<T, 3>& vb = vobj.bvalue();
    vb(0) = scale * aobj.bvalue * x(0);
    vb(1) = scale * aobj.bvalue * x(1);
    vb(2) = scale * aobj.bvalue * x(2);
  }
  void reverse() {
    const Vec<T, 3>& vb = vobj.bvalue();
    aobj.bvalue += scale * Vec3DotCore<T>(x, vb);
  }

  const T scale;
  ADScalar<T>& aobj;
  const Vec<T, 3>& x;
  const Vec<T, 3>& y;
  ADVec<Vec<T, 3>>& vobj;
};

template <typename T>
class ADVec3VecADScalarAxpyExpr {
 public:
  ADVec3VecADScalarAxpyExpr(ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj,
                            const Vec<T, 3>& y, ADVec<Vec<T, 3>>& vobj)
      : scale(1.0), aobj(aobj), xobj(xobj), y(y), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(aobj.value, x, y, v);
  }
  ADVec3VecADScalarAxpyExpr(const T& scale, ADScalar<T>& aobj,
                            ADVec<Vec<T, 3>>& xobj, const Vec<T, 3>& y,
                            ADVec<Vec<T, 3>>& vobj)
      : scale(scale), aobj(aobj), xobj(xobj), y(y), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(scale * aobj.value, x, y, v);
  }
  void forward() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& xb = xobj.bvalue();
    Vec<T, 3>& vb = vobj.bvalue();

    vb(0) = scale * (aobj.bvalue * x(0) + aobj.value * xb(0));
    vb(1) = scale * (aobj.bvalue * x(1) + aobj.value * xb(1));
    vb(2) = scale * (aobj.bvalue * x(2) + aobj.value * xb(2));
  }
  void reverse() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& vb = vobj.bvalue();
    Vec<T, 3>& xb = xobj.bvalue();

    aobj.bvalue += scale * Vec3DotCore<T>(x, vb);
    xb(0) += scale * aobj.value * vb(0);
    xb(1) += scale * aobj.value * vb(1);
    xb(2) += scale * aobj.value * vb(2);
  }

  const T scale;
  ADScalar<T>& aobj;
  ADVec<Vec<T, 3>>& xobj;
  const Vec<T, 3>& y;
  ADVec<Vec<T, 3>>& vobj;
};

template <typename T>
class ADVec3ADVecScalarAxpyExpr {
 public:
  ADVec3ADVecScalarAxpyExpr(const T& alpha, ADVec<Vec<T, 3>>& xobj,
                            ADVec<Vec<T, 3>>& yobj, ADVec<Vec<T, 3>>& vobj)
      : scale(1.0), alpha(alpha), xobj(xobj), yobj(yobj), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(alpha, x, y, v);
  }
  ADVec3ADVecScalarAxpyExpr(const T& scale, const T& alpha,
                            ADVec<Vec<T, 3>>& xobj, ADVec<Vec<T, 3>>& yobj,
                            ADVec<Vec<T, 3>>& vobj)
      : scale(scale), alpha(alpha), xobj(xobj), yobj(yobj), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(scale * alpha, x, y, v);
  }
  void forward() {
    const Vec<T, 3>& xb = xobj.bvalue();
    const Vec<T, 3>& yb = yobj.bvalue();
    Vec<T, 3>& vb = vobj.bvalue();

    vb(0) = scale * (alpha * xb(0)) + yb(0);
    vb(1) = scale * (alpha * xb(1)) + yb(1);
    vb(2) = scale * (alpha * xb(2)) + yb(2);
  }
  void reverse() {
    const Vec<T, 3>& vb = vobj.bvalue();
    Vec<T, 3>& xb = xobj.bvalue();
    Vec<T, 3>& yb = yobj.bvalue();

    xb(0) += scale * alpha * vb(0);
    xb(1) += scale * alpha * vb(1);
    xb(2) += scale * alpha * vb(2);

    yb(0) += vb(0);
    yb(1) += vb(1);
    yb(2) += vb(2);
  }

  const T scale;
  const T& alpha;
  ADVec<Vec<T, 3>>& xobj;
  ADVec<Vec<T, 3>>& yobj;
  ADVec<Vec<T, 3>>& vobj;
};

template <typename T>
class ADVec3AxpyExpr {
 public:
  ADVec3AxpyExpr(ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj,
                 ADVec<Vec<T, 3>>& yobj, ADVec<Vec<T, 3>>& vobj)
      : scale(1.0), aobj(aobj), xobj(xobj), yobj(yobj), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(aobj.value, x, y, v);
  }
  ADVec3AxpyExpr(const T& scale, ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj,
                 ADVec<Vec<T, 3>>& yobj, ADVec<Vec<T, 3>>& vobj)
      : scale(scale), aobj(aobj), xobj(xobj), yobj(yobj), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3AXPYCore(scale * aobj.value, x, y, v);
  }
  void forward() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    const Vec<T, 3>& xb = xobj.bvalue();
    const Vec<T, 3>& yb = yobj.bvalue();
    Vec<T, 3>& vb = vobj.bvalue();

    vb(0) = scale * (aobj.bvalue * x(0) + aobj.value * xb(0)) + yb(0);
    vb(1) = scale * (aobj.bvalue * x(1) + aobj.value * xb(1)) + yb(1);
    vb(2) = scale * (aobj.bvalue * x(2) + aobj.value * xb(2)) + yb(2);
  }
  void reverse() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& vb = vobj.bvalue();

    Vec<T, 3>& xb = xobj.bvalue();
    Vec<T, 3>& yb = yobj.bvalue();

    aobj.bvalue += scale * Vec3DotCore<T>(x, vb);

    xb(0) += scale * aobj.value * vb(0);
    xb(1) += scale * aobj.value * vb(1);
    xb(2) += scale * aobj.value * vb(2);

    yb(0) += vb(0);
    yb(1) += vb(1);
    yb(2) += vb(2);
  }

  const T scale;
  ADScalar<T>& aobj;
  ADVec<Vec<T, 3>>& xobj;
  ADVec<Vec<T, 3>>& yobj;
  ADVec<Vec<T, 3>>& vobj;
};

template <typename T>
inline void Vec3Axpy(const T& alpha, const Vec<T, 3>& x, const Vec<T, 3>& y,
                     Vec<T, 3>& v) {
  Vec3AXPYCore(alpha, x, y, v);
}

template <typename T>
inline void Vec3Axpy(const T& scale, const T& alpha, const Vec<T, 3>& x,
                     const Vec<T, 3>& y, Vec<T, 3>& v) {
  Vec3AXPYCore(scale * alpha, x, y, v);
}

template <typename T>
inline Vec3VecADScalarAxpyExpr<T> Vec3Axpy(ADScalar<T>& aobj,
                                           const Vec<T, 3>& x,
                                           const Vec<T, 3>& y,
                                           ADVec<Vec<T, 3>>& vobj) {
  return Vec3VecADScalarAxpyExpr<T>(aobj, x, y, vobj);
}

template <typename T>
inline Vec3VecADScalarAxpyExpr<T> Vec3Axpy(const T& scale, ADScalar<T>& aobj,
                                           const Vec<T, 3>& x,
                                           const Vec<T, 3>& y,
                                           ADVec<Vec<T, 3>>& vobj) {
  return Vec3VecADScalarAxpyExpr<T>(scale, aobj, x, y, vobj);
}

template <typename T>
inline ADVec3VecADScalarAxpyExpr<T> Vec3Axpy(ADScalar<T>& aobj,
                                             ADVec<Vec<T, 3>>& xobj,
                                             const Vec<T, 3>& y,
                                             ADVec<Vec<T, 3>>& vobj) {
  return ADVec3VecADScalarAxpyExpr<T>(aobj, xobj, y, vobj);
}

template <typename T>
inline ADVec3VecADScalarAxpyExpr<T> Vec3Axpy(const T& scale, ADScalar<T>& aobj,
                                             ADVec<Vec<T, 3>>& xobj,
                                             const Vec<T, 3>& y,
                                             ADVec<Vec<T, 3>>& vobj) {
  return ADVec3VecADScalarAxpyExpr<T>(scale, aobj, xobj, y, vobj);
}

template <typename T>
inline ADVec3ADVecScalarAxpyExpr<T> Vec3Axpy(const T& alpha,
                                             ADVec<Vec<T, 3>>& xobj,
                                             ADVec<Vec<T, 3>>& yobj,
                                             ADVec<Vec<T, 3>>& vobj) {
  return ADVec3ADVecScalarAxpyExpr<T>(alpha, xobj, yobj, vobj);
}

template <typename T>
inline ADVec3ADVecScalarAxpyExpr<T> Vec3Axpy(const T& scale, const T& alpha,
                                             ADVec<Vec<T, 3>>& xobj,
                                             ADVec<Vec<T, 3>>& yobj,
                                             ADVec<Vec<T, 3>>& vobj) {
  return ADVec3ADVecScalarAxpyExpr<T>(scale, alpha, xobj, yobj, vobj);
}

template <typename T>
inline ADVec3AxpyExpr<T> Vec3Axpy(ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj,
                                  ADVec<Vec<T, 3>>& yobj,
                                  ADVec<Vec<T, 3>>& vobj) {
  return ADVec3AxpyExpr<T>(aobj, xobj, yobj, vobj);
}

template <typename T>
inline ADVec3AxpyExpr<T> Vec3Axpy(const T& scale, ADScalar<T>& aobj,
                                  ADVec<Vec<T, 3>>& xobj,
                                  ADVec<Vec<T, 3>>& yobj,
                                  ADVec<Vec<T, 3>>& vobj) {
  return ADVec3AxpyExpr<T>(scale, aobj, xobj, yobj, vobj);
}

/*
  Vec3Dot
  This represents vector-vector dot product (inner product) of two 3 vectors
  alpha = xObj . yObj
  A.K.A.
  alpha = xObj @ yObj (in numpy notation)
*/
template <typename T>
inline void Vec3Dot(const Vec<T, 3>& x,
                    const Vec<T, 3>& y,
                    T& a) {
  a = Vec3DotCore<T>(x, y);
}

template <typename T>
class ADVec3DotVecExpr {
 public:
  ADVec3DotVecExpr(ADVec<Vec<T, 3>>& xObj,
                   const Vec<T, 3>& y,
                   ADScalar<T>& aObj)
      : xObj(xObj), y(y), aObj(aObj) {
    const Vec<T, 3>& x = xObj.value();

    aObj.value = Vec3DotCore<T, Vec<T, 3>>(x, y);
  }

  void forward() {
//            input:
    const Vec<T, 3>& xb = xObj.bvalue();
//            output:
    T& ab = aObj.bvalue;
//            operations:
    ab = Vec3DotCore<T, Vec<T, 3>>(xb, y);
  }

  void reverse() {
//            input:
    const T& ab = aObj.bvalue;
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Vec3ScaleCore(ab, y, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  ADScalar<T>& aObj;
};

template <typename T>
inline ADVec3DotVecExpr<T> Vec3Dot(ADVec<Vec<T, 3>>& x,
                                   const Vec<T, 3>& y,
                                   ADScalar<T>& a) {
  return ADVec3DotVecExpr<T>(x, y, a);
}

template <typename T>
class ADVec3DotADVecExpr {
 public:
  ADVec3DotADVecExpr(ADVec<Vec<T, 3>>& xObj,
                     ADVec<Vec<T, 3>>& yObj,
                     ADScalar<T>& aObj)
      : xObj(xObj), yObj(yObj), aObj(aObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();

    aObj.value = Vec3DotCore<T, Vec<T, 3>>(x, y);
  }

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& yb = yObj.bvalue();
//            output:
    T& ab = aObj.bvalue;
//            operations:
    ab = Vec3DotCore<T, Vec<T, 3>>(xb, y) + Vec3DotCore<T, Vec<T, 3>>(x, yb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const T& ab = aObj.bvalue;
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
//            operations:
    Vec3ScaleCore(ab, y, xb);
    Vec3ScaleCore(ab, x, yb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& yObj;
  ADScalar<T>& aObj;
};

template <typename T>
inline ADVec3DotADVecExpr<T> Vec3Dot(ADVec<Vec<T, 3>>& x,
                                     ADVec<Vec<T, 3>>& y,
                                     ADScalar<T>& a) {
  return ADVec3DotADVecExpr<T>(x, y, a);
}

/*
  Vec3CrossProduct
*/
template <typename T>
class ADVec3CrossProductExpr {
 public:
  ADVec3CrossProductExpr(ADVec<Vec<T, 3>>& xobj, ADVec<Vec<T, 3>>& yobj,
                         ADVec<Vec<T, 3>>& vobj)
      : xobj(xobj), yobj(yobj), vobj(vobj) {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    Vec<T, 3>& v = vobj.value();
    Vec3CrossProductCore(x, y, v);
  }
  void forward() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& xb = xobj.bvalue();
    const Vec<T, 3>& y = yobj.value();
    const Vec<T, 3>& yb = yobj.bvalue();

    Vec<T, 3>& vb = vobj.bvalue();

    Vec3CrossProductCore(xb, y, vb);
    Vec3CrossProductAddCore(x, yb, vb);
  }
  void reverse() {
    const Vec<T, 3>& x = xobj.value();
    const Vec<T, 3>& y = yobj.value();
    const Vec<T, 3>& vb = vobj.bvalue();

    Vec<T, 3>& xb = xobj.bvalue();
    Vec<T, 3>& yb = yobj.bvalue();

    Vec3CrossProductAddCore(y, vb, xb);
    Vec3CrossProductAddCore(vb, x, yb);
  }

  ADVec<Vec<T, 3>>& xobj;
  ADVec<Vec<T, 3>>& yobj;
  ADVec<Vec<T, 3>>& vobj;
};

template <typename T>
inline void Vec3Cross(const Vec<T, 3>& x, const Vec<T, 3>& y, Vec<T, 3>& v) {
  Vec3CrossProductCore(x, y, v);
}

template <typename T>
inline ADVec3CrossProductExpr<T> Vec3Cross(ADVec<Vec<T, 3>>& xobj,
                                           ADVec<Vec<T, 3>>& yobj,
                                           ADVec<Vec<T, 3>>& vobj) {
  return ADVec3CrossProductExpr<T>(xobj, yobj, vobj);
}

/*
  Vec3Normalize
  This represents the normalization of a vector
  vObj = xObj / ||xObj||
*/
template <typename T>
inline void Vec3Normalize(const Vec<T, 3>& x,
                          Vec<T, 3>& v) {
  Vec3ScaleCore(1 / Vec3NormCore<T>(x), x, v);
}

template <typename T>
class ADVec3NormalizeExpr {
 public:
  ADVec3NormalizeExpr(ADVec<Vec<T, 3>>& xObj,
                      ADVec<Vec<T, 3>>& vObj)
      : xObj(xObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    negativeSquaredNormInv = 1 / Vec3DotCore<T>(x, x);
    normInv = sqrt(negativeSquaredNormInv);
    negativeSquaredNormInv *= -1;
    Vec3ScaleCore(normInv, x, v);
  }

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
    const Vec<T, 3>& v = vObj.value();
//            output:
    Vec<T, 3>& vb = vObj.bvalue();
//            operations:
    Vec3AXPBYCore(normInv, xb, negativeSquaredNormInv * Vec3DotCore<T>(xb, x), v, vb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& v = vObj.value();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Vec3AXPBYCore(normInv, vb, negativeSquaredNormInv * Vec3DotCore<T>(vb, x), v, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& vObj;

 private:
  T normInv;
  T negativeSquaredNormInv;
};

template <typename T>
inline ADVec3NormalizeExpr<T> Vec3Normalize(ADVec<Vec<T, 3>>& x,
                                            ADVec<Vec<T, 3>>& v) {
  return ADVec3NormalizeExpr<T>(x, v);
}

/*
  Vec3ScaleSymmetricOuterProduct
  This represents a scaled, 3 vector outer product with itself
  sObj = aObj * xObj . (xObj^T)
*/
template <typename T>
inline void Vec3ScaleSymmetricOuterProduct(const T& a,
                                           const Vec<T, 3>& x,
                                           Mat<T, 3, 3>& S) {
  Vec3OuterProductScaleCore(a, x, x, S);
}

template <typename T>
class ADVec3ScaleSymmetricOuterProductExpr {
 public:
  ADVec3ScaleSymmetricOuterProductExpr(const T& a,
                                       ADVec<Vec<T, 3>>& xObj,
                                       ADMat<Mat<T, 3, 3>>& sObj)
      : a(a), xObj(xObj), sObj(sObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
//            output:
    Vec<T, 3> ax;
    Mat<T, 3, 3>& Sb = sObj.bvalue();
//            operations:
    Vec3ScaleCore<T, Vec<T, 3>>(a, x, ax);
    Vec3OuterProductCore(xb, ax, Sb);
    Vec3OuterProductAddCore(ax, xb, Sb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
//            operations:
    Mat3x3VecMultScaleCore(2 * a, Sb, x, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADMat<Mat<T, 3, 3>>& sObj;
  const T& a;
};

template <typename T>
inline ADVec3ScaleSymmetricOuterProductExpr<T> Vec3ScaleSymmetricOuterProduct(const T& a,
                                                                              ADVec<Vec<T, 3>>& x,
                                                                              ADMat<Mat<T, 3, 3>>& S) {
  return ADVec3ScaleSymmetricOuterProductExpr<T>(a, x, S);
}

template <typename T>
class Vec3ADScaleSymmetricOuterProductExpr {
 public:
  Vec3ADScaleSymmetricOuterProductExpr(ADScalar<T>& aObj,
                                       const Vec<T, 3>& x,
                                       ADMat<Mat<T, 3, 3>>& sObj)
      : aObj(aObj), x(x), sObj(sObj) {
    const T& a = aObj.value;
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void forward() {
//            input:
    const T& ab = aObj.bvalue;
//            output:
    Mat<T, 3, 3>& Sb = sObj.bvalue();
//            operations:
    Vec3OuterProductAddScaleCore(ab, x, x, Sb);
  }

  void reverse() {
//            input:
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
    const T& a = aObj.value;
//            output:
    T& ab = aObj.bvalue;
//            operations:
    ab = Mat3x3InnerProductCore<T>(Sb, x, x);
  }

  const Vec<T, 3>& x;
  ADMat<Mat<T, 3, 3>>& sObj;
  ADScalar<T>& aObj;
};

template <typename T>
inline Vec3ADScaleSymmetricOuterProductExpr<T> Vec3ScaleSymmetricOuterProduct(ADScalar<T>& a,
                                                                              const Vec<T, 3>& x,
                                                                              ADMat<Mat<T, 3, 3>>& S) {
  return Vec3ADScaleSymmetricOuterProductExpr<T>(a, x, S);
}

template <typename T>
class ADVec3ADScaleSymmetricOuterProductExpr {
 public:
  ADVec3ADScaleSymmetricOuterProductExpr(ADScalar<T>& aObj,
                                         ADVec<Vec<T, 3>>& xObj,
                                         ADMat<Mat<T, 3, 3>>& sObj)
      : aObj(aObj), xObj(xObj), sObj(sObj) {
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void forward() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& xb = xObj.bvalue();
    const T& a = aObj.value;
    const T& ab = aObj.bvalue;
//            output:
    Vec<T, 3> ax;
    Mat<T, 3, 3>& Sb = sObj.bvalue();
//            operations:
    Vec3ScaleCore<T, Vec<T, 3>>(a, x, ax);
    Vec3OuterProductCore<Vec<T, 3>, Mat<T, 3, 3>>(xb, ax, Sb);
    Vec3OuterProductAddCore<Vec<T, 3>, Mat<T, 3, 3>>(ax, xb, Sb);
    Vec3OuterProductAddScaleCore<T, Vec<T, 3>, Mat<T, 3, 3>>(ab, x, x, Sb);
  }

  void reverse() {
//            input:
    const Vec<T, 3>& x = xObj.value();
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
    const T& a = aObj.value;
//            output:
    Vec<T, 3>& xb = xObj.bvalue();
    T& ab = aObj.bvalue;
//            operations:
    ab = Mat3x3InnerProductCore<T, Vec<T, 3>, Mat<T, 3, 3>>(Sb, x, x);
    Mat3x3VecMultScaleCore(2 * a, Sb, x, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADMat<Mat<T, 3, 3>>& sObj;
  ADScalar<T>& aObj;
};

template <typename T>
inline ADVec3ADScaleSymmetricOuterProductExpr<T> Vec3ScaleSymmetricOuterProduct(ADScalar<T>& a,
                                                                                ADVec<Vec<T, 3>>& x,
                                                                                ADMat<Mat<T, 3, 3>>& S) {
  return ADVec3ADScaleSymmetricOuterProductExpr<T>(a, x, S);
}

}  // namespace A2D

#endif  // A2D_VEC_OPS_H
