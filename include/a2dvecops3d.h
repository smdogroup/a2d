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
A2D_INLINE_FUNCTION void Vec3Norm(const Vec<T, 3>& x, T& norm) {
  norm = Vec3NormCore<T>(x);
}

template <typename T>
class ADVec3NormExpr {
 public:
  ADVec3NormExpr(ADVec<Vec<T, 3>>& xObj, ADScalar<T>& normObj)
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
A2D_INLINE_FUNCTION ADVec3NormExpr<T> Vec3Norm(ADVec<Vec<T, 3>>& x,
                                               ADScalar<T>& norm) {
  return ADVec3NormExpr<T>(x, norm);
}

/*
  Vec3Scale
  This represents scalar multiplication of a 3 vector
  yObj = aObj * xObj
*/
template <typename T>
A2D_INLINE_FUNCTION void Vec3Scale(const Vec<T, 3>& x, const T& a,
                                   Vec<T, 3>& v) {
  Vec3ScaleCore(a, x, v);
}

template <typename T>
class ADVec3ScaleExpr {
 public:
  ADVec3ScaleExpr(ADVec<Vec<T, 3>>& xObj, const T& a, ADVec<Vec<T, 3>>& vObj)
      : xObj(xObj), vObj(vObj), a(a) {
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
A2D_INLINE_FUNCTION ADVec3ScaleExpr<T> Vec3Scale(ADVec<Vec<T, 3>>& x,
                                                 const T& a,
                                                 ADVec<Vec<T, 3>>& v) {
  return ADVec3ScaleExpr<T>(x, a, v);
}

template <typename T>
class Vec3ADScaleExpr {
 public:
  Vec3ADScaleExpr(const Vec<T, 3>& x, ADScalar<T>& aObj, ADVec<Vec<T, 3>>& vObj)
      : x(x), vObj(vObj), aObj(aObj) {
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
A2D_INLINE_FUNCTION Vec3ADScaleExpr<T> Vec3Scale(const Vec<T, 3>& x,
                                                 ADScalar<T>& a,
                                                 ADVec<Vec<T, 3>>& v) {
  return Vec3ADScaleExpr<T>(x, a, v);
}

template <typename T>
class ADVec3ADScaleExpr {
 public:
  ADVec3ADScaleExpr(ADVec<Vec<T, 3>>& xObj, ADScalar<T>& aObj,
                    ADVec<Vec<T, 3>>& vObj)
      : xObj(xObj), vObj(vObj), aObj(aObj) {
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
A2D_INLINE_FUNCTION ADVec3ADScaleExpr<T> Vec3Scale(ADVec<Vec<T, 3>>& x,
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
A2D_INLINE_FUNCTION void Vec3Axpy(const T& alpha, const Vec<T, 3>& x,
                                  const Vec<T, 3>& y, Vec<T, 3>& v) {
  Vec3AXPYCore(alpha, x, y, v);
}

template <typename T>
A2D_INLINE_FUNCTION void Vec3Axpy(const T& scale, const T& alpha,
                                  const Vec<T, 3>& x, const Vec<T, 3>& y,
                                  Vec<T, 3>& v) {
  Vec3AXPYCore(scale * alpha, x, y, v);
}

template <typename T>
A2D_INLINE_FUNCTION Vec3VecADScalarAxpyExpr<T> Vec3Axpy(
    ADScalar<T>& aobj, const Vec<T, 3>& x, const Vec<T, 3>& y,
    ADVec<Vec<T, 3>>& vobj) {
  return Vec3VecADScalarAxpyExpr<T>(aobj, x, y, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION Vec3VecADScalarAxpyExpr<T> Vec3Axpy(
    const T& scale, ADScalar<T>& aobj, const Vec<T, 3>& x, const Vec<T, 3>& y,
    ADVec<Vec<T, 3>>& vobj) {
  return Vec3VecADScalarAxpyExpr<T>(scale, aobj, x, y, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3VecADScalarAxpyExpr<T> Vec3Axpy(
    ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj, const Vec<T, 3>& y,
    ADVec<Vec<T, 3>>& vobj) {
  return ADVec3VecADScalarAxpyExpr<T>(aobj, xobj, y, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3VecADScalarAxpyExpr<T> Vec3Axpy(
    const T& scale, ADScalar<T>& aobj, ADVec<Vec<T, 3>>& xobj,
    const Vec<T, 3>& y, ADVec<Vec<T, 3>>& vobj) {
  return ADVec3VecADScalarAxpyExpr<T>(scale, aobj, xobj, y, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3ADVecScalarAxpyExpr<T> Vec3Axpy(
    const T& alpha, ADVec<Vec<T, 3>>& xobj, ADVec<Vec<T, 3>>& yobj,
    ADVec<Vec<T, 3>>& vobj) {
  return ADVec3ADVecScalarAxpyExpr<T>(alpha, xobj, yobj, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3ADVecScalarAxpyExpr<T> Vec3Axpy(
    const T& scale, const T& alpha, ADVec<Vec<T, 3>>& xobj,
    ADVec<Vec<T, 3>>& yobj, ADVec<Vec<T, 3>>& vobj) {
  return ADVec3ADVecScalarAxpyExpr<T>(scale, alpha, xobj, yobj, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3AxpyExpr<T> Vec3Axpy(ADScalar<T>& aobj,
                                               ADVec<Vec<T, 3>>& xobj,
                                               ADVec<Vec<T, 3>>& yobj,
                                               ADVec<Vec<T, 3>>& vobj) {
  return ADVec3AxpyExpr<T>(aobj, xobj, yobj, vobj);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3AxpyExpr<T> Vec3Axpy(const T& scale,
                                               ADScalar<T>& aobj,
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
A2D_INLINE_FUNCTION void Vec3Dot(const Vec<T, 3>& x, const Vec<T, 3>& y, T& a) {
  a = Vec3DotCore<T>(x, y);
}

template <typename T>
class ADVec3DotVecExpr {
 public:
  ADVec3DotVecExpr(ADVec<Vec<T, 3>>& xObj, const Vec<T, 3>& y,
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
A2D_INLINE_FUNCTION ADVec3DotVecExpr<T> Vec3Dot(ADVec<Vec<T, 3>>& x,
                                                const Vec<T, 3>& y,
                                                ADScalar<T>& a) {
  return ADVec3DotVecExpr<T>(x, y, a);
}

template <typename T>
class ADVec3DotADVecExpr {
 public:
  ADVec3DotADVecExpr(ADVec<Vec<T, 3>>& xObj, ADVec<Vec<T, 3>>& yObj,
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
A2D_INLINE_FUNCTION ADVec3DotADVecExpr<T> Vec3Dot(ADVec<Vec<T, 3>>& x,
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
A2D_INLINE_FUNCTION void Vec3Cross(const Vec<T, 3>& x, const Vec<T, 3>& y,
                                   Vec<T, 3>& v) {
  Vec3CrossProductCore(x, y, v);
}

template <typename T>
A2D_INLINE_FUNCTION ADVec3CrossProductExpr<T> Vec3Cross(
    ADVec<Vec<T, 3>>& xobj, ADVec<Vec<T, 3>>& yobj, ADVec<Vec<T, 3>>& vobj) {
  return ADVec3CrossProductExpr<T>(xobj, yobj, vobj);
}

/*
  Vec3Normalize
  This represents the normalization of a vector
  vObj = xObj / ||xObj||
*/
template <typename T>
A2D_INLINE_FUNCTION void Vec3Normalize(const Vec<T, 3>& x, Vec<T, 3>& v) {
  Vec3ScaleCore(1 / Vec3NormCore<T>(x), x, v);
}

template <typename T>
class ADVec3NormalizeExpr {
 public:
  ADVec3NormalizeExpr(ADVec<Vec<T, 3>>& xObj, ADVec<Vec<T, 3>>& vObj)
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
    Vec3AXPBYCore(normInv, xb, negativeSquaredNormInv * Vec3DotCore<T>(xb, x),
                  v, vb);
  }

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& v = vObj.value();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3AXPBYCore(normInv, vb, negativeSquaredNormInv * Vec3DotCore<T>(vb, x),
                  v, xb);
  }

  ADVec<Vec<T, 3>>& xObj;
  ADVec<Vec<T, 3>>& vObj;

 private:
  T normInv;
  T negativeSquaredNormInv;
};

template <typename T>
A2D_INLINE_FUNCTION ADVec3NormalizeExpr<T> Vec3Normalize(ADVec<Vec<T, 3>>& x,
                                                         ADVec<Vec<T, 3>>& v) {
  return ADVec3NormalizeExpr<T>(x, v);
}

/*
  Vec3ScaleSymmetricOuterProduct
  This represents a scaled, 3 vector outer product with itself
  sObj = aObj * xObj . (xObj^T)
*/
template <typename T>
A2D_INLINE_FUNCTION void Vec3ScaleSymmetricOuterProduct(const T& a,
                                                        const Vec<T, 3>& x,
                                                        Mat<T, 3, 3>& S) {
  Vec3OuterProductScaleCore(a, x, x, S);
}

template <typename T>
class ADVec3ScaleSymmetricOuterProductExpr {
 public:
  ADVec3ScaleSymmetricOuterProductExpr(const T& a, ADVec<Vec<T, 3>>& xObj,
                                       ADMat<Mat<T, 3, 3>>& sObj)
      : xObj(xObj), sObj(sObj), a(a) {
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
A2D_INLINE_FUNCTION ADVec3ScaleSymmetricOuterProductExpr<T>
Vec3ScaleSymmetricOuterProduct(const T& a, ADVec<Vec<T, 3>>& x,
                               ADMat<Mat<T, 3, 3>>& S) {
  return ADVec3ScaleSymmetricOuterProductExpr<T>(a, x, S);
}

template <typename T>
class Vec3ADScaleSymmetricOuterProductExpr {
 public:
  Vec3ADScaleSymmetricOuterProductExpr(ADScalar<T>& aObj, const Vec<T, 3>& x,
                                       ADMat<Mat<T, 3, 3>>& sObj)
      : x(x), sObj(sObj), aObj(aObj) {
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
A2D_INLINE_FUNCTION Vec3ADScaleSymmetricOuterProductExpr<T>
Vec3ScaleSymmetricOuterProduct(ADScalar<T>& a, const Vec<T, 3>& x,
                               ADMat<Mat<T, 3, 3>>& S) {
  return Vec3ADScaleSymmetricOuterProductExpr<T>(a, x, S);
}

template <typename T>
class ADVec3ADScaleSymmetricOuterProductExpr {
 public:
  ADVec3ADScaleSymmetricOuterProductExpr(ADScalar<T>& aObj,
                                         ADVec<Vec<T, 3>>& xObj,
                                         ADMat<Mat<T, 3, 3>>& sObj)
      : xObj(xObj), sObj(sObj), aObj(aObj) {
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
A2D_INLINE_FUNCTION ADVec3ADScaleSymmetricOuterProductExpr<T>
Vec3ScaleSymmetricOuterProduct(ADScalar<T>& a, ADVec<Vec<T, 3>>& x,
                               ADMat<Mat<T, 3, 3>>& S) {
  return ADVec3ADScaleSymmetricOuterProductExpr<T>(a, x, S);
}

/********************************************************************************************************************
 **                                            A2D Vector Operations **
 ********************************************************************************************************************/

/*
    Vec3Norm
    This represents the L2 norm (or 2 norm) of a 3 vector
    normObj = ||xObj||
*/
template <int N, typename T>
class A2DVec3NormExpr {
 public:
  A2DVec3NormExpr(A2DVec<Vec<T, 3>>& xObj, A2DScalar<T>& normObj)
      : xObj(xObj), normObj(normObj) {
    const Vec<T, 3>& x = xObj.value();
    normObj.value = Vec3NormCore<T>(x);
    Vec3ScaleCore(1 / normObj.value, x, normX);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& n = normObj.value;
    const T& nb = normObj.bvalue;
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3ScaleCore(nb / n, x, xb);
  };

  void hforward() {
    /*
    //            input:
                const Vec<T, 3>& x = xObj.value();
                const T& norm = normObj.value;
    //            main loop: (uses decrements because comparison to 0 is faster)
                for (int i = N - 1; i != 0; i--) {
    //    //                input:
    //                    const Vec<T, 3>& xp = xObj.pvalue(i);
    //    //                output:
    //                    T& normp = normObj.pvalue[i];
    //    //                operations:
    //                    normp = Vec3DotCore<T>(x, xp) / norm;
                    normObj.pvalue[i] = Vec3DotCore<T>(normX, xObj.pvalue(i));
                }
    //            const Vec<T, 3>& xp = xObj.pvalue(0);
    //            normObj.pvalue[0] = Vec3DotCore<T>(x, xp) / norm;
                normObj.pvalue[0] = Vec3DotCore<T>(normX, xObj.pvalue(0));
    */

    //            main loop:
    for (int i = 0; i < N; i++) {
      normObj.pvalue[i] = Vec3DotCore<T>(normX, xObj.pvalue(i));
    }
    //            TODO: there are faster ways to do this instead of repeatedly
    //            calling Vec3DotCore
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& norm = normObj.value;
    const T& nb = normObj.bvalue;
    T invNorm = 1 / norm;
    /*
    //            main loop:
                for (int i = N - 1; i != 0; i--) {
    //                input:
                    const T& nh = normObj.hvalue(i);
                    const Vec<T, 3>& xp = xObj.pvalue(i);
    //                output:
                    Vec<T, 3>& xh = xObj.hvalue(i);
    //                operations:
                    Vec3ScaleCore(-Vec3DotCore<T>(normX, xp), normX, xh);
                    Vec3AddThenScaleCore(nb * invNorm, xp, xh);
                    Vec3ScaleThenAddCore(nh, normX, xh);
                }
    //                input:
                const T& nh = normObj.hvalue(0);
                const Vec<T, 3>& xp = xObj.pvalue(0);
    //                output:
                Vec<T, 3>& xh = xObj.hvalue(0);
    //                operations:
                Vec3ScaleCore(-Vec3DotCore<T>(normX, xp), normX, xh);
                Vec3AddThenScaleCore(nb * invNorm, xp, xh);
                Vec3ScaleThenAddCore(nh, normX, xh);
    */

    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const T& nh = normObj.hvalue[i];
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(-Vec3DotCore<T>(normX, xp), normX, xh);
      Vec3AddThenScaleCore(nb * invNorm, xp, xh);
      Vec3ScaleAndAddCore(nh, normX, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DScalar<T>& normObj;

 private:
  Vec<T, 3> normX;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3NormExpr<N, T> Vec3Norm(A2DVec<Vec<T, 3>>& x,
                                                   A2DScalar<T>& norm) {
  return A2DVec3NormExpr<N, T>(x, norm);
}

/*
    Vec3Scale
    This represents scalar multiplication of a 3 vector
    yObj = aObj * xObj
*/
template <int N, typename T>
class A2DVec3ScaleExpr {
 public:
  A2DVec3ScaleExpr(A2DVec<Vec<T, 3>>& xObj, const T& a, A2DVec<Vec<T, 3>>& vObj)
      : xObj(xObj), a(a), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3ScaleCore(a, vb, xb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3ScaleCore(a, xp, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(a, vh, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& vObj;
  const T& a;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3ScaleExpr<N, T> Vec3Scale(A2DVec<Vec<T, 3>>& x,
                                                     const T& a,
                                                     A2DVec<Vec<T, 3>>& v) {
  return A2DVec3ScaleExpr<N, T>(x, a, v);
}

template <int N, typename T>
class Vec3A2DScaleExpr {
 public:
  Vec3A2DScaleExpr(const Vec<T, 3>& x, A2DScalar<T>& aObj,
                   A2DVec<Vec<T, 3>>& vObj)
      : x(x), aObj(aObj), vObj(vObj) {
    const T& a = aObj.value;
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T, Vec<T, 3>>(vb, x);
  };

  void hforward() {
    //            input:
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3ScaleCore(ap, x, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x);
    }
  };

  const Vec<T, 3>& x;
  A2DVec<Vec<T, 3>>& vObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3A2DScaleExpr<N, T> Vec3Scale(const Vec<T, 3>& x,
                                                     A2DScalar<T>& a,
                                                     A2DVec<Vec<T, 3>>& v) {
  return Vec3A2DScaleExpr<N, T>(x, a, v);
}

template <int N, typename T>
class A2DVec3A2DScaleExpr {
 public:
  A2DVec3A2DScaleExpr(A2DVec<Vec<T, 3>>& xObj, A2DScalar<T>& aObj,
                      A2DVec<Vec<T, 3>>& vObj)
      : xObj(xObj), aObj(aObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    Vec<T, 3>& v = vObj.value();
    Vec3ScaleCore(a, x, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T, Vec<T, 3>>(vb, x);
    Vec3ScaleCore(a, vb, xb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPBYCore(ap, x, a, xp, vp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; i++) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& vh = vObj.hvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x) + Vec3DotCore<T>(vb, xp);
      Vec3AXPBYCore(a, vh, ap, vb, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& vObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3A2DScaleExpr<N, T> Vec3Scale(A2DVec<Vec<T, 3>>& x,
                                                        A2DScalar<T>& a,
                                                        A2DVec<Vec<T, 3>>& v) {
  return A2DVec3A2DScaleExpr<N, T>(x, a, v);
}

/*
    Vec3Axpy
    This operation is a compound scalar multiplication of a vector (Ax)
    and addition of another unscaled vector (py).  So Axpy means:
    "a (scalar) times (vector) x plus (vector) y"
*/
template <int N, typename T>
class Vec3VecA2DScalarAxpyExpr {
 public:
  Vec3VecA2DScalarAxpyExpr(A2DScalar<T>& aObj, const Vec<T, 3>& x,
                           const Vec<T, 3>& y, A2DVec<Vec<T, 3>>& vObj)
      : aObj(aObj), x(x), y(y), vObj(vObj) {
    const T& a = aObj.value;
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T>(vb, x);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3ScaleCore(ap, x, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x);
    }
  };

  A2DScalar<T>& aObj;
  const Vec<T, 3>& x;
  const Vec<T, 3>& y;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3VecA2DScalarAxpyExpr<N, T> Vec3Axpy(
    A2DScalar<T>& a, const Vec<T, 3>& x, const Vec<T, 3>& y,
    A2DVec<Vec<T, 3>>& v) {
  return Vec3VecA2DScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class A2DVec3VecScalarAxpyExpr {
 public:
  A2DVec3VecScalarAxpyExpr(const T& a, A2DVec<Vec<T, 3>>& xObj,
                           const Vec<T, 3>& y, A2DVec<Vec<T, 3>>& vObj)
      : a(a), xObj(xObj), y(y), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3ScaleCore(a, vb, xb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3ScaleCore(a, xp, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(a, vh, xh);
    }
  };

  const T& a;
  A2DVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3VecScalarAxpyExpr<N, T> Vec3Axpy(
    const T& a, A2DVec<Vec<T, 3>>& x, const Vec<T, 3>& y,
    A2DVec<Vec<T, 3>>& v) {
  return A2DVec3VecScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class Vec3A2DVecScalarAxpyExpr {
 public:
  Vec3A2DVecScalarAxpyExpr(const T& a, const Vec<T, 3>& x,
                           A2DVec<Vec<T, 3>>& yObj, A2DVec<Vec<T, 3>>& vObj)
      : a(a), x(x), yObj(yObj), vObj(vObj) {
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    yb(0) = vb(0);
    yb(1) = vb(1);
    yb(2) = vb(2);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      vp(0) = yp(0);
      vp(1) = yp(1);
      vp(2) = yp(2);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      yh(0) = vh(0);
      yh(1) = vh(1);
      yh(2) = vh(2);
    }
  };

  const T& a;
  const Vec<T, 3>& x;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3A2DVecScalarAxpyExpr<N, T> Vec3Axpy(
    const T& a, const Vec<T, 3>& x, A2DVec<Vec<T, 3>>& y,
    A2DVec<Vec<T, 3>>& v) {
  return Vec3A2DVecScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class A2DVec3VecA2DScalarAxpyExpr {
 public:
  A2DVec3VecA2DScalarAxpyExpr(A2DScalar<T>& aObj, A2DVec<Vec<T, 3>>& xObj,
                              const Vec<T, 3>& y, A2DVec<Vec<T, 3>>& vObj)
      : aObj(aObj), xObj(xObj), y(y), vObj(vObj) {
    const T& a = aObj.value;
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T>(vb, x);
    Vec3ScaleCore(a, vb, xb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ap = aObj.pvalue[i];
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPBYCore(ap, x, a, xp, vp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x) + Vec3DotCore<T>(vb, xp);
      Vec3AXPBYCore(a, vh, ap, vb, xh);
    }
  };

  A2DScalar<T>& aObj;
  A2DVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3VecA2DScalarAxpyExpr<N, T> Vec3Axpy(
    A2DScalar<T>& a, A2DVec<Vec<T, 3>>& x, const Vec<T, 3>& y,
    A2DVec<Vec<T, 3>>& v) {
  return A2DVec3VecA2DScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class Vec3A2DVecA2DScalarAxpyExpr {
 public:
  Vec3A2DVecA2DScalarAxpyExpr(A2DScalar<T>& aObj, const Vec<T, 3>& x,
                              A2DVec<Vec<T, 3>>& yObj, A2DVec<Vec<T, 3>>& vObj)
      : aObj(aObj), x(x), yObj(yObj), vObj(vObj) {
    const T& a = aObj.value;
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const T& a = aObj.value;
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T>(vb, x);
    yb(0) = vb(0);
    yb(1) = vb(1);
    yb(2) = vb(2);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ap = aObj.pvalue[i];
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPYCore(ap, x, yp, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x);
      yh(0) = vh(0);
      yh(1) = vh(1);
      yh(2) = vh(2);
    }
  };

  A2DScalar<T>& aObj;
  const Vec<T, 3>& x;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3A2DVecA2DScalarAxpyExpr<N, T> Vec3Axpy(
    A2DScalar<T>& a, const Vec<T, 3>& x, A2DVec<Vec<T, 3>>& y,
    A2DVec<Vec<T, 3>>& v) {
  return Vec3A2DVecA2DScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class A2DVec3A2DVecScalarAxpyExpr {
 public:
  A2DVec3A2DVecScalarAxpyExpr(const T& a, A2DVec<Vec<T, 3>>& xObj,
                              A2DVec<Vec<T, 3>>& yObj, A2DVec<Vec<T, 3>>& vObj)
      : a(a), xObj(xObj), yObj(yObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    Vec3ScaleCore(a, vb, xb);
    yb(0) = vb(0);
    yb(1) = vb(1);
    yb(2) = vb(2);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPYCore(a, xp, yp, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(a, vh, xh);
      yh(0) = vh(0);
      yh(1) = vh(1);
      yh(2) = vh(2);
    }
  };

  const T& a;
  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3A2DVecScalarAxpyExpr<N, T> Vec3Axpy(
    const T& a, A2DVec<Vec<T, 3>>& x, A2DVec<Vec<T, 3>>& y,
    A2DVec<Vec<T, 3>>& v) {
  return A2DVec3A2DVecScalarAxpyExpr<N, T>(a, x, y, v);
}

template <int N, typename T>
class A2DVec3A2DVecA2DScalarAxpyExpr {
 public:
  A2DVec3A2DVecA2DScalarAxpyExpr(A2DScalar<T>& aObj, A2DVec<Vec<T, 3>>& xObj,
                                 A2DVec<Vec<T, 3>>& yObj,
                                 A2DVec<Vec<T, 3>>& vObj)
      : aObj(aObj), xObj(xObj), yObj(yObj), vObj(vObj) {
    const T& a = aObj.value;
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3AXPYCore(a, x, y, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    aObj.bvalue = Vec3DotCore<T>(vb, x);
    Vec3ScaleCore(a, vb, xb);
    yb(0) = vb(0);
    yb(1) = vb(1);
    yb(2) = vb(2);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ap = aObj.pvalue[i];
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPBYCore(ap, x, a, xp, vp);
      Vec3AddInPlaceCore(yp, vp);
      //                Vec3AddCore(vp, yp, vp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const T& a = aObj.value;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Vec3DotCore<T>(vh, x) + Vec3DotCore<T>(vb, xp);
      Vec3AXPBYCore(a, vh, ap, vb, xh);
      yh(0) = vh(0);
      yh(1) = vh(1);
      yh(2) = vh(2);
    }
  };

  A2DScalar<T>& aObj;
  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3A2DVecA2DScalarAxpyExpr<N, T> Vec3Axpy(
    A2DScalar<T>& a, A2DVec<Vec<T, 3>>& x, A2DVec<Vec<T, 3>>& y,
    A2DVec<Vec<T, 3>>& v) {
  return A2DVec3A2DVecA2DScalarAxpyExpr<N, T>(a, x, y, v);
}

/*
    Vec3Dot
    This represents vector-vector dot product (inner product) of two 3 vectors
    aObj = xObj . yObj
    A.K.A.
    aObj = xObj @ yObj (in numpy notation)
*/
template <int N, typename T>
class A2DVec3DotVecExpr {
 public:
  A2DVec3DotVecExpr(A2DVec<Vec<T, 3>>& xObj, const Vec<T, 3>& y,
                    A2DScalar<T>& aObj)
      : xObj(xObj), y(y), aObj(aObj) {
    const Vec<T, 3>& x = xObj.value();
    aObj.value = Vec3DotCore<T>(x, y);
  };

  void reverse() {
    //            input:
    const T& ab = aObj.bvalue;
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3ScaleCore(ab, y, xb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                operations:
      aObj.pvalue[i] = Vec3DotCore<T>(xp, y);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ah = aObj.hvalue[i];
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(ah, y, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3DotVecExpr<N, T> Vec3Dot(A2DVec<Vec<T, 3>>& x,
                                                    const Vec<T, 3>& y,
                                                    A2DScalar<T>& a) {
  return A2DVec3DotVecExpr<N, T>(x, y, a);
}

template <int N, typename T>
class Vec3DotA2DVecExpr {
 public:
  Vec3DotA2DVecExpr(const Vec<T, 3>& x, A2DVec<Vec<T, 3>>& yObj,
                    A2DScalar<T>& aObj)
      : x(x), yObj(yObj), aObj(aObj) {
    const Vec<T, 3>& y = yObj.value();
    aObj.value = Vec3DotCore<T>(x, y);
  };

  void reverse() {
    //            input:
    const T& ab = aObj.bvalue;
    //            output:
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    Vec3ScaleCore(ab, x, yb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                operations:
      aObj.pvalue[i] = Vec3DotCore<T>(x, yp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ah = aObj.hvalue[i];
      //                output:
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      Vec3ScaleCore(ah, x, yh);
    }
  };

  const Vec<T, 3>& x;
  A2DVec<Vec<T, 3>>& yObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3DotA2DVecExpr<N, T> Vec3Dot(const Vec<T, 3>& x,
                                                    A2DVec<Vec<T, 3>>& y,
                                                    A2DScalar<T>& a) {
  return Vec3DotA2DVecExpr<N, T>(x, y, a);
}

template <int N, typename T>
class A2DVec3DotA2DVecExpr {
 public:
  A2DVec3DotA2DVecExpr(A2DVec<Vec<T, 3>>& xObj, A2DVec<Vec<T, 3>>& yObj,
                       A2DScalar<T>& aObj)
      : xObj(xObj), yObj(yObj), aObj(aObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();

    aObj.value = Vec3DotCore<T>(x, y);
  };

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
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                operations:
      aObj.pvalue[i] = Vec3DotCore<T>(xp, y) + Vec3DotCore<T>(x, yp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const T& ab = aObj.bvalue;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ah = aObj.hvalue[i];
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      Vec3AXPBYCore(ah, y, ab, yp, xh);
      Vec3AXPBYCore(ah, x, ab, xp, yh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& yObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3DotA2DVecExpr<N, T> Vec3Dot(A2DVec<Vec<T, 3>>& x,
                                                       A2DVec<Vec<T, 3>>& y,
                                                       A2DScalar<T>& a) {
  return A2DVec3DotA2DVecExpr<N, T>(x, y, a);
}

/*
    Vec3CrossProduct
    This represents vector-vector cross product of two 3 vectors
*/
template <int N, typename T>
class A2DVec3CrossVecExpr {
 public:
  A2DVec3CrossVecExpr(A2DVec<Vec<T, 3>>& xObj, const Vec<T, 3>& y,
                      A2DVec<Vec<T, 3>>& vObj)
      : xObj(xObj), y(y), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3CrossProductCore(x, y, v);
  }

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3CrossProductAddCore(y, vb, xb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3CrossProductCore(xp, y, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Vec3CrossProductCore(y, vh, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  const Vec<T, 3>& y;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3CrossVecExpr<N, T> Vec3Cross(A2DVec<Vec<T, 3>>& x,
                                                        const Vec<T, 3>& y,
                                                        A2DVec<Vec<T, 3>>& v) {
  return A2DVec3CrossVecExpr<N, T>(x, y, v);
}

template <int N, typename T>
class Vec3CrossA2DVecExpr {
 public:
  Vec3CrossA2DVecExpr(const Vec<T, 3>& x, A2DVec<Vec<T, 3>>& yObj,
                      A2DVec<Vec<T, 3>>& vObj)
      : x(x), yObj(yObj), vObj(vObj) {
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3CrossProductCore(x, y, v);
  }

  void reverse() {
    //            input:
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    Vec3CrossProductAddCore(vb, x, yb);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3CrossProductCore(x, yp, vp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      Vec3CrossProductCore(vh, x, yh);
    }
  };

  const Vec<T, 3>& x;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3CrossA2DVecExpr<N, T> Vec3Cross(const Vec<T, 3>& x,
                                                        A2DVec<Vec<T, 3>>& y,
                                                        A2DVec<Vec<T, 3>>& v) {
  return Vec3CrossA2DVecExpr<N, T>(x, y, v);
}

template <int N, typename T>
class A2DVec3CrossA2DVecExpr {
 public:
  A2DVec3CrossA2DVecExpr(A2DVec<Vec<T, 3>>& xObj, A2DVec<Vec<T, 3>>& yObj,
                         A2DVec<Vec<T, 3>>& vObj)
      : xObj(xObj), yObj(yObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    Vec<T, 3>& v = vObj.value();
    Vec3CrossProductCore(x, y, v);
  }

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    Vec<T, 3>& yb = yObj.bvalue();
    //            operations:
    Vec3CrossProductAddCore(y, vb, xb);
    Vec3CrossProductAddCore(vb, x, yb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3CrossProductCore(xp, y, vp);
      Vec3CrossProductAddCore(x, yp, vp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& y = yObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& yp = yObj.pvalue(i);
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      Vec<T, 3>& yh = yObj.hvalue(i);
      //                operations:
      Vec3CrossProductCore(y, vh, xh);
      Vec3CrossProductAddCore(yp, vb, xh);
      Vec3CrossProductCore(vh, x, yh);
      Vec3CrossProductAddCore(vb, xp, yh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& yObj;
  A2DVec<Vec<T, 3>>& vObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3CrossA2DVecExpr<N, T> Vec3Cross(
    A2DVec<Vec<T, 3>>& x, A2DVec<Vec<T, 3>>& y, A2DVec<Vec<T, 3>>& v) {
  return A2DVec3CrossA2DVecExpr<N, T>(x, y, v);
}

/*
    Vec3Normalize
    This represents the normalization of a vector
    vObj = xObj / ||xObj||
*/
template <int N, typename T>
class A2DVec3NormalizeExpr {
 public:
  A2DVec3NormalizeExpr(A2DVec<Vec<T, 3>>& xObj, A2DVec<Vec<T, 3>>& vObj)
      : xObj(xObj), vObj(vObj) {
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();
    normInv = sqrt(1 / Vec3DotCore<T>(x, x));
    Vec3ScaleCore(normInv, x, v);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& v = vObj.value();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Vec3AXPBYCore(normInv, vb, -normInv * normInv * Vec3DotCore<T>(vb, x), v,
                  xb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& v = vObj.value();
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& vp = vObj.pvalue(i);
      //                operations:
      Vec3AXPYCore(-Vec3DotCore<T>(xp, v), v, xp, vp);
      Vec3ScaleCore(normInv, vp, vp);
    }
  };

  void hreverse() {
    //            input:
    const Vec<T, 3>& v = vObj.value();
    const Vec<T, 3>& vb = vObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    //            operations:
    T temp1{Vec3DotCore<T>(vb, v)}, temp2;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const Vec<T, 3>& vh = vObj.hvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      temp2 = Vec3DotCore<T>(xp, v);
      Vec3AXPBYCore(-temp2, vb, -temp1, xp, xh);
      Vec3AXPYCore(
          3 * temp2 * temp1 - Vec3DotCore<T>(vh, x) - Vec3DotCore<T>(vb, xp), v,
          xh, xh);
      Vec3AXPYCore(normInv, xh, vh, xh);
      Vec3ScaleCore(normInv, xh, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DVec<Vec<T, 3>>& vObj;

 private:
  T normInv;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3NormalizeExpr<N, T> Vec3Normalize(
    A2DVec<Vec<T, 3>>& x, A2DVec<Vec<T, 3>>& v) {
  return A2DVec3NormalizeExpr<N, T>(x, v);
}

/*
    Vec3ScaleSymmetricOuterProduct
    This represents a scaled, 3 vector outer product with itself
    sObj = aObj * xObj . (xObj^T)
*/
template <int N, typename T>
class A2DVec3ScaleSymmetricOuterProductExpr {
 public:
  A2DVec3ScaleSymmetricOuterProductExpr(const T& a, A2DVec<Vec<T, 3>>& xObj,
                                        A2DMat<Mat<T, 3, 3>>& sObj)
      : a(a), xObj(xObj), sObj(sObj) {
    const Vec<T, 3>& x = xObj.value();
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    Mat3x3VecMultScaleCore(2 * a, Sb, x, xb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3> ax;
    //            operations:
    Vec3ScaleCore<T, Vec<T, 3>>(a, x, ax);
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Mat<T, 3, 3>& Sp = sObj.pvalue(i);
      //                operations:
      Vec3OuterProductCore<Vec<T, 3>, Mat<T, 3, 3>>(xp, ax, Sp);
      Vec3OuterProductAddCore<Vec<T, 3>, Mat<T, 3, 3>>(ax, xp, Sp);
    }
  };

  void hreverse() {
    //            input:
    const Mat<T, 3, 3>& sb = sObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3> sbxp, shx, x_sx;
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Mat<T, 3, 3>& sh = sObj.hvalue(i);
      const Vec<T, 3>& xp = xObj.pvalue(i);
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Mat3x3VecMultCore(sb, xp, sbxp);
      Mat3x3VecMultCore(sh, x, shx);
      Vec3AddCore(shx, sbxp, x_sx);
      Vec3ScaleCore(2 * a, x_sx, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DMat<Mat<T, 3, 3>>& sObj;
  const T& a;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3ScaleSymmetricOuterProductExpr<N, T>
Vec3ScaleSymmetricOuterProduct(const T& a, A2DVec<Vec<T, 3>>& x,
                               A2DMat<Mat<T, 3, 3>>& S) {
  return A2DVec3ScaleSymmetricOuterProductExpr<N, T>(a, x, S);
}

template <int N, typename T>
class Vec3A2DScaleSymmetricOuterProductExpr {
 public:
  Vec3A2DScaleSymmetricOuterProductExpr(A2DScalar<T>& aObj, const Vec<T, 3>& x,
                                        A2DMat<Mat<T, 3, 3>>& sObj)
      : aObj(aObj), x(x), sObj(sObj) {
    const T& a = aObj.value;
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void reverse() {
    //            input:
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
    //            operations:
    aObj.bvalue = Mat3x3InnerProductCore<T, Vec<T, 3>, Mat<T, 3, 3>>(Sb, x, x);
  };

  void hforward() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const T& ap = aObj.pvalue[i];
      //                output:
      Mat<T, 3, 3>& Sp = sObj.pvalue(i);
      //                operations:
      Vec3OuterProductScaleCore(ap, x, x, Sp);
    }
  };

  void hreverse() {
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Mat<T, 3, 3>& sh = sObj.hvalue(i);
      //                operations:
      aObj.hvalue[i] = Mat3x3InnerProductCore<T>(sh, x, x);
    }
  };

  const Vec<T, 3>& x;
  A2DMat<Mat<T, 3, 3>>& sObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION Vec3A2DScaleSymmetricOuterProductExpr<N, T>
Vec3ScaleSymmetricOuterProduct(A2DScalar<T>& a, const Vec<T, 3>& x,
                               A2DMat<Mat<T, 3, 3>>& S) {
  return Vec3A2DScaleSymmetricOuterProductExpr<N, T>(a, x, S);
}

template <int N, typename T>
class A2DVec3A2DScaleSymmetricOuterProductExpr {
 public:
  A2DVec3A2DScaleSymmetricOuterProductExpr(A2DScalar<T>& aObj,
                                           A2DVec<Vec<T, 3>>& xObj,
                                           A2DMat<Mat<T, 3, 3>>& sObj)
      : aObj(aObj), xObj(xObj), sObj(sObj) {
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    Mat<T, 3, 3>& S = sObj.value();
    Vec3OuterProductScaleCore(a, x, x, S);
  };

  void reverse() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const Mat<T, 3, 3>& Sb = sObj.bvalue();
    const T& a = aObj.value;
    //            output:
    Vec<T, 3>& xb = xObj.bvalue();
    //            operations:
    aObj.bvalue = Mat3x3InnerProductCore<T, Vec<T, 3>, Mat<T, 3, 3>>(Sb, x, x);
    Mat3x3VecMultScaleCore(2 * a, Sb, x, xb);
  };

  void hforward() {
    //            input:
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    //            operations:
    Vec<T, 3> ax;
    Vec3ScaleCore<T, Vec<T, 3>>(a, x, ax);
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Mat<T, 3, 3>& Sp = sObj.pvalue(i);
      //                operations:
      Vec3OuterProductCore<Vec<T, 3>, Mat<T, 3, 3>>(xp, ax, Sp);
      Vec3OuterProductAddCore<Vec<T, 3>, Mat<T, 3, 3>>(ax, xp, Sp);
      Vec3OuterProductAddScaleCore<T, Vec<T, 3>, Mat<T, 3, 3>>(ap, x, x, Sp);
    }
  };

  void hreverse() {
    //            input:
    const Mat<T, 3, 3>& sb = sObj.bvalue();
    const Vec<T, 3>& x = xObj.value();
    const T& a = aObj.value;
    Vec<T, 3> sbxp, shx, sbx, a_sx, x_sx;
    //            operations:
    Mat3x3VecMultCore(sb, x, sbx);
    //            main loop:
    for (int i = 0; i < N; ++i) {
      //                input:
      const Mat<T, 3, 3>& sh = sObj.hvalue(i);
      const Vec<T, 3>& xp = xObj.pvalue(i);
      const T& ap = aObj.pvalue[i];
      //                output:
      Vec<T, 3>& xh = xObj.hvalue(i);
      //                operations:
      Mat3x3VecMultCore(sb, xp, sbxp);
      Mat3x3VecMultCore(sh, x, shx);
      Vec3AXPYCore(2, sbxp, shx, a_sx);
      Vec3AddCore(shx, sbxp, x_sx);
      aObj.hvalue[i] = Vec3DotCore<T>(x, a_sx);
      Vec3AXPBYCore(2 * a, x_sx, 2 * ap, sbx, xh);
    }
  };

  A2DVec<Vec<T, 3>>& xObj;
  A2DMat<Mat<T, 3, 3>>& sObj;
  A2DScalar<T>& aObj;
};

template <int N, typename T>
A2D_INLINE_FUNCTION A2DVec3A2DScaleSymmetricOuterProductExpr<N, T>
Vec3ScaleSymmetricOuterProduct(A2DScalar<T>& a, A2DVec<Vec<T, 3>>& x,
                               A2DMat<Mat<T, 3, 3>>& S) {
  return A2DVec3A2DScaleSymmetricOuterProductExpr<N, T>(a, x, S);
}

}  // namespace A2D

#endif  // A2D_VEC_OPS_H
