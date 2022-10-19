#ifndef A2D_VEC_OPS_H
#define A2D_VEC_OPS_H

#include "a2dobjs.h"
#include "a2dtypes.h"
#include "a2dveccore3d.h"

namespace A2D {

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

}  // namespace A2D

#endif  // A2D_VEC_OPS_H
