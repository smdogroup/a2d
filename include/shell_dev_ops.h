//
// Created by James on 10/28/2022.
//

#ifndef A2D_SHELL_DEV_OPS
#define A2D_SHELL_DEV_OPS

namespace A2D {

/**
 * ScalarMult operation (z = a * b)
 * */
template <int N, typename T>  // TODO: make the ScalarMult operation (z = a * b)
class A2DScalarScalarMultExpr {
// TODO
 public:
  A2DScalarScalarMultExpr(A2DScalar<N, T>& aObj,
                          const T& b,
                          A2DScalar<N, T>& zObj) {

  };
};

template <int N, typename T>
class A2DScalarA2DScalarMultExpr {
 public:
  A2DScalarA2DScalarMultExpr(A2DScalar<N, T>& aObj,
                             A2DScalar<N, T>& bObj,
                             A2DScalar<N, T>& zObj)
      : aObj(aObj), bObj(bObj), zObj(zObj) {
    zObj.value = aObj.value * bObj.value;
  };

  void reverse();
  void hforward();
  void hreverse();

  A2DScalar<N, T>& aObj;
  A2DScalar<N, T>& bObj;
  A2DScalar<N, T>& zObj;
};

/**
 *  Vec3ScaleDiv operation (v = (1/a) * x)
 * */
template <int N, typename T>  // TODO: make the Vec3ScaleDiv operation (v = (1/a) * x)
class A2DVec3A2DScaleDivExpr {
 public:
  A2DVec3A2DScaleDivExpr(A2DVec<N, Vec<T, 3>>& xObj,
                         A2DScalar<N, T>& aObj,
                         A2DVec<N, Vec<T, 3>>& vObj)
      : xObj(xObj), aObj(aObj), vObj(vObj) {
    const T a = aObj.value;
    const Vec<T, 3>& x = xObj.value();
    Vec<T, 3>& v = vObj.value();

    v(0) = x(0) / a;
    v(1) = x(1) / a;
    v(2) = x(2) / a;
  };

  void reverse();
  void hforward();
  void hreverse();

  A2DVec<N, Vec<T, 3>>& xObj;
  A2DScalar<N, T>& aObj;
  A2DVec<N, Vec<T, 3>>& vObj;
};

/**
 *  ScalarAxpay operation (z = a * (x + y))
 * */
template <int N, typename T>  // TODO: make the ScalarAxpay operation (z = a * (x + y))
class ScalarA2DScalarA2DScalarAxpayExpr {
 public:
  ScalarA2DScalarA2DScalarAxpayExpr(const T& a,
                                    A2DScalar<N, T>& xObj,
                                    A2DScalar<N, T>& yObj,
                                    A2DScalar<N, T>& zObj)
      : a(a), xObj(xObj), yObj(yObj), zObj(zObj) {
    zObj.value = a * (xObj.value + yObj.value);
  };

  void reverse();
  void hforward();
  void hreverse();

  const T& a;
  A2DScalar<N, T>& xObj;
  A2DScalar<N, T>& yObj;
  A2DScalar<N, T>& zObj;
};

/**
 *  ScalarAxpby operation (z = a * x + b * y)
 * */
template <int N, typename T>  // TODO: make the ScalarAxpby operation (z = a * x + b * y)
class ScalarA2DScalarScalarA2DScalarAxpbyExpr {
 public:
  ScalarA2DScalarScalarA2DScalarAxpbyExpr(const T& a,
                                          A2DScalar<N, T>& xObj,
                                          const T& b,
                                          A2DScalar<N, T>& yObj,
                                          A2DScalar<N, T>& zObj)
      : a(a), xObj(xObj), b(b), yObj(yObj), zObj(zObj) {
    zObj.value = a * xObj.value + b * yObj.value;
  };

  const T& a;
  A2DScalar<N, T>& xObj;
  const T& b;
  A2DScalar<N, T>& yObj;
  A2DScalar<N, T>& zObj;
};

template <int N, typename T>
class A2DScalarA2DScalarA2DScalarA2DScalarAxpbyExpr {
 public:
  A2DScalarA2DScalarA2DScalarA2DScalarAxpbyExpr(A2DScalar<N, T>& aObj,
                                                A2DScalar<N, T>& xObj,
                                                A2DScalar<N, T>& bObj,
                                                A2DScalar<N, T>& yObj,
                                                A2DScalar<N, T>& zObj)
      : aObj(aObj), xObj(xObj), bObj(bObj), yObj(yObj), zObj(zObj) {
    zObj.value = aObj.value * xObj.value + bObj.value * yObj.value;
  };

  void reverse();
  void hforward();
  void hreverse();

  A2DScalar<N, T>& aObj;
  A2DScalar<N, T>& xObj;
  A2DScalar<N, T>& bObj;
  A2DScalar<N, T>& yObj;
  A2DScalar<N, T>& zObj;
};

/**
 *  MatInnerProduct operation (a = x.A.y)
 * */
template <typename T, int P, int Q>
// TODO: inspect this to see how it is compiled
inline T MatInnerProductCore(const Mat<T, P, Q>& A,
                             const Vec<T, P>& x,
                             const Vec<T, Q>& y) {
  T a{0}, b;
  for (int i = 0; i < P; ++i) {
    b = 0;
    for (int j = 0; j < Q; ++j) {
      b += A(i, j) * y(j);
    }
    a += b * x(i);
  }
  return a;
};

template <int N, typename T, int P, int Q>  // TODO: make the MatInnerProduct operation (a = x.A.y)
class MatA2DVecA2DVecInnerProductExpr {
 public:
  MatA2DVecA2DVecInnerProductExpr(const Mat<T, P, Q>& A,
                                  A2DVec<N, Vec<T, P>>& xObj,
                                  A2DVec<N, Vec<T, Q>>& yObj,
                                  A2DScalar<N, T>& aObj)
      : A(A), xObj(xObj), yObj(yObj), aObj(aObj) {
    const Vec<T, P>& x = xObj.value();
    const Vec<T, Q>& y = yObj.value();
    aObj.value = MatInnerProductCore(A, x, y);
  };

  void reverse();
  void hforward();
  void hreverse();

  const Mat<T, P, Q>& A;
  A2DVec<N, Vec<T, P>>& xObj;
  A2DVec<N, Vec<T, Q>>& yObj;
  A2DScalar<N, T>& aObj;
};

/**
 * Vec5ScalarAssembly operation (v = {x0, x1, x2, x3, x4})
 * */
template <int N, typename T>  // TODO: make the Vec5ScalarAssembly operation (v = {x0, x1, x2, x3, x4})
class A2DScalar5VecAssemblyExpr {
 public:
  A2DScalar5VecAssemblyExpr(A2DVec<N, Vec<T, 5>>& vObj,
                            A2DScalar<N, T>& x0Obj,
                            A2DScalar<N, T>& x1Obj,
                            A2DScalar<N, T>& x2Obj,
                            A2DScalar<N, T>& x3Obj,
                            A2DScalar<N, T>& x4Obj)
      : vObj(vObj), x0Obj(x0Obj), x1Obj(x1Obj), x2Obj(x2Obj), x3Obj(x3Obj), x4Obj(x4Obj) {
    Vec<T, 5>& v = vObj.value();
    v(0) = x0Obj.value;
    v(1) = x1Obj.value;
    v(2) = x2Obj.value;
    v(3) = x3Obj.value;
    v(4) = x4Obj.value;
  };

  void reverse();
  void hforward();
  void hreverse();

  A2DVec<N, Vec<T, 5>>& vObj;
  A2DScalar<N, T>& x0Obj;
  A2DScalar<N, T>& x1Obj;
  A2DScalar<N, T>& x2Obj;
  A2DScalar<N, T>& x3Obj;
  A2DScalar<N, T>& x4Obj;
};


// TODO: figure out how to do this at a later date.
/*template <int N, typename T, int D>  // TODO: make the VecScalarAssembly operation (v = {v1, v2, v3, ..., vD})
class A2DScalarVecAssemblyExpr {
  A2DScalarVecAssemblyExpr(A2DVec<N, Vec<T, D>>& vObj,
                           A2DScalar<N, T> x[D]){

    for (int i = 0; i < D; ++i) {
      v.value()(i) = x[i].value;
    }
  };
};
template <int N, typename T, int D>
A2DScalarVecAssemblyExpr<N, T, D> VecScalarAssembly(A2DVec<N, Vec<T, D>>& v, A2DScalar<N, T> x[D]) {
  return A2DScalarVecAssemblyExpr<N, T, D>(v, x);
};*/

}  // namespace A2D

#endif //A2D_SHELL_DEV_OPS
