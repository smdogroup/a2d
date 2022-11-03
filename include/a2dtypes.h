#ifndef A2D_TYPES_H
#define A2D_TYPES_H

#include "a2dobjs.h"

namespace A2D {

template <class A>
class ADExpression {
 public:
  A2D_INLINE_FUNCTION A& cast() { return static_cast<A&>(*this); }
  A2D_INLINE_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  A2D_INLINE_FUNCTION void forward() { cast().forward(); }
  A2D_INLINE_FUNCTION void reverse() { cast().reverse(); }
};

template <class A>
class A2DExpression {
 public:
  A2D_INLINE_FUNCTION A& cast() { return static_cast<A&>(*this); }
  A2D_INLINE_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  A2D_INLINE_FUNCTION void reverse() { cast().reverse(); }
  A2D_INLINE_FUNCTION void hforward() { cast().hforward(); }
  A2D_INLINE_FUNCTION void hreverse() { cast().hreverse(); }
};

template <class ScalarType>
class ADScalar {
 public:
  A2D_INLINE_FUNCTION ADScalar(ScalarType value = 0.0, ScalarType bvalue = 0.0)
      : value(value), bvalue(bvalue) {}

  ScalarType value;
  ScalarType bvalue;
};

template <class MatType>
class ADMat {
 public:
  A2D_INLINE_FUNCTION ADMat(MatType& A, MatType& Ab) : A(A), Ab(Ab) {}

  A2D_INLINE_FUNCTION MatType& value() { return A; }
  A2D_INLINE_FUNCTION const MatType& value() const { return A; }

  A2D_INLINE_FUNCTION MatType& bvalue() { return Ab; }
  A2D_INLINE_FUNCTION const MatType& bvalue() const { return Ab; }

  MatType& A;   // Matrix
  MatType& Ab;  // Reverse mode derivative value
};

template <class VecType>
class ADVec {
 public:
  A2D_INLINE_FUNCTION ADVec(VecType& V, VecType& Vb) : V(V), Vb(Vb) {}

  A2D_INLINE_FUNCTION VecType& value() { return V; }
  A2D_INLINE_FUNCTION const VecType& value() const { return V; }

  A2D_INLINE_FUNCTION VecType& bvalue() { return Vb; }
  A2D_INLINE_FUNCTION const VecType& bvalue() const { return Vb; }

  VecType& V;   // Vector
  VecType& Vb;  // Reverse mode derivative value
};

template <class ScalarType>
class A2DScalar {
 public:
  A2D_INLINE_FUNCTION A2DScalar(ScalarType value = 0.0, ScalarType bvalue = 0.0,
                                ScalarType pvalue = 0.0,
                                ScalarType hvalue = 0.0)
      : value(value), bvalue(bvalue), pvalue(pvalue), hvalue(hvalue) {}

  ScalarType value;
  ScalarType bvalue;
  ScalarType pvalue;
  ScalarType hvalue;
};

template <class VecType>
class A2DVec {
 public:
  A2D_INLINE_FUNCTION A2DVec() {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x) : x(x) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb)
      : x(x), xb(xb) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb,
                             const VecType& xp)
      : x(x), xb(xb), xp(xp) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb,
                             const VecType& xp, const VecType& xh)
      : x(x), xb(xb), xp(xp), xh(xh) {}

  A2D_INLINE_FUNCTION VecType& value() { return x; }
  A2D_INLINE_FUNCTION const VecType& value() const { return x; }

  A2D_INLINE_FUNCTION VecType& bvalue() { return xb; }
  A2D_INLINE_FUNCTION const VecType& bvalue() const { return xb; }

  A2D_INLINE_FUNCTION VecType& pvalue() { return xp; }
  A2D_INLINE_FUNCTION const VecType& pvalue() const { return xp; }

  A2D_INLINE_FUNCTION VecType& hvalue() { return xh; }
  A2D_INLINE_FUNCTION const VecType& hvalue() const { return *xh; }

  VecType x;
  VecType xb;
  VecType xp;
  VecType xh;
};

template <class MatType>
class A2DMat {
 public:
  A2D_INLINE_FUNCTION A2DMat() {}
  A2D_INLINE_FUNCTION A2DMat(const MatType& A) : A(A) {}
  A2D_INLINE_FUNCTION A2DMat(const MatType& A, const MatType& Ab)
      : A(A), Ab(Ab) {}
  A2D_INLINE_FUNCTION A2DMat(const MatType& A, const MatType& Ab,
                             const MatType& Ap)
      : A(A), Ab(Ab), Ap(Ap) {}
  A2D_INLINE_FUNCTION A2DMat(const MatType& A, const MatType& Ab,
                             const MatType& Ap, const MatType& Ah)
      : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  A2D_INLINE_FUNCTION MatType& value() { return A; }
  A2D_INLINE_FUNCTION const MatType& value() const { return A; }

  A2D_INLINE_FUNCTION void set_bvalue(const MatType& val) { Ab.set(val); }
  A2D_INLINE_FUNCTION void get_bvalue(MatType& val) { Ab.get(val); }
  A2D_INLINE_FUNCTION MatType& bvalue() { return Ab; }
  A2D_INLINE_FUNCTION const MatType& bvalue() const { return Ab; }

  A2D_INLINE_FUNCTION void set_pvalue(const MatType& val) { Ap.set(val); }
  A2D_INLINE_FUNCTION void get_pvalue(MatType& val) { Ap.get(val); }
  A2D_INLINE_FUNCTION MatType& pvalue() { return Ap; }
  A2D_INLINE_FUNCTION const MatType& pvalue() const { return Ap; }

  A2D_INLINE_FUNCTION void set_hvalue(const MatType& val) { Ah.set(val); }
  A2D_INLINE_FUNCTION void get_hvalue(MatType& val) { Ah.get(val); }
  A2D_INLINE_FUNCTION MatType& hvalue() { return Ah; }
  A2D_INLINE_FUNCTION const MatType& hvalue() const { return Ah; }

  MatType A;   // Matrix
  MatType Ab;  // Reverse mode derivative value
  MatType Ap;  // Projected second derivative value
  MatType Ah;  // Reverse mode second derivative
};

}  // namespace A2D

#endif  // A2D_TYPES_H
