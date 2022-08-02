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

template <int N, class ScalarType>
class A2DScalar {
 public:
  A2D_INLINE_FUNCTION A2DScalar(ScalarType value = 0.0, ScalarType bvalue = 0.0)
      : value(value), bvalue(bvalue) {
    for (int i = 0; i < N; i++) {
      pvalue[i] = 0.0;
      hvalue[i] = 0.0;
    }
  }

  ScalarType value;
  ScalarType bvalue;
  ScalarType pvalue[N];
  ScalarType hvalue[N];
};

template <int N, class MatType>
class A2DMat {
 public:
  A2D_INLINE_FUNCTION A2DMat(MatType& A, MatType& Ab) : A(A), Ab(Ab) {}

  A2D_INLINE_FUNCTION MatType& value() { return A; }
  A2D_INLINE_FUNCTION const MatType& value() const { return A; }

  A2D_INLINE_FUNCTION MatType& bvalue() { return Ab; }
  A2D_INLINE_FUNCTION const MatType& bvalue() const { return Ab; }

  A2D_INLINE_FUNCTION MatType& pvalue(const int i) { return Ap[i]; }
  A2D_INLINE_FUNCTION const MatType& pvalue(const int i) const { return Ap[i]; }

  A2D_INLINE_FUNCTION MatType& hvalue(const int i) { return Ah[i]; }
  A2D_INLINE_FUNCTION const MatType& hvalue(const int i) const { return Ah[i]; }

  MatType& A;     // Matrix
  MatType& Ab;    // Reverse mode derivative value
  MatType Ap[N];  // Projected second derivative value
  MatType Ah[N];  // Reverse mode second derivative
};

}  // namespace A2D

#endif  // A2D_TYPES_H
