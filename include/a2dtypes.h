#ifndef A2D_TYPES_H
#define A2D_TYPES_H

namespace A2D {

template <class A>
class ADExpression {
 public:
  A& cast() { return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void forward() { cast().forward(); }
  void reverse() { cast().reverse(); }
};

template <class A>
class A2DExpression {
 public:
  A& cast() { return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void reverse() { cast().reverse(); }
  void hforward() { cast().hforward(); }
  void hreverse() { cast().hreverse(); }
};

template <class ScalarType>
class ADScalar {
 public:
  ADScalar(ScalarType value = 0.0, ScalarType bvalue = 0.0)
      : value(value), bvalue(bvalue) {}

  ScalarType value;
  ScalarType bvalue;
};

template <class MatType>
class ADMat {
 public:
  ADMat(MatType& A, MatType& Ab) : A(A), Ab(Ab) {}

  MatType& value() { return A; }
  const MatType& value() const { return A; }

  MatType& bvalue() { return Ab; }
  const MatType& bvalue() const { return Ab; }

  MatType& A;   // Matrix
  MatType& Ab;  // Reverse mode derivative value
};

template <class VecType>
class ADVec {
 public:
  ADVec(VecType& V, VecType& Vb) : V(V), Vb(Vb) {}

  VecType& value() { return V; }
  const VecType& value() const { return V; }

  VecType& bvalue() { return Vb; }
  const VecType& bvalue() const { return Vb; }

  VecType& V;   // Vector
  VecType& Vb;  // Reverse mode derivative value
};

template <int N, class ScalarType>
class A2DScalar {
 public:
  A2DScalar(ScalarType value = 0.0, ScalarType bvalue = 0.0)
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
  A2DMat(MatType& A, MatType& Ab) : A(A), Ab(Ab) {}

  MatType& value() { return A; }
  const MatType& value() const { return A; }

  MatType& bvalue() { return Ab; }
  const MatType& bvalue() const { return Ab; }

  MatType& pvalue(const int i) { return Ap[i]; }
  const MatType& pvalue(const int i) const { return Ap[i]; }

  MatType& hvalue(const int i) { return Ah[i]; }
  const MatType& hvalue(const int i) const { return Ah[i]; }

  MatType& A;     // Matrix
  MatType& Ab;    // Reverse mode derivative value
  MatType Ap[N];  // Projected second derivative value
  MatType Ah[N];  // Reverse mode second derivative
};

}  // namespace A2D

#endif  // A2D_TYPES_H
