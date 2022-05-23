#ifndef A2D_TYPES_H
#define A2D_TYPES_H

namespace A2D {

template <class A>
class ADExpression {
public:
  A& cast(){ return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void forward(){ cast().forward(); }
  void reverse(){ cast().reverse(); }
};

template <class A>
class A2DExpression {
public:
  A& cast(){ return static_cast<A&>(*this); }
  const A& cast() const { return static_cast<const A&>(*this); }

  void forward(){ cast().forward(); }
  void reverse(){ cast().reverse(); }

  void hforward(){ cast().hforward(); }
  void hproduct(){ cast().hproduct(); }
  void hreverse(){ cast().hreverse(); }
};

template<class ScalarType>
class ADScalar {
 public:
 ADScalar( ScalarType value=0.0, ScalarType bvalue=0.0 ) : value(value), bvalue(bvalue) {}

  ScalarType value;
  ScalarType bvalue;
};

template<class ScalarType>
class A2DScalar {
public:
  A2DScalar( ScalarType value=0.0, ScalarType bvalue=0.0,
             ScalarType pvalue=0.0, ScalarType hvalue=0.0 ) :
  value(value), bvalue(bvalue), pvalue(pvalue), hvalue(hvalue) {}

  ScalarType value;
  ScalarType bvalue;
  ScalarType pvalue;
  ScalarType hvalue;
};


template <class MatType>
class ADMat {
public:
 ADMat( MatType& A, MatType& Ab ) : A(A), Ab(Ab) {}

  MatType& value(){ return A; }
  const MatType& value() const { return A; }

  MatType& bvalue(){ return Ab; }
  const MatType& bvalue() const { return Ab; }

  MatType& A; // Matrix
  MatType& Ab; // Reverse mode derivative value
};

template<class MatType>
class A2DMat {
public:
 A2DMat( MatType& A, MatType& Ab, MatType& Ap, MatType& Ah ) : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  MatType& value(){ return A; }
  const MatType& value() const { return A; }

  MatType& bvalue(){ return Ab; }
  const MatType& bvalue() const { return Ab; }

  MatType& pvalue(){ return Ap; }
  const MatType& pvalue() const { return Ap; }

  MatType& hvalue(){ return Ah; }
  const MatType& hvalue() const { return Ah; }

  MatType& A; // Matrix
  MatType& Ab; // Reverse mode derivative value
  MatType& Ap; // Projected second derivative value
  MatType& Ah; // Reverse mode second derivative
};

} // namespace A2D

#endif // A2D_TYPES_H
