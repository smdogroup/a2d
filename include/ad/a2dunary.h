#ifndef A2D_UNARY_OPS_H
#define A2D_UNARY_OPS_H

#include "a2ddefs.h"

namespace A2D {

/*
  Definitions for memory-less forward and reverse-mode first-order AD

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  DERIVBODY: Body of the derivative evaluation
*/
#define A2D_1ST_UNARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, DERIVBODY)            \
                                                                               \
  template <class A, class Ta, class T, bool CA>                               \
  class OBJNAME : public ADExpr<OBJNAME<A, Ta, T, CA>, T> {                    \
   public:                                                                     \
    using expr_t = typename std::conditional<CA, const ADExpr<A, Ta>,          \
                                             ADExpr<A, Ta>>::type;             \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0) : a(a0.self()), val(0.0), bval(0.0) {} \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      val = (FUNCBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void forward() {                                           \
      a.forward();                                                             \
      bval = (DERIVBODY)*a.bvalue();                                           \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (DERIVBODY)*bval;                                          \
      a.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                            \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
                                                                               \
   private:                                                                    \
    A_t a;                                                                     \
    T val, bval;                                                               \
  };                                                                           \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(const ADExpr<A, Ta>& a) {                      \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, true>(a);                                         \
  }                                                                            \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(ADExpr<A, Ta>& a) {                            \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, false>(a);                                        \
  }

// A2D_1ST_UNARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, DERIVBODY)
A2D_1ST_UNARY_BASIC(ExpExpr, exp, std::exp(a.value()), val)
A2D_1ST_UNARY_BASIC(UnaryPos, operator+, a.value(), T(1.0))
A2D_1ST_UNARY_BASIC(UnaryNeg, operator-, -a.value(), -T(1.0))

/*
  Definitions for second-derivative implementations with zero second derivative

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  DERIVBODY: Body of the derivative evaluation
*/
#define A2D_2ND_UNARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, DERIVBODY)            \
                                                                               \
  template <class A, class Ta, class T, bool CA>                               \
  class OBJNAME : public A2DExpr<OBJNAME<A, Ta, T, CA>, T> {                   \
   public:                                                                     \
    using expr_t = typename std::conditional<CA, const A2DExpr<A, Ta>,         \
                                             A2DExpr<A, Ta>>::type;            \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0)                                        \
        : a(a0.self()), val(0.0), bval(0.0), pval(0.0), hval(0.0), tmp(0.0) {} \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      val = (FUNCBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (DERIVBODY)*bval;                                          \
      a.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION void hforward() {                                          \
      a.hforward();                                                            \
      pval = (DERIVBODY)*a.pvalue();                                           \
    }                                                                          \
    KOKKOS_FUNCTION void hreverse() {                                          \
      a.hvalue() += (DERIVBODY)*hval;                                          \
      a.hreverse();                                                            \
    }                                                                          \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                            \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                            \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                               \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }                   \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                               \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }                   \
                                                                               \
   private:                                                                    \
    A_t a;                                                                     \
    T val, bval, pval, hval, tmp;                                              \
  };                                                                           \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(const A2DExpr<A, Ta>& a) {                     \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, true>(a);                                         \
  }                                                                            \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(A2DExpr<A, Ta>& a) {                           \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, false>(a);                                        \
  }

A2D_2ND_UNARY_BASIC(UnaryPos2, operator+, a.value(), T(1.0))
A2D_2ND_UNARY_BASIC(UnaryNeg2, operator-, -a.value(), -T(1.0))

/*
  Definitions for forward and reverse-mode first-order AD with temporary

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  TEMPBODY: Body of the temporary variable calculation (computed in eval)
  DERIVBODY: Body of the derivative evaluation
*/
#define A2D_1ST_UNARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, DERIVBODY) \
                                                                        \
  template <class A, class Ta, class T, bool CA>                        \
  class OBJNAME : public ADExpr<OBJNAME<A, Ta, T, CA>, T> {             \
   public:                                                              \
    using expr_t = typename std::conditional<CA, const ADExpr<A, Ta>,   \
                                             ADExpr<A, Ta>>::type;      \
    using A_t = typename std::conditional<CA, A, A&>::type;             \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0)                                 \
        : a(a0.self()), val(0.0), bval(0.0), tmp(0.0) {}                \
    KOKKOS_FUNCTION void eval() {                                       \
      a.eval();                                                         \
      val = (FUNCBODY);                                                 \
      tmp = (TEMPBODY);                                                 \
    }                                                                   \
    KOKKOS_FUNCTION void forward() {                                    \
      a.forward();                                                      \
      bval = (DERIVBODY)*a.bvalue();                                    \
    }                                                                   \
    KOKKOS_FUNCTION void reverse() {                                    \
      a.bvalue() += (DERIVBODY)*bval;                                   \
      a.reverse();                                                      \
    }                                                                   \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                     \
    KOKKOS_FUNCTION T& value() { return val; }                          \
    KOKKOS_FUNCTION const T& value() const { return val; }              \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                        \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }            \
                                                                        \
   private:                                                             \
    A_t a;                                                              \
    T val, bval, tmp;                                                   \
  };                                                                    \
  template <class A, class Ta>                                          \
  KOKKOS_FUNCTION auto OPERNAME(const ADExpr<A, Ta>& a) {               \
    using T = typename remove_const_and_refs<Ta>::type;                 \
    return OBJNAME<A, Ta, T, true>(a);                                  \
  }                                                                     \
  template <class A, class Ta>                                          \
  KOKKOS_FUNCTION auto OPERNAME(ADExpr<A, Ta>& a) {                     \
    using T = typename remove_const_and_refs<Ta>::type;                 \
    return OBJNAME<A, Ta, T, false>(a);                                 \
  }

// A2D_1ST_UNARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, DERIVBODY)
A2D_1ST_UNARY(SinExpr, sin, std::sin(a.value()), std::cos(a.value()), tmp)
A2D_1ST_UNARY(CosExpr, cos, std::cos(a.value()), std::sin(a.value()), -tmp)
A2D_1ST_UNARY(SqrtExpr, sqrt, std::sqrt(a.value()), 1.0 / val, 0.5 * tmp)
A2D_1ST_UNARY(LogExpr, log, std::log(a.value()), 1.0 / a.value(), tmp)

/*
  Definitions for forward and reverse-mode first-order AD with temporary

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  TEMPBODY: Body of the temporary variable calculation (computed in eval)
  DERIVBODY: Body of the derivative evaluation
  DERIV2BODY: Body of the second derivative evaluation
*/
#define A2D_2ND_UNARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, DERIVBODY,        \
                      DERIV2BODY)                                              \
                                                                               \
  template <class A, class Ta, class T, bool CA>                               \
  class OBJNAME : public A2DExpr<OBJNAME<A, Ta, T, CA>, T> {                   \
   public:                                                                     \
    using expr_t = typename std::conditional<CA, const A2DExpr<A, Ta>,         \
                                             A2DExpr<A, Ta>>::type;            \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0)                                        \
        : a(a0.self()), val(0.0), bval(0.0), pval(0.0), hval(0.0), tmp(0.0) {} \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      val = (FUNCBODY);                                                        \
      tmp = (TEMPBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (DERIVBODY)*bval;                                          \
      a.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION void hforward() {                                          \
      a.hforward();                                                            \
      pval = (DERIVBODY)*a.pvalue();                                           \
    }                                                                          \
    KOKKOS_FUNCTION void hreverse() {                                          \
      a.hvalue() += (DERIVBODY)*hval + (DERIV2BODY)*bval * a.pvalue();         \
      a.hreverse();                                                            \
    }                                                                          \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                            \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                            \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                               \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }                   \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                               \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }                   \
                                                                               \
   private:                                                                    \
    A_t a;                                                                     \
    T val, bval, pval, hval, tmp;                                              \
  };                                                                           \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(const A2DExpr<A, Ta>& a) {                     \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, true>(a);                                         \
  }                                                                            \
  template <class A, class Ta>                                                 \
  KOKKOS_FUNCTION auto OPERNAME(A2DExpr<A, Ta>& a) {                           \
    using T = typename remove_const_and_refs<Ta>::type;                        \
    return OBJNAME<A, Ta, T, false>(a);                                        \
  }

// A2D_1ST_UNARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, DERIVBODY)
A2D_2ND_UNARY(ExpExpr2, exp, std::exp(a.value()), val, tmp, tmp)
A2D_2ND_UNARY(SinExpr2, sin, std::sin(a.value()), std::cos(a.value()), tmp,
              -val)
A2D_2ND_UNARY(CosExpr2, cos, std::cos(a.value()), std::sin(a.value()), -tmp,
              -val)
A2D_2ND_UNARY(SqrtExpr2, sqrt, std::sqrt(a.value()), 1.0 / val, 0.5 * tmp,
              -0.25 * tmp * tmp * tmp)
A2D_2ND_UNARY(LogExpr2, log, std::log(a.value()), 1.0 / a.value(), tmp,
              -tmp* tmp)

}  // namespace A2D

#endif  // A2D_UNARY_OPS_H