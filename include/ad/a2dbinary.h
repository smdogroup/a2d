#ifndef A2D_BINARY_OPS_H
#define A2D_BINARY_OPS_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

#define A2D_1ST_BINARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, FORWARDBODY,         \
                             AREVBODY, BREVBODY)                               \
                                                                               \
  template <class A, class B, class T, bool CA, bool CB>                       \
  class OBJNAME : public ADExpr<OBJNAME<A, B, T, CA, CB>, T> {                 \
   public:                                                                     \
    using Aexpr_t =                                                            \
        typename std::conditional<CA, const ADExpr<A, T>, ADExpr<A, T>>::type; \
    using Bexpr_t =                                                            \
        typename std::conditional<CB, const ADExpr<B, T>, ADExpr<B, T>>::type; \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    using B_t = typename std::conditional<CB, B, B&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(Aexpr_t& a0, Bexpr_t& b0)                          \
        : a(a0.self()), b(b0.self()), val(0.0), bval(0.0) {}                   \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      b.eval();                                                                \
      val = (FUNCBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void forward() {                                           \
      a.forward();                                                             \
      b.forward();                                                             \
      bval = (FORWARDBODY);                                                    \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (AREVBODY);                                                \
      b.bvalue() += (BREVBODY);                                                \
      a.reverse();                                                             \
      b.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
                                                                               \
   private:                                                                    \
    A_t a;                                                                     \
    B_t b;                                                                     \
    T val, bval;                                                               \
  };                                                                           \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(const ADExpr<A, T>& a, const ADExpr<B, T>& b) {         \
    return OBJNAME<A, B, T, true, true>(a, b);                                 \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(const ADExpr<A, T>& a, ADExpr<B, T>& b) {               \
    return OBJNAME<A, B, T, true, false>(a, b);                                \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(ADExpr<A, T>& a, const ADExpr<B, T>& b) {               \
    return OBJNAME<A, B, T, false, true>(a, b);                                \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(ADExpr<A, T>& a, ADExpr<B, T>& b) {                     \
    return OBJNAME<A, B, T, false, false>(a, b);                               \
  }

// A2D_1ST_BINARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, FORWARDBODY, AREVBODY,
//                      BREVBODY)
A2D_1ST_BINARY_BASIC(AddExpr, operator+, a.value() + b.value(),
                     a.bvalue() + b.bvalue(), bval, bval)
A2D_1ST_BINARY_BASIC(SubExpr, operator-, a.value() - b.value(),
                     a.bvalue() - b.bvalue(), bval, -bval)
A2D_1ST_BINARY_BASIC(MultExpr, operator*, a.value() * b.value(),
                     a.bvalue() * b.value() + a.value() * b.bvalue(),
                     b.value() * bval, a.value() * bval)

#define A2D_1ST_BINARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, FORWARDBODY,     \
                       AREVBODY, BREVBODY)                                     \
                                                                               \
  template <class A, class B, class T, bool CA, bool CB>                       \
  class OBJNAME : public ADExpr<OBJNAME<A, B, T, CA, CB>, T> {                 \
   public:                                                                     \
    using Aexpr_t =                                                            \
        typename std::conditional<CA, const ADExpr<A, T>, ADExpr<A, T>>::type; \
    using Bexpr_t =                                                            \
        typename std::conditional<CB, const ADExpr<B, T>, ADExpr<B, T>>::type; \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    using B_t = typename std::conditional<CB, B, B&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(Aexpr_t& a0, Bexpr_t& b0)                          \
        : a(a0.self()), b(b0.self()), tmp(0.0), val(0.0), bval(0.0) {}         \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      b.eval();                                                                \
      val = (FUNCBODY);                                                        \
      tmp = (TEMPBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void forward() {                                           \
      a.forward();                                                             \
      b.forward();                                                             \
      bval = (FORWARDBODY);                                                    \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (AREVBODY);                                                \
      b.bvalue() += (BREVBODY);                                                \
      a.reverse();                                                             \
      b.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
                                                                               \
   private:                                                                    \
    A_t a;                                                                     \
    B_t b;                                                                     \
    T tmp, val, bval;                                                          \
  };                                                                           \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(const ADExpr<A, T>& a, const ADExpr<B, T>& b) {         \
    return OBJNAME<A, B, T, true, true>(a, b);                                 \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(const ADExpr<A, T>& a, ADExpr<B, T>& b) {               \
    return OBJNAME<A, B, T, true, false>(a, b);                                \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(ADExpr<A, T>& a, const ADExpr<B, T>& b) {               \
    return OBJNAME<A, B, T, false, true>(a, b);                                \
  }                                                                            \
  template <class A, class B, class T>                                         \
  inline auto OPERNAME(ADExpr<A, T>& a, ADExpr<B, T>& b) {                     \
    return OBJNAME<A, B, T, false, false>(a, b);                               \
  }

// A2D_1ST_BINARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, FORWARDBODY,
//                AREVBODY, BREVBODY)
A2D_1ST_BINARY(Divide, operator/, a.value() / b.value(), T(1.0) / b.value(),
               tmp*(a.bvalue() - tmp * a.value() * b.bvalue()), tmp* bval,
               -tmp* tmp* a.value() * bval)

#define A2D_2ND_BINARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, AREVBODY, BREVBODY, \
                             HFORWARDBODY, HAREVBODY, HBREVBODY)              \
                                                                              \
  template <class A, class B, class T, bool CA, bool CB>                      \
  class OBJNAME : public A2DExpr<OBJNAME<A, B, T, CA, CB>, T> {               \
   public:                                                                    \
    using Aexpr_t = typename std::conditional<CA, const A2DExpr<A, T>,        \
                                              A2DExpr<A, T>>::type;           \
    using Bexpr_t = typename std::conditional<CB, const A2DExpr<B, T>,        \
                                              A2DExpr<B, T>>::type;           \
    using A_t = typename std::conditional<CA, A, A&>::type;                   \
    using B_t = typename std::conditional<CB, B, B&>::type;                   \
    KOKKOS_FUNCTION OBJNAME(Aexpr_t& a0, Bexpr_t& b0)                         \
        : a(a0.self()),                                                       \
          b(b0.self()),                                                       \
          val(0.0),                                                           \
          bval(0.0),                                                          \
          pval(0.0),                                                          \
          hval(0.0) {}                                                        \
    KOKKOS_FUNCTION void eval() {                                             \
      a.eval();                                                               \
      b.eval();                                                               \
      val = (FUNCBODY);                                                       \
    }                                                                         \
    KOKKOS_FUNCTION void reverse() {                                          \
      a.bvalue() += (AREVBODY);                                               \
      b.bvalue() += (BREVBODY);                                               \
      a.reverse();                                                            \
      b.reverse();                                                            \
    }                                                                         \
    KOKKOS_FUNCTION void hforward() {                                         \
      a.hforward();                                                           \
      b.hforward();                                                           \
      pval = (HFORWARDBODY);                                                  \
    }                                                                         \
    KOKKOS_FUNCTION void hreverse() {                                         \
      a.hvalue() += (HAREVBODY);                                              \
      b.hvalue() += (HBREVBODY);                                              \
      a.hreverse();                                                           \
      b.hreverse();                                                           \
    }                                                                         \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                           \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                           \
    KOKKOS_FUNCTION T& value() { return val; }                                \
    KOKKOS_FUNCTION const T& value() const { return val; }                    \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                              \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                  \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                              \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }                  \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                              \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }                  \
                                                                              \
   private:                                                                   \
    A_t a;                                                                    \
    B_t b;                                                                    \
    T val, bval, pval, hval;                                                  \
  };                                                                          \
  template <class A, class B, class T>                                        \
  inline auto OPERNAME(const A2DExpr<A, T>& a, const A2DExpr<B, T>& b) {      \
    return OBJNAME<A, B, T, true, true>(a, b);                                \
  }                                                                           \
  template <class A, class B, class T>                                        \
  inline auto OPERNAME(const A2DExpr<A, T>& a, A2DExpr<B, T>& b) {            \
    return OBJNAME<A, B, T, true, false>(a, b);                               \
  }                                                                           \
  template <class A, class B, class T>                                        \
  inline auto OPERNAME(A2DExpr<A, T>& a, const A2DExpr<B, T>& b) {            \
    return OBJNAME<A, B, T, false, true>(a, b);                               \
  }                                                                           \
  template <class A, class B, class T>                                        \
  inline auto OPERNAME(A2DExpr<A, T>& a, A2DExpr<B, T>& b) {                  \
    return OBJNAME<A, B, T, false, false>(a, b);                              \
  }

// A2D_2ND_BINARY_BASIC(OBJNAME, OPERNAME, FUNCBODY, AREVBODY, BREVBODY,
//                      HFORWARDBODY, HAREVBODY, HBREVBODY)
A2D_2ND_BINARY_BASIC(AddExpr2, operator+, a.value() + b.value(), bval, bval,
                     a.pvalue() + b.pvalue(), hval, hval)
A2D_2ND_BINARY_BASIC(SubExpr2, operator-, a.value() - b.value(), bval, -bval,
                     a.pvalue() - b.pvalue(), hval, -hval)
A2D_2ND_BINARY_BASIC(MultExpr2, operator*, a.value() * b.value(),
                     b.value() * bval, a.value() * bval,
                     a.pvalue() * b.value() + a.value() * b.pvalue(),
                     b.value() * hval + bval * b.pvalue(),
                     a.value() * hval + bval * a.pvalue())

#define A2D_2ND_BINARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, AREVBODY,  \
                       BREVBODY, HFORWARDBODY, HAREVBODY, HBREVBODY)     \
                                                                         \
  template <class A, class B, class T, bool CA, bool CB>                 \
  class OBJNAME : public A2DExpr<OBJNAME<A, B, T, CA, CB>, T> {          \
   public:                                                               \
    using Aexpr_t = typename std::conditional<CA, const A2DExpr<A, T>,   \
                                              A2DExpr<A, T>>::type;      \
    using Bexpr_t = typename std::conditional<CB, const A2DExpr<B, T>,   \
                                              A2DExpr<B, T>>::type;      \
    using A_t = typename std::conditional<CA, A, A&>::type;              \
    using B_t = typename std::conditional<CB, B, B&>::type;              \
    KOKKOS_FUNCTION OBJNAME(Aexpr_t& a0, Bexpr_t& b0)                    \
        : a(a0.self()),                                                  \
          b(b0.self()),                                                  \
          tmp(0.0),                                                      \
          val(0.0),                                                      \
          bval(0.0),                                                     \
          pval(0.0),                                                     \
          hval(0.0) {}                                                   \
    KOKKOS_FUNCTION void eval() {                                        \
      a.eval();                                                          \
      b.eval();                                                          \
      val = (FUNCBODY);                                                  \
      tmp = (TEMPBODY);                                                  \
    }                                                                    \
    KOKKOS_FUNCTION void reverse() {                                     \
      a.bvalue() += (AREVBODY);                                          \
      b.bvalue() += (BREVBODY);                                          \
      a.reverse();                                                       \
      b.reverse();                                                       \
    }                                                                    \
    KOKKOS_FUNCTION void hforward() {                                    \
      a.hforward();                                                      \
      b.hforward();                                                      \
      pval = (HFORWARDBODY);                                             \
    }                                                                    \
    KOKKOS_FUNCTION void hreverse() {                                    \
      a.hvalue() += (HAREVBODY);                                         \
      b.hvalue() += (HBREVBODY);                                         \
      a.hreverse();                                                      \
      b.hreverse();                                                      \
    }                                                                    \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                      \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                      \
    KOKKOS_FUNCTION T& value() { return val; }                           \
    KOKKOS_FUNCTION const T& value() const { return val; }               \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                         \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }             \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                         \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }             \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                         \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }             \
                                                                         \
   private:                                                              \
    A_t a;                                                               \
    B_t b;                                                               \
    T tmp, val, bval, pval, hval;                                        \
  };                                                                     \
  template <class A, class B, class T>                                   \
  inline auto OPERNAME(const A2DExpr<A, T>& a, const A2DExpr<B, T>& b) { \
    return OBJNAME<A, B, T, true, true>(a, b);                           \
  }                                                                      \
  template <class A, class B, class T>                                   \
  inline auto OPERNAME(const A2DExpr<A, T>& a, A2DExpr<B, T>& b) {       \
    return OBJNAME<A, B, T, true, false>(a, b);                          \
  }                                                                      \
  template <class A, class B, class T>                                   \
  inline auto OPERNAME(A2DExpr<A, T>& a, const A2DExpr<B, T>& b) {       \
    return OBJNAME<A, B, T, false, true>(a, b);                          \
  }                                                                      \
  template <class A, class B, class T>                                   \
  inline auto OPERNAME(A2DExpr<A, T>& a, A2DExpr<B, T>& b) {             \
    return OBJNAME<A, B, T, false, false>(a, b);                         \
  }

// A2D_2ND_BINARY(OBJNAME, OPERNAME, FUNCBODY, TEMPBODY, AREVBODY,
//                BREVBODY, HFORWARDBODY, HAREVBODY, HBREVBODY)
A2D_2ND_BINARY(Divide2, operator/, a.value() / b.value(), T(1.0) / b.value(),
               tmp* bval, -tmp* tmp* a.value() * bval,
               tmp*(a.pvalue() - tmp * a.value() * b.pvalue()),
               tmp*(hval - tmp * bval * b.pvalue()),
               tmp* tmp * (2.0 * tmp * a.value() * bval * b.pvalue() -
                           a.value() * hval - bval * a.pvalue()))

/*
  Definitions for memory-less forward and reverse-mode first-order AD

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  FORWARD: Body of the forward derivative
  REVERSE: Body of the reverse derivative
*/
#define A2D_1ST_BINARY_LEFT_BASIC(OBJNAME, OPERNAME, FUNCBODY, FORWARD,        \
                                  REVERSE)                                     \
                                                                               \
  template <class A, class B, class T, bool CA>                                \
  class OBJNAME : public ADExpr<OBJNAME<A, B, T, CA>, T> {                     \
   public:                                                                     \
    using expr_t =                                                             \
        typename std::conditional<CA, const ADExpr<A, T>, ADExpr<A, T>>::type; \
    using A_t = typename std::conditional<CA, A, A&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0, const B& b)                            \
        : a(a0.self()), b(b), val(0.0), bval(0.0) {}                           \
    KOKKOS_FUNCTION void eval() {                                              \
      a.eval();                                                                \
      val = (FUNCBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void forward() {                                           \
      a.forward();                                                             \
      bval = (FORWARD);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      a.bvalue() += (REVERSE);                                                 \
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
    const B& b;                                                                \
    T val, bval;                                                               \
  };                                                                           \
  template <class A, class B, class T,                                         \
            std::enable_if_t<is_scalar_type<B>::value, bool> = true>           \
  KOKKOS_FUNCTION auto OPERNAME(const ADExpr<A, T>& a, const B& b) {           \
    return OBJNAME<A, B, T, true>(a, b);                                       \
  }                                                                            \
  template <class A, class B, class T,                                         \
            std::enable_if_t<is_scalar_type<B>::value, bool> = true>           \
  KOKKOS_FUNCTION auto OPERNAME(ADExpr<A, T>& a, const B& b) {                 \
    return OBJNAME<A, B, T, false>(a, b);                                      \
  }

// A2D_1ST_BINARY_LEFT_BASIC(OBJNAME, OPERNAME, FUNCBODY, DERIVBODY)
A2D_1ST_BINARY_LEFT_BASIC(LAddExpr, operator+, a.value() + b, a.bvalue(), bval)
A2D_1ST_BINARY_LEFT_BASIC(LSubExpr, operator-, a.value() - b, a.bvalue(), bval)
A2D_1ST_BINARY_LEFT_BASIC(LMultExpr, operator*, a.value() * b, a.bvalue() * b,
                          b* bval)

/*
  Definitions for memory-less forward and reverse-mode first-order AD

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  REVERSE: Body of the reverse derivative
  HFORWARD: Body of the second-order forward derivative
  HREVERSE: Body of the seond-order reverse derivative
*/
#define A2D_2ND_BINARY_LEFT_BASIC(OBJNAME, OPERNAME, FUNCBODY, REVERSE,    \
                                  HFORWARD, HREVERSE)                      \
                                                                           \
  template <class A, class B, class T, bool CA>                            \
  class OBJNAME : public A2DExpr<OBJNAME<A, B, T, CA>, T> {                \
   public:                                                                 \
    using expr_t = typename std::conditional<CA, const A2DExpr<A, T>,      \
                                             A2DExpr<A, T>>::type;         \
    using A_t = typename std::conditional<CA, A, A&>::type;                \
    KOKKOS_FUNCTION OBJNAME(expr_t& a0, const B& b)                        \
        : a(a0.self()), b(b), val(0.0), bval(0.0), pval(0.0), hval(0.0) {} \
    KOKKOS_FUNCTION void eval() {                                          \
      a.eval();                                                            \
      val = (FUNCBODY);                                                    \
    }                                                                      \
    KOKKOS_FUNCTION void reverse() {                                       \
      a.bvalue() += (REVERSE);                                             \
      a.reverse();                                                         \
    }                                                                      \
    KOKKOS_FUNCTION void hforward() {                                      \
      a.hforward();                                                        \
      pval = (HFORWARD);                                                   \
    }                                                                      \
    KOKKOS_FUNCTION void hreverse() {                                      \
      a.hvalue() += (HREVERSE);                                            \
      a.hreverse();                                                        \
    }                                                                      \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                        \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                        \
    KOKKOS_FUNCTION T& value() { return val; }                             \
    KOKKOS_FUNCTION const T& value() const { return val; }                 \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                           \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }               \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                           \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }               \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                           \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }               \
                                                                           \
   private:                                                                \
    A_t a;                                                                 \
    const B& b;                                                            \
    T val, bval, pval, hval;                                               \
  };                                                                       \
  template <class A, class B, class T,                                     \
            std::enable_if_t<is_scalar_type<B>::value, bool> = true>       \
  KOKKOS_FUNCTION auto OPERNAME(const A2DExpr<A, T>& a, const B& b) {      \
    return OBJNAME<A, B, T, true>(a, b);                                   \
  }                                                                        \
  template <class A, class B, class T,                                     \
            std::enable_if_t<is_scalar_type<B>::value, bool> = true>       \
  KOKKOS_FUNCTION auto OPERNAME(A2DExpr<A, T>& a, const B& b) {            \
    return OBJNAME<A, B, T, false>(a, b);                                  \
  }

// A2D_2ND_BINARY_LEFT_BASIC(OBJNAME, OPERNAME, FUNCBODY, REVERSE,
//                           HFORWARD, HREVERSE)
A2D_2ND_BINARY_LEFT_BASIC(LAddExpr2, operator+, a.value() + b, bval, a.pvalue(),
                          hval)
A2D_2ND_BINARY_LEFT_BASIC(LSubExpr2, operator-, a.value() - b, bval, a.pvalue(),
                          hval)
A2D_2ND_BINARY_LEFT_BASIC(LMultExpr2, operator*, a.value() * b, b* bval,
                          b* a.pvalue(), b* hval)

/*
  Definitions for memory-less forward and reverse-mode first-order AD

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  FORWARD: Body of the forward derivative
  REVERSE: Body of the reverse derivative
*/
#define A2D_1ST_BINARY_RIGHT_BASIC(OBJNAME, OPERNAME, FUNCBODY, FORWARD,       \
                                   REVERSE)                                    \
                                                                               \
  template <class A, class B, class T, bool CB>                                \
  class OBJNAME : public ADExpr<OBJNAME<A, B, T, CB>, T> {                     \
   public:                                                                     \
    using expr_t =                                                             \
        typename std::conditional<CB, const ADExpr<B, T>, ADExpr<B, T>>::type; \
    using B_t = typename std::conditional<CB, B, B&>::type;                    \
    KOKKOS_FUNCTION OBJNAME(const A& a, expr_t& b0)                            \
        : a(a), b(b0.self()), val(0.0), bval(0.0) {}                           \
    KOKKOS_FUNCTION void eval() {                                              \
      b.eval();                                                                \
      val = (FUNCBODY);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void forward() {                                           \
      b.forward();                                                             \
      bval = (FORWARD);                                                        \
    }                                                                          \
    KOKKOS_FUNCTION void reverse() {                                           \
      b.bvalue() += (REVERSE);                                                 \
      b.reverse();                                                             \
    }                                                                          \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                            \
    KOKKOS_FUNCTION T& value() { return val; }                                 \
    KOKKOS_FUNCTION const T& value() const { return val; }                     \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                               \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }                   \
                                                                               \
   private:                                                                    \
    const A& a;                                                                \
    B_t b;                                                                     \
    T val, bval;                                                               \
  };                                                                           \
  template <class A, class B, class T,                                         \
            std::enable_if_t<is_scalar_type<A>::value, bool> = true>           \
  KOKKOS_FUNCTION auto OPERNAME(const A& a, const ADExpr<B, T>& b) {           \
    return OBJNAME<A, B, T, true>(a, b);                                       \
  }                                                                            \
  template <class A, class B, class T,                                         \
            std::enable_if_t<is_scalar_type<A>::value, bool> = true>           \
  KOKKOS_FUNCTION auto OPERNAME(const A& a, ADExpr<B, T>& b) {                 \
    return OBJNAME<A, B, T, false>(a, b);                                      \
  }

// A2D_1ST_BINARY_RIGHT_BASIC(OBJNAME, OPERNAME, FUNCBODY, DERIVBODY)
A2D_1ST_BINARY_RIGHT_BASIC(RAddExpr, operator+, a + b.value(), b.bvalue(), bval)
A2D_1ST_BINARY_RIGHT_BASIC(RSubExpr, operator-, a - b.value(), -b.bvalue(),
                           -bval)
A2D_1ST_BINARY_RIGHT_BASIC(RMultExpr, operator*, a* b.value(), a* b.bvalue(),
                           a* bval)

/*
  Definitions for memory-less forward and reverse-mode first-order AD

  OBJNAME: Expression template name
  OPERNAME: Name of the operator
  FUNCBODY: Body of the function evaluation
  FORWARD: Body of the forward derivative
  REVERSE: Body of the reverse derivative
*/
#define A2D_2ND_BINARY_RIGHT_BASIC(OBJNAME, OPERNAME, FUNCBODY, REVERSE,   \
                                   HFORWARD, HREVERSE)                     \
                                                                           \
  template <class A, class B, class T, bool CB>                            \
  class OBJNAME : public A2DExpr<OBJNAME<A, B, T, CB>, T> {                \
   public:                                                                 \
    using expr_t = typename std::conditional<CB, const A2DExpr<B, T>,      \
                                             A2DExpr<B, T>>::type;         \
    using B_t = typename std::conditional<CB, B, B&>::type;                \
    KOKKOS_FUNCTION OBJNAME(const A& a, expr_t& b0)                        \
        : a(a), b(b0.self()), val(0.0), bval(0.0), pval(0.0), hval(0.0) {} \
    KOKKOS_FUNCTION void eval() {                                          \
      b.eval();                                                            \
      val = (FUNCBODY);                                                    \
    }                                                                      \
    KOKKOS_FUNCTION void reverse() {                                       \
      b.bvalue() += (REVERSE);                                             \
      b.reverse();                                                         \
    }                                                                      \
    KOKKOS_FUNCTION void hforward() {                                      \
      b.hforward();                                                        \
      pval = (HFORWARD);                                                   \
    }                                                                      \
    KOKKOS_FUNCTION void hreverse() {                                      \
      b.hvalue() += (HREVERSE);                                            \
      b.hreverse();                                                        \
    }                                                                      \
    KOKKOS_FUNCTION void bzero() { bval = T(0.0); }                        \
    KOKKOS_FUNCTION void hzero() { hval = T(0.0); }                        \
    KOKKOS_FUNCTION T& value() { return val; }                             \
    KOKKOS_FUNCTION const T& value() const { return val; }                 \
    KOKKOS_FUNCTION T& bvalue() { return bval; }                           \
    KOKKOS_FUNCTION const T& bvalue() const { return bval; }               \
    KOKKOS_FUNCTION T& pvalue() { return pval; }                           \
    KOKKOS_FUNCTION const T& pvalue() const { return pval; }               \
    KOKKOS_FUNCTION T& hvalue() { return hval; }                           \
    KOKKOS_FUNCTION const T& hvalue() const { return hval; }               \
                                                                           \
   private:                                                                \
    const A& a;                                                            \
    B_t b;                                                                 \
    T val, bval, pval, hval;                                               \
  };                                                                       \
  template <class A, class B, class T,                                     \
            std::enable_if_t<is_scalar_type<A>::value, bool> = true>       \
  KOKKOS_FUNCTION auto OPERNAME(const A& a, const A2DExpr<B, T>& b) {      \
    return OBJNAME<A, B, T, true>(a, b);                                   \
  }                                                                        \
  template <class A, class B, class T,                                     \
            std::enable_if_t<is_scalar_type<A>::value, bool> = true>       \
  KOKKOS_FUNCTION auto OPERNAME(const A& a, A2DExpr<B, T>& b) {            \
    return OBJNAME<A, B, T, false>(a, b);                                  \
  }

// A2D_2ND_BINARY_RIGHT_BASIC(OBJNAME, OPERNAME, FUNCBODY, REVERSE,
//                            HFORWARD, HREVERSE)
A2D_2ND_BINARY_RIGHT_BASIC(RAddExpr2, operator+, a + b.value(), bval,
                           b.pvalue(), hval)
A2D_2ND_BINARY_RIGHT_BASIC(RSubExpr2, operator-, a - b.value(), -bval,
                           -b.pvalue(), -hval)
A2D_2ND_BINARY_RIGHT_BASIC(RMultExpr2, operator*, a* b.value(), a* bval,
                           a* b.pvalue(), a* hval)

}  // namespace A2D

#endif  // A2D_BINARY_OPS_H