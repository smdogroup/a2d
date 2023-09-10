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

}  // namespace A2D

#endif  // A2D_BINARY_OPS_H