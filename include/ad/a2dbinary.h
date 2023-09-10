#ifndef A2D_BINARY_OPS_H
#define A2D_BINARY_OPS_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

#define A2D_1ST_BINARY(OBJNAME, OPERNAME, FUNCBODY, FORWARDBODY, AREVBODY,     \
                       BREVBODY)                                               \
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

// A2D_1ST_BINARY(OBJNAME, OPERNAME, FUNCBODY, FORWARDBODY, AREVBODY, BREVBODY)
// A2D_1ST_BINARY(AddExpr, operator+, a.value() + b.value(),
//                a.bvalue() + b.bvalue(), bval, bval)
// A2D_1ST_BINARY(SubExpr, operator-, a.value() - b.value(),
//                a.bvalue() - b.bvalue(), bval, -bval)
// A2D_1ST_BINARY(MultExpr, operator*, a.value() * b.value(),
//                a.bvalue() * b.value() + a.value() * b.bvalue(),
//                b.value() * bval, a.value() * bval)

template <class A, class B, class T, bool CA, bool CB>
class AddExpr : public ADExpr<AddExpr<A, B, T, CA, CB>, T> {
 public:
  using Aexpr_t =
      typename std::conditional<CA, const ADExpr<A, T>, ADExpr<A, T>>::type;
  using Bexpr_t =
      typename std::conditional<CB, const ADExpr<B, T>, ADExpr<B, T>>::type;
  using A_t = typename std::conditional<CA, A, A&>::type;
  using B_t = typename std::conditional<CB, B, B&>::type;
  KOKKOS_FUNCTION AddExpr(Aexpr_t& a0, Bexpr_t& b0)
      : a(a0.self()), b(b0.self()), val(0.0), bval(0.0) {}
  KOKKOS_FUNCTION void eval() {
    a.eval();
    b.eval();
    val = (a.value() + b.value());
  }
  KOKKOS_FUNCTION void forward() {
    a.forward();
    b.forward();
    bval = (a.bvalue() + b.bvalue());
  }
  KOKKOS_FUNCTION void reverse() {
    a.bvalue() += (bval);
    b.bvalue() += (bval);
    a.reverse();
    b.reverse();
  }
  KOKKOS_FUNCTION T& value() { return val; }
  KOKKOS_FUNCTION const T& value() const { return val; }
  KOKKOS_FUNCTION T& bvalue() { return bval; }
  KOKKOS_FUNCTION const T& bvalue() const { return bval; }

 private:
  A_t a;
  B_t b;
  T val, bval;
};
template <class A, class B, class T>
inline auto operator+(const ADExpr<A, T>& a, const ADExpr<B, T>& b) {
  return AddExpr<A, B, T, true, true>(a, b);
}
template <class A, class B, class T>
inline auto operator+(const ADExpr<A, T>& a, ADExpr<B, T>& b) {
  return AddExpr<A, B, T, true, false>(a, b);
}
template <class A, class B, class T>
inline auto operator+(ADExpr<A, T>& a, const ADExpr<B, T>& b) {
  return AddExpr<A, B, T, false, true>(a, b);
}
template <class A, class B, class T>
inline auto operator+(ADExpr<A, T>& a, ADExpr<B, T>& b) {
  return AddExpr<A, B, T, false, false>(a, b);
}

template <class A, class B, class T, bool CA, bool CB>
class MultExpr : public ADExpr<MultExpr<A, B, T, CA, CB>, T> {
 public:
  using Aexpr_t =
      typename std::conditional<CA, const ADExpr<A, T>, ADExpr<A, T>>::type;
  using Bexpr_t =
      typename std::conditional<CB, const ADExpr<B, T>, ADExpr<B, T>>::type;
  using A_t = typename std::conditional<CA, A, A&>::type;
  using B_t = typename std::conditional<CB, B, B&>::type;
  KOKKOS_FUNCTION MultExpr(Aexpr_t& a0, Bexpr_t& b0)
      : a(a0.self()), b(b0.self()), val(0.0), bval(0.0) {}
  KOKKOS_FUNCTION void eval() {
    a.eval();
    b.eval();
    val = (a.value() * b.value());
  }
  KOKKOS_FUNCTION void forward() {
    a.forward();
    b.forward();
    bval = (a.bvalue() * b.value() + a.value() * b.bvalue());
  }
  KOKKOS_FUNCTION void reverse() {
    a.bvalue() += (b.value() * bval);
    b.bvalue() += (a.value() * bval);
    a.reverse();
    b.reverse();
  }
  KOKKOS_FUNCTION T& value() { return val; }
  KOKKOS_FUNCTION const T& value() const { return val; }
  KOKKOS_FUNCTION T& bvalue() { return bval; }
  KOKKOS_FUNCTION const T& bvalue() const { return bval; }

 private:
  A_t a;
  B_t b;
  T val, bval;
};
template <class A, class B, class T>
inline auto operator*(const ADExpr<A, T>& a, const ADExpr<B, T>& b) {
  return MultExpr<A, B, T, true, true>(a, b);
}
template <class A, class B, class T>
inline auto operator*(const ADExpr<A, T>& a, ADExpr<B, T>& b) {
  return MultExpr<A, B, T, true, false>(a, b);
}
template <class A, class B, class T>
inline auto operator*(ADExpr<A, T>& a, const ADExpr<B, T>& b) {
  return MultExpr<A, B, T, false, true>(a, b);
}
template <class A, class B, class T>
inline auto operator*(ADExpr<A, T>& a, ADExpr<B, T>& b) {
  return MultExpr<A, B, T, false, false>(a, b);
}

template <class A, class B, class T, bool CA, bool CB>
class MultiplyExpr2 : public A2DExpr<MultiplyExpr2<A, B, T, CA, CB>, T> {
 public:
  using Aexpr_t =
      typename std::conditional<CA, const A2DExpr<A, T>, A2DExpr<A, T>>::type;
  using Bexpr_t =
      typename std::conditional<CB, const A2DExpr<B, T>, A2DExpr<B, T>>::type;
  using A_t = typename std::conditional<CA, A, A&>::type;
  using B_t = typename std::conditional<CB, B, B&>::type;

  KOKKOS_FUNCTION MultiplyExpr2(Aexpr_t& a0, Bexpr_t& b0)
      : a(a0.self()), b(b0.self()), val(0.0), bval(0.0), pval(0.0), hval(0.0) {}

  // Evaluation and derivatives
  KOKKOS_FUNCTION void eval() {
    a.eval();
    b.eval();
    val = a.value() * b.value();
  }
  KOKKOS_FUNCTION void reverse() {
    a.bvalue() += b.value() * bval;
    b.bvalue() += a.value() * bval;
    a.reverse();
    b.reverse();
  }
  KOKKOS_FUNCTION void hforward() {
    a.hforward();
    b.hforward();
    pval = a.pvalue() * b.value() + a.value() * b.pvalue();
  }
  KOKKOS_FUNCTION void hreverse() {
    a.hvalue() += b.value() * hval + bval * b.pvalue();
    b.hvalue() += a.value() * hval + bval * a.pvalue();
    a.hreverse();
    b.hreverse();
  }
  KOKKOS_FUNCTION void bzero() { bval = T(0.0); }
  KOKKOS_FUNCTION void hzero() { hval = T(0.0); }
  KOKKOS_FUNCTION T& value() { return val; }
  KOKKOS_FUNCTION const T& value() const { return val; }
  KOKKOS_FUNCTION T& bvalue() { return bval; }
  KOKKOS_FUNCTION const T& bvalue() const { return bval; }
  KOKKOS_FUNCTION T& pvalue() { return pval; }
  KOKKOS_FUNCTION const T& pvalue() const { return pval; }
  KOKKOS_FUNCTION T& hvalue() { return hval; }
  KOKKOS_FUNCTION const T& hvalue() const { return hval; }

 private:
  A_t a;
  B_t b;
  T val, bval, pval, hval;
};

template <class A, class B, class T>
inline auto operator*(const A2DExpr<A, T>& a, const A2DExpr<B, T>& b) {
  return MultiplyExpr2<A, B, T, true, true>(a, b);
}

template <class A, class B, class T>
inline auto operator*(const A2DExpr<A, T>& a, A2DExpr<B, T>& b) {
  return MultiplyExpr2<A, B, T, true, false>(a, b);
}

template <class A, class B, class T>
inline auto operator*(A2DExpr<A, T>& a, const A2DExpr<B, T>& b) {
  return MultiplyExpr2<A, B, T, false, true>(a, b);
}

template <class A, class B, class T>
inline auto operator*(A2DExpr<A, T>& a, A2DExpr<B, T>& b) {
  return MultiplyExpr2<A, B, T, false, false>(a, b);
}

}  // namespace A2D

#endif  // A2D_BINARY_OPS_H