#ifndef A2D_SCALAR_OPS_H
#define A2D_SCALAR_OPS_H

#include "a2dbinary.h"
#include "a2denum.h"
#include "a2dobjs.h"
#include "a2dunary.h"

namespace A2D {

template <class A, class T>
class EvalExpr {
 public:
  KOKKOS_FUNCTION EvalExpr(ADExpr<A, T>&& expr, ADObj<T>& out)
      : expr(std::forward<ADExpr<A, T>>(expr)), out(out) {}

  KOKKOS_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(forder == ADorder::FIRST,
                  "EvalExpr only works for first-order AD");
    expr.forward();
    out.bvalue() = expr.bvalue();
  }
  KOKKOS_FUNCTION void reverse() {
    expr.bvalue() = out.bvalue();
    expr.reverse();
  }

 private:
  ADExpr<A, T> expr;
  ADObj<T>& out;
};

template <class A, class T>
auto Eval(ADExpr<A, T>&& expr, ADObj<T>& out) {
  return EvalExpr<A, T>(std::forward<ADExpr<A, T>>(expr), out);
}

/*
template <class A, class T>
class EvalExpr2 {
 public:
  KOKKOS_FUNCTION EvalExpr(A2DExpr<A, T> expr, A2DObj<T>& out)
      : expr(expr), out(out) {}

  KOKKOS_FUNCTION void eval() {
    expr.eval();
    out.value = expr.value();
  }
  KOKKOS_FUNCTION void reverse() {
    expr.bvalue() = out.bvalue();
    expr.reverse();
  }
  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(forder == ADorder::SECOND,
                  "EvalExpr2 only works for second-order AD");
    expr.hforward();
    out.pvalue = expr.pvalue();
  }
  KOKKOS_FUNCTION void hreverse() {
    expr.hvalue() = out.hvalue();
    expr.hreverse();
  }

 private:
  A2DExpr<A, T> expr;
  A2DObj<T> out;
};
*/

// template <class A, class T>
// auto Eval(A2DExpr<A, T> expr, A2DObj<T>& out) {
//   return EvalExpr2<A, T>(expr, out);
// }

/*
template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Log(const T a, T& b) {
  b = std::log(a);
}

template <class atype, class btype>
class LogExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION LogExpr(atype& a, btype& b) : a(a), b(b) {}

  KOKKOS_FUNCTION void eval() { get_data(b) = std::log(get_data(a)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) = GetSeed<seed>::get_data(a) / get_data(a);
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) += GetSeed<seed>::get_data(b) / get_data(a);
  }

  KOKKOS_FUNCTION void hreverse() {
    T inv = 1.0 / get_data(a);
    GetSeed<ADseed::h>::get_data(a) +=
        inv * (GetSeed<ADseed::h>::get_data(b) -
               inv * GetSeed<ADseed::p>::get_data(a) *
                   GetSeed<ADseed::b>::get_data(b));
  }

  atype& a;
  btype& b;
};

template <class atype, class btype>
KOKKOS_FUNCTION auto Log(ADObj<atype>& a, ADObj<btype>& b) {
  return LogExpr<ADObj<atype>, ADObj<btype>>(a, b);
}

template <class atype, class btype>
KOKKOS_FUNCTION auto Log(A2DObj<atype>& a, A2DObj<btype>& b) {
  return LogExpr<A2DObj<atype>, A2DObj<btype>>(a, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Exp(const T a, T& b) {
  b = std::exp(a);
}

template <class atype, class btype>
class ExpExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION ExpExpr(atype& a, btype& b) : a(a), b(b) {}

  KOKKOS_FUNCTION void eval() { get_data(b) = std::exp(get_data(a)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) = GetSeed<seed>::get_data(a) * get_data(b);
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) += GetSeed<seed>::get_data(b) * get_data(b);
  }

  KOKKOS_FUNCTION void hreverse() {
    GetSeed<ADseed::h>::get_data(a) +=
        get_data(b) *
        (GetSeed<ADseed::h>::get_data(b) +
         GetSeed<ADseed::p>::get_data(a) * GetSeed<ADseed::b>::get_data(b));
  }

  atype& a;
  btype& b;
};

template <class atype, class btype>
KOKKOS_FUNCTION auto Exp(ADObj<atype>& a, ADObj<btype>& b) {
  return ExpExpr<ADObj<atype>, ADObj<btype>>(a, b);
}

template <class atype, class btype>
KOKKOS_FUNCTION auto Exp(A2DObj<atype>& a, A2DObj<btype>& b) {
  return ExpExpr<A2DObj<atype>, A2DObj<btype>>(a, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Sin(const T a, T& b) {
  b = std::sin(a);
}

template <class atype, class btype>
class SinExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION SinExpr(atype& a, btype& b) : a(a), b(b) {}

  KOKKOS_FUNCTION void eval() { get_data(b) = std::sin(get_data(a)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) =
        GetSeed<seed>::get_data(a) * std::cos(get_data(a));
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) +=
        GetSeed<seed>::get_data(b) * std::cos(get_data(a));
  }

  KOKKOS_FUNCTION void hreverse() {
    GetSeed<ADseed::h>::get_data(a) +=
        (GetSeed<ADseed::h>::get_data(b) * std::cos(get_data(a)) -
         get_data(b) * GetSeed<ADseed::p>::get_data(a) *
             GetSeed<ADseed::b>::get_data(b));
  }

  atype& a;
  btype& b;
};

template <class atype, class btype>
KOKKOS_FUNCTION auto Sin(ADObj<atype>& a, ADObj<btype>& b) {
  return SinExpr<ADObj<atype>, ADObj<btype>>(a, b);
}

template <class atype, class btype>
KOKKOS_FUNCTION auto Sin(A2DObj<atype>& a, A2DObj<btype>& b) {
  return SinExpr<A2DObj<atype>, A2DObj<btype>>(a, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Cos(const T a, T& b) {
  b = std::cos(a);
}

template <class atype, class btype>
class CosExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION CosExpr(atype& a, btype& b) : a(a), b(b) {}

  KOKKOS_FUNCTION void eval() { get_data(b) = std::cos(get_data(a)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) =
        -GetSeed<seed>::get_data(a) * std::sin(get_data(a));
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) +=
        -GetSeed<seed>::get_data(b) * std::sin(get_data(a));
  }

  KOKKOS_FUNCTION void hreverse() {
    GetSeed<ADseed::h>::get_data(a) +=
        -(GetSeed<ADseed::h>::get_data(b) * std::sin(get_data(a)) +
          get_data(b) * GetSeed<ADseed::b>::get_data(b) *
              GetSeed<ADseed::p>::get_data(a));
  }

  atype& a;
  btype& b;
};

template <class atype, class btype>
KOKKOS_FUNCTION auto Cos(ADObj<atype>& a, ADObj<btype>& b) {
  return CosExpr<ADObj<atype>, ADObj<btype>>(a, b);
}

template <class atype, class btype>
KOKKOS_FUNCTION auto Cos(A2DObj<atype>& a, A2DObj<btype>& b) {
  return CosExpr<A2DObj<atype>, A2DObj<btype>>(a, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Pow(const T a, const T exponent, T& b) {
  b = std::pow(a, exponent);
}

template <class atype, class btype>
class PowExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION PowExpr(atype& a, const T exponent, btype& b)
      : a(a), exponent(exponent), b(b) {}

  KOKKOS_FUNCTION void eval() {
    get_data(b) = std::pow(get_data(a), exponent);
    inv = 1.0 / get_data(a);
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) =
        exponent * inv * get_data(b) * GetSeed<seed>::get_data(a);
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) +=
        exponent * inv * get_data(b) * GetSeed<seed>::get_data(b);
  }

  KOKKOS_FUNCTION void hreverse() {
    GetSeed<ADseed::h>::get_data(a) +=
        get_data(b) *
        (exponent * inv * GetSeed<ADseed::h>::get_data(b) +
         exponent * (exponent - 1.0) * inv * inv *
             GetSeed<ADseed::b>::get_data(b) * GetSeed<ADseed::p>::get_data(a));
  }

  T inv;
  atype& a;
  const T exponent;
  btype& b;
};

template <class atype, typename T, class btype>
KOKKOS_FUNCTION auto Pow(ADObj<atype>& a, const T exponent, ADObj<btype>& b) {
  return PowExpr<ADObj<atype>, ADObj<btype>>(a, exponent, b);
}

template <class atype, typename T, class btype>
KOKKOS_FUNCTION auto Pow(A2DObj<atype>& a, const T exponent, A2DObj<btype>& b) {
  return PowExpr<A2DObj<atype>, A2DObj<btype>>(a, exponent, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Sqrt(const T a, T& b) {
  b = std::sqrt(a);
}

template <class atype, class btype>
class SqrtExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  KOKKOS_FUNCTION SqrtExpr(atype& a, btype& b) : a(a), b(b) {}

  KOKKOS_FUNCTION void eval() {
    get_data(b) = std::sqrt(get_data(a));
    inv = T(1.0) / get_data(b);
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(b) = 0.5 * inv * GetSeed<seed>::get_data(a);
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    GetSeed<seed>::get_data(a) += 0.5 * inv * GetSeed<seed>::get_data(b);
  }

  KOKKOS_FUNCTION void hreverse() {
    GetSeed<ADseed::h>::get_data(a) +=
        (0.5 * inv * GetSeed<ADseed::h>::get_data(b) -
         0.25 * inv * inv * inv * GetSeed<ADseed::b>::get_data(b) *
             GetSeed<ADseed::p>::get_data(a));
  }

  T inv;
  atype& a;
  btype& b;
};

template <class atype, class btype>
KOKKOS_FUNCTION auto Sqrt(ADObj<atype>& a, ADObj<btype>& b) {
  return SqrtExpr<ADObj<atype>, ADObj<btype>>(a, b);
}

template <class atype, class btype>
KOKKOS_FUNCTION auto Sqrt(A2DObj<atype>& a, A2DObj<btype>& b) {
  return SqrtExpr<A2DObj<atype>, A2DObj<btype>>(a, b);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Mult(const T a, const T b, T& c) {
  c = a * b;
}

template <class atype, class btype, class ctype>
class MultExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  // Get the types of the different objects
  static constexpr ADiffType ada = get_diff_type<atype>::diff_type;
  static constexpr ADiffType adb = get_diff_type<btype>::diff_type;

  KOKKOS_FUNCTION MultExpr(atype a, btype b, ctype& c) : a(a), b(b), c(c) {}

  KOKKOS_FUNCTION void eval() { get_data(c) = get_data(a) * get_data(b); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = (GetSeed<seed>::get_data(a) * get_data(b) +
                                    get_data(a) * GetSeed<seed>::get_data(b));
    } else if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = GetSeed<seed>::get_data(a) * get_data(b);
    } else if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = get_data(a) * GetSeed<seed>::get_data(b);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(a) += get_data(b) * GetSeed<seed>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(b) += get_data(a) * GetSeed<seed>::get_data(c);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) +=
          get_data(b) * GetSeed<ADseed::h>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(b) +=
          get_data(a) * GetSeed<ADseed::h>::get_data(c);
    }
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) +=
          GetSeed<ADseed::p>::get_data(b) * GetSeed<ADseed::b>::get_data(c);
      GetSeed<ADseed::h>::get_data(b) +=
          GetSeed<ADseed::p>::get_data(a) * GetSeed<ADseed::b>::get_data(c);
    }
  }

  atype a;
  btype b;
  ctype& c;
};

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(ADObj<atype>& a, ADObj<btype>& b, ADObj<ctype>& c) {
  return MultExpr<ADObj<atype>&, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(const atype a, ADObj<btype>& b, ADObj<ctype>& c) {
  return MultExpr<const atype, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(ADObj<atype>& a, const btype b, ADObj<ctype>& c) {
  return MultExpr<ADObj<atype>&, const btype, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(A2DObj<atype>& a, A2DObj<btype>& b,
                          A2DObj<ctype>& c) {
  return MultExpr<A2DObj<atype>&, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(const atype a, A2DObj<btype>& b, A2DObj<ctype>& c) {
  return MultExpr<const atype, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Mult(A2DObj<atype>& a, const btype b, A2DObj<ctype>& c) {
  return MultExpr<A2DObj<atype>&, const btype, A2DObj<ctype>>(a, b, c);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Divide(const T a, const T b, T& c) {
  c = a / b;
}

template <class atype, class btype, class ctype>
class DivideExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  // Get the types of the different objects
  static constexpr ADiffType ada = get_diff_type<atype>::diff_type;
  static constexpr ADiffType adb = get_diff_type<btype>::diff_type;

  KOKKOS_FUNCTION DivideExpr(atype a, btype b, ctype& c) : a(a), b(b), c(c) {}

  KOKKOS_FUNCTION void eval() {
    inv = 1.0 / get_data(b);
    get_data(c) = inv * get_data(a);
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) =
          inv * (GetSeed<seed>::get_data(a) -
                 inv * get_data(a) * GetSeed<seed>::get_data(b));
    } else if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = inv * GetSeed<seed>::get_data(a);
    } else if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) =
          -inv * inv * get_data(a) * GetSeed<seed>::get_data(b);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(a) += inv * GetSeed<seed>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(b) -=
          inv * inv * get_data(a) * GetSeed<seed>::get_data(c);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) += inv * GetSeed<ADseed::h>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(b) -=
          inv * inv * get_data(a) * GetSeed<ADseed::h>::get_data(c);
      GetSeed<ADseed::h>::get_data(b) += 2.0 * inv * inv * inv * get_data(a) *
                                         GetSeed<ADseed::p>::get_data(b) *
                                         GetSeed<ADseed::b>::get_data(c);
    }
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) -= inv * inv *
                                         GetSeed<ADseed::p>::get_data(b) *
                                         GetSeed<ADseed::b>::get_data(c);
      GetSeed<ADseed::h>::get_data(b) -= inv * inv *
                                         GetSeed<ADseed::p>::get_data(a) *
                                         GetSeed<ADseed::b>::get_data(c);
    }
  }

  T inv;
  atype a;
  btype b;
  ctype& c;
};

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(ADObj<atype>& a, ADObj<btype>& b, ADObj<ctype>& c) {
  return DivideExpr<ADObj<atype>&, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(const atype a, ADObj<btype>& b, ADObj<ctype>& c) {
  return DivideExpr<const atype, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(ADObj<atype>& a, const btype b, ADObj<ctype>& c) {
  return DivideExpr<ADObj<atype>&, const btype, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(A2DObj<atype>& a, A2DObj<btype>& b,
                            A2DObj<ctype>& c) {
  return DivideExpr<A2DObj<atype>&, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(const atype a, A2DObj<btype>& b, A2DObj<ctype>& c) {
  return DivideExpr<const atype, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Divide(A2DObj<atype>& a, const btype b, A2DObj<ctype>& c) {
  return DivideExpr<A2DObj<atype>&, const btype, A2DObj<ctype>>(a, b, c);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Sum(const T a, const T b, T& c) {
  c = a + b;
}

template <class atype, class btype, class ctype>
class SumExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  // Get the types of the different objects
  static constexpr ADiffType ada = get_diff_type<atype>::diff_type;
  static constexpr ADiffType adb = get_diff_type<btype>::diff_type;

  KOKKOS_FUNCTION SumExpr(atype a, btype b, ctype& c) : a(a), b(b), c(c) {}

  KOKKOS_FUNCTION void eval() { get_data(c) = get_data(a) + get_data(b); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) =
          (GetSeed<seed>::get_data(a) + GetSeed<seed>::get_data(b));
    } else if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = GetSeed<seed>::get_data(a);
    } else if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = GetSeed<seed>::get_data(b);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(a) += GetSeed<seed>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(b) += GetSeed<seed>::get_data(c);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) += GetSeed<ADseed::h>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(b) += GetSeed<ADseed::h>::get_data(c);
    }
  }

  atype a;
  btype b;
  ctype& c;
};

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(ADObj<atype>& a, ADObj<btype>& b, ADObj<ctype>& c) {
  return SumExpr<ADObj<atype>&, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const atype a, ADObj<btype>& b, ADObj<ctype>& c) {
  return SumExpr<const atype, ADObj<btype>&, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(ADObj<atype>& a, const btype b, ADObj<ctype>& c) {
  return SumExpr<ADObj<atype>&, const btype, ADObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(A2DObj<atype>& a, A2DObj<btype>& b, A2DObj<ctype>& c) {
  return SumExpr<A2DObj<atype>&, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const atype a, A2DObj<btype>& b, A2DObj<ctype>& c) {
  return SumExpr<const atype, A2DObj<btype>&, A2DObj<ctype>>(a, b, c);
}

template <class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(A2DObj<atype>& a, const btype b, A2DObj<ctype>& c) {
  return SumExpr<A2DObj<atype>&, const btype, A2DObj<ctype>>(a, b, c);
}

template <typename T, std::enable_if_t<is_numeric_type<T>::value, bool> = true>
KOKKOS_FUNCTION void Sum(const T c1, const T a, const T c2, const T b, T& c) {
  c = c1 * a + c2 * b;
}

template <class atype, class btype, class ctype>
class SumScaleExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<btype>::type T;

  // Get the types of the different objects
  static constexpr ADiffType ada = get_diff_type<atype>::diff_type;
  static constexpr ADiffType adb = get_diff_type<btype>::diff_type;

  KOKKOS_FUNCTION SumScaleExpr(const T c1, atype a, const T c2, btype b,
                               ctype& c)
      : c1(c1), a(a), c2(c2), b(b), c(c) {}

  KOKKOS_FUNCTION void eval() {
    get_data(c) = c1 * get_data(a) + c2 * get_data(b);
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (ada == ADiffType::ACTIVE && adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) =
          c1 * GetSeed<seed>::get_data(a) + c2 * GetSeed<seed>::get_data(b);
    } else if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = GetSeed<seed>::get_data(a);
    } else if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(c) = GetSeed<seed>::get_data(b);
    }
  }

  KOKKOS_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(a) += c1 * GetSeed<seed>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(b) += c2 * GetSeed<seed>::get_data(c);
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (ada == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(a) += c1 * GetSeed<ADseed::h>::get_data(c);
    }
    if constexpr (adb == ADiffType::ACTIVE) {
      GetSeed<ADseed::h>::get_data(b) += c2 * GetSeed<ADseed::h>::get_data(c);
    }
  }

  const T c1;
  atype a;
  const T c2;
  btype b;
  ctype& c;
};

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, ADObj<atype>& a, const T c2,
                         ADObj<btype>& b, ADObj<ctype>& c) {
  return SumScaleExpr<ADObj<atype>&, ADObj<btype>&, ADObj<ctype>>(c1, a, c2, b,
                                                                  c);
}

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, const atype a, const T c2, ADObj<btype>& b,
                         ADObj<ctype>& c) {
  return SumScaleExpr<const atype, ADObj<btype>&, ADObj<ctype>>(c1, a, c2, b,
                                                                c);
}

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, ADObj<atype>& a, const T c2, const btype b,
                         ADObj<ctype>& c) {
  return SumScaleExpr<ADObj<atype>&, const btype, ADObj<ctype>>(c1, a, c2, b,
                                                                c);
}

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, A2DObj<atype>& a, const T c2,
                         A2DObj<btype>& b, A2DObj<ctype>& c) {
  return SumScaleExpr<A2DObj<atype>&, A2DObj<btype>&, A2DObj<ctype>>(c1, a, c2,
                                                                     b, c);
}

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, const atype a, const T c2,
                         A2DObj<btype>& b, A2DObj<ctype>& c) {
  return SumScaleExpr<const atype, A2DObj<btype>&, A2DObj<ctype>>(c1, a, c2, b,
                                                                  c);
}

template <typename T, class atype, class btype, class ctype>
KOKKOS_FUNCTION auto Sum(const T c1, A2DObj<atype>& a, const T c2,
                         const btype b, A2DObj<ctype>& c) {
  return SumScaleExpr<A2DObj<atype>&, const btype, A2DObj<ctype>>(c1, a, c2, b,
                                                                  c);
}
*/
namespace Test {

template <typename T>
class ScalarTest : public A2DTest<T, T, T> {
 public:
  using Input = VarTuple<T, T>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return std::string("ScalarOperations"); }

  void get_point(Input& x) { x[0] = 0.35; }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T a, b, c, d, e, f, q, r, s, t, u;
    x.get_values(a);

    u = a * a + a;

    // log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a);

    // exp(sin(cos(a)) * a * cos(-a) * sin(a *
    // sin(sqrt(a) * exp(a) * a)));
    // u = exp(sin(a));  // sin(a * cos(a) * sin(a * sin(a * exp(a) * a)));
    // expr.eval();
    // u = expr.value();

    // Log(a, b);
    // Sin(b, c);
    // Exp(c, d);
    // Pow(d, T(3.5), e);
    // Cos(e, f);
    // Mult(c, f, q);
    // Sum(d, q, r);
    // Sum(T(-1.0), q, T(3.14), r, s);
    // Divide(q, s, t);
    // Sqrt(t, u);

    return MakeVarTuple<T>(u);
  }

  template <class A, class Ty>
  void evaluate(ADExpr<A, Ty> expr) {
    expr.eval();
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> a, u;  // , b, c, d, e, f, q, r, s, t, u;
    x.get_values(a.value());
    auto expr = a * a;

    std::cout << "Address: " << &a << std::endl;

    evaluate(expr);

    // auto stack = MakeStack(Eval(a * a + a, u));
    //         Eval(log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a), u));
    // auto stack = MakeStack(Log(a, b),                       //
    //                        Sin(b, c),                       //
    //                        Exp(c, d),                       //
    //                        Pow(d, T(3.5), e),               //
    //                        Cos(e, f),                       //
    //                        Mult(c, f, q),                   //
    //                        Sum(d, q, r),                    //
    //                        Sum(T(-1.0), q, T(3.14), r, s),  //
    //                        Divide(q, s, t),                 //
    //                        Sqrt(t, u));

    // seed.get_values(u.bvalue());
    // stack.reverse();
    // g.set_values(a.bvalue());

    // //  + a;  // log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a);
    // // auto expr = log(a * a * sqrt(exp(a * sin(a))));
    // // exp(sin(cos(a)) * a * cos(-a) * sin(a * sin(sqrt(a) * exp(a) * a)));

    // // sin(a * cos(a) * sin(a * sin(a * exp(a) * a)));
    // expr.eval();
    // expr.bvalue() = seed[0];
    // expr.reverse();

    // g[0] = a.bvalue();
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> a;
    x.get_values(a.value());
    p.get_values(a.pvalue());

    auto expr = log(a * a * sqrt(exp(a * sin(a))));
    expr.eval();
    expr.bvalue() = seed[0];
    expr.reverse();

    expr.hforward();
    expr.hvalue() = hval[0];
    expr.hreverse();

    h[0] = a.hvalue();

    // A2DObj<T> a, b, c, d, e, f, q, r, s, t, u;
    // x.get_values(a.value());
    // p.get_values(a.pvalue());
    // auto stack = MakeStack(Log(a, b),                       //
    //                        Sin(b, c),                       //
    //                        Exp(c, d),                       //
    //                        Pow(d, T(3.5), e),               //
    //                        Cos(e, f),                       //
    //                        Mult(c, f, q),                   //
    //                        Sum(d, q, r),                    //
    //                        Sum(T(-1.0), q, T(3.14), r, s),  //
    //                        Divide(q, s, t),                 //
    //                        Sqrt(t, u));
    // seed.get_values(u.bvalue());
    // hval.get_values(u.hvalue());
    // stack.reverse();
    // stack.hforward();
    // stack.hreverse();
    // h.set_values(a.hvalue());
  }
};

bool ScalarTestAll(bool component, bool write_output) {
  using Tc = std::complex<double>;

  bool passed = true;
  ScalarTest<Tc> test1;
  passed = passed && Run(test1, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SCALAR_OPS_H