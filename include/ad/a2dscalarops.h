#ifndef A2D_SCALAR_OPS_H
#define A2D_SCALAR_OPS_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T>
KOKKOS_FUNCTION void Log(const T a, T& b) {
  b = std::log(a);
}

template <typename T, ADorder order>
class LogExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Log(ADObj<T>& a, ADObj<T>& b) {
  return LogExpr<T, ADorder::FIRST>(a, b);
}

template <typename T>
KOKKOS_FUNCTION auto Log(A2DObj<T>& a, A2DObj<T>& b) {
  return LogExpr<T, ADorder::SECOND>(a, b);
}

template <typename T>
KOKKOS_FUNCTION void Exp(const T a, T& b) {
  b = std::exp(a);
}

template <typename T, ADorder order>
class ExpExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Exp(ADObj<T>& a, ADObj<T>& b) {
  return ExpExpr<T, ADorder::FIRST>(a, b);
}

template <typename T>
KOKKOS_FUNCTION auto Exp(A2DObj<T>& a, A2DObj<T>& b) {
  return ExpExpr<T, ADorder::SECOND>(a, b);
}

template <typename T>
KOKKOS_FUNCTION void Sin(const T a, T& b) {
  b = std::sin(a);
}

template <typename T, ADorder order>
class SinExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Sin(ADObj<T>& a, ADObj<T>& b) {
  return SinExpr<T, ADorder::FIRST>(a, b);
}

template <typename T>
KOKKOS_FUNCTION auto Sin(A2DObj<T>& a, A2DObj<T>& b) {
  return SinExpr<T, ADorder::SECOND>(a, b);
}

template <typename T>
KOKKOS_FUNCTION void Cos(const T a, T& b) {
  b = std::cos(a);
}

template <typename T, ADorder order>
class CosExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Cos(ADObj<T>& a, ADObj<T>& b) {
  return CosExpr<T, ADorder::FIRST>(a, b);
}

template <typename T>
KOKKOS_FUNCTION auto Cos(A2DObj<T>& a, A2DObj<T>& b) {
  return CosExpr<T, ADorder::SECOND>(a, b);
}

template <typename T>
KOKKOS_FUNCTION void Pow(const T a, const T exponent, T& b) {
  b = std::pow(a, exponent);
}

template <typename T, ADorder order>
class PowExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Pow(ADObj<T>& a, const T exponent, ADObj<T>& b) {
  return PowExpr<T, ADorder::FIRST>(a, exponent, b);
}

template <typename T>
KOKKOS_FUNCTION auto Pow(A2DObj<T>& a, const T exponent, A2DObj<T>& b) {
  return PowExpr<T, ADorder::SECOND>(a, exponent, b);
}

template <typename T>
KOKKOS_FUNCTION void Sqrt(const T a, T& b) {
  b = std::sqrt(a);
}

template <typename T, ADorder order>
class SqrtExpr {
 public:
  using atype = ADScalarType<ADiffType::ACTIVE, order, T>;
  using btype = ADScalarType<ADiffType::ACTIVE, order, T>;

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
  const T exponent;
  btype& b;
};

template <typename T>
KOKKOS_FUNCTION auto Sqrt(ADObj<T>& a, ADObj<T>& b) {
  return SqrtExpr<T, ADorder::FIRST>(a, b);
}

template <typename T>
KOKKOS_FUNCTION auto Sqrt(A2DObj<T>& a, A2DObj<T>& b) {
  return SqrtExpr<T, ADorder::SECOND>(a, b);
}

template <typename T>
KOKKOS_FUNCTION void Mult(const T a, const T b, T& c) {
  c = a * b;
}

template <typename T, ADorder order, ADiffType ada, ADiffType adb>
class MultExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using btype = ADScalarInputType<adb, order, T>;
  using ctype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Mult(ADObj<T>& a, ADObj<T>& b, ADObj<T>& c) {
  return MultExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(a, b,
                                                                           c);
}

template <typename T>
KOKKOS_FUNCTION auto Mult(const T a, ADObj<T>& b, ADObj<T>& c) {
  return MultExpr<T, ADorder::FIRST, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Mult(ADObj<T>& a, const T b, ADObj<T>& c) {
  return MultExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Mult(A2DObj<T>& a, A2DObj<T>& b, A2DObj<T>& c) {
  return MultExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Mult(const T a, A2DObj<T>& b, A2DObj<T>& c) {
  return MultExpr<T, ADorder::SECOND, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Mult(A2DObj<T>& a, const T b, A2DObj<T>& c) {
  return MultExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION void Divide(const T a, const T b, T& c) {
  c = a / b;
}

template <typename T, ADorder order, ADiffType ada, ADiffType adb>
class DivideExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using btype = ADScalarInputType<adb, order, T>;
  using ctype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Divide(ADObj<T>& a, ADObj<T>& b, ADObj<T>& c) {
  return DivideExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Divide(const T a, ADObj<T>& b, ADObj<T>& c) {
  return DivideExpr<T, ADorder::FIRST, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Divide(ADObj<T>& a, const T b, ADObj<T>& c) {
  return DivideExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Divide(A2DObj<T>& a, A2DObj<T>& b, A2DObj<T>& c) {
  return DivideExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Divide(const T a, A2DObj<T>& b, A2DObj<T>& c) {
  return DivideExpr<T, ADorder::SECOND, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Divide(A2DObj<T>& a, const T b, A2DObj<T>& c) {
  return DivideExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION void Sum(const T a, const T b, T& c) {
  c = a + b;
}

template <typename T, ADorder order, ADiffType ada, ADiffType adb>
class SumExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using btype = ADScalarInputType<adb, order, T>;
  using ctype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Sum(ADObj<T>& a, ADObj<T>& b, ADObj<T>& c) {
  return SumExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(a, b,
                                                                          c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T a, ADObj<T>& b, ADObj<T>& c) {
  return SumExpr<T, ADorder::FIRST, ADiffType::PASSIVE, ADiffType::ACTIVE>(a, b,
                                                                           c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(ADObj<T>& a, const T b, ADObj<T>& c) {
  return SumExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::PASSIVE>(a, b,
                                                                           c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(A2DObj<T>& a, A2DObj<T>& b, A2DObj<T>& c) {
  return SumExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::ACTIVE>(a, b,
                                                                           c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T a, A2DObj<T>& b, A2DObj<T>& c) {
  return SumExpr<T, ADorder::SECOND, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(A2DObj<T>& a, const T b, A2DObj<T>& c) {
  return SumExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      a, b, c);
}

template <typename T>
KOKKOS_FUNCTION void Sum(const T c1, const T a, const T c2, const T b, T& c) {
  c = c1 * a + c2 * b;
}

template <typename T, ADorder order, ADiffType ada, ADiffType adb>
class SumScaleExpr {
 public:
  using atype = ADScalarInputType<ada, order, T>;
  using btype = ADScalarInputType<adb, order, T>;
  using ctype = ADScalarType<ADiffType::ACTIVE, order, T>;

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

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, ADObj<T>& a, const T c2, ADObj<T>& b,
                         ADObj<T>& c) {
  return SumScaleExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      c1, a, c2, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, const T a, const T c2, ADObj<T>& b,
                         ADObj<T>& c) {
  return SumScaleExpr<T, ADorder::FIRST, ADiffType::PASSIVE, ADiffType::ACTIVE>(
      c1, a, c2, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, ADObj<T>& a, const T c2, const T b,
                         ADObj<T>& c) {
  return SumScaleExpr<T, ADorder::FIRST, ADiffType::ACTIVE, ADiffType::PASSIVE>(
      c1, a, c2, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, A2DObj<T>& a, const T c2, A2DObj<T>& b,
                         A2DObj<T>& c) {
  return SumScaleExpr<T, ADorder::SECOND, ADiffType::ACTIVE, ADiffType::ACTIVE>(
      c1, a, c2, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, const T a, const T c2, A2DObj<T>& b,
                         A2DObj<T>& c) {
  return SumScaleExpr<T, ADorder::SECOND, ADiffType::PASSIVE,
                      ADiffType::ACTIVE>(c1, a, c2, b, c);
}

template <typename T>
KOKKOS_FUNCTION auto Sum(const T c1, A2DObj<T>& a, const T c2, const T b,
                         A2DObj<T>& c) {
  return SumScaleExpr<T, ADorder::SECOND, ADiffType::ACTIVE,
                      ADiffType::PASSIVE>(c1, a, c2, b, c);
}

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
    Log(a, b);
    Sin(b, c);
    Exp(c, d);
    Pow(d, T(3.5), e);
    Cos(e, f);
    Mult(c, f, q);
    Sum(d, q, r);
    Sum(T(-1.0), q, T(3.14), r, s);
    Divide(q, s, t);
    Sqrt(t, u);

    return MakeVarTuple<T>(u);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> a, b, c, d, e, f, q, r, s, t, u;
    x.get_values(a.value());
    auto stack = MakeStack(Log(a, b),                       //
                           Sin(b, c),                       //
                           Exp(c, d),                       //
                           Pow(d, T(3.5), e),               //
                           Cos(e, f),                       //
                           Mult(c, f, q),                   //
                           Sum(d, q, r),                    //
                           Sum(T(-1.0), q, T(3.14), r, s),  //
                           Divide(q, s, t),                 //
                           Sqrt(t, u));

    seed.get_values(u.bvalue());
    stack.reverse();
    g.set_values(a.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> a, b, c, d, e, f, q, r, s, t, u;
    x.get_values(a.value());
    p.get_values(a.pvalue());
    auto stack = MakeStack(Log(a, b),                       //
                           Sin(b, c),                       //
                           Exp(c, d),                       //
                           Pow(d, T(3.5), e),               //
                           Cos(e, f),                       //
                           Mult(c, f, q),                   //
                           Sum(d, q, r),                    //
                           Sum(T(-1.0), q, T(3.14), r, s),  //
                           Divide(q, s, t),                 //
                           Sqrt(t, u));
    seed.get_values(u.bvalue());
    hval.get_values(u.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(a.hvalue());
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