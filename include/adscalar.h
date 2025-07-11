#ifndef A2D_ADSCALAR_H
#define A2D_ADSCALAR_H

#include <complex>
#include <type_traits>

#include "a2ddefs.h"

namespace A2D {

template <class T, int N>
class ADScalar {
 public:
  using type = T;

  A2D_FUNCTION ADScalar() {}

  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION ADScalar(const R value) : value(value), deriv{0.0} {}

  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION ADScalar(const R value, const T d[]) : value(value) {
    for (int i = 0; i < N; i++) {
      deriv[i] = d[i];
    }
  }

  A2D_FUNCTION ADScalar(const ADScalar<T, N> &r) : value(r.value) {
    for (int i = 0; i < N; i++) {
      deriv[i] = r.deriv[i];
    }
  }

  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline ADScalar<T, N> &operator=(const R &r) {
    value = r;
    for (int i = 0; i < N; i++) {
      deriv[i] = 0.0;
    }
    return *this;
  }

  // Comparison operators
  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline bool operator<(const R &rhs) const {
    return value < rhs;
  }
  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline bool operator<=(const R &rhs) const {
    return value <= rhs;
  }
  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline bool operator>(const R &rhs) const {
    return value > rhs;
  }
  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline bool operator>=(const R &rhs) const {
    return value >= rhs;
  }
  template <typename R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline bool operator!=(const R &rhs) const {
    return value != rhs;
  }

  template <typename X, int M>
  A2D_FUNCTION inline bool operator<(const ADScalar<X, M> &rhs) const {
    return value < rhs.value;
  }
  template <typename X, int M>
  A2D_FUNCTION inline bool operator<=(const ADScalar<X, M> &rhs) const {
    return value <= rhs.value;
  }
  template <typename X, int M>
  A2D_FUNCTION inline bool operator>(const ADScalar<X, M> &rhs) const {
    return value > rhs.value;
  }
  template <typename X, int M>
  A2D_FUNCTION inline bool operator>=(const ADScalar<X, M> &rhs) const {
    return value >= rhs.value;
  }

  // Operator +=, -=, *=, /=
  A2D_FUNCTION inline ADScalar<T, N> &operator+=(const ADScalar<T, N> &r) {
    value += r.value;
    for (int i = 0; i < N; i++) {
      deriv[i] += r.deriv[i];
    }
    return *this;
  }
  template <class R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline ADScalar<T, N> &operator+=(const R &r) {
    value += r;
    return *this;
  }
  A2D_FUNCTION inline ADScalar<T, N> &operator-=(const ADScalar<T, N> &r) {
    value -= r.value;
    for (int i = 0; i < N; i++) {
      deriv[i] -= r.deriv[i];
    }
    return *this;
  }
  template <class R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline ADScalar<T, N> &operator-=(const R &r) {
    value -= r;
    return *this;
  }
  A2D_FUNCTION inline ADScalar<T, N> &operator*=(const ADScalar<T, N> &r) {
    for (int i = 0; i < N; i++) {
      deriv[i] = r.value * deriv[i] + value * r.deriv[i];
    }
    value *= r.value;
    return *this;
  }
  template <class R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline ADScalar<T, N> &operator*=(const R &r) {
    value *= r;
    for (int i = 0; i < N; i++) {
      deriv[i] = r * deriv[i];
    }
    return *this;
  }

  A2D_FUNCTION inline ADScalar<T, N> &operator/=(const ADScalar<T, N> &r) {
    T inv = 1.0 / r.value;
    T inv2 = value * inv * inv;
    value *= inv;
    for (int i = 0; i < N; i++) {
      deriv[i] = inv * deriv[i] - inv2 * r.deriv[i];
    }
    return *this;
  }
  template <class R, typename = std::enable_if_t<is_scalar_type<R>::value>>
  A2D_FUNCTION inline ADScalar<T, N> operator/=(const R &r) {
    T inv = 1.0 / r;
    value *= inv;
    for (int i = 0; i < N; i++) {
      deriv[i] = inv * deriv[i];
    }
    return *this;
  }

  A2D_FUNCTION inline ADScalar<T, N> operator-() const {
    T negderivs[N];
    for (int i = 0; i < N; i++) {
      negderivs[i] = -deriv[i];
    }
    return ADScalar<T, N>(-value,
                          negderivs);  // Return by value, not by reference
  }

  //  private:
  T value;
  T deriv[N];
};

// Addition
template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> operator+(const ADScalar<X, M> &l,
                                             const ADScalar<X, M> &r) {
  ADScalar<X, M> out(l.value + r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = l.deriv[i] + r.deriv[i];
  }
  return out;
}
template <class X, int M, class L,
          typename = std::enable_if_t<is_scalar_type<L>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator+(const L &l,
                                             const ADScalar<X, M> &r) {
  return ADScalar<X, M>(r.value + l, r.deriv);
}
template <class X, int M, class R,
          typename = std::enable_if_t<is_scalar_type<R>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator+(const ADScalar<X, M> &l,
                                             const R &r) {
  return ADScalar<X, M>(l.value + r, l.deriv);
}

// Subtraction
template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> operator-(const ADScalar<X, M> &l,
                                             const ADScalar<X, M> &r) {
  ADScalar<X, M> out(l.value - r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = l.deriv[i] - r.deriv[i];
  }
  return out;
}
template <class X, int M, class L,
          typename = std::enable_if_t<is_scalar_type<L>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator-(const L &l,
                                             const ADScalar<X, M> &r) {
  ADScalar<X, M> out(l - r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = -r.deriv[i];
  }
  return out;
}
template <class X, int M, class R,
          typename = std::enable_if_t<is_scalar_type<R>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator-(const ADScalar<X, M> &l,
                                             const R &r) {
  return ADScalar<X, M>(l.value - r, l.deriv);
}

// Multiplication
template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> operator*(const ADScalar<X, M> &l,
                                             const ADScalar<X, M> &r) {
  ADScalar<X, M> out(l.value * r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = r.value * l.deriv[i] + r.deriv[i] * l.value;
  }
  return out;
}
template <class X, int M, class L,
          typename = std::enable_if_t<is_scalar_type<L>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator*(const L &l,
                                             const ADScalar<X, M> &r) {
  ADScalar<X, M> out(l * r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = r.deriv[i] * l;
  }
  return out;
}
template <class X, int M, class R,
          typename = std::enable_if_t<is_scalar_type<R>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator*(const ADScalar<X, M> &l,
                                             const R &r) {
  ADScalar<X, M> out(l.value * r);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = l.deriv[i] * r;
  }
  return out;
}

// Division
template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> operator/(const ADScalar<X, M> &l,
                                             const ADScalar<X, M> &r) {
  X inv = 1.0 / r.value;
  X inv2 = l.value * inv * inv;
  ADScalar<X, M> out(inv * l.value);

  for (int i = 0; i < M; i++) {
    out.deriv[i] = inv * l.deriv[i] - inv2 * r.deriv[i];
  }
  return out;
}
template <class X, int M, class L,
          typename = std::enable_if_t<is_scalar_type<L>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator/(const L &l,
                                             const ADScalar<X, M> &r) {
  X inv = 1.0 / r.value;
  X inv2 = l * inv * inv;
  ADScalar<X, M> out(inv * l);

  for (int i = 0; i < M; i++) {
    out.deriv[i] = -inv2 * r.deriv[i];
  }
  return out;
}
template <class X, int M, class R,
          typename = std::enable_if_t<is_scalar_type<R>::value>>
A2D_FUNCTION inline ADScalar<X, M> operator/(const ADScalar<X, M> &l,
                                             const R &r) {
  X inv = 1.0 / r;
  ADScalar<X, M> out(inv * l.value);

  for (int i = 0; i < M; i++) {
    out.deriv[i] = inv * l.deriv[i];
  }
  return out;
}

// sign function
// template <class X, int M>
// A2D_FUNCTION inline ADScalar<X, M> fsgn(const ADScalar<X, M> &r) {
//   X sign = 1.0;
//   if (r.value < 0.0) {
//     sign = -1.0;
//   }
//   // device compatible fsgn
//   ADScalar<X, M> out(::fsgn(r.value));
//   for (int i = 0; i < M; i++) {
//     out.deriv[i] = sign * r.deriv[i];
//   }
//   return out;
// }

// fabs, sqrt
template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> fabs(const ADScalar<X, M> &r) {
  X scalar = 1.0;
  if (r.value < 0.0) {
    scalar = -1.0;
  }
  // device compatible fabs
  ADScalar<X, M> out(::fabs(r.value));
  for (int i = 0; i < M; i++) {
    out.deriv[i] = scalar * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> sqrt(const ADScalar<X, M> &r) {
  // device compatible sqrt
  X value = ::sqrt(r.value);
  ADScalar<X, M> out(value);
  X inv = 0.5 / value;
  for (int i = 0; i < M; i++) {
    out.deriv[i] = inv * r.deriv[i];
  }
  return out;
}

template <class X, int M, class R,
          typename = std::enable_if_t<is_scalar_type<R>::value>>
A2D_FUNCTION inline ADScalar<X, M> pow(const ADScalar<X, M> &r,
                                       const R &exponent) {
  // device compatible pow
  X value = ::pow(r.value, exponent);
  ADScalar<X, M> out(value);
  X inv = exponent * value / r.value;
  for (int i = 0; i < M; i++) {
    out.deriv[i] = inv * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> exp(const ADScalar<X, M> &r) {
  // device compatible exp
  X value = ::exp(r.value);
  ADScalar<X, M> out(value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = value * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> log(const ADScalar<X, M> &r) {
  // device compatible log
  ADScalar<X, M> out(::log(r.value));
  X inv = 1.0 / r.value;
  for (int i = 0; i < M; i++) {
    out.deriv[i] = inv * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> sin(const ADScalar<X, M> &r) {
  // device compatible sin, cos
  ADScalar<X, M> out(::sin(r.value));
  X d = ::cos(r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = d * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> cos(const ADScalar<X, M> &r) {
  // device compatible sin, cos
  ADScalar<X, M> out(::cos(r.value));
  X d = -::sin(r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = d * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> atan(const ADScalar<X, M> &r) {
  // device compatible sin, cos
  ADScalar<X, M> out(::atan(r.value));
  X d = 1.0 / (1.0 + r.value * r.value);  // 1/(1+x^2)
  for (int i = 0; i < M; i++) {
    out.deriv[i] = d * r.deriv[i];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> atan2(const ADScalar<X, M> &y,
                                         const ADScalar<X, M> &x) {
  /** atan2(y,x) => theta */
  ADScalar<X, M> out(::atan2(y.value, x.value));
  X denom = x.value * x.value + y.value * y.value;
  X dx = -y.value / denom;
  X dy = x.value / denom;

  for (int i = 0; i < M; i++) {
    out.deriv[i] = dx * x.deriv[0] + dy * y.deriv[0];
  }
  return out;
}

template <class X, int M>
A2D_FUNCTION inline ADScalar<X, M> tanh(const ADScalar<X, M> &r) {
  // for smooth sign function essentially
  ADScalar<X, M> out(::tanh(r.value));
  X d = 1.0 / ::cosh(r.value) / ::cosh(r.value);
  for (int i = 0; i < M; i++) {
    out.deriv[i] = d * r.deriv[i];
  }
  return out;
}

// template <int N>
// struct __get_object_numeric_type<ADScalar<double,N>> {
//   using type = ADScalar<double,N>;
// };

// template <int N>
// struct __get_object_numeric_type<ADScalar<A2D_complex_t<double>,N>> {
//   using type = ADScalar<A2D_complex_t<double>,N>;
// };

}  // namespace A2D

#endif  // A2D_ADSCALAR_H
