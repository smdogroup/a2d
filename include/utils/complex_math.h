#ifndef A2D_COMPLEX_MATH
#define A2D_COMPLEX_MATH

#include <complex>

namespace A2D {

// The following two operator overloadings compute the product of a complex
// number and an integer
template <typename T, typename I,
          std::enable_if_t<std::is_integral<I>::value, bool> = true>
constexpr std::complex<T> operator*(const std::complex<T>& lhs, const I& rhs) {
  return std::complex<T>(lhs.real() * rhs, lhs.imag() * rhs);
}

template <typename T, typename I,
          std::enable_if_t<std::is_integral<I>::value, bool> = true>
constexpr std::complex<T> operator*(const I& lhs, const std::complex<T>& rhs) {
  return std::complex<T>(rhs.real() * lhs, rhs.imag() * lhs);
}

}  // namespace A2D

#endif  // A2D_COMPLEX_MATH