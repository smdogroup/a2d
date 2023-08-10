#ifndef A2D_VEC_H
#define A2D_VEC_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T, int N>
class Vec {
 public:
  typedef T type;

  A2D_INLINE_FUNCTION Vec() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  A2D_INLINE_FUNCTION Vec(const T* vals) {
    for (int i = 0; i < N; i++) {
      x[i] = vals[i];
    }
  }
  template <class VecType>
  A2D_INLINE_FUNCTION Vec(const VecType& vec) {
    for (int i = 0; i < N; i++) {
      x[i] = vec(i);
    }
  }
  A2D_INLINE_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      x[i] = 0.0;
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION T& operator()(const IdxType i) {
    return x[i];
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType i) const {
    return x[i];
  }

  T* data() { return x; }

  T x[N];
};

template <class VecType>
class ADVec {
 public:
  A2D_INLINE_FUNCTION ADVec(VecType& V, VecType& Vb) : V(V), Vb(Vb) {}

  A2D_INLINE_FUNCTION VecType& value() { return V; }
  A2D_INLINE_FUNCTION const VecType& value() const { return V; }

  A2D_INLINE_FUNCTION VecType& bvalue() { return Vb; }
  A2D_INLINE_FUNCTION const VecType& bvalue() const { return Vb; }

  VecType& V;   // Vector
  VecType& Vb;  // Reverse mode derivative value
};

template <class VecType>
class A2DVec {
 public:
  A2D_INLINE_FUNCTION A2DVec() {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x) : x(x) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb)
      : x(x), xb(xb) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb,
                             const VecType& xp)
      : x(x), xb(xb), xp(xp) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& x, const VecType& xb,
                             const VecType& xp, const VecType& xh)
      : x(x), xb(xb), xp(xp), xh(xh) {}

  A2D_INLINE_FUNCTION VecType& value() { return x; }
  A2D_INLINE_FUNCTION const VecType& value() const { return x; }

  A2D_INLINE_FUNCTION VecType& bvalue() { return xb; }
  A2D_INLINE_FUNCTION const VecType& bvalue() const { return xb; }

  A2D_INLINE_FUNCTION VecType& pvalue() { return xp; }
  A2D_INLINE_FUNCTION const VecType& pvalue() const { return xp; }

  A2D_INLINE_FUNCTION VecType& hvalue() { return xh; }
  A2D_INLINE_FUNCTION const VecType& hvalue() const { return *xh; }

  VecType x;
  VecType xb;
  VecType xp;
  VecType xh;
};

/**
 * @brief Select type based on whether the vector is passive or active (can be
 * differentiated)
 *
 * @tparam adiff_type passive or active
 * @tparam VecType the numeric type of the vector
 */
template <ADiffType adiff_type, class VecType>
using ADVecType = typename std::conditional<adiff_type == ADiffType::ACTIVE,
                                            ADVec<VecType>, VecType>::type;

template <ADiffType adiff_type, class VecType>
using A2DVecType = typename std::conditional<adiff_type == ADiffType::ACTIVE,
                                             A2DVec<VecType>, VecType>::type;

}  // namespace A2D
#endif  // A2D_VEC_H