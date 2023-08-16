#ifndef A2D_VEC_H
#define A2D_VEC_H

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T, int N>
class Vec {
 public:
  typedef T type;

  static const index_t num_components = N;

  A2D_INLINE_FUNCTION Vec() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  A2D_INLINE_FUNCTION Vec(const T* vals) {
    for (int i = 0; i < N; i++) {
      V[i] = vals[i];
    }
  }
  template <class VecType>
  A2D_INLINE_FUNCTION Vec(const VecType& vec) {
    for (int i = 0; i < N; i++) {
      V[i] = vec(i);
    }
  }
  A2D_INLINE_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION T& operator()(const IdxType i) {
    return V[i];
  }
  template <class IdxType>
  A2D_INLINE_FUNCTION const T& operator()(const IdxType i) const {
    return V[i];
  }

  T* data() { return V; }

  // private:

  template <typename I>
  A2D_INLINE_FUNCTION T& operator[](const I i) {
    return V[i];
  }
  template <typename I>
  A2D_INLINE_FUNCTION const T& operator[](const I i) const {
    return V[i];
  }

  T V[N];
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
  A2D_INLINE_FUNCTION A2DVec(const VecType& V) : V(V) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& V, const VecType& Vb)
      : V(V), Vb(Vb) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& V, const VecType& Vb,
                             const VecType& Vp)
      : V(V), Vb(Vb), Vp(Vp) {}
  A2D_INLINE_FUNCTION A2DVec(const VecType& V, const VecType& Vb,
                             const VecType& Vp, const VecType& Vh)
      : V(V), Vb(Vb), Vp(Vp), Vh(Vh) {}

  A2D_INLINE_FUNCTION VecType& value() { return V; }
  A2D_INLINE_FUNCTION const VecType& value() const { return V; }

  A2D_INLINE_FUNCTION VecType& bvalue() { return Vb; }
  A2D_INLINE_FUNCTION const VecType& bvalue() const { return Vb; }

  A2D_INLINE_FUNCTION VecType& pvalue() { return Vp; }
  A2D_INLINE_FUNCTION const VecType& pvalue() const { return Vp; }

  A2D_INLINE_FUNCTION VecType& hvalue() { return Vh; }
  A2D_INLINE_FUNCTION const VecType& hvalue() const { return *Vh; }

  VecType V;
  VecType Vb;
  VecType Vp;
  VecType Vh;
};

/**
 * @brief Get data pointers from Mat/ADMat/A2DMat objects
 */
template <typename T, int n>
A2D_INLINE_FUNCTION T* get_data(Vec<T, n>& vec) {
  return vec.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_data(ADVec<Vec<T, n>>& vec) {
  return vec.V.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_data(A2DVec<Vec<T, n>>& vec) {
  return vec.V.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_bdata(ADVec<Vec<T, n>>& vec) {
  return vec.Vb.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_bdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vb.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_pdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vp.V;
}

template <typename T, int n>
A2D_INLINE_FUNCTION T* get_hdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vh.V;
}

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