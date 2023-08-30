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

  KOKKOS_FUNCTION Vec() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  KOKKOS_FUNCTION Vec(const T* vals) {
    for (int i = 0; i < N; i++) {
      V[i] = vals[i];
    }
  }
  template <class VecType>
  KOKKOS_FUNCTION Vec(const VecType& vec) {
    for (int i = 0; i < N; i++) {
      V[i] = vec(i);
    }
  }
  KOKKOS_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  template <class IdxType>
  KOKKOS_FUNCTION T& operator()(const IdxType i) {
    return V[i];
  }
  template <class IdxType>
  KOKKOS_FUNCTION const T& operator()(const IdxType i) const {
    return V[i];
  }

  T* data() { return V; }

  // private:

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return V[i];
  }
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return V[i];
  }

  T V[N];
};

template <class VecType>
class ADVec {
 public:
  KOKKOS_FUNCTION ADVec(VecType& V, VecType& Vb) : V(V), Vb(Vb) {}

  KOKKOS_FUNCTION VecType& value() { return V; }
  KOKKOS_FUNCTION const VecType& value() const { return V; }

  KOKKOS_FUNCTION VecType& bvalue() { return Vb; }
  KOKKOS_FUNCTION const VecType& bvalue() const { return Vb; }

  VecType& V;   // Vector
  VecType& Vb;  // Reverse mode derivative value
};

template <class VecType>
class A2DVec {
 public:
  KOKKOS_FUNCTION A2DVec() {}
  KOKKOS_FUNCTION A2DVec(const VecType& V) : V(V) {}
  KOKKOS_FUNCTION A2DVec(const VecType& V, const VecType& Vb) : V(V), Vb(Vb) {}
  KOKKOS_FUNCTION A2DVec(const VecType& V, const VecType& Vb, const VecType& Vp)
      : V(V), Vb(Vb), Vp(Vp) {}
  KOKKOS_FUNCTION A2DVec(const VecType& V, const VecType& Vb, const VecType& Vp,
                         const VecType& Vh)
      : V(V), Vb(Vb), Vp(Vp), Vh(Vh) {}

  KOKKOS_FUNCTION VecType& value() { return V; }
  KOKKOS_FUNCTION const VecType& value() const { return V; }

  KOKKOS_FUNCTION VecType& bvalue() { return Vb; }
  KOKKOS_FUNCTION const VecType& bvalue() const { return Vb; }

  KOKKOS_FUNCTION VecType& pvalue() { return Vp; }
  KOKKOS_FUNCTION const VecType& pvalue() const { return Vp; }

  KOKKOS_FUNCTION VecType& hvalue() { return Vh; }
  KOKKOS_FUNCTION const VecType& hvalue() const { return *Vh; }

  VecType V;
  VecType Vb;
  VecType Vp;
  VecType Vh;
};

/**
 * @brief Get data pointers from Mat/ADMat/A2DMat objects
 */
template <typename T, int n>
KOKKOS_FUNCTION T* get_data(Vec<T, n>& vec) {
  return vec.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(ADVec<Vec<T, n>>& vec) {
  return vec.V.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_data(A2DVec<Vec<T, n>>& vec) {
  return vec.V.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_bdata(ADVec<Vec<T, n>>& vec) {
  return vec.Vb.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_bdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vb.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_pdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vp.V;
}

template <typename T, int n>
KOKKOS_FUNCTION T* get_hdata(A2DVec<Vec<T, n>>& vec) {
  return vec.Vh.V;
}

/**
 * @brief Select type based on whether the vector is passive or active (can be
 * differentiated)
 *
 * @tparam adiff_type passive or active
 * @tparam VecType the numeric type of the vector
 */
template <ADiffType adiff_type, ADorder order, class VecType>
using ADVecType = typename std::conditional<
    adiff_type == ADiffType::ACTIVE,
    typename std::conditional<order == ADorder::FIRST, ADVec<VecType>,
                              A2DVec<VecType>>::type,
    const VecType>::type;

}  // namespace A2D
#endif  // A2D_VEC_H