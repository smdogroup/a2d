#ifndef A2D_VEC_H
#define A2D_VEC_H

#include "a2ddefs.h"

namespace A2D {

template <typename T, int N>
class Vec {
 public:
  typedef T type;
  static const ADObjType obj_type = ADObjType::VECTOR;
  static const index_t ncomp = N;

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

  T* get_data() { return V; }
  const T* get_data() const { return V; }

  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I i) {
    return V[i];
  }
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I i) const {
    return V[i];
  }

 private:
  T V[N];
};

}  // namespace A2D
#endif  // A2D_VEC_H