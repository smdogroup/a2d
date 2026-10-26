#ifndef A2D_VEC_H
#define A2D_VEC_H

#include "../a2ddefs.h"

namespace A2D {

template <typename T, int N>
class Vec {
 public:
  typedef T type;
  static const ADObjType obj_type = ADObjType::VECTOR;
  static const int ncomp = N;

  A2D_FUNCTION explicit Vec() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  template <typename T2>
  A2D_FUNCTION explicit Vec(const T2* vals) {
    for (int i = 0; i < N; i++) {
      V[i] = vals[i];
    }
  }
  template <typename T2>
  A2D_FUNCTION explicit Vec(const Vec<T2, N>& src) {
    for (int i = 0; i < N; i++) {
      V[i] = src(i);
    }
  }
  A2D_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      V[i] = 0.0;
    }
  }
  template <typename T2>
  A2D_FUNCTION void copy(const Vec<T2, N>& vec) {
    for (int i = 0; i < N; i++) {
      V[i] = vec(i);
    }
  }
  template <class IdxType>
  A2D_FUNCTION T& operator()(const IdxType i) {
    return V[i];
  }
  template <class IdxType>
  A2D_FUNCTION const T& operator()(const IdxType i) const {
    return V[i];
  }

  A2D_FUNCTION T* get_data() { return V; }
  A2D_FUNCTION const T* get_data() const { return V; }

  template <typename I>
  A2D_FUNCTION T& operator[](const I i) {
    return V[i];
  }
  template <typename I>
  A2D_FUNCTION const T& operator[](const I i) const {
    return V[i];
  }

 private:
  T V[N];
};

template <typename T, int N>
class VecRef {
 public:
  typedef T type;
  static const ADObjType obj_type = ADObjType::VECTOR;
  static const int ncomp = N;

  A2D_FUNCTION explicit VecRef(T* ptr) : ptr(ptr) {}
  A2D_FUNCTION void zero() {
    for (int i = 0; i < N; i++) {
      ptr[i] = 0.0;
    }
  }
  template <typename T2>
  A2D_FUNCTION void copy(const Vec<T2, N>& vec) {
    for (int i = 0; i < N; i++) {
      ptr[i] = vec(i);
    }
  }
  template <class IdxType>
  A2D_FUNCTION T& operator()(const IdxType i) {
    return ptr[i];
  }
  template <class IdxType>
  A2D_FUNCTION const T& operator()(const IdxType i) const {
    return ptr[i];
  }

  A2D_FUNCTION T* get_data() { return ptr; }
  A2D_FUNCTION const T* get_data() const { return ptr; }

  template <typename I>
  A2D_FUNCTION T& operator[](const I i) {
    return ptr[i];
  }
  template <typename I>
  A2D_FUNCTION const T& operator[](const I i) const {
    return ptr[i];
  }

 private:
  T* ptr;
};

}  // namespace A2D

#endif  // A2D_VEC_H
