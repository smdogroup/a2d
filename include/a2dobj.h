#ifndef A2D_OBJECTS_H
#define A2D_OBJECTS_H

#include "a2ddefs.h"
#include "a2dmat.h"
#include "a2dvec.h"
#include "adscalar.h"

namespace A2D {

/*
  Expression template for first derivatives
*/
template <class A, typename T>
class ADExpr {
 public:
  A2D_FUNCTION A& self() { return static_cast<A&>(*this); }
  A2D_FUNCTION const A& self() const { return static_cast<const A&>(*this); }

  // Evaluation and derivatives
  A2D_FUNCTION void eval() { self().eval(); }
  A2D_FUNCTION void forward() { self().forward(); }
  A2D_FUNCTION void reverse() { self().reverse(); }

  // Access the values and derivatives
  A2D_FUNCTION T& value() { return self().value(); }
  A2D_FUNCTION const T& value() const { return self().value(); }
  A2D_FUNCTION T& bvalue() { return self().bvalue(); }
  A2D_FUNCTION const T& bvalue() const { return self().bvalue(); }
};

/*
  The base first-order AD class
*/
template <class T>
class ADObj : public ADExpr<ADObj<T>, T> {
 public:
  using type = object_scalar_t<T>;
  static constexpr ADObjType obj_type = a2d_obj_type_v<T>;

  // If this is a scalar, non-reference type - initialize values to zero
  A2D_FUNCTION ADObj() {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      A = Ab = type(0.0);
    }
  }

  // If this is a scalar, non-reference type - initialize values to zero
  A2D_FUNCTION ADObj(const T& A) : A(A) {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      Ab = type(0.0);
    }
  }

  // Initialize with both values
  template <typename U = T,
            std::enable_if_t<std::is_reference<U>::value, int> = 0>
  A2D_FUNCTION ADObj(T& A, T& Ab) : A(A), Ab(Ab) {}

  template <typename U = T,
            std::enable_if_t<!std::is_reference<U>::value, int> = 0>
  A2D_FUNCTION ADObj(const T& A, const T& Ab) : A(A), Ab(Ab) {}

  // Evaluation and derivatives
  A2D_FUNCTION void eval() {}
  A2D_FUNCTION void forward() {}
  A2D_FUNCTION void reverse() {}
  A2D_FUNCTION void bzero() {
    if constexpr (obj_type == ADObjType::SCALAR) {
      Ab = type(0.0);
    } else {
      Ab.zero();
    }
  }

  // Extract the data from the underlying objects
  A2D_FUNCTION T& value() { return A; }
  A2D_FUNCTION const T& value() const { return A; }
  A2D_FUNCTION T& bvalue() { return Ab; }
  A2D_FUNCTION const T& bvalue() const { return Ab; }

  template <typename I, typename U = T,
            std::enable_if_t<is_a2d_vector_v<U>, int> = 0>
  A2D_FUNCTION ADObj<type&> operator[](const I i) {
    return ADObj<type&>(A[i], Ab[i]);
  }

  template <typename I, typename U = T,
            std::enable_if_t<is_a2d_vector_v<U>, int> = 0>
  A2D_FUNCTION ADObj<type&> operator()(const I i) {
    return ADObj<type&>(A[i], Ab[i]);
  }

  template <
      typename I, typename U = T,
      std::enable_if_t<is_a2d_matrix_v<U> || is_a2d_sym_matrix_v<U>, int> = 0>
  A2D_FUNCTION ADObj<type&> operator()(const I i, const I j) {
    return ADObj<type&>(A(i, j), Ab(i, j));
  }

 private:
  T A;   // Object
  T Ab;  // Reverse mode derivative value
};

/*
  Expression template for second derivatives
*/
template <class A, typename T>
class A2DExpr {
 public:
  A2D_FUNCTION A& self() { return static_cast<A&>(*this); }
  A2D_FUNCTION const A& self() const { return static_cast<const A&>(*this); }

  // Evaluation and derivatives
  A2D_FUNCTION void eval() { self().eval(); }
  A2D_FUNCTION void reverse() { self().reverse(); }
  A2D_FUNCTION void hforward() { self().hforward(); }
  A2D_FUNCTION void hreverse() { self().hreverse(); }
  A2D_FUNCTION void bzero() { self().bzero(); }
  A2D_FUNCTION void hzero() { self().hzero(); }

  // Access the values and derivatives
  A2D_FUNCTION T& value() { return self().value(); }
  A2D_FUNCTION const T& value() const { return self().value(); }
  A2D_FUNCTION T& bvalue() { return self().bvalue(); }
  A2D_FUNCTION const T& bvalue() const { return self().bvalue(); }
  A2D_FUNCTION T& pvalue() { return self().pvalue(); }
  A2D_FUNCTION const T& pvalue() const { return self().pvalue(); }
  A2D_FUNCTION T& hvalue() { return self().hvalue(); }
  A2D_FUNCTION const T& hvalue() const { return self().hvalue(); }
};

/*
  The base second-order AD class
*/
template <class T>
class A2DObj : public A2DExpr<A2DObj<T>, T> {
 public:
  using type = object_scalar_t<T>;
  static constexpr ADObjType obj_type = a2d_obj_type_v<T>;

  A2D_FUNCTION A2DObj() {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      A = Ab = Ap = Ah = type(0.0);
    }
  }
  A2D_FUNCTION A2DObj(const T& A) : A(A) {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      Ab = Ap = Ah = type(0.0);
    }
  }
  A2D_FUNCTION A2DObj(const T& A, const T& Ab) : A(A), Ab(Ab) {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      Ap = Ah = type(0.0);
    }
  }
  A2D_FUNCTION A2DObj(const T& A, const T& Ab, const T& Ap)
      : A(A), Ab(Ab), Ap(Ap) {
    if constexpr (is_numeric_v<T> && !std::is_reference<T>::value) {
      Ah = type(0.0);
    }
  }
  template <typename U = T,
            std::enable_if_t<std::is_reference<U>::value, int> = 0>
  A2D_FUNCTION A2DObj(T& A, T& Ab, T& Ap, T& Ah)
      : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  template <typename U = T,
            std::enable_if_t<!std::is_reference<U>::value, int> = 0>
  A2D_FUNCTION A2DObj(const T& A, const T& Ab, const T& Ap, const T& Ah)
      : A(A), Ab(Ab), Ap(Ap), Ah(Ah) {}

  // Evaluation and derivatives
  A2D_FUNCTION void eval() {}
  A2D_FUNCTION void reverse() {}
  A2D_FUNCTION void hforward() {}
  A2D_FUNCTION void hreverse() {}

  A2D_FUNCTION void bzero() {
    if constexpr (obj_type == ADObjType::SCALAR) {
      Ab = type(0.0);
    } else {
      Ab.zero();
    }
  }
  A2D_FUNCTION void hzero() {
    if constexpr (obj_type == ADObjType::SCALAR) {
      Ah = type(0.0);
    } else {
      Ah.zero();
    }
  }
  A2D_FUNCTION T& value() { return A; }
  A2D_FUNCTION const T& value() const { return A; }
  A2D_FUNCTION T& bvalue() { return Ab; }
  A2D_FUNCTION const T& bvalue() const { return Ab; }
  A2D_FUNCTION T& pvalue() { return Ap; }
  A2D_FUNCTION const T& pvalue() const { return Ap; }
  A2D_FUNCTION T& hvalue() { return Ah; }
  A2D_FUNCTION const T& hvalue() const { return Ah; }

  template <typename I, typename U = T,
            std::enable_if_t<is_a2d_vector_v<U>, int> = 0>
  A2D_FUNCTION A2DObj<type&> operator[](const I i) {
    return A2DObj<type&>(A[i], Ab[i], Ap[i], Ah[i]);
  }

  template <typename I, typename U = T,
            std::enable_if_t<is_a2d_vector_v<U>, int> = 0>
  A2D_FUNCTION A2DObj<type&> operator()(const I i) {
    return A2DObj<type&>(A[i], Ab[i], Ap[i], Ah[i]);
  }

  template <
      typename I, typename U = T,
      std::enable_if_t<is_a2d_matrix_v<U> || is_a2d_sym_matrix_v<U>, int> = 0>
  A2D_FUNCTION A2DObj<type&> operator()(const I i, const I j) {
    return A2DObj<type&>(A(i, j), Ab(i, j), Ap(i, j), Ah(i, j));
  }

 private:
  T A;   // Object
  T Ab;  // Reverse mode derivative value
  T Ap;  // Projected second derivative value
  T Ah;  // Reverse mode second derivative
};

// Extract underlying values from AD wrappers
template <typename T>
A2D_FUNCTION decltype(auto) get_value(T& obj) {
  return (obj);
}

template <typename T>
A2D_FUNCTION decltype(auto) get_value(ADObj<T>& obj) {
  return obj.value();
}

template <typename T>
A2D_FUNCTION decltype(auto) get_value(A2DObj<T>& obj) {
  return obj.value();
}

// Extract data from the underlying value
template <typename T>
A2D_FUNCTION decltype(auto) get_data(T& obj) {
  auto&& value = get_value(obj);

  if constexpr (is_numeric_v<decltype(value)>) {
    return (value);
  } else {
    return value.get_data();
  }
}

/**
 * @brief Get objects and pointers to seed data (bvalue(), pvalue(), hvalue)
 */
template <ADseed seed>
class GetSeed {
 public:
  template <typename T>
  static A2D_FUNCTION T& get_obj(ADObj<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    return value.bvalue();
  }

  template <typename T>
  static A2D_FUNCTION T& get_obj(A2DObj<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      return value.bvalue();
    } else if constexpr (seed == ADseed::p) {
      return value.pvalue();
    } else {  // seed == ADseed::h
      return value.hvalue();
    }
  }

  template <typename T>
  static A2D_FUNCTION decltype(auto) get_data(ADObj<T>& value) {
    static_assert(seed == ADseed::b, "Incompatible seed type for ADObj");
    if constexpr (is_numeric_v<T>) {
      return (value.bvalue());
    } else {
      return value.bvalue().get_data();
    }
  }

  template <typename T>
  static A2D_FUNCTION decltype(auto) get_data(A2DObj<T>& value) {
    static_assert(seed == ADseed::b or seed == ADseed::p or seed == ADseed::h,
                  "Incompatible seed type for A2DObj");
    if constexpr (seed == ADseed::b) {
      if constexpr (is_numeric_v<T>) {
        return (value.bvalue());
      } else {
        return value.bvalue().get_data();
      }
    } else if constexpr (seed == ADseed::p) {
      if constexpr (is_numeric_v<T>) {
        return (value.pvalue());
      } else {
        return value.pvalue().get_data();
      }
    } else {  // seed == ADseed::h
      if constexpr (is_numeric_v<T>) {
        return (value.hvalue());
      } else {
        return value.hvalue().get_data();
      }
    }
  }
};

}  // namespace A2D

#endif  // A2D_OBJECTS_H
