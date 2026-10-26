#ifndef A2D_VAR_TUPLE_H
#define A2D_VAR_TUPLE_H

#include <type_traits>

#include "a2ddefs.h"
#include "a2dmat.h"
#include "a2dobj.h"
#include "a2dvec.h"
#include "adscalar.h"

namespace A2D {

template <class Var, typename Enable = void>
struct VarTraits;  // Unsupported types produce a compile error

// Scalar
template <typename T>
struct VarTraits<T, std::enable_if_t<is_numeric_v<T>>> {
  static constexpr int ncomp = 1;
  using ScalarType = T;
  using Ref = T&;
  using ConstRef = const T&;
};

// Vector
template <typename T, int N>
struct VarTraits<Vec<T, N>> {
  static constexpr int ncomp = N;
  using ScalarType = T;
  using Ref = VecRef<T, N>;
  using ConstRef = VecRef<const T, N>;
};

// Matrix
template <typename T, int M, int N>
struct VarTraits<Mat<T, M, N>> {
  static constexpr int ncomp = M * N;
  using ScalarType = T;
  using Ref = MatRef<T, M, N>;
  using ConstRef = MatRef<const T, M, N>;
};

// Symmetric matrix
template <typename T, int N>
struct VarTraits<SymMat<T, N>> {
  static constexpr int ncomp = N * (N + 1) / 2;
  using ScalarType = T;
  using Ref = SymMatRef<T, N>;
  using ConstRef = SymMatRef<const T, N>;
};

// Existing references
template <typename T, int N>
struct VarTraits<VecRef<T, N>> : VarTraits<Vec<T, N>> {};

template <typename T, int M, int N>
struct VarTraits<MatRef<T, M, N>> : VarTraits<Mat<T, M, N>> {};

template <typename T, int N>
struct VarTraits<SymMatRef<T, N>> : VarTraits<SymMat<T, N>> {};

// Constant references
template <typename Var>
struct VarTraits<const Var> : VarTraits<Var> {
  using Ref = typename VarTraits<Var>::ConstRef;
};

// Compute the offset based on the variable list
template <int I, typename First, typename... Remain>
A2D_FUNCTION constexpr int var_offset_impl() {
  if constexpr (I == 0) {
    return 0;
  } else {
    return VarTraits<First>::ncomp + var_offset_impl<I - 1, Remain...>();
  }
}

template <int I, class... Vars>
A2D_FUNCTION constexpr int var_offset() {
  static_assert(I < sizeof...(Vars));
  return var_offset_impl<I, Vars...>();
}

/*
 * Extract the type of the variable at the given index
 */
template <int I, typename... Vars>
struct var_element;

template <int I, typename First, typename... Remain>
struct var_element<I, First, Remain...> : var_element<I - 1, Remain...> {};

template <typename First, typename... Remain>
struct var_element<0, First, Remain...> {
  using type = First;
};

template <int I, typename... Vars>
using var_element_t = typename var_element<I, Vars...>::type;

template <typename T, class... Vars>
class VarTuple {
 public:
  static_assert((std::is_same_v<T, typename VarTraits<Vars>::ScalarType> &&
                 ...),
                "All variables must use the same scalar type");
  static constexpr int ncomp = (int{0} + ... + VarTraits<Vars>::ncomp);

  // Get the number of components
  int get_ncomp() const { return ncomp; }

  template <typename I>
  A2D_FUNCTION T& operator[](const I i) {
    return values[i];
  }

  template <typename I>
  A2D_FUNCTION const T& operator[](const I i) const {
    return values[i];
  }

  // Zero the values
  A2D_FUNCTION void zero() {
    for (int i = 0; i < ncomp; i++) {
      values[i] = 0.0;
    }
  }

  // Set a random set of values on an interval
  A2D_FUNCTION void set_rand(const T low = T(-1.0), const T high = T(1.0)) {
    for (int i = 0; i < ncomp; i++) {
      values[i] = low + (high - low) * (static_cast<double>(rand()) / RAND_MAX);
    }
  }

  // Set values into the VarTuple based on the object values
  A2D_FUNCTION void set_values(const Vars&... vars) {
    set_values_impl<0, Vars...>(vars...);
  }

  // Get the object values from the VarTuple
  A2D_FUNCTION void get_values(Vars&... vars) const {
    get_values_impl<0, Vars...>(vars...);
  }

  T values[ncomp];

 private:
  template <int I, class Var, class... Remain>
  A2D_FUNCTION void set_values_impl(const Var& var, const Remain&... remain) {
    set_value<I, Var>(var);
    if constexpr (sizeof...(Remain) > 0) {
      set_values_impl<I + 1, Remain...>(remain...);
    }
  }

  template <int I, class Var>
  A2D_FUNCTION void set_value(const Var& var) {
    constexpr int offset = var_offset<I, Vars...>();
    constexpr int ncomp = VarTraits<Var>::ncomp;

    if constexpr (is_numeric_v<Var>) {
      values[offset] = var;
    } else {
      for (int i = 0; i < ncomp; i++) {
        values[offset + i] = var[i];
      }
    }
  }

  template <int I, class Var, class... Remain>
  A2D_FUNCTION void get_values_impl(Var& var, Remain&... remain) const {
    get_value<I, Var>(var);
    if constexpr (sizeof...(Remain) > 0) {
      get_values_impl<I + 1, Remain...>(remain...);
    }
  }

  template <int I, class Var>
  A2D_FUNCTION void get_value(Var& var) const {
    constexpr int offset = var_offset<I, Vars...>();
    constexpr int ncomp = VarTraits<Var>::ncomp;

    if constexpr (is_numeric_v<Var>) {
      var = values[offset];
    } else {
      for (int i = 0; i < ncomp; i++) {
        var[i] = values[offset + i];
      }
    }
  }
};

template <typename T, class... Vars>
A2D_FUNCTION VarTuple<T, Vars...> MakeVarTuple(const Vars&... vars) {
  VarTuple<T, Vars...> t;
  t.set_values(vars...);
  return t;
}

template <int I, typename T, class... Vars>
A2D_FUNCTION decltype(auto) get(VarTuple<T, Vars...>& t) {
  static_assert(I < sizeof...(Vars));

  using Var = var_element_t<I, Vars...>;
  using Ref = typename VarTraits<Var>::Ref;

  constexpr int offset = var_offset<I, Vars...>();

  if constexpr (is_numeric_v<Var>) {
    return (t.values[offset]);  // T&
  } else {
    return Ref(t.values + offset);
  }
}

template <int I, typename T, class... Vars>
A2D_FUNCTION decltype(auto) get(const VarTuple<T, Vars...>& t) {
  static_assert(I < sizeof...(Vars));

  using Var = var_element_t<I, Vars...>;
  using ConstRef = typename VarTraits<Var>::ConstRef;

  constexpr int offset = var_offset<I, Vars...>();

  if constexpr (is_numeric_v<Var>) {
    return (t.values[offset]);  // T&
  } else {
    return ConstRef(t.values + offset);
  }
}

}  // namespace A2D

#endif  // A2D_VAR_TUPLE_H
