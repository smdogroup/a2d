#ifndef A2D_VAR_TUPLE_H
#define A2D_VAR_TUPLE_H

#include <tuple>
#include <type_traits>

#include "a2denum.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T>
struct __is_complex : public std::false_type {};

template <typename T>
struct __is_complex<std::complex<T>> : public std::true_type {};

template <typename T>
struct __is_scalar_type {
  static const bool value =
      std::is_arithmetic<T>::value || __is_complex<T>::value;
};

struct __basic_arithmetic_type {
  static const index_t ncomp = 1;
};

template <class... Vars>
struct __count_var_components;

template <>
struct __count_var_components<> {
  static const index_t ncomp = 0;
};

template <class First, class... Remain>
struct __count_var_components<First, Remain...> {
  static const index_t ncomp =
      std::conditional<__is_scalar_type<First>::value, __basic_arithmetic_type,
                       First>::type::ncomp +
      __count_var_components<Remain...>::ncomp;
};

template <typename T, class... Vars>
class VarTupleBase {
 public:
  /**
   * @brief Number of components in the tuple of variables
   */
  static constexpr index_t ncomp = __count_var_components<Vars...>::ncomp;

  /**
   * @brief Get the number of components
   */
  KOKKOS_FUNCTION index_t get_num_components() const { return ncomp; }

 protected:
  template <typename I, index_t index, class TupleObj, class First,
            class... Remain>
  KOKKOS_FUNCTION T& get_value_(TupleObj& var, const I comp) {
    if constexpr (__is_scalar_type<First>::value) {
      if (comp == 0) {
        return std::get<index>(var);
      } else if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var);
      } else {
        return get_value_<I, index + 1, TupleObj, Remain...>(var, comp - 1);
      }
    } else {
      if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var)[comp];
      } else if (comp < First::ncomp) {
        return std::get<index>(var)[comp];
      } else {
        return get_value_<I, index + 1, TupleObj, Remain...>(
            var, comp - First::ncomp);
      }
    }
  }

  template <typename I, index_t index, class TupleObj, class First,
            class... Remain>
  KOKKOS_FUNCTION const T& get_value_const_(const TupleObj& var,
                                            const I comp) const {
    if constexpr (__is_scalar_type<First>::value) {
      if (comp == 0) {
        return std::get<index>(var);
      } else if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var);
      } else {
        return get_value_const_<I, index + 1, TupleObj, Remain...>(var,
                                                                   comp - 1);
      }
    } else {
      if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var)[comp];
      } else if (comp < First::ncomp) {
        return std::get<index>(var)[comp];
      } else {
        return get_value_const_<I, index + 1, TupleObj, Remain...>(
            var, comp - First::ncomp);
      }
    }
  }

  template <index_t index, class TupleObj, class First, class... Remain>
  KOKKOS_FUNCTION void set_values_(TupleObj& var, const First& f,
                                   const Remain&... r) {
    if constexpr (__is_scalar_type<First>::value) {
      std::get<index>(var) = f;
    } else {
      First& val = std::get<index>(var);
      for (index_t i = 0; i < First::ncomp; i++) {
        val[i] = f[i];
      }
    }
    if constexpr (sizeof...(Remain) > 0) {
      set_values_<index + 1, TupleObj, Remain...>(var, r...);
    }
  }

  template <index_t index, class TupleObj, class First, class... Remain>
  KOKKOS_FUNCTION void get_values_(const TupleObj& var, First& f,
                                   Remain&... r) const {
    if constexpr (__is_scalar_type<First>::value) {
      f = std::get<index>(var);
    } else {
      const First& val = std::get<index>(var);
      for (index_t i = 0; i < First::ncomp; i++) {
        f[i] = val[i];
      }
    }
    if constexpr (sizeof...(Remain) > 0) {
      get_values_<index + 1, TupleObj, Remain...>(var, r...);
    }
  }

  template <index_t index, class TupleObj, class First, class... Remain>
  KOKKOS_FUNCTION void zero_(TupleObj& var) {
    if constexpr (__is_scalar_type<First>::value) {
      std::get<index>(var) = T(0.0);
    } else {
      std::get<index>(var).zero();
    }
    if constexpr (sizeof...(Remain) > 0) {
      zero_<index + 1, TupleObj, Remain...>(var);
    }
  }

  template <index_t index, class TupleObj, class First, class... Remain>
  KOKKOS_FUNCTION void set_rand_(TupleObj& var, const T low, const T high) {
    if constexpr (__is_scalar_type<First>::value) {
      std::get<index>(var) =
          low + (high - low) * (static_cast<double>(rand()) / RAND_MAX);
    } else {
      First& val = std::get<index>(var);
      for (index_t i = 0; i < First::ncomp; i++) {
        val[i] = low + (high - low) * (static_cast<double>(rand()) / RAND_MAX);
      }
    }
    if constexpr (sizeof...(Remain) > 0) {
      set_rand_<index + 1, TupleObj, Remain...>(var, low, high);
    }
  }
};

template <typename T, class... Vars>
class VarTuple : public VarTupleBase<T, Vars...> {
 public:
  using VarTupleObj = std::tuple<Vars...>;

  KOKKOS_FUNCTION VarTuple() {}
  KOKKOS_FUNCTION VarTuple(const Vars&... s) {
    this->template set_values_<0, VarTupleObj, Vars...>(var, s...);
  }

  /// @brief Access a reference to an object indexed by this class
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
    return this->template get_value_<I, 0, VarTupleObj, Vars...>(var, comp);
  }

  /// @brief Access a reference to an object indexed by this class
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
    return this->template get_value_const_<I, 0, VarTupleObj, Vars...>(var,
                                                                       comp);
  }

  /// @brief Zero the components of the tuple
  KOKKOS_FUNCTION void zero() {
    this->template zero_<0, VarTupleObj, Vars...>(var);
  }

  /// @brief Set a random set of values on an interval
  KOKKOS_FUNCTION void set_rand(T low = T(-1.0), T high = T(1.0)) {
    this->template set_rand_<0, VarTupleObj, Vars...>(var, low, high);
  }

  /// @brief Set values into the tuple from a list of objects
  KOKKOS_FUNCTION void set_values(const Vars&... s) {
    this->template set_values_<0, VarTupleObj, Vars...>(var, s...);
  }

  /// @brief Place the values from the tuple into list of objects
  KOKKOS_FUNCTION void get_values(Vars&... s) const {
    this->template get_values_<0, VarTupleObj, Vars...>(var, s...);
  }

 private:
  VarTupleObj var;
};

template <typename T, class... Vars>
KOKKOS_FUNCTION auto MakeVarTuple(Vars&... s) {
  return VarTuple<T, Vars...>(s...);
}

template <typename T, class... Vars>
class TieTuple : public VarTupleBase<T, Vars...> {
 public:
  using VarTupleObj = std::tuple<Vars&...>;

  KOKKOS_FUNCTION TieTuple(Vars&... s) : var(s...) {}

  /// @brief Access a reference to an object indexed by this class
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
    return this->template get_value_<I, 0, VarTupleObj, Vars...>(var, comp);
  }

  /// @brief Access a reference to an object indexed by this class
  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
    return this->template get_value_const_<I, 0, VarTupleObj, Vars...>(var,
                                                                       comp);
  }

  /// @brief Zero the components of the tuple
  KOKKOS_FUNCTION void zero() {
    this->template zero_<0, VarTupleObj, Vars...>(var);
  }

  /// @brief Set a random set of values on an interval
  KOKKOS_FUNCTION void set_rand(T low = T(-1.0), T high = T(1.0)) {
    this->template set_rand_<0, VarTupleObj, Vars...>(var, low, high);
  }

  /// @brief Set values into the tuple from a list of objects
  KOKKOS_FUNCTION void set_values(const Vars&... s) {
    this->template set_values_<0, VarTupleObj, Vars...>(var, s...);
  }

  /// @brief Place the values from the tuple into list of objects
  KOKKOS_FUNCTION void get_values(Vars&... s) const {
    this->template get_values_<0, VarTupleObj, Vars...>(var, s...);
  }

 private:
  VarTupleObj var;
};

template <typename T, class... Vars>
KOKKOS_FUNCTION auto MakeTieTuple(Vars&... s) {
  return TieTuple<T, Vars...>(s...);
}

template <typename T, ADseed seed, class... Vars>
KOKKOS_FUNCTION auto MakeTieTuple(Vars&... s) {
  return MakeTieTuple<T>(GetSeed<seed>::get_obj(s)...);
}

}  // namespace A2D

#endif  // A2D_VAR_TUPLE_H