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
  static const index_t num_components = 1;
};

template <class... Vars>
struct __count_var_components;

template <>
struct __count_var_components<> {
  static const index_t num_components = 0;
};

template <class First, class... Remain>
struct __count_var_components<First, Remain...> {
  static const index_t num_components =
      std::conditional<__is_scalar_type<First>::value, __basic_arithmetic_type,
                       First>::type::num_components +
      __count_var_components<Remain...>::num_components;
};

template <typename T, class... Vars>
class VarTuple {
 public:
  typedef std::tuple<Vars...> VarTupleObj;

  VarTuple() {}
  VarTuple(const Vars&... s) { set_values_<0, Vars...>(s...); }

  /**
   * @brief Number of components in the tuple of variables
   */
  static constexpr index_t num_components =
      __count_var_components<Vars...>::num_components;

  /**
   * @brief Get the number of components
   */
  index_t get_num_components() const { return num_components; }

  /**
   * @brief Access a reference to an object indexed by this class
   */
  template <typename I>
  T& operator[](const I comp) {
    return get_value<I, 0, Vars...>(comp);
  }

  /**
   * @brief Access a reference to an object indexed by this class
   */
  template <typename I>
  const T& operator[](const I comp) const {
    return get_value<I, 0, Vars...>(comp);
  }

  /**
   * @brief Set a random set of values on an interval
   */
  void set_rand(T low = T(-1.0), T high = T(1.0)) {
    for (index_t i = 0; i < num_components; i++) {
      (*this)[i] =
          low + (high - low) * (static_cast<double>(rand()) / RAND_MAX);
    }
  }

  /**
   * @brief Set values into the tuple from a list of objects
   *
   * @param s The list of const references to objects
   */
  void set_values(const Vars&... s) { set_values_<0, Vars...>(s...); }

  /**
   * @brief Place the values from the tuple into list of objects
   *
   * @param s The list of references to objects
   */
  void get_values(Vars&... s) const { get_values_<0, Vars...>(s...); }

 private:
  VarTupleObj var;

  template <typename I, index_t index, class First, class... Remain>
  T& get_value(const I comp) {
    if constexpr (__is_scalar_type<First>::value) {
      if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var);
      } else if (comp == 0) {
        return std::get<index>(var);
      } else {
        return get_value<I, index + 1, Remain...>(comp - 1);
      }
    } else {
      if constexpr (sizeof...(Remain) == 0) {
        return std::get<index>(var)[comp];
      } else if (comp < First::num_components) {
        return std::get<index>(var)[comp];
      } else {
        return get_value<I, index + 1, Remain...>(comp - First::num_components);
      }
    }
  }

  template <index_t index, class First, class... Remain>
  void set_values_(const First& f, const Remain&... r) {
    if constexpr (__is_scalar_type<First>::value) {
      std::get<index>(var) = f;
    } else {
      First& val = std::get<index>(var);
      for (index_t i = 0; i < First::num_components; i++) {
        val[i] = f[i];
      }
    }
    if constexpr (sizeof...(Remain) > 0) {
      set_values_<index + 1, Remain...>(r...);
    }
  }

  template <index_t index, class First, class... Remain>
  void get_values_(First& f, Remain&... r) const {
    if constexpr (__is_scalar_type<First>::value) {
      f = std::get<index>(var);
    } else {
      const First& val = std::get<index>(var);
      for (index_t i = 0; i < First::num_components; i++) {
        f[i] = val[i];
      }
    }
    if constexpr (sizeof...(Remain) > 0) {
      get_values_<index + 1, Remain...>(r...);
    }
  }
};

template <typename T, class... Vars>
VarTuple<T, Vars...> MakeVarTuple(Vars&... s) {
  return VarTuple<T, Vars...>(s...);
}

}  // namespace A2D

#endif  // A2D_VAR_TUPLE_H