#ifndef A2D_STACK_H
#define A2D_STACK_H

namespace A2D {

template <typename T>
struct tag {
  using type = T;
};

template <typename... Ts>
struct select_last {
  // Use a fold-expression to fold the comma operator over the parameter pack.
  using type = typename decltype((tag<Ts>{}, ...))::type;
};

template <class... Operations>
class Stack {
 public:
  using StackTuple = std::tuple<Operations&...>;
  static constexpr index_t num_ops = sizeof...(Operations);

  Stack(Operations&... s) : stack(std::tie(s...)) {}

  // First-order AD
  void forward() { forward_<0>(); }
  void reverse() { reverse_<num_ops - 1>(); }

  // Second-order AD
  void hforward() { hforward_<0>(); }
  void hreverse() { hreverse_<num_ops - 1>(); }

 private:
  StackTuple stack;

  template <index_t index>
  void forward_() {
    std::get<index>(stack).forward();
    if constexpr (index < num_ops - 1) {
      forward_<index + 1>();
    }
  }

  template <index_t index>
  void reverse_() {
    std::get<index>(stack).reverse();
    if constexpr (index) {
      reverse_<index - 1>();
    }
  }

  template <index_t index>
  void hforward_() {
    std::get<index>(stack).hforward();
    if constexpr (index < num_ops - 1) {
      hforward_<index + 1>();
    }
  }

  template <index_t index>
  void hreverse_() {
    std::get<index>(stack).hreverse();
    if constexpr (index) {
      hreverse_<index - 1>();
    }
  }
};

template <class... Operations>
Stack<Operations...> MakeStack(Operations&... s) {
  return Stack<Operations...>(s...);
}

}  // namespace A2D
#endif  // A2D_STACK_H