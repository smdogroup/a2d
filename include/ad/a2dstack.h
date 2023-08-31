#ifndef A2D_STACK_H
#define A2D_STACK_H

namespace A2D {

template <class... Operations>
class OperationRefStack {
 public:
  using StackTuple = std::tuple<Operations&...>;
  static constexpr index_t num_ops = sizeof...(Operations);

  OperationRefStack(Operations&... s) : stack(std::tie(s...)) {}

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
    std::get<index>(stack).template forward<ADorder::FIRST>();
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
    std::get<index>(stack).template forward<ADorder::SECOND>();
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
class OperationStack {
 public:
  using StackTuple = std::tuple<Operations...>;
  static constexpr index_t num_ops = sizeof...(Operations);

  OperationStack(Operations&&... s) : stack(std::forward<Operations>(s)...) {}

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
    std::get<index>(stack).template forward<ADorder::FIRST>();
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
    std::get<index>(stack).template forward<ADorder::SECOND>();
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
OperationRefStack<Operations...> MakeStack(Operations&... s) {
  return OperationRefStack<Operations...>(s...);
}

template <class... Operations>
OperationStack<Operations...> MakeStack(Operations&&... s) {
  return OperationStack<Operations...>(std::forward<Operations>(s)...);
}

}  // namespace A2D

#endif  // A2D_STACK_H
