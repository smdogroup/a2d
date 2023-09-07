#ifndef A2D_STACK_H
#define A2D_STACK_H

#include "a2denum.h"

namespace A2D {

template <class... Operations>
class OperationStack {
 public:
  using StackTuple = std::tuple<Operations...>;
  static constexpr index_t num_ops = sizeof...(Operations);

  KOKKOS_FUNCTION OperationStack(Operations &&...s)
      : stack(std::forward<Operations>(s)...) {
    eval_<0>();
  }

  // First-order AD
  KOKKOS_FUNCTION void forward() { forward_<0>(); }
  KOKKOS_FUNCTION void reverse() { reverse_<num_ops - 1>(); }

  // Second-order AD
  KOKKOS_FUNCTION void hforward() { hforward_<0>(); }
  KOKKOS_FUNCTION void hreverse() { hreverse_<num_ops - 1>(); }

  // Apply Hessian-vector products to extract derivatives
  template <typename T, index_t Ninput, index_t Noutput, class Input,
            class Output, class Intermediate, class Jacobian>
  KOKKOS_FUNCTION void hextract(Intermediate &inter, Input &p, Output &Jp,
                                Jacobian &jac) {
    for (index_t i = 0; i < Ninput; i++) {
      // Zero all the intermeidate values. This inter object must include the
      // input values, and all values included.
      inter.zero();
      p.zero();
      Jp.zero();

      p[i] = T(1.0);

      // Forward sweep
      hforward();

      // Reverse sweep
      hreverse();

      // Extract the number of columns
      for (index_t j = 0; j < Noutput; j++) {
        jac(j, i) = Jp[j];
      }
    }
  }

 private:
  StackTuple stack;

  template <index_t index>
  KOKKOS_FUNCTION void eval_() {
    std::get<index>(stack).eval();
    if constexpr (index < num_ops - 1) {
      eval_<index + 1>();
    }
  }

  template <index_t index>
  KOKKOS_FUNCTION void forward_() {
    std::get<index>(stack).template forward<ADorder::FIRST>();
    if constexpr (index < num_ops - 1) {
      forward_<index + 1>();
    }
  }

  template <index_t index>
  KOKKOS_FUNCTION void reverse_() {
    std::get<index>(stack).reverse();
    if constexpr (index) {
      reverse_<index - 1>();
    }
  }

  template <index_t index>
  KOKKOS_FUNCTION void hforward_() {
    std::get<index>(stack).template forward<ADorder::SECOND>();
    if constexpr (index < num_ops - 1) {
      hforward_<index + 1>();
    }
  }

  template <index_t index>
  KOKKOS_FUNCTION void hreverse_() {
    std::get<index>(stack).hreverse();
    if constexpr (index) {
      hreverse_<index - 1>();
    }
  }
};

template <class... Operations>
KOKKOS_FUNCTION OperationStack<Operations...> MakeStack(Operations &&...s) {
  return OperationStack<Operations...>(std::forward<Operations>(s)...);
}

}  // namespace A2D

#endif  // A2D_STACK_H
