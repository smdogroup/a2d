#ifndef A2D_STACK_H
#define A2D_STACK_H

#ifdef __CUDACC__
#include <cuda/std/tuple>
#else
#include <tuple>
#endif

#include "a2ddefs.h"

namespace A2D {

template <class... Operations>
class OperationStack {
 public:
  static constexpr int num_ops = sizeof...(Operations);
#ifdef __CUDACC__
  using StackTuple = cuda::std::tuple<Operations...>;
#else
  using StackTuple = std::tuple<Operations...>;
#endif

  template <class... Args>
  A2D_FUNCTION OperationStack(Args&&... s) : stack(std::forward<Args>(s)...) {
    eval_<0>();
  }

  // First-order AD
  A2D_FUNCTION void bzero() { bzero_<0>(); }
  A2D_FUNCTION void forward() { forward_<0>(); }
  A2D_FUNCTION void reverse() { reverse_<num_ops - 1>(); }

  // Second-order AD
  A2D_FUNCTION void hzero() { hzero_<0>(); }
  A2D_FUNCTION void hforward() { hforward_<0>(); }
  A2D_FUNCTION void hreverse() { hreverse_<num_ops - 1>(); }

  // Perform a Hessian-vector product
  A2D_FUNCTION void hproduct() {
    reverse();
    hforward();
    hreverse();
  }

  // Apply Hessian-vector products to extract derivatives
  template <class Input, class Output, class Jacobian>
  A2D_FUNCTION void hextract(Input& p, Output& Jp, Jacobian& jac) {
    reverse();

    for (int i = 0; i < Input::ncomp; i++) {
      // Zero all the intermeidate values. This inter object must include the
      // input values, and all values included.
      p.zero();
      Jp.zero();
      hzero();

      p[i] = 1.0;

      // Forward sweep
      hforward();

      // Reverse sweep
      hreverse();

      // Extract the number of columns
      for (int j = 0; j < Output::ncomp; j++) {
        jac(j, i) = Jp[j];
      }
    }
  }

 private:
  StackTuple stack;

  template <int index>
  A2D_FUNCTION inline decltype(auto) get(StackTuple& stack) {
#ifdef __CUDACC__
    return cuda::std::get<index>(stack);
#else
    return std::get<index>(stack);
#endif
  }

  template <int index>
  A2D_FUNCTION void eval_() {
    get<index>(stack).eval();
    if constexpr (index < num_ops - 1) {
      eval_<index + 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void bzero_() {
    get<index>(stack).bzero();
    if constexpr (index < num_ops - 1) {
      bzero_<index + 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void forward_() {
    get<index>(stack).template forward<ADorder::FIRST>();
    if constexpr (index < num_ops - 1) {
      forward_<index + 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void reverse_() {
    get<index>(stack).reverse();
    if constexpr (index) {
      reverse_<index - 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void hzero_() {
    get<index>(stack).hzero();
    if constexpr (index < num_ops - 1) {
      hzero_<index + 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void hforward_() {
    get<index>(stack).template forward<ADorder::SECOND>();
    if constexpr (index < num_ops - 1) {
      hforward_<index + 1>();
    }
  }

  template <int index>
  A2D_FUNCTION void hreverse_() {
    get<index>(stack).hreverse();
    if constexpr (index) {
      hreverse_<index - 1>();
    }
  }
};

/**
 * @brief Make an operations stack for automatic differentiation
 *
 * Operations are evaluated immediately on construction
 *
 * @tparam Operations Template parameter list deduced from context
 * @param s The operator objects
 * @return The list of operations
 */
template <class... Operations>
A2D_FUNCTION auto MakeStack(Operations&&... s) {
  return OperationStack<Operations...>(std::forward<Operations>(s)...);
}

}  // namespace A2D

#endif  // A2D_STACK_H
