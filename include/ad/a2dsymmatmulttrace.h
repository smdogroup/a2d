#ifndef A2D_SYM_MAT_TRACE_H
#define A2D_SYM_MAT_TRACE_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "ad/core/a2dsymmatmulttracecore.h"

namespace A2D {

template <typename T, int N>
KOKKOS_FUNCTION void SymMatMultTrace(SymMat<T, N>& S, SymMat<T, N>& E,
                                     T& output) {
  output = SymMatMultTraceCore<T, N>(get_data(S), get_data(E));
}

template <typename T, int N, ADorder order>
class SymMatMultTraceExpr {
 public:
  using Stype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION SymMatMultTraceExpr(Stype& S, Stype& E, dtype& out)
      : S(S), E(E), out(out) {}

  KOKKOS_FUNCTION void eval() {
    get_data(out) = SymMatMultTraceCore<T, N>(get_data(S), get_data(E));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    GetSeed<seed>::get_data(out) =
        SymMatMultTraceCore<T, N>(GetSeed<seed>::get_data(S), get_data(E)) +
        SymMatMultTraceCore<T, N>(get_data(S), GetSeed<seed>::get_data(E));
  }

  KOKKOS_FUNCTION void reverse() {
    SymMatMultTraceReverseCore<T, N>(out.bvalue(), get_data(S),
                                     GetSeed<ADseed::b>::get_data(E));
    SymMatMultTraceReverseCore<T, N>(out.bvalue(), get_data(E),
                                     GetSeed<ADseed::b>::get_data(S));
  }

  KOKKOS_FUNCTION void hreverse() {
    SymMatMultTraceReverseCore<T, N>(out.bvalue(),
                                     GetSeed<ADseed::p>::get_data(S),
                                     GetSeed<ADseed::h>::get_data(E));
    SymMatMultTraceReverseCore<T, N>(out.hvalue(), get_data(S),
                                     GetSeed<ADseed::h>::get_data(E));

    SymMatMultTraceReverseCore<T, N>(out.bvalue(),
                                     GetSeed<ADseed::p>::get_data(E),
                                     GetSeed<ADseed::h>::get_data(S));
    SymMatMultTraceReverseCore<T, N>(out.hvalue(), get_data(E),
                                     GetSeed<ADseed::h>::get_data(S));
  }

  Stype &S, &E;
  dtype& out;
};

template <typename T, int N>
KOKKOS_FUNCTION auto SymMatMultTrace(ADObj<SymMat<T, N>>& S,
                                     ADObj<SymMat<T, N>>& E, ADObj<T>& output) {
  return SymMatMultTraceExpr<T, N, ADorder::FIRST>(S, E, output);
}

template <typename T, int N>
KOKKOS_FUNCTION auto SymMatMultTrace(A2DObj<SymMat<T, N>>& S,
                                     A2DObj<SymMat<T, N>>& E,
                                     A2DObj<T>& output) {
  return SymMatMultTraceExpr<T, N, ADorder::SECOND>(S, E, output);
}

namespace Test {

template <typename T, int N>
class SymMatMultTraceTest : public A2DTest<T, T, SymMat<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatMultTrace<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T output;
    SymMat<T, N> S, E;
    x.get_values(S, E);
    SymMatMultTrace(S, E, output);
    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> output;
    ADObj<SymMat<T, N>> S, E;

    x.get_values(S.value(), E.value());
    auto stack = MakeStack(SymMatMultTrace(S, E, output));
    seed.get_values(output.bvalue());
    stack.reverse();
    g.set_values(S.bvalue(), E.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> output;
    A2DObj<SymMat<T, N>> S, E;

    x.get_values(S.value(), E.value());
    p.get_values(S.pvalue(), E.pvalue());
    auto stack = MakeStack(SymMatMultTrace(S, E, output));
    seed.get_values(output.bvalue());
    hval.get_values(output.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(S.hvalue(), E.hvalue());
  }
};

bool SymMatMultTraceTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  SymMatMultTraceTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatMultTraceTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  SymMatMultTraceTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYM_MAT_TRACE_H