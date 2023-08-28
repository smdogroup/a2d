#ifndef A2D_SYM_MAT_TRACE_H
#define A2D_SYM_MAT_TRACE_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "a2dstack.h"
#include "a2dtest.h"
#include "ad/core/a2dsymtracecore.h"

namespace A2D {

template <typename T, int N>
KOKKOS_FUNCTION void SymMatTrace(SymMat<T, N>& S, SymMat<T, N>& E, T& output) {
  output = SymMatTraceCore<T, N>(get_data(S), get_data(E));
}

template <typename T, int N, ADorder order>
class SymMatTraceExpr {
 public:
  using Stype = ADMatType<ADiffType::ACTIVE, order, SymMat<T, N>>;
  using Scalar = ADScalarType<ADiffType::ACTIVE, order, T>;

  KOKKOS_FUNCTION SymMatTraceExpr(Stype& S, Stype& E, Scalar& out)
      : S(S), E(E), out(out) {
    get_data(out) = SymMatTraceCore<T, N>(get_data(S), get_data(E));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    GetSeed<seed>::get_data(out) =
        SymMatTraceCore<T, N>(GetSeed<seed>::get_data(S), get_data(E)) +
        SymMatTraceCore<T, N>(get_data(S), GetSeed<seed>::get_data(E));
  }

  KOKKOS_FUNCTION void reverse() {
    SymMatTraceReverseCore<T, N>(out.bvalue, get_data(S),
                                 GetSeed<ADseed::b>::get_data(E));
    SymMatTraceReverseCore<T, N>(out.bvalue, get_data(E),
                                 GetSeed<ADseed::b>::get_data(S));
  }

  KOKKOS_FUNCTION void hreverse() {
    SymMatTraceReverseCore<T, N>(out.bvalue, GetSeed<ADseed::p>::get_data(S),
                                 GetSeed<ADseed::h>::get_data(E));
    SymMatTraceReverseCore<T, N>(out.hvalue, get_data(S),
                                 GetSeed<ADseed::h>::get_data(E));

    SymMatTraceReverseCore<T, N>(out.bvalue, GetSeed<ADseed::p>::get_data(E),
                                 GetSeed<ADseed::h>::get_data(S));
    SymMatTraceReverseCore<T, N>(out.hvalue, get_data(E),
                                 GetSeed<ADseed::h>::get_data(S));
  }

  Stype &S, &E;
  Scalar& out;
};

template <typename T, int N>
KOKKOS_FUNCTION auto SymMatTrace(ADMat<SymMat<T, N>>& S, ADMat<SymMat<T, N>>& E,
                                 ADScalar<T>& output) {
  return SymMatTraceExpr<T, N, ADorder::FIRST>(S, E, output);
}

template <typename T, int N>
KOKKOS_FUNCTION auto SymMatTrace(A2DMat<SymMat<T, N>>& S,
                                 A2DMat<SymMat<T, N>>& E,
                                 A2DScalar<T>& output) {
  return SymMatTraceExpr<T, N, ADorder::SECOND>(S, E, output);
}

namespace Test {

template <typename T, int N>
class SymMatTraceTest : public A2DTest<T, T, SymMat<T, N>, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "SymMatTrace<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T output;
    SymMat<T, N> S, E;
    x.get_values(S, E);
    SymMatTrace(S, E, output);
    return MakeVarTuple<T>(output);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADScalar<T> output;
    SymMat<T, N> S0, Sb, E0, Eb;
    ADMat<SymMat<T, N>> S(S0, Sb), E(E0, Eb);

    x.get_values(S0, E0);
    auto op = SymMatTrace(S, E, output);
    auto stack = MakeStack(op);
    seed.get_values(output.bvalue);
    stack.reverse();
    g.set_values(Sb, Eb);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DScalar<T> output;
    A2DMat<SymMat<T, N>> S, E;
    x.get_values(S.value(), E.value());
    p.get_values(S.pvalue(), E.pvalue());

    auto op = SymMatTrace(S, E, output);
    auto stack = MakeStack(op);

    seed.get_values(output.bvalue);
    hval.get_values(output.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(S.hvalue(), E.hvalue());
  }
};

bool SymMatTraceTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  SymMatTraceTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  SymMatTraceTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);
  SymMatTraceTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SYM_MAT_TRACE_H