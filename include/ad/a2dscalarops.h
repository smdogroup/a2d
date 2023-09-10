#ifndef A2D_SCALAR_OPS_H
#define A2D_SCALAR_OPS_H

#include "a2dbinary.h"
#include "a2denum.h"
#include "a2dobjs.h"
#include "a2dunary.h"

namespace A2D {

template <class Expr, class T>
class EvalExpr {
 public:
  KOKKOS_FUNCTION EvalExpr(Expr&& expr, ADObj<T>& out)
      : expr(std::forward<Expr>(expr)), out(out) {}

  KOKKOS_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(forder == ADorder::FIRST,
                  "EvalExpr only works for first-order AD");
    expr.forward();
    out.bvalue() = expr.bvalue();
  }
  KOKKOS_FUNCTION void reverse() {
    expr.bvalue() = out.bvalue();
    expr.reverse();
  }

 private:
  Expr expr;
  ADObj<T>& out;
};

template <class Expr, class T>
auto Eval(Expr&& expr, ADObj<T>& out) {
  return EvalExpr<Expr, T>(std::forward<Expr>(expr), out);
}

template <class Expr, class T>
class EvalExpr2 {
 public:
  KOKKOS_FUNCTION EvalExpr2(Expr&& expr, A2DObj<T>& out)
      : expr(std::forward<Expr>(expr)), out(out) {}

  KOKKOS_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  KOKKOS_FUNCTION void reverse() {
    expr.bvalue() += out.bvalue();
    expr.reverse();
  }
  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(forder == ADorder::SECOND,
                  "EvalExpr2 only works for second-order AD");
    expr.hforward();
    out.pvalue() = expr.pvalue();
  }
  KOKKOS_FUNCTION void hreverse() {
    expr.hvalue() += out.hvalue();
    expr.hreverse();
  }

 private:
  Expr expr;
  A2DObj<T>& out;
};

template <class Expr, class T>
auto Eval(Expr&& expr, A2DObj<T>& out) {
  return EvalExpr2<Expr, T>(std::forward<Expr>(expr), out);
}

namespace Test {
template <typename T>
class ScalarTest : public A2DTest<T, T, T, T> {
 public:
  using Input = VarTuple<T, T, T>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() { return std::string("ScalarOperations"); }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T a, b, f;
    x.get_values(a, b);

    f = log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a) - a / b;

    return MakeVarTuple<T>(f);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> a, b, f;
    x.get_values(a.value(), b.value());

    auto stack = MakeStack(Eval(
        log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a) - a / b, f));

    seed.get_values(f.bvalue());
    stack.reverse();
    g.set_values(a.bvalue(), b.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> a, b, f;
    x.get_values(a.value(), b.value());
    p.get_values(a.pvalue(), b.pvalue());

    auto stack = MakeStack(Eval(
        log(a * a * sqrt(exp(a * sin(a) + a)) + a * a * a * a) - a / b, f));

    seed.get_values(f.bvalue());
    hval.get_values(f.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(a.hvalue(), b.hvalue());
  }
};

bool ScalarTestAll(bool component, bool write_output) {
  using Tc = std::complex<double>;

  bool passed = true;
  ScalarTest<Tc> test1;
  passed = passed && Run(test1, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SCALAR_OPS_H