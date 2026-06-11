#ifndef A2D_SCALAR_OPS_H
#define A2D_SCALAR_OPS_H

#include "../a2ddefs.h"
#include "a2dbinary.h"
#include "a2dtest.h"
#include "a2dunary.h"

namespace A2D {

template <class Expr, class T>
class EvalExpr {
 public:
  A2D_FUNCTION EvalExpr(Expr&& expr, ADObj<T>& out)
      : expr(a2d_forward<Expr>(expr)), out(out) {}

  A2D_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  A2D_FUNCTION void bzero() {
    out.bzero();
    expr.bzero();
  }
  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(forder == ADorder::FIRST,
                  "EvalExpr only works for first-order AD");
    expr.forward();
    out.bvalue() = expr.bvalue();
  }
  A2D_FUNCTION void reverse() {
    expr.bvalue() = out.bvalue();
    expr.reverse();
  }

 private:
  Expr expr;
  ADObj<T>& out;
};

template <class Expr, class T>
class EvalExprRef {
 public:
  A2D_FUNCTION EvalExprRef(Expr&& expr, ADObj<T&> out)
      : expr(a2d_forward<Expr>(expr)), out(out) {}

  A2D_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  A2D_FUNCTION void bzero() {
    out.bzero();
    expr.bzero();
  }
  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(forder == ADorder::FIRST,
                  "EvalExprRef only works for first-order AD");
    expr.forward();
    out.bvalue() = expr.bvalue();
  }
  A2D_FUNCTION void reverse() {
    expr.bvalue() = out.bvalue();
    expr.reverse();
  }

 private:
  Expr expr;
  ADObj<T&> out;
};

template <class Expr, class T,
          std::enable_if_t<!std::is_reference<T>::value, bool> = true>
A2D_FUNCTION auto Eval(Expr&& expr, ADObj<T>& out) {
  return EvalExpr<Expr, T>(a2d_forward<Expr>(expr), out);
}

template <class Expr, class T,
          std::enable_if_t<!std::is_reference<T>::value, bool> = true>
A2D_FUNCTION auto Eval(Expr&& expr, ADObj<T&> out) {
  return EvalExprRef<Expr, T>(a2d_forward<Expr>(expr), out);
}

template <class Expr, class T>
class EvalExpr2 {
 public:
  A2D_FUNCTION EvalExpr2(Expr&& expr, A2DObj<T>& out)
      : expr(a2d_forward<Expr>(expr)), out(out) {}

  A2D_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  A2D_FUNCTION void bzero() {
    out.bzero();
    expr.bzero();
  }
  A2D_FUNCTION void reverse() {
    expr.bvalue() += out.bvalue();
    expr.reverse();
  }
  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(forder == ADorder::SECOND,
                  "EvalExpr2 only works for second-order AD");
    expr.hforward();
    out.pvalue() = expr.pvalue();
  }
  A2D_FUNCTION void hzero() {
    out.hzero();
    expr.hzero();
  }
  A2D_FUNCTION void hreverse() {
    expr.hvalue() += out.hvalue();
    expr.hreverse();
  }

 private:
  Expr expr;
  A2DObj<T>& out;
};

template <class Expr, class T>
class EvalExprRef2 {
 public:
  A2D_FUNCTION EvalExprRef2(Expr&& expr, A2DObj<T&> out)
      : expr(a2d_forward<Expr>(expr)), out(out) {}

  A2D_FUNCTION void eval() {
    expr.eval();
    out.value() = expr.value();
  }
  A2D_FUNCTION void bzero() {
    out.bzero();
    expr.bzero();
  }
  A2D_FUNCTION void reverse() {
    expr.bvalue() += out.bvalue();
    expr.reverse();
  }
  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(forder == ADorder::SECOND,
                  "EvalExprRef2 only works for second-order AD");
    expr.hforward();
    out.pvalue() = expr.pvalue();
  }
  A2D_FUNCTION void hzero() {
    out.hzero();
    expr.hzero();
  }
  A2D_FUNCTION void hreverse() {
    expr.hvalue() += out.hvalue();
    expr.hreverse();
  }

 private:
  Expr expr;
  A2DObj<T&> out;
};

template <class Expr, class T,
          std::enable_if_t<!std::is_reference<T>::value, bool> = true>
A2D_FUNCTION auto Eval(Expr&& expr, A2DObj<T>& out) {
  return EvalExpr2<Expr, T>(a2d_forward<Expr>(expr), out);
}

template <class Expr, class T,
          std::enable_if_t<!std::is_reference<T>::value, bool> = true>
A2D_FUNCTION auto Eval(Expr&& expr, A2DObj<T&> out) {
  return EvalExprRef2<Expr, T>(a2d_forward<Expr>(expr), out);
}

}  // namespace A2D

#endif  // A2D_SCALAR_OPS_H
