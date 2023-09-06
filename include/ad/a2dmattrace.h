#ifndef A2D_MAT_TRACE_H
#define A2D_MAT_TRACE_H

#include <type_traits>

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T, int M>
KOKKOS_FUNCTION T MatTraceCore(const T* A) {
  T trace = T(0.0);
  for (int i = 0; i < M; i++) {
    trace += A[0];
    A += M + 1;
  }
  return trace;
}

template <typename T, int M>
KOKKOS_FUNCTION T SymMatTraceCore(const T* S) {
  T trace = T(0.0);
  for (int i = 0; i < M; i++, S++) {
    S += i;
    trace += S[0];
  }
  return trace;
}

template <typename T, int M>
KOKKOS_FUNCTION void MatAddDiagCore(const T diag, T* A) {
  for (int i = 0; i < M; i++) {
    A[0] += diag;
    A += M + 1;
  }
}

template <typename T, int M>
KOKKOS_FUNCTION void SymMatAddDiagCore(const T diag, T* S) {
  for (int i = 0; i < M; i++, S++) {
    S += i;
    S[0] += diag;
  }
}

template <typename T, int M>
KOKKOS_FUNCTION void MatTrace(Mat<T, M, M>& A, T& trace) {
  trace = MatTraceCore<T, M>(get_data(A));
}

template <typename T, int M>
KOKKOS_FUNCTION void MatTrace(SymMat<T, M>& S, T& trace) {
  trace = SymMatTraceCore<T, M>(get_data(S));
}

template <typename T, int M, ADorder order, ADiffType adA>
class MatTraceExpr {
 private:
  using Atype = ADMatType<adA, order, Mat<T, M, M>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  KOKKOS_FUNCTION MatTraceExpr(Atype& A, dtype& tr) : A(A), tr(tr) {}

  KOKKOS_FUNCTION void eval() {
    get_data(tr) = MatTraceCore<T, M>(get_data(A));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(tr) =
          MatTraceCore<T, M>(GetSeed<seed>::get_data(A));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      MatAddDiagCore<T, M>(GetSeed<ADseed::b>::get_data(tr),
                           GetSeed<ADseed::b>::get_data(A));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      MatAddDiagCore<T, M>(GetSeed<ADseed::h>::get_data(tr),
                           GetSeed<ADseed::h>::get_data(A));
    }
  }

 private:
  Atype& A;
  dtype& tr;
};

template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(ADObj<Mat<T, M, M>>& A, ADObj<T>& tr) {
  return MatTraceExpr<T, M, ADorder::FIRST, ADiffType::ACTIVE>(A, tr);
}

template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(A2DObj<Mat<T, M, M>>& A, A2DObj<T>& tr) {
  return MatTraceExpr<T, M, ADorder::SECOND, ADiffType::ACTIVE>(A, tr);
}

template <typename T, int M, ADorder order, ADiffType adA>
class SymTraceExpr {
 private:
  using Stype = ADMatType<adA, order, SymMat<T, M>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  KOKKOS_FUNCTION SymTraceExpr(Stype& S, dtype& tr) : S(S), tr(tr) {}

  KOKKOS_FUNCTION void eval() {
    get_data(tr) = SymMatTraceCore<T, M>(get_data(S));
  }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    if constexpr (adA == ADiffType::ACTIVE) {
      GetSeed<seed>::get_data(tr) =
          SymMatTraceCore<T, M>(GetSeed<seed>::get_data(S));
    }
  }

  KOKKOS_FUNCTION void reverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      SymMatAddDiagCore<T, M>(GetSeed<ADseed::b>::get_data(tr),
                              GetSeed<ADseed::b>::get_data(S));
    }
  }

  KOKKOS_FUNCTION void hreverse() {
    if constexpr (adA == ADiffType::ACTIVE) {
      SymMatAddDiagCore<T, M>(GetSeed<ADseed::h>::get_data(tr),
                              GetSeed<ADseed::h>::get_data(S));
    }
  }

 private:
  Stype& S;
  dtype& tr;
};

template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(ADObj<SymMat<T, M>>& S, ADObj<T>& tr) {
  return SymTraceExpr<T, M, ADorder::FIRST, ADiffType::ACTIVE>(S, tr);
}

template <typename T, int M>
KOKKOS_FUNCTION auto MatTrace(A2DObj<SymMat<T, M>>& S, A2DObj<T>& tr) {
  return SymTraceExpr<T, M, ADorder::SECOND, ADiffType::ACTIVE>(S, tr);
}

namespace Test {

template <typename T, int N>
class MatTraceTest : public A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatTrace<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T trace;
    Mat<T, N, N> A;
    x.get_values(A);
    MatTrace(A, trace);
    return MakeVarTuple<T>(trace);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> trace;
    ADObj<Mat<T, N, N>> A;

    x.get_values(A.value());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> trace;
    A2DObj<Mat<T, N, N>> A;
    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    hval.get_values(trace.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

template <typename T, int N>
class SymTraceTest : public A2DTest<T, T, SymMat<T, N>> {
 public:
  using Input = VarTuple<T, SymMat<T, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatTrace<" << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T trace;
    SymMat<T, N> A;
    x.get_values(A);
    MatTrace(A, trace);
    return MakeVarTuple<T>(trace);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> trace;
    ADObj<SymMat<T, N>> A;

    x.get_values(A.value());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> trace;
    A2DObj<SymMat<T, N>> A;

    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatTrace(A, trace));
    seed.get_values(trace.bvalue());
    hval.get_values(trace.hvalue());
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

bool MatTraceTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatTraceTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatTraceTest<Tc, 4> test2;
  passed = passed && Run(test2, component, write_output);

  SymTraceTest<Tc, 4> test3;
  passed = passed && Run(test3, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_MAT_TRACE_H