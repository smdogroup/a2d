#ifndef A2D_MAT_DET_H
#define A2D_MAT_DET_H

#include "a2denum.h"
#include "a2dmat.h"
#include "a2dobjs.h"
#include "ad/core/a2dmatdetcore.h"

namespace A2D {

template <typename T, int N>
KOKKOS_FUNCTION void MatDet(Mat<T, N, N>& A, T& det) {
  det = MatDetCore<T, N>(get_data(A));
}

template <typename T, int N, ADorder order>
class MatDetExpr {
 private:
  using Atype = ADMatType<ADiffType::ACTIVE, order, Mat<T, N, N>>;
  using dtype = ADScalarType<ADiffType::ACTIVE, order, T>;

 public:
  KOKKOS_FUNCTION MatDetExpr(Atype& A, dtype& det) : A(A), det(det) {}

  KOKKOS_FUNCTION void eval() { get_data(det) = MatDetCore<T, N>(get_data(A)); }

  template <ADorder forder>
  KOKKOS_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    GetSeed<seed>::get_data(det) =
        MatDetForwardCore<T, N>(get_data(A), GetSeed<seed>::get_data(A));
  }

  KOKKOS_FUNCTION void reverse() {
    MatDetReverseCore<T, N>(GetSeed<ADseed::b>::get_data(det), get_data(A),
                            GetSeed<ADseed::b>::get_data(A));
  }

  KOKKOS_FUNCTION void hreverse() {
    static_assert(order == ADorder::SECOND,
                  "hreverse() can be called for only second order objects.");

    MatDetHReverseCore<T, N>(GetSeed<ADseed::b>::get_data(det),
                             GetSeed<ADseed::h>::get_data(det), get_data(A),
                             GetSeed<ADseed::p>::get_data(A),
                             GetSeed<ADseed::h>::get_data(A));
  }

  Atype& A;
  dtype& det;
};

template <typename T, int N>
KOKKOS_FUNCTION auto MatDet(ADMat<Mat<T, N, N>>& A, ADScalar<T>& det) {
  return MatDetExpr<T, N, ADorder::FIRST>(A, det);
}

template <typename T, int N>
KOKKOS_FUNCTION auto MatDet(A2DMat<Mat<T, N, N>>& A, A2DScalar<T>& det) {
  return MatDetExpr<T, N, ADorder::SECOND>(A, det);
}

namespace Test {

template <typename T, int N>
class MatDetTest : public A2DTest<T, T, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, T>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatDet<" << N << "," << N << ">";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    T det;
    Mat<T, N, N> A;
    x.get_values(A);
    MatDet(A, det);
    return MakeVarTuple<T>(det);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADScalar<T> det;
    Mat<T, N, N> A0, Ab;
    ADMat<Mat<T, N, N>> A(A0, Ab);

    x.get_values(A0);
    auto op = MatDet(A, det);
    auto stack = MakeStack(op);
    seed.get_values(det.bvalue);
    stack.reverse();
    g.set_values(Ab);
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DScalar<T> det;
    A2DMat<Mat<T, N, N>> A;
    x.get_values(A.value());
    p.get_values(A.pvalue());

    auto op = MatDet(A, det);
    auto stack = MakeStack(op);

    seed.get_values(det.bvalue);
    hval.get_values(det.hvalue);
    stack.reverse();
    stack.hforward();
    stack.hreverse();
    h.set_values(A.hvalue());
  }
};

bool MatDetTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  bool passed = true;
  MatDetTest<Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatDetTest<Tc, 3> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  //  A2D_MAT_DET_H