#ifndef A2D_MAT_DET_H
#define A2D_MAT_DET_H

#include "a2ddefs.h"
#include "a2dmat.h"
#include "a2dmatinv.h"
#include "ad/core/a2dmatdetcore.h"

namespace A2D {

template <typename T, int N>
KOKKOS_FUNCTION void MatDet(const Mat<T, N, N>& A, T& det) {
  det = MatDetCore<T, N>(get_data(A));
}

template <class Atype, class dtype>
class MatDetExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<dtype>::type T;

  // Extract the dimensions of the matrix
  static constexpr int N = get_matrix_rows<Atype>::size;
  static constexpr int M = get_matrix_columns<Atype>::size;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<dtype>::order;

  // Assert that the matrix is square
  static_assert(N == M, "Matrix must be square");

  // Make sure that the order is correct
  static_assert(get_diff_order<Atype>::order == order,
                "ADorder does not match");

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

template <class Atype, class dtype>
KOKKOS_FUNCTION auto MatDet(ADObj<Atype>& A, ADObj<dtype>& det) {
  return MatDetExpr<ADObj<Atype>, ADObj<dtype>>(A, det);
}

template <class Atype, class dtype>
KOKKOS_FUNCTION auto MatDet(A2DObj<Atype>& A, A2DObj<dtype>& det) {
  return MatDetExpr<A2DObj<Atype>, A2DObj<dtype>>(A, det);
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
    Mat<T, N, N> A, Ainv;
    x.get_values(A);
    MatDet(A, det);
    MatInv(A, Ainv);
    return MakeVarTuple<T>(det);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<T> det;
    ADObj<Mat<T, N, N>> A, Ainv;

    x.get_values(A.value());
    auto stack = MakeStack(MatDet(A, det), MatInv(A, Ainv));
    seed.get_values(det.bvalue());
    stack.reverse();
    g.set_values(A.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<T> det;
    A2DObj<Mat<T, N, N>> A, Ainv;
    x.get_values(A.value());
    p.get_values(A.pvalue());
    auto stack = MakeStack(MatDet(A, det), MatInv(A, Ainv));
    seed.get_values(det.bvalue());
    hval.get_values(det.hvalue());
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