#ifndef A2D_SHELL_ASSEMBLE_FRAME_H
#define A2D_SHELL_ASSEMBLE_FRAME_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dstack.h"
#include "../a2dtest.h"
#include "../a2dvec.h"

namespace A2D {

template <typename T>
A2D_FUNCTION void ShellAssembleFrameCore(const T Axi[], const T n[],
                                         T frame[]) {
  // Axi usually 3x2 matrix, n is length 3 vec => assembled to 3x3 frame matrix
  frame[0] = Axi[0];
  frame[1] = Axi[1];
  frame[2] = n[0];

  frame[3] = Axi[2];
  frame[4] = Axi[3];
  frame[5] = n[1];

  frame[6] = Axi[4];
  frame[7] = Axi[5];
  frame[8] = n[2];
}

template <typename T>
A2D_FUNCTION void ShellAssembleFrameReverseCore(const T frameb[], T Axib[],
                                                T nb[]) {
  // backprop sensitivities from frame to Axi, n
  Axib[0] += frameb[0];
  Axib[1] += frameb[1];
  nb[0] += frameb[2];

  Axib[2] += frameb[3];
  Axib[3] += frameb[4];
  nb[1] += frameb[5];

  Axib[4] += frameb[6];
  Axib[5] += frameb[7];
  nb[2] += frameb[8];
}

template <typename T>
A2D_FUNCTION void ShellAssembleFrame(const Mat<T, 3, 2> &Axi,
                                     const Vec<T, 3> &n, Mat<T, 3, 3> &frame) {
  ShellAssembleFrameCore<T>(get_data(Axi), get_data(n), get_data(frame));
}

template <class Axitype, class ntype, class frametype>
class ShellAssembleFrameExpr {
 public:
  // extract numeric type to use
  typedef typename get_object_numeric_type<frametype>::type T;

  // Get the sizes of the matrices, vectors
  static constexpr int N = get_matrix_rows<Axitype>::size;
  static constexpr int M = get_matrix_columns<Axitype>::size;
  static constexpr int K = get_vec_size<ntype>::size;
  static constexpr int L = get_matrix_rows<frametype>::size;
  static constexpr int P = get_matrix_columns<frametype>::size;

  // Get the types of the matrices and scalars
  static constexpr ADiffType adAxi = get_diff_type<Axitype>::diff_type;
  static constexpr ADiffType adn = get_diff_type<ntype>::diff_type;
  static constexpr ADiffType adframe = get_diff_type<frametype>::diff_type;

  // assert all are same type
  //   static_assert(((get_a2d_object_type<Axitype>::value ==
  //               get_a2d_object_type<ntype>::value) &&
  //              (get_a2d_object_type<ntype>::value ==
  //               get_a2d_object_type<frametype>::value)),
  //             "Inputs are not all of the same type");

  // assert correct sizes
  static_assert((N == 3) && (M == 2) && (K == 3) && (L == 3) && (P == 3),
                "AssembleFrame (3,2);(3) => (3,3)");

  // get the differentiation order
  static constexpr ADorder order = get_diff_order<frametype>::order;

  A2D_FUNCTION ShellAssembleFrameExpr(Axitype &Axi, ntype &n, frametype &frame)
      : Axi(Axi), n(n), frame(frame) {}

  A2D_FUNCTION void eval() {
    ShellAssembleFrameCore<T>(get_data(Axi), get_data(n), get_data(frame));
  }

  A2D_FUNCTION void bzero() { frame.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    if constexpr (adAxi == ADiffType::ACTIVE && adn == ADiffType::ACTIVE) {
      ShellAssembleFrameCore<T>(GetSeed<seed>::get_data(Axi),
                                GetSeed<seed>::get_data(n),
                                GetSeed<seed>::get_data(frame));
    } else if constexpr (adn == ADiffType::ACTIVE) {
      Mat<T, N, M> Axi_void;
      ShellAssembleFrameCore<T>(get_data(Axi_void), GetSeed<seed>::get_data(n),
                                GetSeed<seed>::get_data(frame));
    } else if constexpr (adAxi == ADiffType::ACTIVE) {
      Vec<T, K> n_void;
      ShellAssembleFrameCore<T>(GetSeed<seed>::get_data(Axi), get_data(n_void),
                                GetSeed<seed>::get_data(frame));
    }
  }

  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    if constexpr (adn == ADiffType::ACTIVE) {
      Mat<T, N, M> Axi_void;
      ShellAssembleFrameReverseCore<T>(GetSeed<seed>::get_data(frame),
                                       get_data(Axi_void),
                                       GetSeed<seed>::get_data(n));
    }
    if constexpr (adAxi == ADiffType::ACTIVE) {
      Vec<T, K> n_void;
      ShellAssembleFrameReverseCore<T>(GetSeed<seed>::get_data(frame),
                                       GetSeed<seed>::get_data(Axi),
                                       get_data(n_void));
    }
  }

  A2D_FUNCTION void hzero() { frame.hzero(); }

  A2D_FUNCTION void hreverse() {
    constexpr ADseed seed = ADseed::h;
    if constexpr (adn == ADiffType::ACTIVE) {
      Mat<T, N, M> Axi_void;
      ShellAssembleFrameReverseCore<T>(GetSeed<seed>::get_data(frame),
                                       get_data(Axi_void),
                                       GetSeed<seed>::get_data(n));
    }
    if constexpr (adAxi == ADiffType::ACTIVE) {
      Vec<T, K> n_void;
      ShellAssembleFrameReverseCore<T>(GetSeed<seed>::get_data(frame),
                                       GetSeed<seed>::get_data(Axi),
                                       get_data(n_void));
    }
  }

  Axitype &Axi;
  ntype &n;
  frametype &frame;
};  // end of ShellAssembleFrameExpr class definition

// Full active variants
template <class Axitype, class ntype, class frametype>
A2D_FUNCTION auto ShellAssembleFrame(ADObj<Axitype> &Axi, ADObj<ntype> &n,
                                     ADObj<frametype> &frame) {
  return ShellAssembleFrameExpr<ADObj<Axitype>, ADObj<ntype>, ADObj<frametype>>(
      Axi, n, frame);
}

template <class Axitype, class ntype, class frametype>
A2D_FUNCTION auto ShellAssembleFrame(Axitype &Axi, ntype &n,
                                     ADObj<frametype> &frame) {
  return ShellAssembleFrameExpr<Axitype, ntype, ADObj<frametype>>(Axi, n,
                                                                  frame);
}

template <class Axitype, class ntype, class frametype>
A2D_FUNCTION auto ShellAssembleFrame(A2DObj<Axitype> &Axi, A2DObj<ntype> &n,
                                     A2DObj<frametype> &frame) {
  return ShellAssembleFrameExpr<A2DObj<Axitype>, A2DObj<ntype>,
                                A2DObj<frametype>>(Axi, n, frame);
}

template <class Axitype, class ntype, class frametype>
A2D_FUNCTION auto ShellAssembleFrame(Axitype &Axi, ntype &n,
                                     A2DObj<frametype> &frame) {
  return ShellAssembleFrameExpr<Axitype, ntype, A2DObj<frametype>>(Axi, n,
                                                                   frame);
}

namespace Test {

template <typename T>
class ShellAssembleFrameTest
    : public A2DTest<T, Mat<T, 3, 3>, Mat<T, 3, 2>, Vec<T, 3>> {
 public:
  using Input = VarTuple<T, Mat<T, 3, 2>, Vec<T, 3>>;
  using Output = VarTuple<T, Mat<T, 3, 3>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "ShellAssembleFrame>";
    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input &x) {
    Mat<T, 3, 2> Axi;
    Vec<T, 3> n;
    Mat<T, 3, 3> frame;

    x.get_values(Axi, n);
    ShellAssembleFrame(Axi, n, frame);
    return MakeVarTuple<T>(frame);
  }

  // Compute the derivative
  void deriv(const Output &seed, const Input &x, Input &g) {
    ADObj<Mat<T, 3, 2>> Axi;
    ADObj<Vec<T, 3>> n;
    ADObj<Mat<T, 3, 3>> frame;

    x.get_values(Axi.value(), n.value());
    auto stack = MakeStack(ShellAssembleFrame(Axi, n, frame));
    seed.get_values(frame.bvalue());
    stack.reverse();
    g.set_values(Axi.bvalue(), n.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output &seed, const Output &hval, const Input &x,
             const Input &p, Input &h) {
    A2DObj<Mat<T, 3, 2>> Axi;
    A2DObj<Vec<T, 3>> n;
    A2DObj<Mat<T, 3, 3>> frame;

    x.get_values(Axi.value(), n.value());
    p.get_values(Axi.pvalue(), n.pvalue());
    auto stack = MakeStack(ShellAssembleFrame(Axi, n, frame));
    seed.get_values(frame.bvalue());
    hval.get_values(frame.hvalue());
    stack.hproduct();
    h.set_values(Axi.hvalue(), n.hvalue());
  }
};

bool ShellAssembleFrameTestAll(bool component = false,
                               bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  ShellAssembleFrameTest<Tc> test;
  bool passed = Run(test, component, write_output);
  return passed;
}

}  // namespace Test

}  // namespace A2D
#endif  // A2D_SHELL_ASSEMBLE_FRAME_H
