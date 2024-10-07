#ifndef A2D_SHELL_STRAIN_H
#define A2D_SHELL_STRAIN_H

#include <type_traits>

#include "../../a2ddefs.h"
#include "../a2dmat.h"
#include "../a2dstack.h"
#include "../a2dtest.h"
// #include "../core/a2dgreenstraincore.h"

namespace A2D {

enum class ShellStrainType { LINEAR, NONLINEAR };

template <typename T>
A2D_FUNCTION void LinearShellStrainCore(const T u0x[], const T u1x[],
                                const T e0ty[], const T et[], T e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0]; // e11
    e[1] = e0ty[3]; // e22
    e[2] = 2.0 * e0ty[1]; // e12

    // Compute the bending strain
    e[3] = u1x[0]; // k11
    e[4] = u1x[4]; // k22
    e[5] = u1x[1] + u1x[3]; // k12

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4]; // e23, transverse shear
    e[7] = 2.0 * e0ty[2]; // e13, transverse shear
    e[8] = et[0]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void LinearShellStrainForwardCore(const T u0x[], const T u1x[],
                                const T e0ty[], const T et[], T e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0]; // e11
    e[1] = e0ty[3]; // e22
    e[2] = 2.0 * e0ty[1]; // e12

    // Compute the bending strain
    e[3] = u1x[0]; // k11
    e[4] = u1x[4]; // k22
    e[5] = u1x[1] + u1x[3]; // k12

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4]; // e23, transverse shear
    e[7] = 2.0 * e0ty[2]; // e13, transverse shear
    e[8] = et[0]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void LinearShellStrainReverseCore(const T eb[],
    T u0xb[], T u1xb[], T e0tyb[], T etb[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e0tyb[0] += eb[0]; // e1
    e0tyb[3] += eb[1]; // e22
    e0tyb[1] += 2.0 * eb[2]; // e12

    // Compute the bending strain
    u1xb[0] += eb[3]; // k11
    u1xb[4] += eb[4]; // k22
    u1xb[1] += eb[5]; // k12
    u1xb[3] += eb[5]; // k12

    // Add the components of the shear strain
    e0tyb[4] += 2.0 * eb[6]; // e23, transverse shear
    e0tyb[2] += 2.0 * eb[7]; // e13, transverse shear
    etb[0] += eb[8]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void NonlinearShellStrainCore(const T u0x[], const T u1x[],
                                const T e0ty[], const T et[], T e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0]; // e11
    e[1] = e0ty[3]; // e22
    e[2] = 2.0 * e0ty[1]; // e12

    // Compute the bending strain (here's where nonlinearity comes in)
    e[3] = u1x[0] + (u0x[0] * u1x[0] + u0x[3] * u1x[3] + u0x[6] * u1x[6]); // k11
    e[4] = u1x[4] + (u0x[1] * u1x[1] + u0x[4] * u1x[4] + u0x[7] * u1x[7]); // k22
    e[5] = u1x[1] + u1x[3] +
           (u0x[0] * u1x[1] + u0x[3] * u1x[4] + u0x[6] * u1x[7] +
            u1x[0] * u0x[1] + u1x[3] * u0x[4] + u1x[6] * u0x[7]); // k12

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4]; // e23, transverse shear
    e[7] = 2.0 * e0ty[2]; // e13, transverse shear
    e[8] = et[0]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void NonlinearShellStrainForwardCore(const T u0x[], const T u1x[],
                                const T e0ty[], const T et[], T e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0]; // e11
    e[1] = e0ty[3]; // e22
    e[2] = 2.0 * e0ty[1]; // e12

    // Compute the bending strain (here's where nonlinearity comes in)
    e[3] = u1x[0] + (u0x[0] * u1x[0] + u0x[3] * u1x[3] + u0x[6] * u1x[6]); // k11
    e[4] = u1x[4] + (u0x[1] * u1x[1] + u0x[4] * u1x[4] + u0x[7] * u1x[7]); // k22
    e[5] = u1x[1] + u1x[3] +
           (u0x[0] * u1x[1] + u0x[3] * u1x[4] + u0x[6] * u1x[7] +
            u1x[0] * u0x[1] + u1x[3] * u0x[4] + u1x[6] * u0x[7]); // k12

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4]; // e23, transverse shear
    e[7] = 2.0 * e0ty[2]; // e13, transverse shear
    e[8] = et[0]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void NonlinearShellStrainReverseCore(
  const T eb[],
  const T u0x[], const T u1x[], const T e0ty[], const T et[], 
  T u0xb[], T u1xb[], T e0tyb[], T etb[]) {

    // This is just 1st order backprop version
    // Evaluate the in-plane strains from the tying strain expressions
    // ----------------------
    e0tyb[0] += eb[0]; // e1
    e0tyb[3] += eb[1]; // e22
    e0tyb[1] += 2.0 * eb[2]; // e12

    // Compute the bending strain
    // --------------------------
    // k11 computation
    u1xb[0] += eb[3];
    u0xb[0] += u1x[0] * eb[3]; u1xb[0] += u0x[0] * eb[3];
    u0xb[3] += u1x[3] * eb[3]; u1xb[3] += u0x[3] * eb[3];
    u0xb[6] += u1x[6] * eb[3]; u1xb[6] += u0x[6] * eb[3];
    // k22 computation
    u1xb[4] += eb[4];
    u0xb[1] += u1x[1] * eb[4]; u1xb[1] += u0x[1] * eb[4];
    u0xb[4] += u1x[4] * eb[4]; u1xb[4] += u0x[4] * eb[4];
    u0xb[7] += u1x[7] * eb[4]; u1xb[7] += u0x[7] * eb[4];
    // k12 computation
    u1xb[1] += eb[5]; u1xb[3] += eb[5];
    u0xb[0] += u1x[1] * eb[5]; u1xb[0] += u0x[1] * eb[5];
    u0xb[1] += u1x[0] * eb[5]; u1xb[1] += u0x[0] * eb[5];
    u0xb[3] += u1x[4] * eb[5]; u1xb[3] += u0x[4] * eb[5];
    u0xb[4] += u1x[3] * eb[5]; u1xb[4] += u0x[3] * eb[5];
    u0xb[6] += u1x[7] * eb[5]; u1xb[6] += u0x[7] * eb[5];
    u0xb[7] += u1x[6] * eb[5]; u1xb[7] += u0x[6] * eb[5];

    // Add the components of the shear strain
    // --------------------------------------
    e0tyb[4] += 2.0 * eb[6]; // e23, transverse shear
    e0tyb[2] += 2.0 * eb[7]; // e13, transverse shear
    etb[0] += eb[8]; // e12 (drill strain)
}

template <typename T>
A2D_FUNCTION void NonlinearShellStrainHessianReverseCore(
  const T eh[], const T eb[],
  const T u0x[], const T u1x[], const T e0ty[], const T et[], 
  const T u0xp[], const T u1xp[], const T e0typ[], const T etp[], 
  T u0xh[], T u1xh[], T e0tyh[], T eth[]) {
    // This is the 2nd order backprop version
    // Tip: use alt+shift+leftclick and drag to multi-line select and make this code easier.
    // Evaluate the in-plane strains from the tying strain expressions
    // ----------------------
    e0tyh[0] += eh[0]; // e1
    e0tyh[3] += eh[1]; // e22
    e0tyh[1] += 2.0 * eh[2]; // e12

    // Compute the bending strain
    // --------------------------
    // k11 computation
    u1xh[0] += eh[3];
    //   nonlinear input * h terms
    u0xh[0] += u1x[0] * eh[3]; u1xh[0] += u0x[0] * eh[3];
    u0xh[3] += u1x[3] * eh[3]; u1xh[3] += u0x[3] * eh[3];
    u0xh[6] += u1x[6] * eh[3]; u1xh[6] += u0x[6] * eh[3];
    //   nonlinear bar * ptest terms
    u0xh[0] += u1xp[0] * eb[3]; u1xh[0] += u0xp[0] * eb[3];
    u0xh[3] += u1xp[3] * eb[3]; u1xh[3] += u0xp[3] * eb[3];
    u0xh[6] += u1xp[6] * eb[3]; u1xh[6] += u0xp[6] * eb[3];
    
    // k22 computation
    u1xh[4] += eh[4];
    //   nonlinear input * h terms
    u0xh[1] += u1x[1] * eh[4]; u1xh[1] += u0x[1] * eh[4];
    u0xh[4] += u1x[4] * eh[4]; u1xh[4] += u0x[4] * eh[4];
    u0xh[7] += u1x[7] * eh[4]; u1xh[7] += u0x[7] * eh[4];
    //   nonlinear bar * ptest terms
    u0xh[1] += u1xp[1] * eb[4]; u1xh[1] += u0xp[1] * eb[4];
    u0xh[4] += u1xp[4] * eb[4]; u1xh[4] += u0xp[4] * eb[4];
    u0xh[7] += u1xp[7] * eb[4]; u1xh[7] += u0xp[7] * eb[4];

    // k12 computation
    u1xh[1] += eh[5]; u1xh[3] += eh[5];
    //   nonlinear input * h terms
    u0xh[0] += u1x[1] * eh[5]; u1xh[0] += u0x[1] * eh[5];
    u0xh[1] += u1x[0] * eh[5]; u1xh[1] += u0x[0] * eh[5];
    u0xh[3] += u1x[4] * eh[5]; u1xh[3] += u0x[4] * eh[5];
    u0xh[4] += u1x[3] * eh[5]; u1xh[4] += u0x[3] * eh[5];
    u0xh[6] += u1x[7] * eh[5]; u1xh[6] += u0x[7] * eh[5];
    u0xh[7] += u1x[6] * eh[5]; u1xh[7] += u0x[6] * eh[5];
    //   nonlinear bar * ptest terms
    u0xh[0] += u1xp[1] * eb[5]; u1xh[0] += u0xp[1] * eb[5];
    u0xh[1] += u1xp[0] * eb[5]; u1xh[1] += u0xp[0] * eb[5];
    u0xh[3] += u1xp[4] * eb[5]; u1xh[3] += u0xp[4] * eb[5];
    u0xh[4] += u1xp[3] * eb[5]; u1xh[4] += u0xp[3] * eb[5];
    u0xh[6] += u1xp[7] * eb[5]; u1xh[6] += u0xp[7] * eb[5];
    u0xh[7] += u1xp[6] * eb[5]; u1xh[7] += u0xp[6] * eb[5];

    // Add the components of the shear strain
    // --------------------------------------
    e0tyh[4] += 2.0 * eh[6]; // e23, transverse shear
    e0tyh[2] += 2.0 * eh[7]; // e13, transverse shear
    eth[0] += eh[8]; // e12 (drill strain)
}

template <ShellStrainType straintype, typename T>
A2D_FUNCTION void ShellStrain(const Mat<T,3,3> &u0x, const Mat<T,3,3> &u1x,
                                const SymMat<T,3> &e0ty, const Vec<T,1> et, Vec<T,9> &e) {
  if constexpr (straintype == ShellStrainType::LINEAR) {
    LinearShellStrainCore<T>(get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et), get_data(e));
  } else {
    NonlinearShellStrainCore<T>(get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et), get_data(e));
  }
}

template <ShellStrainType straintype, class u0xtype, class u1xtype, class e0tytype, class ettype, class etype>
class ShellStrainExpr {
 public:
  // Extract the numeric type to use
  typedef typename get_object_numeric_type<etype>::type T;

  // Extract the dimensions of the underlying matrix
  static constexpr int u0x_rows = get_matrix_rows<u0xtype>::size;
  static constexpr int u0x_cols = get_matrix_columns<u0xtype>::size;
  static constexpr int u1x_rows = get_matrix_rows<u1xtype>::size;
  static constexpr int u1x_cols = get_matrix_columns<u1xtype>::size;
  static constexpr int e0ty_size = get_symmatrix_size<e0tytype>::size;
  static constexpr int et_size = get_vec_size<ettype>::size;
  static constexpr int e_size = get_vec_size<etype>::size;

  // make sure the correct sizes
  static_assert(
    (u0x_rows == 3) && (u0x_cols == 3) && (u1x_rows == 3) && (u1x_cols == 3) && (e0ty_size == 3) && (et_size == 1) && (e_size == 9),
    "Shell Strain Expression does not have right size.."
  );

  static constexpr ADiffType adu0x = get_diff_type<u0xtype>::diff_type;
  static constexpr ADiffType adu1x = get_diff_type<u1xtype>::diff_type;
  static constexpr ADiffType ade0ty = get_diff_type<e0tytype>::diff_type;
  static constexpr ADiffType adet = get_diff_type<ettype>::diff_type;

  // Get the differentiation order from the output
  static constexpr ADorder order = get_diff_order<etype>::order;

  // Make sure that the order matches
  static_assert(get_diff_order<u0xtype>::order == order,
                "ADorder does not match");

  A2D_FUNCTION ShellStrainExpr(u0xtype& u0x, u1xtype& u1x, e0tytype& e0ty, ettype& et, etype& e) : 
        u0x(u0x), u1x(u1x), e0ty(e0ty), et(et), e(e) {}

  A2D_FUNCTION void eval() {
    if constexpr (straintype == ShellStrainType::LINEAR) {
      LinearShellStrainCore<T>(get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et), get_data(e));
    } else {
      NonlinearShellStrainCore<T>(get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et), get_data(e));
    }
  }

  A2D_FUNCTION void bzero() { e.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    static_assert(
        !(order == ADorder::FIRST and forder == ADorder::SECOND),
        "Can't perform second order forward with first order objects");
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;

    // need more statements here? (maybe some with only some pvalues transferred forward at a time? see matSum expr)
    if constexpr (straintype == ShellStrainType::LINEAR) {
      LinearShellStrainForwardCore<T>(GetSeed<seed>::get_data(u0x), GetSeed<seed>::get_data(u1x), GetSeed<seed>::get_data(e0ty), 
            GetSeed<seed>::get_data(et), GetSeed<seed>::get_data(e));
    } else {
      NonlinearShellStrainForwardCore<T>(GetSeed<seed>::get_data(u0x), GetSeed<seed>::get_data(u1x), 
            GetSeed<seed>::get_data(e0ty), GetSeed<seed>::get_data(et), GetSeed<seed>::get_data(e));
    }
  }

  A2D_FUNCTION void reverse() {
    constexpr ADseed seed = ADseed::b;
    // need more conditions on which ADseeds are active here
    if constexpr (straintype == ShellStrainType::LINEAR) {
      LinearShellStrainReverseCore<T>(
            GetSeed<seed>::get_data(e),
            GetSeed<seed>::get_data(u0x), GetSeed<seed>::get_data(u1x), 
            GetSeed<seed>::get_data(e0ty), GetSeed<seed>::get_data(et));
    } else {
      NonlinearShellStrainReverseCore<T>(
            GetSeed<seed>::get_data(e),
            get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et),
            GetSeed<seed>::get_data(u0x), GetSeed<seed>::get_data(u1x), 
            GetSeed<seed>::get_data(e0ty), GetSeed<seed>::get_data(et));
    }
  }

  A2D_FUNCTION void hzero() { e.hzero(); }

  A2D_FUNCTION void hreverse() {
    // need more conditions on which ADseeds are active here
    constexpr ADseed seed = ADseed::h;
    if constexpr (straintype == ShellStrainType::LINEAR) {
      LinearShellStrainReverseCore<T>(
            GetSeed<seed>::get_data(e),
            GetSeed<seed>::get_data(u0x), GetSeed<seed>::get_data(u1x), 
            GetSeed<seed>::get_data(e0ty), GetSeed<seed>::get_data(et));
    } else {
      NonlinearShellStrainHessianReverseCore<T>(
          GetSeed<ADseed::h>::get_data(e), GetSeed<ADseed::b>::get_data(e),
          get_data(u0x), get_data(u1x), get_data(e0ty), get_data(et),
          GetSeed<ADseed::p>::get_data(u0x), GetSeed<ADseed::p>::get_data(u1x), 
          GetSeed<ADseed::p>::get_data(e0ty), GetSeed<ADseed::p>::get_data(et),
          GetSeed<ADseed::h>::get_data(u0x), GetSeed<ADseed::h>::get_data(u1x), 
          GetSeed<ADseed::h>::get_data(e0ty), GetSeed<ADseed::h>::get_data(et));
    }
  }

  u0xtype& u0x;
  u1xtype& u1x;
  e0tytype& e0ty;
  ettype& et;
  etype& e;
};

// template <ShellStrainType straintype, typename T>
// A2D_FUNCTION void ShellStrain(const Mat<T,3,3> &u0x, const Mat<T,3,3> &u1x,
//                                 const SymMat<T,3> &e0ty, const T &et, Vec<T,9> &e) {

template <ShellStrainType straintype, class u0xtype, class u1xtype, class e0tytype, class ettype, class etype>
A2D_FUNCTION auto ShellStrain(ADObj<u0xtype> &u0x, ADObj<u1xtype> &u1x, ADObj<e0tytype> &e0ty, ADObj<ettype> &et, ADObj<etype> &e) {
  return ShellStrainExpr<straintype, ADObj<u0xtype>, ADObj<u1xtype>, ADObj<e0tytype>, ADObj<ettype>, ADObj<etype>>(u0x, u1x, e0ty, et, e);
}

template <ShellStrainType straintype, class u0xtype, class u1xtype, class e0tytype, class ettype, class etype>
A2D_FUNCTION auto ShellStrain(A2DObj<u0xtype> &u0x, A2DObj<u1xtype> &u1x, A2DObj<e0tytype> &e0ty, A2DObj<ettype> &et, A2DObj<etype> &e) {
  return ShellStrainExpr<straintype, A2DObj<u0xtype>, A2DObj<u1xtype>, A2DObj<e0tytype>, A2DObj<ettype>, A2DObj<etype>>(u0x, u1x, e0ty, et, e);
}

namespace Test {

// template <ShellStrainType straintype, typename T>
// A2D_FUNCTION void ShellStrain(const Mat<T,3,3> &u0x, const Mat<T,3,3> &u1x,
//                                 const SymMat<T,3> &e0ty, const T &et, Vec<T,9> &e)

template <ShellStrainType straintype, typename T>
class ShellStrainTest : public A2DTest<T, Vec<T,9>, Mat<T,3,3>, Mat<T,3,3>, SymMat<T,3>, Vec<T,1>> {
 public:
  using Input = VarTuple<T, Mat<T,3,3>, Mat<T,3,3>, SymMat<T,3>, Vec<T,1>>;
  using Output = VarTuple<T, Vec<T,9>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "ShellStrain<";
    if (straintype == ShellStrainType::LINEAR) {
      s << "LINEAR>";
    } else {
      s << "NONLINEAR>";
    }

    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T,3,3> u0x, u1x;
    SymMat<T,3> e0ty;
    Vec<T,1> et;
    Vec<T,9> e;
    x.get_values(u0x, u1x, e0ty, et);
    ShellStrain<straintype>(u0x, u1x, e0ty, et, e);
    return MakeVarTuple<T>(e);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T,3,3>> u0x, u1x;
    ADObj<SymMat<T,3>>e0ty;
    ADObj<Vec<T,1>> et;
    ADObj<Vec<T,9>> e;
    x.get_values(u0x.value(), u1x.value(), e0ty.value(), et.value());
    auto stack = MakeStack(ShellStrain<straintype>(u0x, u1x, e0ty, et, e));
    seed.get_values(e.bvalue());
    stack.reverse();
    g.set_values(u0x.bvalue(), u1x.bvalue(), e0ty.bvalue(), et.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T,3,3>> u0x, u1x;
    A2DObj<SymMat<T,3>>e0ty;
    A2DObj<Vec<T,1>> et;
    A2DObj<Vec<T,9>> e;
    x.get_values(u0x.value(), u1x.value(), e0ty.value(), et.value());
    p.get_values(u0x.pvalue(), u1x.pvalue(), e0ty.pvalue(), et.pvalue());
    auto stack = MakeStack(ShellStrain<straintype>(u0x, u1x, e0ty, et, e));
    seed.get_values(e.bvalue());
    hval.get_values(e.hvalue());
    stack.hproduct();
    h.set_values(u0x.hvalue(), u1x.hvalue(), e0ty.hvalue(), et.hvalue());
  }
};

bool ShellStrainTestAll(bool component = false, bool write_output = true) {
  using Tc = std::complex<double>;

  ShellStrainTest<ShellStrainType::LINEAR, Tc> test1;
  bool passed = Run(test1, component, write_output);

  ShellStrainTest<ShellStrainType::NONLINEAR, Tc> test2;
  passed = passed && Run(test2, component, write_output);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_SHELL_STRAIN_H