#include <gtest/gtest.h>

#include "a2dcore.h"

using namespace A2D;

template <GreenStrainType etype, typename T, int N>
class MatGreenStrainTest : public Test::A2DTest<T, SymMat<T, N>, Mat<T, N, N>> {
 public:
  using Input = VarTuple<T, Mat<T, N, N>>;
  using Output = VarTuple<T, SymMat<T, N>>;

  // Assemble a string to describe the test
  std::string name() {
    std::stringstream s;
    s << "MatGreenStrain<";
    if (etype == GreenStrainType::LINEAR) {
      s << "LINEAR," << N << ">";
    } else {
      s << "NONLINEAR," << N << ">";
    }

    return s.str();
  }

  // Evaluate the matrix-matrix product
  Output eval(const Input& x) {
    Mat<T, N, N> Ux;
    SymMat<T, N> E;
    x.get_values(Ux);
    MatGreenStrain<etype>(Ux, E);
    return MakeVarTuple<T>(E);
  }

  // Compute the derivative
  void deriv(const Output& seed, const Input& x, Input& g) {
    ADObj<Mat<T, N, N>> Ux;
    ADObj<SymMat<T, N>> E;

    x.get_values(Ux.value());
    auto stack = MakeStack(MatGreenStrain<etype>(Ux, E));
    seed.get_values(E.bvalue());
    stack.reverse();
    g.set_values(Ux.bvalue());
  }

  // Compute the second-derivative
  void hprod(const Output& seed, const Output& hval, const Input& x,
             const Input& p, Input& h) {
    A2DObj<Mat<T, N, N>> Ux;
    A2DObj<SymMat<T, N>> E;

    x.get_values(Ux.value());
    p.get_values(Ux.pvalue());
    auto stack = MakeStack(MatGreenStrain<etype>(Ux, E));
    seed.get_values(E.bvalue());
    hval.get_values(E.hvalue());
    stack.hproduct();
    h.set_values(Ux.hvalue());
  }
};

inline bool MatGreenStrainTestAll(bool component = false,
                                  bool write_output = true) {
  using Tc = A2D_complex_t<double>;

  bool passed = true;
  MatGreenStrainTest<GreenStrainType::LINEAR, Tc, 2> test1;
  passed = passed && Run(test1, component, write_output);
  MatGreenStrainTest<GreenStrainType::NONLINEAR, Tc, 2> test2;
  passed = passed && Run(test2, component, write_output);

  MatGreenStrainTest<GreenStrainType::LINEAR, Tc, 3> test3;
  passed = passed && Run(test3, component, write_output);
  MatGreenStrainTest<GreenStrainType::NONLINEAR, Tc, 3> test4;
  passed = passed && Run(test4, component, write_output);

  return passed;
}

TEST(test_greenstrain, MatGreenStrainTestAll) {
  EXPECT_TRUE(MatGreenStrainTestAll());
}