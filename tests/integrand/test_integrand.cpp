#include <complex>
#include <iostream>
#include <vector>

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/integrand_elasticity.h"
// #include "multiphysics/integrand_heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_poisson.h"
#include "multiphysics/integrand_test.h"

using namespace A2D;

template <typename T>
using TopoIntegrand = TopoElasticityIntegrand<T, 3, GreenStrainType::LINEAR>;

template <typename T>
using NonlinearTopoIntegrand =
    TopoElasticityIntegrand<T, 3, GreenStrainType::NONLINEAR>;

template <typename T>
using TopoBodyIntegrand = TopoBodyForceIntegrand<T, 3>;

// IntegrandTopoLinearElasticity<T, 3>;

template <typename T>
using PoissonIntegrand = Poisson<T, 3>;

// template <typename T>
// using MixedPoissonIntegrand = MixedPoisson<T, 3>;

// template <typename T>
// using HeatIntegrand = HeatConduction<T, 3>;

bool TestIntegrands(bool component, bool write_output) {
  bool passed = true;

  using Tc = std::complex<double>;

  TopoIntegrand<Tc> integrand1(Tc(0.7), Tc(0.3), Tc(5.0));
  A2D::Test::A2DIntegrandAnalysisTest<TopoIntegrand, Tc> test1(integrand1);
  passed = passed && A2D::Test::Run(test1, component, write_output);

  NonlinearTopoIntegrand<Tc> integrand2(Tc(0.7), Tc(0.3), Tc(5.0));
  A2D::Test::A2DIntegrandAnalysisTest<NonlinearTopoIntegrand, Tc> test2(
      integrand2);
  passed = passed && A2D::Test::Run(test2, component, write_output);

  Tc tx[3] = {1.1, -1.2, -0.8};
  TopoBodyIntegrand<Tc> integrand3(Tc(5.0), tx);
  A2D::Test::A2DIntegrandAnalysisTest<TopoBodyIntegrand, Tc> test3(integrand3);
  passed = passed && A2D::Test::Run(test3, component, write_output);

  PoissonIntegrand<Tc> integrand4;
  A2D::Test::A2DIntegrandAnalysisTest<PoissonIntegrand, Tc> test4(integrand4);
  passed = passed && A2D::Test::Run(test4, component, write_output);

  // MixedPoissonIntegrand<Tc> integrand3;
  // A2D::Test::A2DIntegrandAnalysisTest<MixedPoissonIntegrand, Tc> test3(
  //     integrand3);
  // passed = passed && A2D::Test::Run(test3, component, write_output);

  // HeatIntegrand<Tc> integrand4(Tc(0.7), Tc(0.3), Tc(5.0));
  // A2D::Test::A2DIntegrandAnalysisTest<HeatIntegrand, Tc> test4(integrand4);
  // passed = passed && A2D::Test::Run(test4, component, write_output);

  return passed;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

  // Check for the write_output flag
  for (int i = 0; i < argc; i++) {
    std::string str(argv[i]);
    if (str.compare("--write_output") == 0) {
      write_output = true;
    }
    if (str.compare("--component") == 0) {
      component = true;
    }
  }

  typedef std::function<bool(bool, bool)> TestFunc;
  std::vector<TestFunc> tests;

  tests.push_back(TestIntegrands);

  bool passed = true;
  for (int i = 0; i < tests.size(); i++) {
    bool test_passed = tests[i](component, write_output);

    // If the test fails, perform a component test and print the results
    // to screen
    if (!test_passed) {
      bool comp = true;
      bool write = true;
      tests[i](comp, write);
    }

    passed = test_passed && passed;
  }

  // Convert to a fail flag
  int fail = 0;
  if (!passed) {
    fail = 1;
  }

  Kokkos::finalize();
  return fail;
}
