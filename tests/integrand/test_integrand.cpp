#include <complex>
#include <iostream>
#include <vector>

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/integrand_poisson.h"
// #include "multiphysics/integrand_heat_conduction.h"
#include "multiphysics/integrand_test.h"

using namespace A2D;

template <typename T>
using TopoIntegrand = IntegrandTopoLinearElasticity<T, 3>;

template <typename T>
using PoissonIntegrand = Poisson<T, 3>;

// template <typename T>
// using HeatIntegrand = HeatConduction<T, 3>;

bool TestIntegrands(bool component, bool write_output) {
  bool passed = true;

  using Tc = std::complex<double>;

  TopoIntegrand<Tc> integrand1(Tc(0.7), Tc(0.3), Tc(5.0));
  A2D::Test::A2DIntegrandAnalysisTest<TopoIntegrand, Tc> test1(integrand1);
  passed = passed && A2D::Test::Run(test1, component, write_output);

  PoissonIntegrand<Tc> integrand2;
  A2D::Test::A2DIntegrandAnalysisTest<PoissonIntegrand, Tc> test2(integrand2);
  passed = passed && A2D::Test::Run(test2, component, write_output);

  // HeatIntegrand<Tc> integrand2(Tc(0.7), Tc(0.3), Tc(5.0));
  // A2D::Test::A2DIntegrandAnalysisTest<HeatIntegrand, Tc> test2(integrand2);
  // passed = passed && A2D::Test::Run(test2, component, write_output);

  return passed;
}

int main() {
  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

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
  return fail;
}
