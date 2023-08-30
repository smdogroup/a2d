#include <vector>

#include "ad/a2dgemm.h"
#include "ad/a2dgreenstrain.h"
#include "ad/a2disotropic.h"
#include "ad/a2dmatdet.h"
#include "ad/a2dmatinv.h"
#include "ad/a2dmatsum.h"
#include "ad/a2dmattrace.h"
#include "ad/a2dsymrk.h"
#include "ad/a2dsymtrace.h"

int main() {
  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

  typedef std::function<bool(bool, bool)> TestFunc;
  std::vector<TestFunc> tests;

  tests.push_back(A2D::Test::MatMatMultTestAll);
  tests.push_back(A2D::Test::MatDetTestAll);
  tests.push_back(A2D::Test::MatInvTestAll);
  tests.push_back(A2D::Test::MatTraceTestAll);
  tests.push_back(A2D::Test::MatGreenStrainTestAll);
  tests.push_back(A2D::Test::SymMatTraceTestAll);
  // tests.push_back(A2D::Test::SymIsotropicTestAll);
  tests.push_back(A2D::Test::MatSumTestAll);
  tests.push_back(A2D::Test::SymMatRKTestAll);

  bool passed = true;
  for (int i = 0; i < tests.size(); i++) {
    bool test_passed = tests[i](component, write_output);

    // If the test fails, perform a component test and print the results to
    // screen
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