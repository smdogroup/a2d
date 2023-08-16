#ifndef A2D_TEST_H
#define A2D_TEST_H

#include <iomanip>

#include "a2dvartuple.h"

namespace A2D {

namespace Test {

template <typename T, class... Inputs>
class A2DTest {
 public:
  /**
   * @brief Construct the Test
   */
  A2DTest(double dh = 1e-50, double rtol = 1e-10, double atol = 1e-30)
      : dh(dh), rtol(rtol), atol(atol) {}

  /**
   * @brief Get the complex-step step size
   */
  double get_step_size() const { return dh; }

  /**
   * @brief Get the point to perform the test at - defaults to random
   *
   * @param x Variable tuple for the inputs
   */

  virtual void get_point(VarTuple<T, Inputs...>& x) { x.set_rand(); }

  /**
   * @brief Evaluate the outputs as a function of the inputs
   *
   * @param x Variable tuple for the input
   * @return The value of the parameter
   */
  virtual T eval(const VarTuple<T, Inputs...>& x) = 0;

  /**
   * @brief Compute the derivative of the output as a function of the inputs
   *
   * @param x Variable tuple for the input
   * @param g Derivative of the output w.r.t. the input
   */
  virtual void deriv(const VarTuple<T, Inputs...>& x,
                     VarTuple<T, Inputs...>& g) = 0;

  /**
   * @brief Check whether two values are close to one another
   *
   * @param test_value The test value (computed using AD)
   * @param ref_value The reference value (usually complex-step)
   * @return true If the test passes
   * @return false If the test fails
   */
  bool is_close(const T test_value, const T ref_value) {
    T abs_err = fabs(std::real(test_value - ref_value));
    T combo = atol + rtol * fabs(std::real(ref_value));
    if (std::real(abs_err) > std::real(combo)) {
      return false;
    }
    return true;
  }

  /**
   * @brief Write a brief, one-line summary of the result
   *
   * @param out Where to write the result
   * @param test_value The test value (computed using AD)
   * @param ref_value The reference value (usually complex-step)
   */
  void write_result(std::string str, std::ostream& out, const T test_value,
                    const T ref_value) {
    // Find the result of the test
    bool passed = is_close(test_value, ref_value);

    // Compute the relative error
    T abs_err = fabs(std::real(test_value - ref_value));
    T rel_err = fabs(std::real((test_value - ref_value) / ref_value));

    if (passed) {
      out << str << " PASSED.";
    } else {
      out << str << " FAILED.";
    }
    out << std::setprecision(9) << " AD: " << std::setw(12)
        << std::real(test_value) << " CS: " << std::setw(12)
        << std::real(ref_value) << " Rel Err: " << std::setw(12)
        << std::real(rel_err) << " Abs Err: " << std::setw(12)
        << std::real(abs_err) << std::endl;
  }

 private:
  double dh;    // Complex-step size
  double rtol;  // relative tolerance
  double atol;  // absolute tolerance
};

// Perform a complex step test
template <typename T, class... Inputs>
bool RunADTest(A2DTest<std::complex<T>, Inputs...>& test) {
  // Declare all of the variables needed
  VarTuple<std::complex<T>, Inputs...> x, g, x1, p;

  // Get the starting point
  test.get_point(x);

  // Evaluate the function and its derivatives
  test.deriv(x, g);

  // Set a random direction for the test
  p.set_rand();

  // Set x1 = x + dh * p1
  double dh = test.get_step_size();
  for (index_t i = 0; i < x.get_num_components(); i++) {
    x1[i] = std::complex<double>(std::real(x[i]), std::real(dh * p[i]));
  }

  // Compute the complex-step result: fd = p^{T} * df/dx
  T fd = std::imag(test.eval(x1)) / dh;

  // Compute the solution from the AD
  T ans = 0.0;
  for (index_t i = 0; i < x.get_num_components(); i++) {
    ans += std::real(g[i] * p[i]);
  }

  bool passed = test.is_close(ans, fd);

  test.write_result("First-order", std::cout, ans, fd);

  return passed;
}

}  // namespace Test

}  // namespace A2D

#endif  // A2D_TEST_H