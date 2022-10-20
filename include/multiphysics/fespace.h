#ifndef A2D_FE_SOLUTION_SPACE_H
#define A2D_FE_SOLUTION_SPACE_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"

namespace A2D {

template <typename T, A2D::index_t D>
class L1ScalarSpace {
 public:
  L1ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const A2D::index_t components = 1;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const A2D::index_t ndof = 1;

  // Set the input seed
  void set_seed(const A2D::index_t seed) { u = 1.0; }
  T get_value(const A2D::index_t seed) { return u; }

  T& get_value() { return u; }
  const T& get_value() const { return u; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, L1ScalarSpace<T, D>& s) {
    s.u = u;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, L1ScalarSpace<T, D>& s) {
    s.u = u;
  }

 private:
  T u;
};

template <typename T, A2D::index_t D>
class H1ScalarSpace {
 public:
  H1ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const A2D::index_t components = 1;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const A2D::index_t ndof = 1 + D;

  // Set the seed value
  void set_seed(const A2D::index_t seed) {
    if (seed == 0) {
      u = 1.0;
    } else {
      grad(seed - 1) = 1.0;
    }
  }
  T get_value(const A2D::index_t seed) {
    if (seed == 0) {
      return u;
    } else {
      return grad(seed - 1);
    }
  }

  T& get_value() { return u; }
  A2D::Vec<T, D>& get_grad() { return grad; }

  const T& get_value() const { return u; }
  const A2D::Vec<T, D>& get_grad() const { return grad; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }

 private:
  T u;
  A2D::Vec<T, D> grad;
};

template <typename T, A2D::index_t C, A2D::index_t D>
class H1Space {
 public:
  H1Space() {}

  // Number of solution components
  static const A2D::index_t components = C;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const A2D::index_t ndof = (D + 1) * C;

  // Set the seed value
  void set_seed(const A2D::index_t seed) {
    if (seed < C) {
      u(seed) = 1.0;
    } else {
      grad((seed - C) / D, (seed - C) % D) = 1.0;
    }
  }
  T get_value(const A2D::index_t seed) {
    if (seed < C) {
      return u(seed);
    } else {
      return grad((seed - C) / D, (seed - C) % D);
    }
  }

  A2D::Vec<T, C>& get_value() { return u; }
  A2D::Mat<T, C, D>& get_grad() { return grad; }

  const A2D::Vec<T, C>& get_value() const { return u; }
  const A2D::Mat<T, C, D>& get_grad() const { return grad; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }

 private:
  A2D::Vec<T, C> u;
  A2D::Mat<T, C, D> grad;
};

template <typename T>
class Hdiv2DSpace {
 public:
  Hdiv2DSpace() { div = 0.0; }
  // Number of solution components
  static const A2D::index_t components = 2;

  // Spatial dimension
  static const A2D::index_t dim = 2;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const A2D::index_t ndof = 3;

  void set_seed(const A2D::index_t seed) {
    if (seed < 2) {
      u(seed) = 1.0;
    } else {
      div = 1.0;
    }
  }
  T get_value(const A2D::index_t seed) {
    if (seed < 2) {
      return u(seed);
    } else {
      return div;
    }
  }

  A2D::Vec<T, 2>& get_value() { return u; }
  T& get_div() { return div; }

  const A2D::Vec<T, 2>& get_value() const { return u; }
  const T& get_div() const { return div; }

  void transform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                 const A2D::Mat<T, 2, 2>& Jinv, Hdiv2DSpace<T>& s) {
    s.u = u;
    s.div = div;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                  const A2D::Mat<T, 2, 2>& Jinv, Hdiv2DSpace<T>& s) {
    s.u = u;
    s.div = div;
  }

 private:
  A2D::Vec<T, 2> u;
  T div;
};

template <typename T>
class Hcurl2DSpace {
 public:
  Hcurl2DSpace() { curl = 0.0; }

  // Number of solution components
  static const A2D::index_t components = 2;

  // Spatial dimension
  static const A2D::index_t dim = 2;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const A2D::index_t ndof = 3;

  void set_seed(const A2D::index_t seed) {
    if (seed < 2) {
      u(seed) = 1.0;
    } else {
      curl = 1.0;
    }
  }
  T get_value(const A2D::index_t seed) {
    if (seed < 2) {
      return u(seed);
    } else {
      return curl;
    }
  }

  A2D::Vec<T, 2>& get_value() { return u; }
  T& get_curl() { return curl; }

  const A2D::Vec<T, 2>& get_value() const { return u; }
  const T& get_curl() const { return curl; }

  void transform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                 const A2D::Mat<T, 2, 2>& Jinv, Hcurl2DSpace<T>& s) {
    s.u = u;
    s.curl = curl;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                  const A2D::Mat<T, 2, 2>& Jinv, Hcurl2DSpace<T>& s) {
    s.u = u;
    s.curl = curl;
  }

 private:
  A2D::Vec<T, 2> u;
  T curl;
};

template <class... Spaces>
struct __count_solution_ndof;

template <>
struct __count_solution_ndof<> {
  static const A2D::index_t ndof = 0;
};

template <class First, class... Remain>
struct __count_solution_ndof<First, Remain...> {
  static const A2D::index_t ndof =
      First::ndof + __count_solution_ndof<Remain...>::ndof;
};

/*
  A collection of spaces used in the finite-element problem
*/
template <typename T, A2D::index_t D, class... Spaces>
class FESpace {
 public:
  typedef std::tuple<Spaces...> SolutionSpace;
  static constexpr A2D::index_t nspaces =
      std::tuple_size<std::tuple<Spaces...>>();

  /*
    Count up the total number of degrees of freedom
  */
  static constexpr A2D::index_t ndof = __count_solution_ndof<Spaces...>::ndof;

  /*
    Set a seed on a given solution field
  */
  void set_seed(const A2D::index_t seed) { set_seed_<0, Spaces...>(seed); }

  /*
    Get a solution value based on the index
  */
  T get_value(const A2D::index_t seed) {
    return get_value_<0, Spaces...>(seed);
  }

  /*
    Extract the specified solution space from the FESpace object
  */
  template <A2D::index_t index>
  typename std::tuple_element<index, SolutionSpace>::type& get() noexcept {
    return std::get<index>(u);
  }

  template <A2D::index_t index>
  const typename std::tuple_element<index, SolutionSpace>::type& get()
      const noexcept {
    return std::get<index>(u);
  }

  /*
    Transform the solution space from the reference element to the physical
    element
  */
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    transform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  /*
    Perform the reverse of the transform - transfer the derivative from the
    physical element to the reference element
  */
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    rtransform_<0, Spaces...>(detJ, J, Jinv, s);
  }

 private:
  // Solution space tuple object
  SolutionSpace u;

  template <A2D::index_t index>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {}

  template <A2D::index_t index, class First, class... Remain>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    std::get<index>(u).transform(detJ, J, Jinv, std::get<index>(s.u));
    transform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <A2D::index_t index>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
  }

  template <A2D::index_t index, class First, class... Remain>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    std::get<index>(u).rtransform(detJ, J, Jinv, std::get<index>(s.u));
    rtransform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <A2D::index_t index>
  void set_seed_(const A2D::index_t seed) {}

  template <A2D::index_t index, class First, class... Remain>
  void set_seed_(const A2D::index_t seed) {
    if (seed < First::ndof) {
      std::get<index>(u).set_seed(seed);
    } else {
      set_seed_<index + 1, Remain...>(seed - First::ndof);
    }
  }

  template <A2D::index_t index>
  T get_value_(const A2D::index_t seed) {
    return 0.0;
  }

  template <A2D::index_t index, class First, class... Remain>
  T get_value_(const A2D::index_t seed) {
    if (seed < First::ndof) {
      return std::get<index>(u).get_value(seed);
    } else {
      return get_value_<index + 1, Remain...>(seed - First::ndof);
    }
  }
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_SPACE_H