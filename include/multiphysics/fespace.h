#ifndef A2D_FE_SOLUTION_SPACE_H
#define A2D_FE_SOLUTION_SPACE_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"

namespace A2D {

/*
  Solution space classes:

  The finite-element function space classes store the solution and derivative
  information needed to compute a weak form of the governing equations.

  The classes are named after their corresponding function space.

  Special scalar versions of the function spaces are also implemented (for H1
  and L2) so that values are scalars (not vectors) and derivatives are vectors
  (not matrices).

  The space class is also responsible for defining an ordering of the solution
  and derivative information. This ordering exists to simplify setting and
  extracting values from the objects. This is used when computing the Jacobian
  matrix.

  Each finite element space object must implement the following functionality:

  0. The constructor must initialize all solution entries to zero.

  1. static const ncomp: Number of component outputs treating solution and
  derivative information as independent.

  2. set_value(comp, val): Sets the value in the given component index, all
  other values are unmodified. The comp value must satisfy 0 <= comp < ncomps.
  Behavior undefined otherwise.

  3. get_value(comp): Returns the value in the given component index. The comp
  value must satisfy 0 <= comp < ncomps. Behavior undefined otherwise.

  4. The space object must implement getters for each of the specific values
  that it computes such as:

  get_value(), get_grad(), get_div(), get_curl() etc.

  5. sref.transform(detJ, J, Jinv, s): This function transforms the current
  sref object from the reference element to the physical element.

  6. s.rtransform(detJ, J, Jinv, sref): Transform the derivatives from the
  physical element back to the reference element.
*/

template <typename T, A2D::index_t D>
class L2ScalarSpace {
 public:
  L2ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const A2D::index_t ncomp = 1;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) { u = val; }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const { return u; }

  // Get the scalar solution value
  T& get_value() { return u; }
  const T& get_value() const { return u; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, L2ScalarSpace<T, D>& s) const {
    s.u = u;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, L2ScalarSpace<T, D>& s) const {
    s.u = u;
  }

 private:
  T u;
};

template <typename T, A2D::index_t C, A2D::index_t D>
class L2Space {
 public:
  L2Space() { u.zero(); }

  // Number of solution components
  static const A2D::index_t ncomp = C;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) { u(comp) = val; }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const { return u(comp); }

  // Get the scalar solution value
  A2D::Vec<T, C>& get_value() { return u; }
  const A2D::Vec<T, C>& get_value() const { return u; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, L2Space<T, C, D>& s) const {
    s.u = u;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, L2Space<T, C, D>& s) const {
    s.u = u;
  }

 private:
  A2D::Vec<T, C> u;
};

template <typename T, A2D::index_t D>
class H1ScalarSpace {
 public:
  H1ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const A2D::index_t ncomp = 1 + D;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) {
    if (comp == 0) {
      u = val;
    } else {
      grad(comp - 1) = val;
    }
  }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const {
    if (comp == 0) {
      return u;
    } else {
      return grad(comp - 1);
    }
  }

  // Get the solution value
  T& get_value() { return u; }
  const T& get_value() const { return u; }

  // Get the gradient of the solution
  A2D::Vec<T, D>& get_grad() { return grad; }
  const A2D::Vec<T, D>& get_grad() const { return grad; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) const {
    s.u = u;
    s.grad = grad;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) const {
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
  static const A2D::index_t ncomp = (D + 1) * C;

  // Spatial dimension
  static const A2D::index_t dim = D;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) {
    if (comp % (D + 1) == 0) {
      u(comp / (D + 1)) = val;
    } else {
      grad(comp / (D + 1), (comp % (D + 1)) - 1) = val;
    }
  }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const {
    if (comp % (D + 1) == 0) {
      return u(comp / (D + 1));
    } else {
      return grad(comp / (D + 1), (comp % (D + 1)) - 1);
    }
  }

  // Get the value of the solution
  A2D::Vec<T, C>& get_value() { return u; }
  const A2D::Vec<T, C>& get_value() const { return u; }

  // Get the gradient of the solution
  A2D::Mat<T, C, D>& get_grad() { return grad; }
  const A2D::Mat<T, C, D>& get_grad() const { return grad; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) const {
    s.u = u;

    // s.grad = grad * Jinv
    MatMatMult(grad, Jinv, s.grad);
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) const {
    // dot{s.grad} = dot{grad} Jinv =>
    // tr(bar{s.grad}^{T} * dot{s.grad})
    // = tr(bar{s.grad}^{T} * dot{grad} * Jinv)
    // = tr(Jinv * bar{s.grad}^{T} * dot{grad})
    // = tr((s.grad * Jinv^{T})^{T} * dot{grad})
    s.u = u;

    MatMatMult(s.grad, transpose(Jinv), grad);
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
  static const A2D::index_t ncomp = 3;

  // Spatial dimension
  static const A2D::index_t dim = 2;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) {
    if (comp < 2) {
      u(comp) = val;
    } else {
      div = val;
    }
  }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const {
    if (comp < 2) {
      return u(comp);
    } else {
      return div;
    }
  }

  // Get the solution value
  A2D::Vec<T, 2>& get_value() { return u; }
  const A2D::Vec<T, 2>& get_value() const { return u; }

  // Get the value of the divergence
  T& get_div() { return div; }
  const T& get_div() const { return div; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, dim, dim>& J,
                 const A2D::Mat<T, dim, dim>& Jinv, Hdiv2DSpace<T>& s) const {
    s.u = u;
    s.div = div;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, dim, dim>& J,
                  const A2D::Mat<T, dim, dim>& Jinv, Hdiv2DSpace<T>& s) const {
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
  static const A2D::index_t ncomp = 3;

  // Spatial dimension
  static const A2D::index_t dim = 2;

  // Set the input value based on the component index
  void set_value(const A2D::index_t comp, T val) {
    if (comp < 2) {
      u(comp) = val;
    } else {
      curl = val;
    }
  }

  // Get the value of the specified component
  T get_value(const A2D::index_t comp) const {
    if (comp < 2) {
      return u(comp);
    } else {
      return curl;
    }
  }

  // Get the solution value
  A2D::Vec<T, 2>& get_value() { return u; }
  const A2D::Vec<T, 2>& get_value() const { return u; }

  // Get the curl (2D so it's a scalar value)
  T& get_curl() { return curl; }
  const T& get_curl() const { return curl; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const A2D::Mat<T, dim, dim>& J,
                 const A2D::Mat<T, dim, dim>& Jinv, Hcurl2DSpace<T>& s) const {
    s.u = u;
    s.curl = curl;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const A2D::Mat<T, dim, dim>& J,
                  const A2D::Mat<T, dim, dim>& Jinv, Hcurl2DSpace<T>& s) const {
    s.u = u;
    s.curl = curl;
  }

 private:
  A2D::Vec<T, dim> u;
  T curl;
};

template <class... Spaces>
struct __count_space_components;

template <>
struct __count_space_components<> {
  static const A2D::index_t ncomp = 0;
};

template <class First, class... Remain>
struct __count_space_components<First, Remain...> {
  static const A2D::index_t ncomp =
      First::ncomp + __count_space_components<Remain...>::ncomp;
};

/*
  FESpace:

  To form a full finite-element space, we use a variadic template class FESpace,
  FESpace takes the scalar type, the physical dimension of the problem 2D, 3D
  etc. and a list of the solution space objects. For instance, this class might
  look like this when declared:

  FESpace<double, 3, H1Space, L2Space, Hdiv3DSpace>

  FESpace implements a static member ncomp as well as set_value(comp, val),
  get_value(comp) in an analogous manner to the solution spaces. However, the
  implementation in FESpace treats these as a concatenated list across all
  solution spaces.

  FESpace implements a get<index>() function that returns a reference to the
  index-th type in the FESpace. When writing to or reading from FESpace, you can
  either use the set_value(comp, val) get_value(comp) functions or the
  get<index>() that returns the solution space directly. Both can be useful.

  T The typename
  D The dimension of the problem
  Spaces The varidic template list of finite element space types
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
  static constexpr A2D::index_t ncomp =
      __count_space_components<Spaces...>::ncomp;

  /*
    Set a value on a given solution field
  */
  void set_value(const A2D::index_t comp, T val) {
    set_value_<0, Spaces...>(comp, val);
  }

  /*
    Get a solution value based on the index
  */
  T get_value(const A2D::index_t comp) const {
    return get_value_<0, Spaces...>(comp);
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
                 const A2D::Mat<T, D, D>& Jinv,
                 FESpace<T, D, Spaces...>& s) const {
    transform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  /*
    Perform the reverse of the transform - transfer the derivative from the
    physical element to the reference element
  */
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv,
                  FESpace<T, D, Spaces...>& s) const {
    rtransform_<0, Spaces...>(detJ, J, Jinv, s);
  }

 private:
  // Solution space tuple object
  SolutionSpace u;

  template <A2D::index_t index>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv,
                  FESpace<T, D, Spaces...>& s) const {}

  template <A2D::index_t index, class First, class... Remain>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv,
                  FESpace<T, D, Spaces...>& s) const {
    std::get<index>(u).transform(detJ, J, Jinv, std::get<index>(s.u));
    transform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <A2D::index_t index>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv,
                   FESpace<T, D, Spaces...>& s) const {}

  template <A2D::index_t index, class First, class... Remain>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv,
                   FESpace<T, D, Spaces...>& s) const {
    std::get<index>(u).rtransform(detJ, J, Jinv, std::get<index>(s.u));
    rtransform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <A2D::index_t index>
  void set_value_(const A2D::index_t comp, T val) {}

  template <A2D::index_t index, class First, class... Remain>
  void set_value_(const A2D::index_t comp, T val) {
    if (comp < First::ncomp) {
      std::get<index>(u).set_value(comp, val);
    } else {
      set_value_<index + 1, Remain...>(comp - First::ncomp, val);
    }
  }

  template <A2D::index_t index>
  T get_value_(const A2D::index_t comp) const {
    return 0.0;
  }

  template <A2D::index_t index, class First, class... Remain>
  T get_value_(const A2D::index_t comp) const {
    if (comp < First::ncomp) {
      return std::get<index>(u).get_value(comp);
    } else {
      return get_value_<index + 1, Remain...>(comp - First::ncomp);
    }
  }
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_SPACE_H