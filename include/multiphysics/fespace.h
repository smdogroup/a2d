#ifndef A2D_FE_SOLUTION_SPACE_H
#define A2D_FE_SOLUTION_SPACE_H

#include <type_traits>

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

  1. static const index_t ncomp: Number of component outputs treating solution
  and derivative information as independent.

  2. static const index_t dim: Number of spatial dimensions

  3. zero(): Zero all component values for this space

  3. get_value(comp): Returns a reference or constant reference to the value in
  the given component index. The comp value must satisfy 0 <= comp < ncomps.
  Behavior undefined otherwise.

  3. The space object must implement getters for each of the specific values
  that it computes such as:

  get_value(), get_grad(), get_div(), get_curl()

  must be defined according to the actual finite element space.

  4. sref.transform(detJ, J, Jinv, s): This function transforms the current
  sref object from the reference element to the physical element.

  5. s.rtransform(detJ, J, Jinv, sref): Transform the derivatives from the
  physical element back to the reference element.
*/

/**
 * @brief The L2 space, which has (scalar or vector) values only, no derivatives
 *
 * @tparam C dimension of the solution variable associated with this space (1 or
 *           more), note: if C == 1, solution value u has type T, otherwise it
 *           has type Vec<T, C>.
 * @tparam D geometrical dimension (usually 2 or 3)
 */
template <typename T, index_t C, index_t D>
class L2Space {
 public:
  using VarType = typename std::conditional<C == 1, T, Vec<T, C>>::type;

  L2Space() {
    if constexpr (C == 1) {
      u = 0.0;
    } else {
      u.zero();
    }
  }

  // Number of solution components
  static const index_t ncomp = C;

  // Spatial dimension
  static const index_t dim = D;

  // Zero the solution
  void zero() {
    if constexpr (C == 1) {
      u = 0.0;
    } else {
      u.zero();
    }
  }

  // Get the value of the specified component
  T& get_value(const index_t comp) {
    if constexpr (C == 1) {
      return u;
    } else {
      return u(comp);
    }
  }

  const T& get_value(const index_t comp) const {
    if constexpr (C == 1) {
      return u;
    } else {
      return u(comp);
    }
  }

  // Get the scalar solution value
  VarType& get_value() { return u; }
  const VarType& get_value() const { return u; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const Mat<T, D, D>& J, const Mat<T, D, D>& Jinv,
                 L2Space<T, C, D>& s) const {
    s.u = u;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const Mat<T, D, D>& J,
                  const Mat<T, D, D>& Jinv, L2Space<T, C, D>& s) const {
    s.u = u;
  }

 private:
  VarType u;
};

template <typename T, index_t C, index_t D>
class H1Space {
 public:
  using VarType = typename std::conditional<C == 1, T, Vec<T, C>>::type;
  using GradType =
      typename std::conditional<C == 1, Vec<T, D>, Mat<T, C, D>>::type;

  H1Space() { zero(); }

  // Number of solution components
  static const index_t ncomp = (D + 1) * C;

  // Spatial dimension
  static const index_t dim = D;

  // Zero the solution
  void zero() {
    if constexpr (C == 1) {
      u = T(0);
    } else {
      u.zero();
    }
    grad.zero();
  }

  // Get the value of the specified component
  T& get_value(const index_t comp) {
    if constexpr (C == 1) {
      if (comp == 0) {
        return u;
      } else {
        return grad(comp - 1);
      }
    } else {
      if (comp % (D + 1) == 0) {
        return u(comp / (D + 1));
      } else {
        return grad(comp / (D + 1), (comp % (D + 1)) - 1);
      }
    }
  }

  const T& get_value(const index_t comp) const {
    if constexpr (C == 1) {
      if (comp == 0) {
        return u;
      } else {
        return grad(comp - 1);
      }
    } else {
      if (comp % (D + 1) == 0) {
        return u(comp / (D + 1));
      } else {
        return grad(comp / (D + 1), (comp % (D + 1)) - 1);
      }
    }
  }

  // Get the value of the solution
  VarType& get_value() { return u; }
  const VarType& get_value() const { return u; }

  // Get the gradient of the solution
  GradType& get_grad() { return grad; }
  const GradType& get_grad() const { return grad; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const Mat<T, D, D>& J, const Mat<T, D, D>& Jinv,
                 H1Space<T, C, D>& s) const {
    s.u = u;

    // s.grad = grad * Jinv
    if constexpr (C == 1) {
      for (index_t i = 0; i < dim; i++) {
        s.grad(i) = 0.0;
        for (index_t j = 0; j < dim; j++) {
          s.grad(i) += grad(j) * Jinv(j, i);
        }
      }
    } else {
      MatMatMult(grad, Jinv, s.grad);
    }
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const Mat<T, D, D>& J,
                  const Mat<T, D, D>& Jinv, H1Space<T, C, D>& s) {
    // dot{s.grad} = dot{grad} Jinv =>
    // tr(bar{s.grad}^{T} * dot{s.grad})
    // = tr(bar{s.grad}^{T} * dot{grad} * Jinv)
    // = tr(Jinv * bar{s.grad}^{T} * dot{grad})
    // = tr((s.grad * Jinv^{T})^{T} * dot{grad})
    s.u = u;

    // s.grad = grad * Jinv^{T}
    if constexpr (C == 1) {
      for (index_t i = 0; i < dim; i++) {
        s.grad(i) = 0.0;
        for (index_t j = 0; j < dim; j++) {
          s.grad(i) += Jinv(i, j) * grad(j);
        }
      }
    } else {
      for (index_t i = 0; i < C; i++) {
        for (index_t j = 0; j < dim; j++) {
          s.grad(i, j) = 0.0;

          for (index_t k = 0; k < dim; k++) {
            s.grad(i, j) += grad(i, k) * Jinv(j, k);
          }
        }
      }
      // MatMatMult<T, false, true>(grad, Jinv, s.grad);
    }
  }

 private:
  VarType u;
  GradType grad;
};

template <typename T, index_t D>
class HdivSpace {
 public:
  HdivSpace() { div = 0.0; }

  // Spatial dimension
  static const index_t dim = D;

  // Number of solution components
  static const index_t ncomp = 1 + D;

  // Zero the solution
  void zero() {
    u.zero();
    div = 0.0;
  }

  // Get the value of the specified component
  T& get_value(const index_t comp) {
    if (comp < D) {
      return u(comp);
    } else {
      return div;
    }
  }
  const T& get_value(const index_t comp) const {
    if (comp < D) {
      return u(comp);
    } else {
      return div;
    }
  }

  // Get the solution value
  Vec<T, D>& get_value() { return u; }
  const Vec<T, D>& get_value() const { return u; }

  // Get the value of the divergence
  T& get_div() { return div; }
  const T& get_div() const { return div; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const Mat<T, dim, dim>& J,
                 const Mat<T, dim, dim>& Jinv, HdivSpace<T, D>& s) const {
    if (D == 2) {
      s.u(0) = (J(0, 0) * u(0) + J(0, 1) * u(1)) / detJ;
      s.u(1) = (J(1, 0) * u(0) + J(1, 1) * u(1)) / detJ;
    } else if (D == 3) {
      s.u(0) = (J(0, 0) * u(0) + J(0, 1) * u(1) + J(0, 2) * u(2)) / detJ;
      s.u(1) = (J(1, 0) * u(0) + J(1, 1) * u(1) + J(1, 2) * u(2)) / detJ;
      s.u(2) = (J(2, 0) * u(0) + J(2, 1) * u(1) + J(2, 2) * u(2)) / detJ;
    }
    s.div = div / detJ;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const Mat<T, dim, dim>& J,
                  const Mat<T, dim, dim>& Jinv, HdivSpace<T, D>& s) const {
    if (D == 2) {
      s.u(0) = (J(0, 0) * u(0) + J(1, 0) * u(1)) / detJ;
      s.u(1) = (J(0, 1) * u(0) + J(1, 1) * u(1)) / detJ;
    } else if (D == 3) {
      s.u(0) = (J(0, 0) * u(0) + J(1, 0) * u(1) + J(2, 0) * u(2)) / detJ;
      s.u(1) = (J(0, 1) * u(0) + J(1, 1) * u(1) + J(2, 1) * u(2)) / detJ;
      s.u(2) = (J(0, 2) * u(0) + J(1, 2) * u(1) + J(2, 2) * u(2)) / detJ;
    }
    s.div = div / detJ;
  }

 private:
  Vec<T, D> u;
  T div;
};

template <typename T>
class Hcurl2DSpace {
 public:
  Hcurl2DSpace() { curl = 0.0; }

  // Number of solution components
  static const index_t ncomp = 3;

  // Spatial dimension
  static const index_t dim = 2;

  // Zero the solution
  void zero() {
    u.zero();
    curl = 0.0;
  }

  // Get the value of the specified component
  T& get_value(const index_t comp) {
    if (comp < 2) {
      return u(comp);
    } else {
      return curl;
    }
  }
  const T& get_value(const index_t comp) const {
    if (comp < 2) {
      return u(comp);
    } else {
      return curl;
    }
  }

  // Get the solution value
  Vec<T, 2>& get_value() { return u; }
  const Vec<T, 2>& get_value() const { return u; }

  // Get the curl (2D so it's a scalar value)
  T& get_curl() { return curl; }
  const T& get_curl() const { return curl; }

  // Transform the values from the reference to the physical space
  void transform(const T& detJ, const Mat<T, dim, dim>& J,
                 const Mat<T, dim, dim>& Jinv, Hcurl2DSpace<T>& s) const {
    s.u = u;
    s.curl = curl;
  }

  // Transform derivatives from the physical to the refernece space
  void rtransform(const T& detJ, const Mat<T, dim, dim>& J,
                  const Mat<T, dim, dim>& Jinv, Hcurl2DSpace<T>& s) const {
    s.u = u;
    s.curl = curl;
  }

 private:
  Vec<T, dim> u;
  T curl;
};

template <class... Spaces>
struct __count_space_components;

template <>
struct __count_space_components<> {
  static const index_t ncomp = 0;
};

template <class First, class... Remain>
struct __count_space_components<First, Remain...> {
  static const index_t ncomp =
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
template <typename T, index_t D, class... Spaces>
class FESpace {
 public:
  typedef std::tuple<Spaces...> SolutionSpace;
  static constexpr index_t nspaces = std::tuple_size<std::tuple<Spaces...>>();

  /*
    Count up the total number of degrees of freedom
  */
  static constexpr index_t ncomp = __count_space_components<Spaces...>::ncomp;

  /*
    Zero all the values in all the spaces
  */
  void zero() { zero_<0, Spaces...>(); }

  /*
    Get a solution value based on the index
  */
  T& operator[](const index_t comp) { return get_value_<0, Spaces...>(comp); }
  const T& operator[](const index_t comp) const {
    return get_value_<0, Spaces...>(comp);
  }

  /*
    Extract the specified solution space from the FESpace object
  */
  template <index_t index>
  typename std::tuple_element<index, SolutionSpace>::type& get() noexcept {
    return std::get<index>(u);
  }
  template <index_t index>
  const typename std::tuple_element<index, SolutionSpace>::type& get()
      const noexcept {
    return std::get<index>(u);
  }

  /*
    Transform the solution space from the reference element to the physical
    element
  */
  void transform(const T& detJ, const Mat<T, D, D>& J, const Mat<T, D, D>& Jinv,
                 FESpace<T, D, Spaces...>& s) const {
    transform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  /*
    Perform the reverse of the transform - transfer the derivative from the
    physical element to the reference element
  */
  void rtransform(const T& detJ, const Mat<T, D, D>& J,
                  const Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    rtransform_<0, Spaces...>(detJ, J, Jinv, s);
  }

 private:
  // Solution space tuple object
  SolutionSpace u;

  template <index_t index>
  void transform_(const T& detJ, const Mat<T, D, D>& J,
                  const Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) const {
  }

  template <index_t index, class First, class... Remain>
  void transform_(const T& detJ, const Mat<T, D, D>& J,
                  const Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) const {
    std::get<index>(u).transform(detJ, J, Jinv, std::get<index>(s.u));
    transform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <index_t index>
  void rtransform_(const T& detJ, const Mat<T, D, D>& J,
                   const Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {}

  template <index_t index, class First, class... Remain>
  void rtransform_(const T& detJ, const Mat<T, D, D>& J,
                   const Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    std::get<index>(u).rtransform(detJ, J, Jinv, std::get<index>(s.u));
    rtransform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <index_t index, class First, class... Remain>
  T& get_value_(const index_t comp) {
    if constexpr (sizeof...(Remain) == 0) {
      return std::get<index>(u).get_value(comp);
    } else if (comp < First::ncomp) {
      return std::get<index>(u).get_value(comp);
    } else {
      return get_value_<index + 1, Remain...>(comp - First::ncomp);
    }
  }

  template <index_t index, class First, class... Remain>
  const T& get_value_(const index_t comp) const {
    if constexpr (sizeof...(Remain) == 0) {
      return std::get<index>(u).get_value(comp);
    } else if (comp < First::ncomp) {
      return std::get<index>(u).get_value(comp);
    } else {
      return get_value_<index + 1, Remain...>(comp - First::ncomp);
    }
  }

  template <index_t index, class First, class... Remain>
  void zero_() {
    std::get<index>(u).zero();
    if constexpr (sizeof...(Remain) > 0) {
      zero_<index + 1, Remain...>();
    }
  }
};

/**
 * @brief A collection of finite-element space objects associated with the
 * quadrature points in an element
 *
 * @tparam Quadrature The Quadrature object
 * @tparam FiniteElementSpace The FESpace object
 */
template <class Quadrature, class FiniteElementSpace>
class QptSpace {
 public:
  FiniteElementSpace& get(const index_t index) { return space[index]; }
  const FiniteElementSpace& get(const index_t index) const {
    return space[index];
  }

 private:
  FiniteElementSpace space[Quadrature::num_quad_points];
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_SPACE_H