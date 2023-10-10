#ifndef A2D_FE_SOLUTION_SPACE_H
#define A2D_FE_SOLUTION_SPACE_H

#include <type_traits>

#include "a2dcore.h"

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

  5. s.btransform(detJ, J, Jinv, sref): Transform the derivatives from the
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

  KOKKOS_FUNCTION L2Space() {
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
  KOKKOS_FUNCTION void zero() {
    if constexpr (C == 1) {
      u = 0.0;
    } else {
      u.zero();
    }
  }

  // Get the value of the specified component
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
    if constexpr (C == 1) {
      return u;
    } else {
      return u(comp);
    }
  }

  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
    if constexpr (C == 1) {
      return u;
    } else {
      return u(comp);
    }
  }

  // Get the scalar solution value
  KOKKOS_FUNCTION VarType& get_value() { return u; }
  KOKKOS_FUNCTION const VarType& get_value() const { return u; }

  // Transform the values from the reference to the physical space
  KOKKOS_FUNCTION void transform(const T& detJ, const Mat<T, D, D>& J,
                                 const Mat<T, D, D>& Jinv,
                                 L2Space<T, C, D>& s) const {
    s.u = u;
  }

  // Transform derivatives from the physical to the refernece space
  KOKKOS_FUNCTION void btransform(const T& detJ, const Mat<T, D, D>& J,
                                  const Mat<T, D, D>& Jinv,
                                  L2Space<T, C, D>& s) const {
    s.u = u;
  }

  template <ADorder forder, typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void forward_transform(const L2Space<T, C, D>& inp,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         L2Space<T, C, D>& outp) const {
    if constexpr (C == 1) {
      outp.u = inp.u;
    } else {
      VecCopyCore<T, C>(get_data(inp.u), get_data(outp.u));
    }
  }

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void reverse_transform(const L2Space<T, C, D>& outb,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         L2Space<T, C, D>& inb) const {
    if constexpr (C == 1) {
      inb.u += outb.u;
    } else {
      VecAddCore<T, C>(get_data(outb.u), get_data(inb.u));
    }
  }

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void hreverse_transform(const L2Space<T, C, D>& outh,
                                          const L2Space<T, C, D>& outb,
                                          const L2Space<T, C, D>& inp,
                                          dtype& detJ, JType& J, JinvType& Jinv,
                                          L2Space<T, C, D>& inh) const {
    if constexpr (C == 1) {
      inh.u += outh.u;
    } else {
      VecAddCore<T, C>(get_data(outh.u), get_data(inh.u));
    }
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

  KOKKOS_FUNCTION H1Space() { zero(); }

  // Number of solution components
  static const index_t ncomp = (D + 1) * C;

  // Spatial dimension
  static const index_t dim = D;

  // Zero the solution
  KOKKOS_FUNCTION void zero() {
    if constexpr (C == 1) {
      u = T(0);
    } else {
      u.zero();
    }
    grad.zero();
  }

  // Get the value of the specified component
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
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

  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
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
  KOKKOS_FUNCTION VarType& get_value() { return u; }
  KOKKOS_FUNCTION const VarType& get_value() const { return u; }

  // Get the gradient of the solution
  KOKKOS_FUNCTION GradType& get_grad() { return grad; }
  KOKKOS_FUNCTION const GradType& get_grad() const { return grad; }

  // Transform the values from the reference to the physical space
  KOKKOS_FUNCTION void transform(const T& detJ, const Mat<T, D, D>& J,
                                 const Mat<T, D, D>& Jinv,
                                 H1Space<T, C, D>& s) const {
    s.u = u;

    // s.grad = grad * Jinv
    if constexpr (C == 1) {
      MatVecCore<T, D, D, MatOp::TRANSPOSE>(get_data(Jinv), get_data(grad),
                                            get_data(s.grad));
    } else {
      MatMatMultCore<T, C, D, D, D, C, D>(get_data(grad), get_data(Jinv),
                                          get_data(s.grad));
    }
  }

  // Transform derivatives from the physical to the reference space
  KOKKOS_FUNCTION void btransform(const T& detJ, const Mat<T, D, D>& J,
                                  const Mat<T, D, D>& Jinv,
                                  H1Space<T, C, D>& s) const {
    // dot{s.grad} = dot{grad} Jinv =>
    // tr(bar{s.grad}^{T} * dot{s.grad})
    // = tr(bar{s.grad}^{T} * dot{grad} * Jinv)
    // = tr(Jinv * bar{s.grad}^{T} * dot{grad})
    // = tr((s.grad * Jinv^{T})^{T} * dot{grad})
    s.u = u;

    // s.grad = grad * Jinv^{T}
    if constexpr (C == 1) {
      MatVecCore<T, D, D>(get_data(Jinv), get_data(grad), get_data(s.grad));
    } else {
      MatMatMultCore<T, C, D, D, D, C, D, MatOp::NORMAL, MatOp::TRANSPOSE>(
          get_data(grad), get_data(Jinv), get_data(s.grad));
    }
  }

  template <ADorder forder, typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void forward_transform(const H1Space<T, C, D>& inp,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         H1Space<T, C, D>& outp) const {
    constexpr ADseed seed = conditional_value<ADseed, forder == ADorder::FIRST,
                                              ADseed::b, ADseed::p>::value;
    const bool additive = true;

    if constexpr (C == 1) {
      outp.u = inp.u;
      MatVecCore<T, D, D, MatOp::TRANSPOSE>(get_data(Jinv), get_data(inp.grad),
                                            get_data(outp.grad));
      MatVecCore<T, D, D, MatOp::TRANSPOSE, additive>(
          GetSeed<seed>::get_data(Jinv), get_data(grad), get_data(outp.grad));
    } else {
      VecCopyCore<T, C>(get_data(inp.u), get_data(outp.u));
      MatMatMultCore<T, C, D, D, D, C, D>(get_data(inp.grad), get_data(Jinv),
                                          get_data(outp.grad));
      MatMatMultCore<T, C, D, D, D, C, D, MatOp::NORMAL, MatOp::NORMAL,
                     additive>(get_data(grad), GetSeed<seed>::get_data(Jinv),
                               get_data(outp.grad));
    }
  }

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void reverse_transform(const H1Space<T, C, D>& outb,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         H1Space<T, C, D>& inb) const {
    const bool additive = true;

    if constexpr (C == 1) {
      inb.u += outb.u;
      MatVecCore<T, D, D, MatOp::NORMAL, additive>(
          get_data(Jinv), get_data(outb.grad), get_data(inb.grad));
      VecOuterCore<T, D, D, additive>(get_data(grad), get_data(outb.grad),
                                      GetSeed<ADseed::b>::get_data(Jinv));
    } else {
      VecAddCore<T, C>(get_data(outb.u), get_data(inb.u));
      MatMatMultCore<T, C, D, D, D, C, D, MatOp::NORMAL, MatOp::TRANSPOSE,
                     additive>(get_data(outb.grad), get_data(Jinv),
                               get_data(inb.grad));
      MatMatMultCore<T, C, D, C, D, D, D, MatOp::TRANSPOSE, MatOp::NORMAL,
                     additive>(get_data(grad), get_data(outb.grad),
                               GetSeed<ADseed::b>::get_data(Jinv));
    }
  }

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void hreverse_transform(const H1Space<T, C, D>& outh,
                                          const H1Space<T, C, D>& outb,
                                          const H1Space<T, C, D>& inp,
                                          dtype& detJ, JType& J, JinvType& Jinv,
                                          H1Space<T, C, D>& inh) const {
    const bool additive = true;

    if constexpr (C == 1) {
      inh.u += outh.u;
      MatVecCore<T, D, D, MatOp::NORMAL, additive>(
          get_data(Jinv), get_data(outh.grad), get_data(inh.grad));
      MatVecCore<T, D, D, MatOp::NORMAL, additive>(
          GetSeed<ADseed::p>::get_data(Jinv), get_data(outb.grad),
          get_data(inh.grad));

      VecOuterCore<T, D, D, additive>(get_data(grad), get_data(outh.grad),
                                      GetSeed<ADseed::h>::get_data(Jinv));
      VecOuterCore<T, D, D, additive>(get_data(inp.grad), get_data(outb.grad),
                                      GetSeed<ADseed::h>::get_data(Jinv));
    } else {
      VecAddCore<T, C>(get_data(outh.u), get_data(inh.u));
      MatMatMultCore<T, C, D, D, D, C, D, MatOp::NORMAL, MatOp::TRANSPOSE,
                     additive>(get_data(outh.grad), get_data(Jinv),
                               get_data(inh.grad));
      MatMatMultCore<T, C, D, D, D, C, D, MatOp::NORMAL, MatOp::TRANSPOSE,
                     additive>(get_data(outb.grad),
                               GetSeed<ADseed::p>::get_data(Jinv),
                               get_data(inh.grad));

      MatMatMultCore<T, C, D, C, D, D, D, MatOp::TRANSPOSE, MatOp::NORMAL,
                     additive>(get_data(grad), get_data(outh.grad),
                               GetSeed<ADseed::h>::get_data(Jinv));
      MatMatMultCore<T, C, D, C, D, D, D, MatOp::TRANSPOSE, MatOp::NORMAL,
                     additive>(get_data(inp.grad), get_data(outb.grad),
                               GetSeed<ADseed::h>::get_data(Jinv));
    }
  }

 private:
  VarType u;
  GradType grad;
};

template <typename T, index_t D>
class HdivSpace {
 public:
  KOKKOS_FUNCTION HdivSpace() { div = 0.0; }

  using VarType = Vec<T, D>;
  using DivType = T;

  // Spatial dimension
  static const index_t dim = D;

  // Number of solution components
  static const index_t ncomp = 1 + D;

  // Zero the solution
  KOKKOS_FUNCTION void zero() {
    u.zero();
    div = 0.0;
  }

  // Get the value of the specified component
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
    if (comp < D) {
      return u(comp);
    } else {
      return div;
    }
  }

  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
    if (comp < D) {
      return u(comp);
    } else {
      return div;
    }
  }

  // Get the solution value
  KOKKOS_FUNCTION Vec<T, D>& get_value() { return u; }
  KOKKOS_FUNCTION const Vec<T, D>& get_value() const { return u; }

  // Get the value of the divergence
  KOKKOS_FUNCTION T& get_div() { return div; }
  KOKKOS_FUNCTION const T& get_div() const { return div; }

  // Transform the values from the reference to the physical space
  KOKKOS_FUNCTION void transform(const T& detJ, const Mat<T, dim, dim>& J,
                                 const Mat<T, dim, dim>& Jinv,
                                 HdivSpace<T, D>& s) const {
    T inv = 1.0 / detJ;
    if (D == 2) {
      s.u(0) = inv * (J(0, 0) * u(0) + J(0, 1) * u(1));
      s.u(1) = inv * (J(1, 0) * u(0) + J(1, 1) * u(1));
    } else if (D == 3) {
      s.u(0) = inv * (J(0, 0) * u(0) + J(0, 1) * u(1) + J(0, 2) * u(2));
      s.u(1) = inv * (J(1, 0) * u(0) + J(1, 1) * u(1) + J(1, 2) * u(2));
      s.u(2) = inv * (J(2, 0) * u(0) + J(2, 1) * u(1) + J(2, 2) * u(2));
    }
    s.div = inv * div;
  }

  // Transform derivatives from the physical to the refernece space
  KOKKOS_FUNCTION void btransform(const T& detJ, const Mat<T, dim, dim>& J,
                                  const Mat<T, dim, dim>& Jinv,
                                  HdivSpace<T, D>& s) const {
    T inv = 1.0 / detJ;
    if (D == 2) {
      s.u(0) = inv * (J(0, 0) * u(0) + J(1, 0) * u(1));
      s.u(1) = inv * (J(0, 1) * u(0) + J(1, 1) * u(1));
    } else if (D == 3) {
      s.u(0) = inv * (J(0, 0) * u(0) + J(1, 0) * u(1) + J(2, 0) * u(2));
      s.u(1) = inv * (J(0, 1) * u(0) + J(1, 1) * u(1) + J(2, 1) * u(2));
      s.u(2) = inv * (J(0, 2) * u(0) + J(1, 2) * u(1) + J(2, 2) * u(2));
    }
    s.div = inv * div;
  }

  template <ADorder forder, typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void forward_transform(const HdivSpace<T, D>& inp,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         HdivSpace<T, D>& outp) const {}

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void reverse_transform(const HdivSpace<T, D>& outb,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         HdivSpace<T, D>& inb) const {
    // const bool additive = true;
    // const Mat<T, D, D>& J = J0.value();

    // T inv = 1.0 / detJ.value();
    // T binv = s.div * div;
    // if (D == 2) {
    //   binv += (s.u(0) * (J(0, 0) * u(0) + J(0, 1) * u(1)) +
    //            s.u(1) * (J(1, 0) * u(0) + J(1, 1) * u(1)));
    // } else if (D == 3) {
    //   binv += (s.u(0) * (J(0, 0) * u(0) + J(0, 1) * u(1) + J(0, 2) * u(2)) +
    //            s.u(1) * (J(1, 0) * u(0) + J(1, 1) * u(1) + J(1, 2) * u(2)) +
    //            s.u(2) * (J(2, 0) * u(0) + J(2, 1) * u(1) + J(2, 2) * u(2)));
    // }
    // detJ.bvalue() -= inv * inv * binv;
    // VecOuterCore<T, D, D, additive>(inv, get_data(s.u), get_data(u),
    //                                 GetSeed<ADseed::b>::get_data(J0));
  }

  template <typename dtype, class JType, class JinvType>
  KOKKOS_FUNCTION void hreverse_transform(const HdivSpace<T, D>& outh,
                                          const HdivSpace<T, D>& outb,
                                          const HdivSpace<T, D>& inp,
                                          dtype& detJ, JType& J, JinvType& Jinv,
                                          HdivSpace<T, D>& inh) const {}

 private:
  VarType u;
  DivType div;
};

template <typename T>
class Hcurl2DSpace {
 public:
  KOKKOS_FUNCTION Hcurl2DSpace() { curl = 0.0; }

  // Number of solution components
  static const index_t ncomp = 3;

  // Spatial dimension
  static const index_t dim = 2;

  // Zero the solution
  KOKKOS_FUNCTION void zero() {
    u.zero();
    curl = 0.0;
  }

  // Get the value of the specified component
  template <typename I>
  KOKKOS_FUNCTION T& operator[](const I comp) {
    if (comp < 2) {
      return u(comp);
    } else {
      return curl;
    }
  }

  template <typename I>
  KOKKOS_FUNCTION const T& operator[](const I comp) const {
    if (comp < 2) {
      return u(comp);
    } else {
      return curl;
    }
  }

  // Get the solution value
  KOKKOS_FUNCTION Vec<T, 2>& get_value() { return u; }
  KOKKOS_FUNCTION const Vec<T, 2>& get_value() const { return u; }

  // Get the curl (2D so it's a scalar value)
  KOKKOS_FUNCTION T& get_curl() { return curl; }
  KOKKOS_FUNCTION const T& get_curl() const { return curl; }

  // Transform the values from the reference to the physical space
  KOKKOS_FUNCTION void transform(const T& detJ, const Mat<T, dim, dim>& J,
                                 const Mat<T, dim, dim>& Jinv,
                                 Hcurl2DSpace<T>& s) const {
    s.u = u;
    s.curl = curl;
  }

  // Transform derivatives from the physical to the refernece space
  KOKKOS_FUNCTION void btransform(const T& detJ, const Mat<T, dim, dim>& J,
                                  const Mat<T, dim, dim>& Jinv,
                                  Hcurl2DSpace<T>& s) const {
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

template <index_t index, class First, class... Spaces>
struct get_fespace_type;

template <index_t index, class First, class... Remain>
struct get_fespace_type : get_fespace_type<index - 1, Remain...> {};

template <class First, class... Remain>
struct get_fespace_type<0, First, Remain...> {
  typedef First type;
};

/*
  FESpace:

  To form a full finite-element space, we use a variadic template class FESpace,
  FESpace takes the scalar type, the physical dimension of the problem 2D, 3D
  etc. and a list of the solution space objects. For instance, this class might
  look like this when declared:

  FESpace<double, 3, H1Space, L2Space, HdivSpace>

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
  // Tuple type for storing the solution space object
  typedef std::tuple<Spaces...> SolutionSpace;

  // Numeric type used for the space
  typedef T type;

  // Type of AD object
  static const ADObjType obj_type = ADObjType::VECTOR;

  // Number of spaces
  static constexpr index_t nspaces = sizeof...(Spaces);

  /*
    The dimension of the problem
  */
  static constexpr index_t dim = D;

  /*
    Count up the total number of degrees of freedom
  */
  static constexpr index_t ncomp = __count_space_components<Spaces...>::ncomp;

  /*
    Zero all the values in all the spaces
  */
  KOKKOS_FUNCTION void zero() {
    if constexpr (sizeof...(Spaces) > 0) {
      zero_<0, Spaces...>();
    }
  }

  /*
    Get a solution value based on the index
  */
  KOKKOS_FUNCTION T& operator[](const index_t comp) {
    return get_value_<0, Spaces...>(comp);
  }
  KOKKOS_FUNCTION const T& operator[](const index_t comp) const {
    return get_value_<0, Spaces...>(comp);
  }

  /*
    Copy the values from another FESpace object
  */
  KOKKOS_FUNCTION void copy(const FESpace<T, D, Spaces...>& src) {
    if constexpr (ncomp > 0) {
      for (index_t i = 0; i < ncomp; i++) {
        (*this)[i] = src[i];
      }
    }
  }

  /*
    Extract the specified solution space from the FESpace object
  */
  template <index_t index>
  KOKKOS_FUNCTION typename std::tuple_element<index, SolutionSpace>::type&
  get() noexcept {
    return std::get<index>(u);
  }
  template <index_t index>
  KOKKOS_FUNCTION const typename std::tuple_element<index, SolutionSpace>::type&
  get() const noexcept {
    return std::get<index>(u);
  }

  /*
    Transform the solution space from the reference element to the physical
    element
  */
  KOKKOS_FUNCTION void transform(const T& detJ, const Mat<T, D, D>& J,
                                 const Mat<T, D, D>& Jinv,
                                 FESpace<T, D, Spaces...>& s) const {
    transform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  /*
    Perform the reverse of the transform - transfer the derivative from the
    physical element to the reference element
  */
  KOKKOS_FUNCTION void btransform(const T& detJ, const Mat<T, D, D>& J,
                                  const Mat<T, D, D>& Jinv,
                                  FESpace<T, D, Spaces...>& s) const {
    btransform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  template <ADorder forder, typename dtype, class JType, class JinvType,
            std::enable_if_t<
                (is_scalar_type<typename remove_a2dobj<dtype>::type>::value &&
                 D == get_matrix_rows<JType>::size &&
                 D == get_matrix_columns<JType>::size &&
                 D == get_matrix_rows<JinvType>::size &&
                 D == get_matrix_columns<JinvType>::size),
                bool> = true>
  KOKKOS_FUNCTION void forward_transform(const FESpace<T, D, Spaces...>& inp,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         FESpace<T, D, Spaces...>& outp) const {
    forward_transform_<0, forder, dtype, JType, JinvType, Spaces...>(
        inp, detJ, J, Jinv, outp);
  }

  template <typename dtype, class JType, class JinvType,
            std::enable_if_t<
                (is_scalar_type<typename remove_a2dobj<dtype>::type>::value &&
                 D == get_matrix_rows<JType>::size &&
                 D == get_matrix_columns<JType>::size &&
                 D == get_matrix_rows<JinvType>::size &&
                 D == get_matrix_columns<JinvType>::size),
                bool> = true>
  KOKKOS_FUNCTION void reverse_transform(const FESpace<T, D, Spaces...>& outb,
                                         dtype& detJ, JType& J, JinvType& Jinv,
                                         FESpace<T, D, Spaces...>& inb) const {
    reverse_transform_<0, dtype, JType, JinvType, Spaces...>(outb, detJ, J,
                                                             Jinv, inb);
  }

  template <typename dtype, class JType, class JinvType,
            std::enable_if_t<
                (is_scalar_type<typename remove_a2dobj<dtype>::type>::value &&
                 D == get_matrix_rows<JType>::size &&
                 D == get_matrix_columns<JType>::size &&
                 D == get_matrix_rows<JinvType>::size &&
                 D == get_matrix_columns<JinvType>::size),
                bool> = true>
  KOKKOS_FUNCTION void hreverse_transform(const FESpace<T, D, Spaces...>& outh,
                                          const FESpace<T, D, Spaces...>& outb,
                                          const FESpace<T, D, Spaces...>& inp,
                                          dtype& detJ, JType& J, JinvType& Jinv,
                                          FESpace<T, D, Spaces...>& inh) const {
    hreverse_transform_<0, dtype, JType, JinvType, Spaces...>(
        outh, outb, inp, detJ, J, Jinv, inh);
  }

 private:
  // Solution space tuple object
  SolutionSpace u;

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION void transform_(const T& detJ, const Mat<T, D, D>& J,
                                  const Mat<T, D, D>& Jinv,
                                  FESpace<T, D, Spaces...>& s) const {
    std::get<index>(u).transform(detJ, J, Jinv, std::get<index>(s.u));
    if constexpr (sizeof...(Remain) > 0) {
      transform_<index + 1, Remain...>(detJ, J, Jinv, s);
    }
  }

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION void btransform_(const T& detJ, const Mat<T, D, D>& J,
                                   const Mat<T, D, D>& Jinv,
                                   FESpace<T, D, Spaces...>& s) const {
    std::get<index>(u).btransform(detJ, J, Jinv, std::get<index>(s.u));
    if constexpr (sizeof...(Remain) > 0) {
      btransform_<index + 1, Remain...>(detJ, J, Jinv, s);
    }
  }

  template <index_t index, ADorder forder, typename dtype, class JType,
            class JinvType, class First, class... Remain>
  KOKKOS_FUNCTION void forward_transform_(
      const FESpace<T, D, Spaces...>& inp, dtype& detJ, JType& J,
      JinvType& Jinv, FESpace<T, D, Spaces...>& outp) const {
    std::get<index>(u).template forward_transform<forder>(
        std::get<index>(inp.u), detJ, J, Jinv, std::get<index>(outp.u));
    if constexpr (sizeof...(Remain) > 0) {
      forward_transform_<index + 1, forder, dtype, JType, JinvType, Remain...>(
          inp, detJ, J, Jinv, outp);
    }
  }

  template <index_t index, typename dtype, class JType, class JinvType,
            class First, class... Remain>
  KOKKOS_FUNCTION void reverse_transform_(const FESpace<T, D, Spaces...>& outb,
                                          dtype& detJ, JType& J, JinvType& Jinv,
                                          FESpace<T, D, Spaces...>& inb) const {
    std::get<index>(u).reverse_transform(std::get<index>(outb.u), detJ, J, Jinv,
                                         std::get<index>(inb.u));
    if constexpr (sizeof...(Remain) > 0) {
      reverse_transform_<index + 1, dtype, JType, JinvType, Remain...>(
          outb, detJ, J, Jinv, inb);
    }
  }

  template <index_t index, typename dtype, class JType, class JinvType,
            class First, class... Remain>
  KOKKOS_FUNCTION void hreverse_transform_(
      const FESpace<T, D, Spaces...>& outh,
      const FESpace<T, D, Spaces...>& outb, const FESpace<T, D, Spaces...>& inp,
      dtype& detJ, JType& J, JinvType& Jinv,
      FESpace<T, D, Spaces...>& inh) const {
    std::get<index>(u).hreverse_transform(
        std::get<index>(outh.u), std::get<index>(outb.u),
        std::get<index>(inp.u), detJ, J, Jinv, std::get<index>(inh.u));
    if constexpr (sizeof...(Remain) > 0) {
      hreverse_transform_<index + 1, dtype, JType, JinvType, Remain...>(
          outh, outb, inp, detJ, J, Jinv, inh);
    }
  }

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION T& get_value_(const index_t comp) {
    if constexpr (sizeof...(Remain) == 0) {
      return std::get<index>(u)[comp];
    } else if (comp < First::ncomp) {
      return std::get<index>(u)[comp];
    } else {
      return get_value_<index + 1, Remain...>(comp - First::ncomp);
    }
  }

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION const T& get_value_(const index_t comp) const {
    if constexpr (sizeof...(Remain) == 0) {
      return std::get<index>(u)[comp];
    } else if (comp < First::ncomp) {
      return std::get<index>(u)[comp];
    } else {
      return get_value_<index + 1, Remain...>(comp - First::ncomp);
    }
  }

  template <index_t index, class First, class... Remain>
  KOKKOS_FUNCTION void zero_() {
    std::get<index>(u).zero();
    if constexpr (sizeof...(Remain) > 0) {
      zero_<index + 1, Remain...>();
    }
  }
};

/**
 * @brief Unpack the value from the FESpace object
 *
 * @tparam index Index of the HdivSpace, L2Space or H1Space object in FESpace
 * @tparam T Deduced scalar type
 * @tparam D Deduced spatial dimension
 * @tparam Spaces Deduced space parameters
 * @param obj Input FESpace object
 * @return VarType&, ADObj<VarType&> or A2DObj<VarType&>
 */
template <index_t index, typename T, index_t D, class... Spaces>
const typename get_fespace_type<index, Spaces...>::type::VarType& get_value(
    const FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_value();
}

template <index_t index, typename T, index_t D, class... Spaces>
typename get_fespace_type<index, Spaces...>::type::VarType& get_value(
    FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_value();
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_value(ADObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::VarType;
  return ADObj<type&>(obj.value().template get<index>().get_value(),
                      obj.bvalue().template get<index>().get_value());
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_value(A2DObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::VarType;
  return A2DObj<type&>(obj.value().template get<index>().get_value(),
                       obj.bvalue().template get<index>().get_value(),
                       obj.pvalue().template get<index>().get_value(),
                       obj.hvalue().template get<index>().get_value());
}

/**
 * @brief Unpack the gradient from the FESpace object
 *
 * @tparam index Index of the H1Space object in FESpace
 * @tparam T Deduced scalar type
 * @tparam D Deduced spatial dimension
 * @tparam Spaces Deduced space parameters
 * @param obj Input FESpace object
 * @return GradType&, ADObj<GradType&> or A2DObj<GradType&>
 */
template <index_t index, typename T, index_t D, class... Spaces>
const typename get_fespace_type<index, Spaces...>::type::GradType& get_grad(
    const FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_grad();
}

template <index_t index, typename T, index_t D, class... Spaces>
typename get_fespace_type<index, Spaces...>::type::GradType& get_grad(
    FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_grad();
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_grad(ADObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::GradType;
  return ADObj<type&>(obj.value().template get<index>().get_grad(),
                      obj.bvalue().template get<index>().get_grad());
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_grad(A2DObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::GradType;
  return A2DObj<type&>(obj.value().template get<index>().get_grad(),
                       obj.bvalue().template get<index>().get_grad(),
                       obj.pvalue().template get<index>().get_grad(),
                       obj.hvalue().template get<index>().get_grad());
}

/**
 * @brief Unpack the divergence from the FESpace object
 *
 * @tparam index Index of the HdivSpace object in FESpace
 * @tparam T Deduced scalar type
 * @tparam D Deduced spatial dimension
 * @tparam Spaces Deduced space parameters
 * @param obj Input FESpace object
 * @return DivType&, ADObj<DivType&> or A2DObj<DivType&>
 */
template <index_t index, typename T, index_t D, class... Spaces>
const typename get_fespace_type<index, Spaces...>::type::DivType& get_div(
    const FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_div();
}

template <index_t index, typename T, index_t D, class... Spaces>
typename get_fespace_type<index, Spaces...>::type::DivType& get_div(
    FESpace<T, D, Spaces...>& obj) {
  return obj.template get<index>().get_div();
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_div(ADObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::DivType;
  return ADObj<type&>(obj.value().template get<index>().get_div(),
                      obj.bvalue().template get<index>().get_div());
}

template <index_t index, typename T, index_t D, class... Spaces>
auto get_div(A2DObj<FESpace<T, D, Spaces...>>& obj) {
  using type = typename get_fespace_type<index, Spaces...>::type::DivType;
  return A2DObj<type&>(obj.value().template get<index>().get_div(),
                       obj.bvalue().template get<index>().get_div(),
                       obj.pvalue().template get<index>().get_div(),
                       obj.hvalue().template get<index>().get_div());
}

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
  KOKKOS_FUNCTION FiniteElementSpace& get(const index_t index) {
    return space[index];
  }
  KOKKOS_FUNCTION const FiniteElementSpace& get(const index_t index) const {
    return space[index];
  }

 private:
  FiniteElementSpace space[Quadrature::num_quad_points];
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_SPACE_H