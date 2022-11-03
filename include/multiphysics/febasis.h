#ifndef A2D_FE_BASIS_H
#define A2D_FE_BASIS_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"
#include "multiphysics/fespace.h"

namespace A2D {

/*
  The Basis class type is a class with all const static member data and
  static member functions.

  The basis classes provide:

  1. static const ndof: The number of degrees of freedom for the basis function.

  2. static const ncomp: The number of components in the function space output.

  3. static function interp: This function evaluates a function space object at
  a quadrature point. This function computes both values and derivatives that
  may be used in the weak form of the governing equations. This can be written
  as

  u = N(pt) * dof

  where u may be the solution as well as the derivatives of the solution, the
  divergence or curl. Here pt is a quadrature point in the reference element
  space.

  4 .static function add: This function adds the residual contribution from the
  corresponding function space object to the degrees of freedom. This can be
  written as

  dof += N(pt)^T * u

  where the dimension of u and N are the same as above.

  5. stride: Indicates that the basis functions are arranged in a special
  structure.

  6. basis: This function provides the full interpolation matrix that takes the
  degrees of freedom and computes the solution. To save space, the structure of
  the matrix may be leveraged for computational efficiency. For instance when
  the static member stride = 2, then the

  dof = [ u1    v1    u2     v2    u3    v3    u4    v4  ]

  u   = [ N1    0     N2     0     N3    0     N4    0   ]
  u,x = [ N1,x  0     N2,x   0     N3,x  0     N4,x  0   ]
  u,y = [ N1,y  0     N2,y   0     N3,y  0     N4,y  0   ]
  u,z = [ N1,z  0     N2,z   0     N3,y  0     N4,z  0   ]
  u   = [ 0     N1    0      N2   0      N3    0     N4  ]
  u,x = [ 0     N1,x  0      N2,x 0      N3,x  0     N4,x]
  u,y = [ 0     N1,y  0      N2,y 0      N3,y  0     N4,y]
  u,z = [ 0     N1,z  0      N2,z 0      N3,y  0     N4,z]

  Note this may be stored as a 4 x 4 matrix, rather than an 8 x 8 matrix.

  In this case, a call to the function basis returns:

  N =
  [ N1    N2    N3    N4
    N1,x  N2,x  N3,x  N4,x
    N1,x  N2,x  N3,x  N4,x
    N1,x  N2,x  N3,x  N4,x ]

  For an H(div) space in 2D, this would look something like this:

  dof = [ dof1   dof2   dof3   dof4   dof5]

  u   = [ N1     N2     N3     N4     N5  ]
  v   = [ N6     N7     N8     N9     N10 ]
  div = [ N11    N12    N13    N14    N15 ]

  where N1 through N15 are each separate functions.

  In this case, stride = 1.

  Note: The entries in the matrix must be consistent in ordering with the
  set_seed()/get_value() calls from the function space objects.
*/

/*
  Lagrange basis for a triangle
*/
template <typename T>
class LagrangeTri0Scalar {
 public:
  static const A2D::index_t ndof = 1;
  static const A2D::index_t ncomp = A2D::L2ScalarSpace<T, 2>::ncomp;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType sol,
                     A2D::L2ScalarSpace<T, 2>& out) {
    T& u = out.get_value();
    u = sol[0];
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::L2ScalarSpace<T, 2>& in,
                  SolnType res) {
    const T& u = in.get_value();
    res[0] += u;
  }

  // Set the matrix stride
  static const A2D::index_t stride = 1;

  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType N) {
    N[0] = 1.0;
  }
};

/*
  Lagrange basis for a triangle, TODO: finish the implementation
*/
template <typename T, A2D::index_t C>
class LagrangeTri0 {
 public:
  static const A2D::index_t ndof = C;
  static const A2D::index_t ncomp = A2D::L2Space<T, C, 2>::ncomp;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType sol,
                     A2D::L2Space<T, C, 2>& out) {
    A2D::Vec<T, 2>& u = out.get_value();
    u(0) = sol[0];
    u(1) = sol[1];
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::L2Space<T, C, 2>& in,
                  SolnType res) {
    const A2D::Vec<T, 2>& u = in.get_value();
    res[0] += u(0);
    res[1] += u(1);
  }

  // Set the matrix stride
  static const A2D::index_t stride = 1;

  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType N) {
    N[0] = 1.0;
  }
};

/*
  Lagrange basis for a triangle
*/
template <typename T, A2D::index_t C>
class LagrangeTri1 {
 public:
  static const A2D::index_t ndof = 3 * C;
  static const A2D::index_t ncomp = A2D::H1Space<T, C, 2>::ncomp;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType sol,
                     A2D::H1Space<T, C, 2>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    A2D::Vec<T, C>& u = out.get_value();
    A2D::Mat<T, C, 2>& grad = out.get_grad();

    for (A2D::index_t i = 0; i < C; i++) {
      u(i) = N[0] * sol[i] + N[1] * sol[C + i] + N[2] * sol[2 * C + i];
      grad(i, 0) = sol[C + i] - sol[i];
      grad(i, 1) = sol[2 * C + i] - sol[i];
    }
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::H1Space<T, C, 2>& in,
                  SolnType res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    const A2D::Vec<T, C>& u = in.get_value();
    const A2D::Mat<T, C, 2>& grad = in.get_grad();

    for (A2D::index_t i = 0; i < C; i++) {
      res[i] += N[0] * u(i) - grad(i, 0) - grad(i, 1);
      res[C + i] += N[1] * u(i) + grad(i, 0);
      res[2 * C + i] += N[2] * u(i) + grad(i, 1);
    }
  }

  // Set the matrix stride
  static const A2D::index_t stride = C;

  // Compute the full matrix of basis functions
  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType N) {
    double pt[2];
    Quadrature::get_point(n, pt);

    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    N[3] = -1.0;
    N[4] = 1.0;
    N[5] = 0.0;

    N[6] = -1.0;
    N[7] = 0.0;
    N[8] = 1.0;
  }
};

/*
  Raviart-Thomas element for H(div) in 2D
*/
template <typename T>
class RT2DTri1 {
 public:
  static const A2D::index_t ndof = 3;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType sol,
                     A2D::Hdiv2DSpace<T>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    A2D::Vec<T, 2>& u = out.get_value();
    T& div = out.get_div();

    u(0) = pt[0] * sol[0] + (pt[0] - 1.0) * sol[1] + pt[0] * sol[2];
    u(1) = pt[1] * sol[0] + pt[1] * sol[1] + (pt[1] - 1.0) * sol[2];
    div = 2.0 * (sol[0] + sol[1] + sol[2]);
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::Hdiv2DSpace<T>& in, SolnType res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    const A2D::Vec<T, 2>& u = in.get_value();
    const T& div = in.get_div();

    res[0] += pt[0] * u(0) + pt[1] * u(1) + 2.0 * div;
    res[1] += (pt[0] - 1.0) * u(0) + pt[1] * u(1) + 2.0 * div;
    res[2] += pt[0] * u(0) + (pt[1] - 1.0) * u(1) + 2.0 * div;
  }

  // Set the matrix stride
  static const A2D::index_t stride = 1;

  // Set the matrix

  // Set the matrix size
  static const A2D::index_t basis_size = 9;

  // Compute the full matrix of basis functions
  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType N) {
    double pt[2];
    Quadrature::get_point(n, pt);

    N[0] = pt[0];
    N[1] = (pt[0] - 1.0);
    N[2] = pt[0];

    N[3] = pt[1];
    N[4] = pt[1];
    N[5] = (pt[1] - 1.0);

    N[6] = 2.0;
    N[7] = 2.0;
    N[8] = 2.0;
  }
};

/*
  Template pattern for computing the number of dof in a basis
*/
template <class... Basis>
struct __count_basis_dof;

template <>
struct __count_basis_dof<> {
  static const A2D::index_t ndof = 0;
};

template <class First, class... Remain>
struct __count_basis_dof<First, Remain...> {
  static const A2D::index_t ndof =
      First::ndof + __count_basis_dof<Remain...>::ndof;
};

/*
  Template pattern for computing the size of the interpolation matrices
*/
template <class... Basis>
struct __count_basis_size;

template <>
struct __count_basis_size<> {
  static const A2D::index_t basis_size = 0;
};

template <class First, class... Remain>
struct __count_basis_size<First, Remain...> {
  static const A2D::index_t basis_size =
      First::basis_size + __count_basis_size<Remain...>::basis_size;
};

/*
  The finite element basis class.

  This class stores a collection of basis function objects
*/
template <typename T, class... Basis>
class FEBasis {
 public:
  typedef std::tuple<Basis...> BasisSpace;

  /**
   * @brief Number of basis function objects
   */
  static constexpr A2D::index_t nbasis =
      std::tuple_size<std::tuple<Basis...>>();

  /**
   * @brief Count up the total number of degrees of freedom for this set of
   * basis functions
   */
  static constexpr A2D::index_t ndof = __count_basis_dof<Basis...>::ndof;

  /**
   * @brief Count up the total basis size required to store all of the
   * interpolation matrices for this set of basis functions
   */
  static constexpr A2D::index_t basis_size =
      __count_basis_size<Basis...>::basis_size;

  /**
   * @brief Get the number of degrees of freedom associated with the given basis
   */
  template <A2D::index_t index>
  static A2D::index_t get_ndof() {
    return std::tuple_element<index, BasisSpace>::type::ndof;
  }

  /**
   * @brief Get the cumulative number of degrees of associated with the basis
   * before the given basis
   */
  template <A2D::index_t index>
  static A2D::index_t get_dof_offset() {
    return get_dof_offset_<0, index, Basis...>();
  }

  /**
   * @brief Interpolate using the basis functions at the quadrature point
   *
   * @param pt The interpolation point index
   * @param dof The degree of freedom object
   * @param s The output finite element space object
   */
  template <class Quadrature, class FEDof, class FiniteElementSpace>
  static void interp(A2D::index_t pt, const FEDof& dof, FiniteElementSpace& s) {
    interp_<Quadrature, FEDof, FiniteElementSpace, 0, Basis...>(pt, dof, s);
  }

  /**
   * @brief Add values to the degree of freedom using the interpolation
   * functions
   *
   * @param pt The interpolation point index
   * @param s The finite element space objec
   * @param dof The degree of freedom object that values are added to
   */
  template <class Quadrature, class FiniteElementSpace, class FEDof>
  static void add(A2D::index_t pt, const FiniteElementSpace& s, FEDof& dof) {
    add_<Quadrature, FiniteElementSpace, FEDof, 0, Basis...>(pt, s, dof);
  }

 private:
  template <class Quadrature, class FEDof, class FiniteElementSpace,
            A2D::index_t index, class First, class... Remain>
  static void interp_(A2D::index_t pt, const FEDof& dof,
                      FiniteElementSpace& s) {
    // Interpolate
    First::template interp<Quadrature>(pt, dof.template get<index>(),
                                       s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, FEDof, FiniteElementSpace, index + 1, Remain...>(
        pt, dof, s);
  }

  template <class Quadrature, class FEDof, class FiniteElementSpace,
            A2D::index_t index>
  static void interp_(A2D::index_t pt, const FEDof& dof,
                      FiniteElementSpace& s) {}

  template <class Quadrature, class FiniteElementSpace, class FEDof,
            A2D::index_t index, class First, class... Remain>
  static void add_(A2D::index_t pt, const FiniteElementSpace& s, FEDof& dof) {
    // Add the interpolation
    First::template add<Quadrature>(pt, s.template get<index>(),
                                    dof.template get<index>());

    // Do the next solution space, if any...
    add_<Quadrature, FiniteElementSpace, FEDof, index + 1, Remain...>(pt, s,
                                                                      dof);
  }

  template <class Quadrature, class FiniteElementSpace, class FEDof,
            A2D::index_t index>
  static void add_(A2D::index_t pt, const FiniteElementSpace& s, FEDof& dof) {}

  // Get the offset recursively
  template <A2D::index_t r, A2D::index_t index, class First, class... Remain>
  static A2D::index_t get_dof_offset_() {
    if (r == index) {
      return 0;
    }
    return First::ndof + get_dof_offset_<r + 1, index, Remain...>();
  }

  template <A2D::index_t r, A2D::index_t index>
  static A2D::index_t get_dof_offset_() {
    return 0;
  }
};

}  // namespace A2D

#endif  // A2D_FE_BASIS_H