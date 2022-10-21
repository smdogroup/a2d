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
class LagrangeTri0 {
 public:
  static const A2D::index_t ndof = 1;
  static const A2D::index_t ncomp = A2D::L1ScalarSpace<T, 2>::ncomp;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType& sol,
                     A2D::L1ScalarSpace<T, 2>& out) {
    T& u = out.get_value();
    u = sol[0];
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::L1ScalarSpace<T, 2>& in,
                  SolnType& res) {
    const T& u = in.get_value();
    res[0] += u;
  }

  // Set the matrix stride
  static const A2D::index_t stride = 1;

  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType& N) {
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
  static void interp(A2D::index_t n, const SolnType& sol,
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
                  SolnType& res) {
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
  static const A2D::index_t stride = 1;

  template <class Quadrature, class BasisType>
  static void basis(A2D::index_t n, BasisType& N) {}
};

/*
  Raviart-Thomas element for H(div) in 2D
*/
template <typename T>
class RT2DTri1 {
 public:
  static const A2D::index_t ndof = 3;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType& sol,
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
  static void add(A2D::index_t n, const A2D::Hdiv2DSpace<T>& in,
                  SolnType& res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    const A2D::Vec<T, 2>& u = in.get_value();
    const T& div = in.get_div();

    res[0] += pt[0] * u(0) + pt[1] * u(1) + 2.0 * div;
    res[1] += (pt[0] - 1.0) * u(0) + pt[1] * u(1) + 2.0 * div;
    res[2] += pt[0] * u(0) + (pt[1] - 1.0) * u(1) + 2.0 * div;
  }
};

}  // namespace A2D

#endif  // A2D_FE_BASIS_H