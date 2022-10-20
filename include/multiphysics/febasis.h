#ifndef A2D_FE_BASIS_H
#define A2D_FE_BASIS_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"
#include "multiphysics/fespace.h"

namespace A2D {

/*
  Lagrange basis for a triangle
*/
template <typename T>
class LagrangeTri0 {
 public:
  static const int ndof = 1;

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
};

/*
  Lagrange basis for a triangle
*/
template <typename T, A2D::index_t C>
class LagrangeTri1 {
 public:
  static const A2D::index_t ndof = 3 * C;

  template <class Quadrature, class SolnType>
  static void interp(A2D::index_t n, const SolnType& sol,
                     A2D::H1ScalarSpace<T, 2>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    T& u = out.get_value();
    A2D::Vec<T, 2>& grad = out.get_grad();

    u = N[0] * sol[0] + N[1] * sol[1] + N[2] * sol[2];
    grad(0) = sol[1] - sol[0];
    grad(1) = sol[2] - sol[0];
  }

  template <class Quadrature, class SolnType>
  static void add(A2D::index_t n, const A2D::H1ScalarSpace<T, 2>& in,
                  SolnType& res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    const T& u = in.get_value();
    const A2D::Vec<T, 2>& grad = in.get_grad();

    res[0] += N[0] * u - grad(0) - grad(1);
    res[1] += N[1] * u + grad(0);
    res[2] += N[2] * u + grad(1);
  }

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
};

/*
  Raviart-Thomas element for H(div) in 2D
*/
template <typename T>
class RT2DTri1 {
 public:
  static const int ndof = 3;

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