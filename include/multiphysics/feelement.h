#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_matrix.h"

namespace A2D {

/*
  In-place element vector implementation
*/
template <typename T, class FiniteElementSpace, class... Basis>
class ElementVector_InPlace {
 public:
  ElementVector_InPlace(A2D::ElementMesh<Basis...>& mesh,
                        A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {}

  A2D::SolutionVector<T>& get_vector() { return vec; }
  void set_solution() {}
  void get_residual() {}

  template <class Quadrature>
  void interp(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    interp_<Quadrature, 0, Basis...>(elem, pt, s);
  }

  template <class Quadrature>
  void add(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s) {
    add_<Quadrature, 0, Basis...>(elem, pt, s);
  }

 private:
  template <class Quadrature, A2D::index_t index, class First, class... Remain>
  void interp_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    for (A2D::index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const A2D::index_t dof = mesh.get_global_dof(elem, index, i);
      values[i] = sign * vec[dof];
    }

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, A2D::index_t index>
  void interp_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {}

  template <class Quadrature, A2D::index_t index, class First, class... Remain>
  void add_(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    std::fill(values, values + First::ndof, T(0.0));

    // Add the interpolation
    First::template add<Quadrature>(pt, s.template get<index>(), values);

    for (A2D::index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const A2D::index_t dof = mesh.get_global_dof(elem, index, i);
      vec[dof] += sign * values[i];
    }

    // Do the next solution space, if any...
    add_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, A2D::index_t index>
  void add_(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s) {}

  A2D::ElementMesh<Basis...>& mesh;
  A2D::SolutionVector<T>& vec;
};

/*
template <typename T, class FiniteElementSpace, class... Basis>
class ElementVector_Kokkos {
 public:
  static const A2D::index_t ndof_per_element = dof_in_basis<Basis...>;

  ElementVector_Kokkos(ElementMesh<Basis...>& mesh, SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {
    A2D::index_t nelems = mesh.get_num_elements();
  }

  SolutionVector<T>& get_vector() { return vec; }

  void set_solution() {
    // Go set values into element_dof and transfer them to the device
    for (A2D::index_t elem = 0; elem < num_elements; elem++) {
      for (A2D::index_t i = 0; i < First::ndof; i++) {
        const int sign = mesh.get_global_dof_sign(elem, index, i);
        const A2D::index_t dof = mesh.get_global_dof(elem, index, i);
        element_dof[elem, i] = sign * vec[dof];
      }
    }
  }

  void get_residual() {
    for (A2D::index_t elem = 0; elem < num_elements; elem++) {
      for (A2D::index_t i = 0; i < First::ndof; i++) {
        const int sign = mesh.get_global_dof_sign(elem, index, i);
        const A2D::index_t dof = mesh.get_global_dof(elem, index, i);
        vec[dof] += sign * element_dof[elem, i];
      }
    }
  }

  template <class Quadrature>
  void interp(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    interp_<Quadrature, 0, Basis...>(elem, pt, s);
  }

 private:
  template <class Quadrature, A2D::index_t index>
  void interp_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {}

  template <class Quadrature, A2D::index_t index, class First, class... Remain>
  void interp_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    const values = subview(element_dof, elem, Kokkos::ALL);

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, A2D::index_t index>
  void add_(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s) {}

  template <class Quadrature, A2D::index_t index, class First, class... Remain>
  void add_(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s) {
    values = subview(element_dof, elem, Kokkos::ALL);

    // Add the interpolation
    First::template add<Quadrature>(pt, s.template get<index>(), values);

    // Do the next solution space, if any...
    add_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  Kokkos::View<T* [ndof_per_element]> element_dof;
};
*/

/*
template <typename T, class FiniteElementSpace, class... Basis>
class ElementMatrix_A2D {
 public:
  ElementMatrix_A2D(A2D::ElementMesh<Basis...>& mesh,
                    A2D::SolutionMatrix<T>& vec)
      : mesh(mesh), vec(vec) {}

  template <class Quadrature>
  void outer(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    outer_<Quadrature, 0, Basis...>(elem, pt, s);
  }

 private:
  template <class Quadrature, A2D::index_t index, class First, class... Remain>
  void outer_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    for (A2D::index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const A2D::index_t dof = mesh.get_global_dof(elem, index, i);
      values[i] = sign * vec[dof];
    }

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    outer_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, A2D::index_t index>
  void outer_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s) {}

  A2D::ElementMesh<Basis...>& mesh;
  A2D::SolutionMatrix<T>& mat;
};
*/

template <typename T, class Quadrature, class PDE, class GeoBasis,
          class... Basis>
class FiniteElement {
 public:
  FiniteElement(A2D::ElementMesh<GeoBasis>& geomesh,
                A2D::SolutionVector<T>& nodes, A2D::ElementMesh<Basis...>& mesh,
                A2D::SolutionVector<T>& solvec, A2D::SolutionVector<T>& resvec)
      : mesh(mesh),
        geomesh(geomesh),
        geo(geomesh, nodes),
        sol(mesh, solvec),
        res(mesh, resvec) {}

  /*
    Copy values to the solution vector
  */
  void set_solution(SolutionVector<T>& sol) {}

  void add_residual() {
    const A2D::index_t num_elements = mesh.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace s;
        sref.transform(detJ, J, Jinv, s);

        // Compute the coefficients for the weak form of the PDE
        typename PDE::FiniteElementSpace coef;
        double weight = Quadrature::get_weight(j);
        PDE::eval_weak_coef(weight * detJ, s, coef);

        // Compute the coefficients in the reference element
        typename PDE::FiniteElementSpace cref;
        coef.rtransform(detJ, J, Jinv, cref);

        // Add the contributions back to the residual
        res.template add<Quadrature>(i, j, cref);
      }
    }
  }

  /*
    Compute y = y + J * x
  */
  void add_jacobian_vector_product(SolutionVector<T>& xvec,
                                   SolutionVector<T>& yvec) {
    A2D::ElementVector_InPlace<T, typename PDE::FiniteElementSpace, Basis...> x(
        mesh, xvec);
    A2D::ElementVector_InPlace<T, typename PDE::FiniteElementSpace, Basis...> y(
        mesh, yvec);

    const A2D::index_t num_elements = mesh.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace xref;
        x.template interp<Quadrature>(i, j, xref);

        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace s, p;
        sref.transform(detJ, J, Jinv, s);
        xref.transform(detJ, J, Jinv, p);

        // Compute the coefficients for the weak form of the PDE
        typename PDE::FiniteElementSpace coef;
        double weight = Quadrature::get_weight(j);
        PDE::eval_weak_jacobian_vec_product(weight * detJ, s, p, coef);

        // Compute the coefficients in the reference element
        typename PDE::FiniteElementSpace cref;
        coef.rtransform(detJ, J, Jinv, cref);

        // Add the contributions back to the residual
        y.template add<Quadrature>(i, j, cref);
      }
    }
  }

  void add_jacobian() {
    const A2D::index_t ndof = PDE::FiniteElementSpace::ndof;
    const A2D::index_t num_elements = mesh.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        typename PDE::FiniteElementSpace cref[ndof];
        for (A2D::index_t k = 0; k < ndof; k++) {
          typename PDE::FiniteElementSpace pref;
          pref.set_seed(k);

          // Transform to the local coordinate system
          typename PDE::FiniteElementSpace s, p;
          sref.transform(detJ, J, Jinv, s);
          pref.transform(detJ, J, Jinv, p);

          // Compute the coefficients for the weak form of the PDE
          typename PDE::FiniteElementSpace coef;
          double weight = Quadrature::get_weight(j);
          PDE::eval_weak_jacobian_vec_product(weight * detJ, s, p, coef);

          // Compute the coefficients in the reference element
          coef.rtransform(detJ, J, Jinv, cref[k]);
        }
      }
    }
  }

 private:
  // The element mesh
  A2D::ElementMesh<Basis...>& mesh;
  A2D::ElementMesh<GeoBasis>& geomesh;

  // Element-wise views of the solution and residual vector
  A2D::ElementVector_InPlace<T, typename PDE::FiniteElementGeometry, GeoBasis>
      geo;
  A2D::ElementVector_InPlace<T, typename PDE::FiniteElementSpace, Basis...> sol;
  A2D::ElementVector_InPlace<T, typename PDE::FiniteElementSpace, Basis...> res;
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H