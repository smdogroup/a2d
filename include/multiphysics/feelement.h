#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_matrix.h"

namespace A2D {

/*
  In-place element vector implementation
*/
template <typename T, class FiniteElementSpace, class Basis>
class ElementVector_InPlace {
 public:
  ElementVector_InPlace(A2D::ElementMesh<Basis>& mesh,
                        A2D::SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {}

  // Required DOF container object (different for each implementation)
  class FEDof {
   public:
    FEDof() {}

    // Variables for all the basis functions
    T dof[Basis::ndof];

    // Zero the values
    void zero() { std::fill(dof, dof + Basis::ndof, T(0.0)); }

    // Get the dof values
    template <A2D::index_t index>
    T* get() {
      return &dof[Basis::template get_dof_offset<index>()];
    }
    template <A2D::index_t index>
    const T* get() const {
      return &dof[Basis::template get_dof_offset<index>()];
    }
  };

  // Get values for this element from the vector
  void get_element_values(A2D::index_t elem, FEDof& dof) {
    get_element_values_<Basis::nbasis - 1>(elem, dof);
  }

  // Add values for this element to the vector
  void add_element_values(A2D::index_t elem, const FEDof& dof) {
    add_element_values_<Basis::nbasis - 1>(elem, dof);
  }

  A2D::SolutionVector<T>& get_vector() { return vec; }
  void set_solution() {}
  void get_residual() {}

 private:
  template <A2D::index_t basis>
  void get_element_values_(A2D::index_t elem, FEDof& dof) {
    T* values = dof.template get<basis>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      values[i] = sign * vec[dof_index];
    }
    get_element_values_<basis - 1>(elem, dof);
  }

  template <>
  void get_element_values_<0>(A2D::index_t elem, FEDof& dof) {
    T* values = dof.template get<0>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<0>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, 0, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, 0, i);
      values[i] = sign * vec[dof_index];
    }
  }

  template <A2D::index_t basis>
  void add_element_values_(A2D::index_t elem, const FEDof& dof) {
    const T* values = dof.template get<basis>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<basis>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, basis, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, basis, i);
      vec[dof_index] += sign * values[i];
    }
    add_element_values_<basis - 1>(elem, dof);
  }

  template <>
  void add_element_values_<0>(A2D::index_t elem, const FEDof& dof) {
    const T* values = dof.template get<0>();
    for (A2D::index_t i = 0; i < Basis::template get_ndof<0>(); i++) {
      const int sign = mesh.get_global_dof_sign(elem, 0, i);
      const A2D::index_t dof_index = mesh.get_global_dof(elem, 0, i);
      vec[dof_index] += sign * values[i];
    }
  }

  A2D::ElementMesh<Basis>& mesh;
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

  template <class Quadrature, A2D::index_t index, class First, class...
Remain> void interp_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace&
s) { const values = subview(element_dof, elem, Kokkos::ALL);

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, A2D::index_t index>
  void add_(A2D::index_t elem, A2D::index_t pt, const FiniteElementSpace& s)
{}

  template <class Quadrature, A2D::index_t index, class First, class...
Remain> void add_(A2D::index_t elem, A2D::index_t pt, const
FiniteElementSpace& s) { values = subview(element_dof, elem, Kokkos::ALL);

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
  template <class Quadrature, A2D::index_t index, class First, class...
Remain> void outer_(A2D::index_t elem, A2D::index_t pt, FiniteElementSpace& s)
{
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

template <typename T, class Quadrature, class PDE, class GeoBasis, class Basis>
class FiniteElement {
 public:
  typedef A2D::ElementVector_InPlace<T, typename PDE::FiniteElementGeometry,
                                     GeoBasis>
      GeoElemVec;
  typedef A2D::ElementVector_InPlace<T, typename PDE::FiniteElementSpace, Basis>
      ElemVec;

  FiniteElement(A2D::ElementMesh<GeoBasis>& geomesh,
                A2D::SolutionVector<T>& nodes, A2D::ElementMesh<Basis>& mesh,
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
      // Get the geometry values
      typename GeoElemVec::FEDof geo_dof;
      geo.get_element_values(i, geo_dof);

      // Get the degrees of freedom for the element
      typename ElemVec::FEDof dof;
      sol.get_element_values(i, dof);

      // Set up values for the residual
      typename ElemVec::FEDof element_res;
      element_res.zero();

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        GeoBasis::template interp<Quadrature>(j, geo_dof, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        Basis::template interp<Quadrature>(j, dof, sref);

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
        Basis::template add<Quadrature>(j, cref, element_res);
      }

      res.add_element_values(i, element_res);
    }
  }

  /*
    Compute y = y + J * x
  */
  void add_jacobian_vector_product(SolutionVector<T>& xvec,
                                   SolutionVector<T>& yvec) {
    ElemVec x(mesh, xvec);
    ElemVec y(mesh, yvec);

    const A2D::index_t num_elements = mesh.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the geometry values
      typename GeoElemVec::FEDof geo_dof;
      geo.get_element_values(i, geo_dof);

      // Get the degrees of freedom for the element
      typename ElemVec::FEDof dof, xdof;
      sol.get_element_values(i, dof);
      x.get_element_values(i, xdof);

      // Set up values for the residual
      typename ElemVec::FEDof ydof;
      ydof.zero();

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        GeoBasis::template interp<Quadrature>(j, geo_dof, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        Basis::template interp<Quadrature>(j, dof, sref);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace xref;
        Basis::template interp<Quadrature>(i, xdof, xref);

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
        Basis::template add<Quadrature>(j, cref, ydof);
      }

      y.add_element_values(i, ydof);
    }
  }

  void add_jacobian() {
    const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;
    const A2D::index_t num_elements = mesh.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the geometry values
      typename GeoElemVec::FEDof geo_dof;
      geo.get_element_values(i, geo_dof);

      // Get the degrees of freedom for the element
      typename ElemVec::FEDof dof;
      sol.get_element_values(i, dof);

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        GeoBasis::template interp<Quadrature>(j, geo_dof, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        Basis::template interp<Quadrature>(j, dof, sref);

        typename PDE::FiniteElementSpace cref[ncomp];
        for (A2D::index_t k = 0; k < ncomp; k++) {
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
  A2D::ElementMesh<Basis> mesh;
  A2D::ElementMesh<GeoBasis>& geomesh;

  // Element-wise views of the solution and residual vector
  GeoElemVec geo;
  ElemVec sol, res;
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H