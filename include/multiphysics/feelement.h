#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include <type_traits>

#include "multiphysics/febase.h"
#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_matrix.h"

namespace A2D {

template <typename T, class PDE, class Quadrature, class DataBasis,
          class GeoBasis, class Basis, bool use_parallel_elemvec = false>
class FiniteElement : public ElementBase<T> {
 public:
  using DataElemVec =
      typename std::conditional<use_parallel_elemvec,
                                A2D::ElementVector_Parallel<T, DataBasis>,
                                A2D::ElementVector_Serial<T, DataBasis>>::type;
  using GeoElemVec =
      typename std::conditional<use_parallel_elemvec,
                                A2D::ElementVector_Parallel<T, GeoBasis>,
                                A2D::ElementVector_Serial<T, GeoBasis>>::type;
  using ElemVec =
      typename std::conditional<use_parallel_elemvec,
                                A2D::ElementVector_Parallel<T, Basis>,
                                A2D::ElementVector_Serial<T, Basis>>::type;

  FiniteElement(DataElemVec& elem_data, GeoElemVec& elem_geo, ElemVec& elem_sol,
                ElemVec& elem_res)
      : elem_data(elem_data),
        elem_geo(elem_geo),
        elem_sol(elem_sol),
        elem_res(elem_res) {}

  void set_geo(SolutionVector<T>& X) {}

  void set_solution(SolutionVector<T>& U) {}

  // TODO: global_res unused!
  void add_residual(SolutionVector<T>& global_res) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);

      // Get the geometry values
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);

      // Get the degrees of freedom for the element
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);

      // Get residual for the element
      typename ElemVec::FEDof res_dof(i, elem_res);
      elem_res.get_element_values(i, res_dof);

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

        // Interpolate the data to the quadrature point
        typename PDE::DataSpace qdata;
        DataBasis::template interp<Quadrature>(j, data_dof, qdata);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        Basis::template interp<Quadrature>(j, sol_dof, sref);

        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace s;
        sref.transform(detJ, J, Jinv, s);

        // Compute the coefficients for the weak form of the PDE
        typename PDE::FiniteElementSpace coef;
        double weight = Quadrature::get_weight(j);

        PDE::weak_coef(weight * detJ, qdata, s, coef);

        // Compute the coefficients in the reference element
        typename PDE::FiniteElementSpace cref;
        coef.rtransform(detJ, J, Jinv, cref);

        // Add the contributions back to the residual
        Basis::template add<Quadrature>(j, cref, res_dof);
      }

      elem_res.add_element_values(i, res_dof);
    }

    // We may call this to add to global dof at once in parallel
    elem_res.add_values();

    return;
  }

  /*
    Compute y = y + J * x
  */
  // void add_jacobian_vector_product(SolutionVector<T>& xvec,
  //                                  SolutionVector<T>& yvec) {
  //   ElemVec x(mesh, xvec);
  //   ElemVec y(mesh, yvec);

  //   const A2D::index_t num_elements = mesh.get_num_elements();
  //   const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

  //   for (A2D::index_t i = 0; i < num_elements; i++) {
  //     // Get the geometry values
  //     typename GeoElemVec::FEDof geo_dof;
  //     elem_geo.get_element_values(i, geo_dof);

  //     // Get the degrees of freedom for the element
  //     typename ElemVec::FEDof dof, xdof;
  //     elem_res.get_element_values(i, dof);
  //     x.get_element_values(i, xdof);

  //     // Set up values for the residual
  //     typename ElemVec::FEDof ydof;

  //     for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
  //       // Extract the Jacobian of the element transformation
  //       typename PDE::FiniteElementGeometry gk;
  //       GeoBasis::template interp<Quadrature>(j, geo_dof, gk);
  //       A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template
  //       get<0>()).get_grad();

  //       // Compute the inverse of the transformation
  //       A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
  //       A2D::MatInverse(J, Jinv);

  //       // Compute the determinant of the Jacobian matrix
  //       T detJ;
  //       A2D::MatDet(J, detJ);

  //       // Interpolate the solution vector using the basis
  //       typename PDE::FiniteElementSpace sref;
  //       Basis::template interp<Quadrature>(j, dof, sref);

  //       // Interpolate the solution vector using the basis
  //       typename PDE::FiniteElementSpace xref;
  //       Basis::template interp<Quadrature>(i, xdof, xref);

  //       // Transform to the local coordinate system
  //       typename PDE::FiniteElementSpace s, p;
  //       sref.transform(detJ, J, Jinv, s);
  //       xref.transform(detJ, J, Jinv, p);

  //       // Compute the coefficients for the weak form of the PDE
  //       typename PDE::FiniteElementSpace coef;
  //       double weight = Quadrature::get_weight(j);
  //       PDE::eval_weak_jacobian_vec_product(weight * detJ, s, p, coef);

  //       // Compute the coefficients in the reference element
  //       typename PDE::FiniteElementSpace cref;
  //       coef.rtransform(detJ, J, Jinv, cref);

  //       // Add the contributions back to the residual
  //       Basis::template add<Quadrature>(j, cref, ydof);
  //     }

  //     y.add_element_values(i, ydof);
  //   }
  // }

  // void add_dof_set(Kokkos::UnorderedMap<COO<I>, void>& node_set) = 0;

  // void add_jacobian() {
  //   const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;
  //   const A2D::index_t num_elements = mesh.get_num_elements();
  //   const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

  //   for (A2D::index_t i = 0; i < num_elements; i++) {
  //     // Get the data for the element
  //     typename DataElemVec::FEDof data_dof(i, data);
  //     data.get_element_values(i, data_dof);

  //     // Get the geometry values
  //     typename GeoElemVec::FEDof geo_dof(i, geo);
  //     geo.get_element_values(i, geo_dof);

  //     // Get the degrees of freedom for the element
  //     typename ElemVec::FEDof dof(i, sol);
  //     sol.get_element_values(i, dof);

  //     // Set up values for the element matrix
  //     typename ElemMat::FEDof element_res(i, res);

  //     for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
  //       // Extract the Jacobian of the element transformation
  //       typename PDE::FiniteElementGeometry gk;
  //       GeoBasis::template interp<Quadrature>(j, geo_dof, gk);
  //       A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template
  //       get<0>()).get_grad();

  //       // Compute the inverse of the transformation
  //       A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
  //       A2D::MatInverse(J, Jinv);

  //       // Compute the determinant of the Jacobian matrix
  //       T detJ;
  //       A2D::MatDet(J, detJ);

  //       // Interpolate the data to the quadrature point
  //       typename PDE::DataSpace qdata;
  //       DataBasis::template interp<Quadrature>(j, data_dof, qdata);

  //       // Interpolate the solution vector using the basis
  //       typename PDE::FiniteElementSpace sref;
  //       Basis::template interp<Quadrature>(j, dof, sref);

  //       // Transform to the local coordinate system
  //       typename PDE::FiniteElementSpace s;
  //       sref.transform(detJ, J, Jinv, s);

  //       // Initialize the Jacobian-vector product functor
  //       double weight = Quadrature::get_weight(j);
  //       typename PDE::JacVecProduct jvp(weight * detJ, qdata, s);

  //       // Fill in the entries of the Jacobian matrix
  //       A2D::Mat<T, ncomp, ncomp> jac;

  //       for (A2D::index_t k = 0; k < ncomp; k++) {
  //         // Set the value into the matrix
  //         typename PDE::FiniteElementSpace pref, Jp;
  //         pref.set_value(k, T(1.0));
  //         pref.transform(detJ, J, Jinv, p);

  //         // Compute the Jacobian-vector product
  //         jvp(p, Jp);

  //         // Transform to back to the reference element
  //         Jp.rtransform(detJ, J, Jinv, pref);

  //         for (A2D::index_t m = 0; m < ncomp; m++) {
  //           jac(m, k) = Jp.get_value(m);
  //         }
  //       }

  //       // Add the results of the matrix-vector product to the
  //       Basis::template add_outer<Quadrature>(j, jac, element_mat);
  //     }
  //   }
  // }

 private:
  DataElemVec& elem_data;
  GeoElemVec& elem_geo;
  ElemVec& elem_sol;
  ElemVec& elem_res;
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H