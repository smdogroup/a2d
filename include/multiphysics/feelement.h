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
class FiniteElement {  // : public ElementBase<T> {
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

  FiniteElement() {}

  void add_residual(DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemVec& elem_res) {
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
  }

  /*
    Compute y = y + J * x
  */
  void add_jacobian_vector_product(DataElemVec& elem_data, GeoElemVec& elem_geo,
                                   ElemVec& elem_sol, ElemVec& elem_xvec,
                                   ElemVec& elem_yvec) {
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

      // Set up the values for the input vector
      typename ElemVec::FEDof x_dof(i, elem_xvec);
      elem_xvec.get_element_values(i, x_dof);

      // Set up values for the output vector
      typename ElemVec::FEDof y_dof(i, elem_yvec);

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

        // Initialize the Jacobian-vector product functor
        double weight = Quadrature::get_weight(j);
        typename PDE::JacVecProduct jvp(weight * detJ, qdata, s);

        // Interpolate the x-dof
        typename PDE::FiniteElementSpace xref, x;
        Basis::template interp<Quadrature>(j, x_dof, xref);

        xref.transform(detJ, J, Jinv, x);

        // Compute the Jacobian-vector product
        typename PDE::FiniteElementSpace y, yref;
        jvp(x, y);

        // Transform to back to the reference element
        y.rtransform(detJ, J, Jinv, yref);

        // Add the contributions back to the residual
        Basis::template add<Quadrature>(j, yref, y_dof);
      }

      elem_yvec.add_element_values(i, y_dof);
    }
  }

  template <class ElemMat>
  void add_jacobian(DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemMat& elem_mat) {
    const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;
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

      // Set up values for the element matrix
      typename ElemMat::FEMat element_mat(i, elem_mat);
      // A2D::Mat<T, Basis::ndof, Basis::ndof> element_mat;

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

        // Initialize the Jacobian-vector product functor
        double weight = Quadrature::get_weight(j);
        typename PDE::JacVecProduct jvp(weight * detJ, qdata, s);

        // Fill in the entries of the Jacobian matrix
        typename PDE::QMatType jac;

        // Temporary vectors
        typename PDE::FiniteElementSpace pref, p, Jp;

        for (A2D::index_t k = 0; k < ncomp; k++) {
          // Set the value into the matrix
          pref.zero();
          pref[k] = T(1.0);
          pref.transform(detJ, J, Jinv, p);

          // Compute the Jacobian-vector product
          jvp(p, Jp);

          // Transform to back to the reference element
          Jp.rtransform(detJ, J, Jinv, pref);

          for (A2D::index_t m = 0; m < ncomp; m++) {
            jac(m, k) = pref[m];
          }
        }

        // Add the results of the matrix-vector product to the
        Basis::template add_outer<Quadrature>(j, jac, element_mat);
      }

      elem_mat.add_element_values(i, element_mat);
    }
  }
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H