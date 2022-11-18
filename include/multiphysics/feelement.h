#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>

#include "multiphysics/febase.h"
#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_matrix.h"

namespace A2D {

template <typename T, class PDE>
void TestPDEImplementation(double dh = 1e-7) {
  typename PDE::DataSpace data;
  typename PDE::FiniteElementGeometry geo;
  typename PDE::FiniteElementSpace s, sref;
  typename PDE::FiniteElementSpace p, pref;
  typename PDE::FiniteElementSpace coef, cref, cref0;
  typename PDE::FiniteElementSpace Jp, Jpref;

  // Generate random data
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> distr(-1.0, 1.0);

  // Set random values for the data
  if constexpr (PDE::DataSpace::ncomp > 0) {
    for (index_t i = 0; i < PDE::DataSpace::ncomp; i++) {
      data[i] = distr(gen);
    }
  }

  // Set random values for the geometry
  for (index_t i = 0; i < PDE::FiniteElementGeometry::ncomp; i++) {
    geo[i] = distr(gen);
  }

  // Set the random values
  for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
    sref[i] = distr(gen);
    pref[i] = distr(gen);
  }

  A2D::Mat<T, PDE::dim, PDE::dim>& J = (geo.template get<0>()).get_grad();
  A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
  A2D::MatInverse(J, Jinv);

  T detJ;
  A2D::MatDet(J, detJ);

  // Compute the coefficients
  sref.transform(detJ, J, Jinv, s);
  PDE::weak_coef(detJ, data, geo, s, coef);
  coef.rtransform(detJ, J, Jinv, cref0);

  if constexpr (std::is_same<T, std::complex<double>>::value) {
    for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
      sref[i] = sref[i] + dh * pref[i] * std::complex<double>(0.0, 1.0);
    }
  } else {
    for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
      sref[i] = sref[i] + dh * pref[i];
    }
  }

  // Compute the coefficients
  sref.transform(detJ, J, Jinv, s);
  PDE::weak_coef(detJ, data, geo, s, coef);
  coef.rtransform(detJ, J, Jinv, cref);

  // Compute the Jacobian-vector product
  typename PDE::JacVecProduct jvp(detJ, data, geo, s);

  // Compute the Jacobian-vector product
  pref.transform(detJ, J, Jinv, p);
  jvp(p, Jp);
  Jp.rtransform(detJ, J, Jinv, Jpref);

  // Compute the finite-difference value
  typename PDE::FiniteElementSpace fd;

  if constexpr (std::is_same<T, std::complex<double>>::value) {
    for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
      fd[i] = std::imag(cref[i]) / dh;
    }
  } else {
    for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
      fd[i] = (cref[i] - cref0[i]) / dh;
    }
  }

  for (index_t i = 0; i < PDE::FiniteElementSpace::ncomp; i++) {
    std::cout << "fd[" << std::setw(2) << i << "]: " << std::setw(12)
              << std::real(fd[i]) << " Jpref[" << std::setw(2) << i
              << "]: " << std::setw(12) << std::real(Jpref[i]) << " err["
              << std::setw(2) << i << "]: " << std::setw(12)
              << std::real((fd[i] - Jpref[i]) / fd[i]) << std::endl;
  }
}

template <typename T, class PDE, class Quadrature, class DataBasis,
          class GeoBasis, class Basis, bool use_parallel_elemvec = false>
class FiniteElement {
 public:
  /*
    This wraps an FEBasis class to make it look like a quadrature class
  */
  template <class MyBasis>
  class DofPoints {
   public:
    static index_t get_num_points() { return MyBasis::ndof; }
    static void get_point(const index_t n, double pt[]) {
      MyBasis::get_dof_point(n, pt);
    }
  };

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

  void get_near_nullspace(const index_t null_index, DataElemVec& elem_data,
                          GeoElemVec& elem_geo, ElemVec& elem_null) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);

      // Get the geometry values
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);

      for (A2D::index_t dof = 0; dof < Basis::ndof; dof++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry geo;
        GeoBasis::template interp<DofPoints<Basis>>(dof, geo_dof, geo);

        // Interpolate the data to the quadrature point
        typename PDE::DataSpace qdata;
        DataBasis::template interp<DofPoints<Basis>>(dof, data_dof, qdata);

        // // Set the solution space
        // PDE::near_nullspace(null_index, qdata)
      }
    }
  }

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
        typename PDE::FiniteElementGeometry geo;
        GeoBasis::template interp<Quadrature>(j, geo_dof, geo);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (geo.template get<0>()).get_grad();

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

        PDE::weak_coef(weight * detJ, qdata, geo, s, coef);

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
        typename PDE::FiniteElementGeometry geo;
        GeoBasis::template interp<Quadrature>(j, geo_dof, geo);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (geo.template get<0>()).get_grad();

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
        typename PDE::JacVecProduct jvp(weight * detJ, qdata, geo, s);

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
        typename PDE::FiniteElementGeometry geo;
        GeoBasis::template interp<Quadrature>(j, geo_dof, geo);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (geo.template get<0>()).get_grad();

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
        typename PDE::JacVecProduct jvp(weight * detJ, qdata, geo, s);

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