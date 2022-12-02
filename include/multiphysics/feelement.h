#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>

#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "sparse/sparse_matrix.h"
#include "utils/a2dprofiler.h"

namespace A2D {

/**
 * @brief Get the coordinates associated with the locations of each of the
 * degrees of freedom
 *
 * This code returns the x, y and z coordinates in separate solution-size
 * vectors
 *
 * @tparam T The solution type
 * @tparam GeoBasis Geometry basis type
 * @tparam Basis Solution basis type
 */
template <typename T, class PDE, class GeoBasis, class Basis>
class DOFCoordinates {
 public:
  // Set the quadrature point
  class DOFQuadrature {
   public:
    static const bool is_tensor_product = false;
    static const index_t num_quad_points = Basis::ndof;
    static index_t get_num_points() { return num_quad_points; }
    static void get_point(const index_t dof, double pt[]) {
      Basis::get_dof_point(dof, pt);
    }
  };

  // Quadrature point object for the geometry
  using QGeoSpace =
      QptSpace<DOFQuadrature, typename PDE::FiniteElementGeometry>;

  template <class GeoElemVec, class ElemVec>
  void get_dof_coordinates(GeoElemVec& elem_geo, ElemVec& elem_x,
                           ElemVec& elem_y, ElemVec& elem_z) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof x_dof(i, elem_x);
      typename ElemVec::FEDof y_dof(i, elem_y);
      typename ElemVec::FEDof z_dof(i, elem_z);

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < Basis::ndof; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementGeometry& gref = geo.get(j);
        A2D::Vec<T, PDE::dim>& X = gref.template get<0>().get_value();

        x_dof[j] = X(0);
        y_dof[j] = X(1);
        z_dof[j] = X(2);
      }

      elem_x.set_element_values(i, x_dof);
      elem_y.set_element_values(i, y_dof);
      elem_z.set_element_values(i, z_dof);
    }
  }
};

template <typename T, class PDE, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class FiniteElement {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename PDE::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace = QptSpace<Quadrature, typename PDE::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace = QptSpace<Quadrature, typename PDE::FiniteElementSpace>;

  FiniteElement() {}

  /**
   * @brief Compute the value of an integral functional over the finite-elemen
   * domain
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The integral over the element
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec>
  T integrate(PDE& pde, DataElemVec& elem_data, GeoElemVec& elem_geo,
              ElemVec& elem_sol) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    T value = 0.0;

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDE
        double weight = Quadrature::get_weight(j);
        typename PDE::FiniteElementSpace coef;
        value += pde.integrand(weight * detJ, data.get(j), gref, s);
      }
    }

    return value;
  }

  /**
   * @brief Compute the maximum value of a quantity at any quadrature point
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The maximum value over the domain
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec>
  T max(PDE& pde, DataElemVec& elem_data, GeoElemVec& elem_geo,
        ElemVec& elem_sol) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    // Need to make this more robust?
    T max_value = -1e20;

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDE
        double weight = Quadrature::get_weight(j);
        typename PDE::FiniteElementSpace coef;
        T value = pde.max(data.get(j), gref, s);
        if (std::real(value) > std::real(max_value)) {
          max_value = value;
        }
      }
    }

    return max_value;
  }

  /**
   * @brief Compute the derivative of the functional value with respect to the
   * data
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @tparam DataDerivElemVec Element vector class for the derivative
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_deriv Output derivative of the functional w.r.t. data
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec,
            class DataDerivElemVec>
  void add_data_derivative(PDE& pde, DataElemVec& elem_data,
                           GeoElemVec& elem_geo, ElemVec& elem_sol,
                           DataDerivElemVec& elem_deriv) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Allocate space for the derivative
      QDataSpace deriv;

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDE
        double weight = Quadrature::get_weight(j);
        pde.data_derivative(weight * detJ, data.get(j), gref, s, deriv.get(j));
      }

      // Add the derivative of the data back to the data space
      typename DataDerivElemVec::FEDof deriv_dof(i, elem_deriv);
      DataBasis::template add(deriv, deriv_dof);
      elem_deriv.add_element_values(i, deriv_dof);
    }
  }

  /**
   * @brief Add the derivative of the adjoint-vector product with respect to
   * data
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution
   * @tparam ElemAdjVec Element vector class for the adjoint
   * @tparam DataDerivElemVec Element vector class for the derivative
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_adj Element adjoint vector
   * @param elem_deriv Data derivative
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec,
            class ElemAdjVec, class DataDerivElemVec>
  void add_adjoint_residual_data_derivative(PDE& pde, DataElemVec& elem_data,
                                            GeoElemVec& elem_geo,
                                            ElemVec& elem_sol,
                                            ElemAdjVec& elem_adj,
                                            DataDerivElemVec& elem_deriv) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Set up the values for the adjoint vector input
      typename ElemAdjVec::FEDof adj_dof(i, elem_adj);
      elem_adj.get_element_values(i, adj_dof);
      QSpace adj;
      Basis::template interp(adj_dof, adj);

      // Allocate space for the derivative
      QDataSpace deriv;

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementSpace& aref = adj.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s, a;
        transform.transform(sref, s);
        transform.transform(aref, a);

        // Allocate the Jacobian-adjoint product functor
        double weight = Quadrature::get_weight(j);
        typename PDE::AdjVecProduct ajp(pde, weight * detJ, data.get(j), gref,
                                        s);

        ajp(a, deriv.get(j));
      }

      // Add the derivative of the data back to the data space
      typename DataDerivElemVec::FEDof deriv_dof(i, elem_deriv);
      DataBasis::template add(deriv, deriv_dof);
      elem_deriv.add_element_values(i, deriv_dof);
    }
  }

  /**
   * @brief Add the residuals for the finite-element problem
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution
   * @tparam ElemResVec Element vector class for the residual
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_res Element residual vector
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec,
            class ElemResVec>
  void add_residual(PDE& pde, DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemResVec& elem_res) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Allocate space for the residual values at each quadrature point
      QSpace res;

      // Compute the weak coefficients at all quadrature points
      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDE
        double weight = Quadrature::get_weight(j);
        typename PDE::FiniteElementSpace coef;
        pde.weak(weight * detJ, data.get(j), gref, s, coef);

        // Transform the coefficents back to the reference element
        typename PDE::FiniteElementSpace& cref = res.get(j);
        transform.rtransform(coef, cref);
      }

      // Add the residual from the quadrature points back to the finite-element
      // mesh
      typename ElemVec::FEDof res_dof(i, elem_res);
      Basis::template add(res, res_dof);
      elem_res.add_element_values(i, res_dof);
    }
  }

  /**
   * @brief Compute the matrix-free Jacobian-vector product y = y + J * x
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param pde Instance of the PDE
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_xvec Element solution vector for storing x-components
   * @param elem_yvec Output element solution vector storing y-components
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void add_jacobian_vector_product(PDE& pde, DataElemVec& elem_data,
                                   GeoElemVec& elem_geo, ElemVec& elem_sol,
                                   ElemVec& elem_xvec, ElemVec& elem_yvec) {
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Set up the values for the input vector
      typename ElemVec::FEDof x_dof(i, elem_xvec);
      elem_xvec.get_element_values(i, x_dof);
      QSpace xsol;
      Basis::template interp(x_dof, xsol);

      // Allocate space for the output vector
      QSpace ysol;

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementSpace& xref = xsol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace x, s;
        transform.transform(sref, s);
        transform.transform(xref, x);

        // Allocate the Jacobian-vector product functor
        double weight = Quadrature::get_weight(j);
        typename PDE::JacVecProduct jvp(pde, weight * detJ, data.get(j), gref,
                                        s);

        // Compute the Jacobian-vector product
        typename PDE::FiniteElementSpace y;
        jvp(x, y);

        // Transform to back to the reference element
        typename PDE::FiniteElementSpace& yref = ysol.get(j);
        transform.rtransform(y, yref);
      }

      // Add the values from the quadrature points back into the finite-element
      // problem
      typename ElemVec::FEDof y_dof(i, elem_yvec);
      Basis::template add(ysol, y_dof);
      elem_yvec.add_element_values(i, y_dof);
    }
  }

  /**
   * @brief Assemble element Jacobian matrices based on the data, geometry and
   * solution vectors.
   *
   * WARNING: This is intended only for the lowest order elements, ie. p = 1. It
   * scales O(p^9) so it is unsuitable for high-order elements!
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @tparam ElemMat The element matrix
   * @param pde The PDE instance
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_mat Element matrix output
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec, class ElemMat>
  void add_jacobian(PDE& pde, DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemMat& elem_mat) {
    Timer timer("FiniteElement::add_jacobian()");
    const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      // Initialize the element matrix
      typename ElemMat::FEMat element_mat(i, elem_mat);

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDE::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDE::FiniteElementSpace s;
        transform.transform(sref, s);

        // // Allocate the Jacobian-vector product functor
        // double weight = Quadrature::get_weight(j);
        // typename PDE::JacVecProduct jvp(pde, weight * detJ, data.get(j),
        // gref, s);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDE::QMatType jac;

        // Temporary vectors
        typename PDE::FiniteElementSpace pref, p, Jp;

        for (A2D::index_t k = 0; k < ncomp; k++) {
          // Set the value into the matrix
          pref.zero();
          pref[k] = T(1.0);
          transform.transform(pref, p);

          // Allocate the Jacobian-vector product functor
          double weight = Quadrature::get_weight(j);
          typename PDE::JacVecProduct jvp(pde, weight * detJ, data.get(j), gref,
                                          s);

          // Compute the Jacobian-vector product
          jvp(p, Jp);

          // Transform to back to the reference element
          transform.rtransform(Jp, pref);

          for (A2D::index_t m = 0; m < ncomp; m++) {
            jac(m, k) = pref[m];
          }
        }

        // Add the results of the outer product
        Basis::template add_outer<Quadrature>(j, jac, element_mat);
      }

      elem_mat.add_element_values(i, element_mat);
    }
  }

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void add_geo_derivative(PDE& pde, DataElemVec& elem_data,
                          GeoElemVec& elem_geo, ElemVec& elem_sol,
                          GeoElemVec& geo_deriv) {}

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void add_adjoint_residual_geo_derivative(PDE& pde, DataElemVec& elem_data,
                                           GeoElemVec& elem_geo,
                                           ElemVec& elem_sol, ElemVec& elem_adj,
                                           GeoElemVec& geo_deriv) {}
};

template <typename T, class PDE, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class MatrixFree {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename PDE::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace = QptSpace<Quadrature, typename PDE::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace = QptSpace<Quadrature, typename PDE::FiniteElementSpace>;

  // Quadrature point view of the Jacobian-matrices
  using QMatSpace = QptSpace<Quadrature, typename PDE::QMatType>;

  MatrixFree() {}

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void initialize(PDE& pde, DataElemVec& elem_data, GeoElemVec& elem_geo,
                  ElemVec& elem_sol) {
    Timer timer("MatrixFree::initialize()");
    // Re-size the vector as needed
    const A2D::index_t num_elements = elem_geo.get_num_elements();
    if (qmat.size() != num_elements) {
      qmat.resize(num_elements);
    }

    // Number of components at the quadrature point
    const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;

    // Get the number of quadrature points
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Get the data for the element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      elem_data.get_element_values(i, data_dof);
      QDataSpace data;
      DataBasis::template interp(data_dof, data);

      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      elem_geo.get_element_values(i, geo_dof);
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      elem_sol.get_element_values(i, sol_dof);
      QSpace sol;
      Basis::template interp(sol_dof, sol);

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace& sref = sol.get(j);
        typename PDE::FiniteElementGeometry& gref = geo.get(j);

        // Get the Jacobian transformation
        A2D::Mat<T, PDE::dim, PDE::dim>& J = gref.template get<0>().get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Transform the solution the physical element
        typename PDE::FiniteElementSpace x, s;
        sref.transform(detJ, J, Jinv, s);

        // // Allocate the Jacobian-vector product functor
        // double weight = Quadrature::get_weight(j);
        // typename PDE::JacVecProduct jvp(pde, weight * detJ, data.get(j),
        // gref, s);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDE::QMatType& jac = qmat[i].get(j);

        // Temporary vectors
        typename PDE::FiniteElementSpace pref, p, Jp;

        for (A2D::index_t k = 0; k < ncomp; k++) {
          // Set the value into the matrix
          pref.zero();
          pref[k] = T(1.0);
          pref.transform(detJ, J, Jinv, p);

          // Allocate the Jacobian-vector product functor
          double weight = Quadrature::get_weight(j);
          typename PDE::JacVecProduct jvp(pde, weight * detJ, data.get(j), gref,
                                          s);

          // Compute the Jacobian-vector product
          jvp(p, Jp);

          // Transform to back to the reference element
          Jp.rtransform(detJ, J, Jinv, pref);

          for (A2D::index_t m = 0; m < ncomp; m++) {
            jac(m, k) = pref[m];
          }
        }
      }
    }
  }

  template <class ElemVec>
  void add_jacobian_vector_product(ElemVec& elem_xvec, ElemVec& elem_yvec) {
    Timer timer("MatrixFree::add_jacobian_vector_product");
    const A2D::index_t num_elements = qmat.size();
    const A2D::index_t num_quadrature_points = Quadrature::get_num_points();
    const A2D::index_t ncomp = PDE::FiniteElementSpace::ncomp;

    for (A2D::index_t i = 0; i < num_elements; i++) {
      // Kokkos::parallel_for(
      //     num_elements, A2D_LAMBDA(index_t i) {
      // Set up the values for the input vector
      typename ElemVec::FEDof x_dof(i, elem_xvec);
      elem_xvec.get_element_values(i, x_dof);
      QSpace xsol;
      Basis::template interp(x_dof, xsol);

      // Allocate space for the output vector
      QSpace ysol;

      for (A2D::index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace& yref = ysol.get(j);
        typename PDE::FiniteElementSpace& xref = xsol.get(j);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDE::QMatType& jac = qmat[i].get(j);

        // Matrix-vector product at the quadrature point
        yref.zero();
        for (A2D::index_t ii = 0; ii < ncomp; ii++) {
          for (A2D::index_t jj = 0; jj < ncomp; jj++) {
            yref[ii] += jac(ii, jj) * xref[jj];
          }
        }
      }

      // Add to the output-vector for the element
      typename ElemVec::FEDof y_dof(i, elem_yvec);
      Basis::template add(ysol, y_dof);
      elem_yvec.add_element_values(i, y_dof);
    }
  }

 private:
  std::vector<QMatSpace> qmat;
};

/**
 * @brief Test the implementation of the PDE to check if the derivatives are
 * consistent with the weak form.
 *
 * @tparam T Solution type
 * @tparam PDE Type of PDE object to test
 * @param pde Instance of the PDE object to test
 * @param dh Finite-difference or complex-step step size
 */
template <typename T, class PDE>
void TestPDEImplementation(PDE& pde, double dh = 1e-7) {
  Timer timer("TestPDEImplementation()");
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
  pde.weak(detJ, data, geo, s, coef);
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
  pde.weak(detJ, data, geo, s, coef);
  coef.rtransform(detJ, J, Jinv, cref);

  // Compute the Jacobian-vector product
  typename PDE::JacVecProduct jvp(pde, detJ, data, geo, s);

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

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H