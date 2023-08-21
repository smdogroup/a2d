#ifndef A2D_FE_ELEMENT_H
#define A2D_FE_ELEMENT_H

#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>

#include "multiphysics/feelementmat.h"
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
template <typename T, class PDEIntegrand, class GeoBasis, class Basis>
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
      QptSpace<DOFQuadrature, typename PDEIntegrand::FiniteElementGeometry>;

  template <ElemVecType evtype, class GeoElemVec, class ElemVec>
  void get_dof_coordinates(ElementVectorBase<evtype, GeoElemVec>& elem_geo,
                           ElementVectorBase<evtype, ElemVec>& elem_x,
                           ElementVectorBase<evtype, ElemVec>& elem_y,
                           ElementVectorBase<evtype, ElemVec>& elem_z) {
    const index_t num_elements = elem_geo.get_num_elements();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_geo.get_values();
    }
    for (index_t i = 0; i < num_elements; i++) {
      // Get the geometry values and interpolate them at all quadrature points
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      if constexpr (evtype == ElemVecType::Serial) {
        elem_geo.get_element_values(i, geo_dof);
      }
      QGeoSpace geo;
      GeoBasis::template interp(geo_dof, geo);

      // Get the degrees of freedom for the element and interpolate the solution
      typename ElemVec::FEDof x_dof(i, elem_x);
      typename ElemVec::FEDof y_dof(i, elem_y);
      typename ElemVec::FEDof z_dof(i, elem_z);

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < Basis::ndof; j++) {
        // Get the solution/geometry in the reference domain
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);
        Vec<T, PDEIntegrand::dim>& X = gref.template get<0>().get_value();

        x_dof[j] = X(0);
        y_dof[j] = X(1);
        z_dof[j] = X(2);
      }

      if constexpr (evtype == ElemVecType::Serial) {
        elem_x.set_element_values(i, x_dof);
        elem_y.set_element_values(i, y_dof);
        elem_z.set_element_values(i, z_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_x.set_values();
      elem_y.set_values();
      elem_z.set_values();
    }
  }
};

template <typename T, class PDEIntegrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class FiniteElement {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename PDEIntegrand::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace =
      QptSpace<Quadrature, typename PDEIntegrand::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace =
      QptSpace<Quadrature, typename PDEIntegrand::FiniteElementSpace>;

  FiniteElement() {}

  /**
   * @brief Compute the value of an integral functional over the finite-element
   * domain
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The integral over the element
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec>
  T integrate(PDEIntegrand& pde,
              ElementVectorBase<evtype, DataElemVec>& elem_data,
              ElementVectorBase<evtype, GeoElemVec>& elem_geo,
              ElementVectorBase<evtype, ElemVec>& elem_sol) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    T value = 0.0;

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDEIntegrand
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::FiniteElementSpace coef;
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
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The maximum value over the domain
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec>
  T max(PDEIntegrand& pde, ElementVectorBase<evtype, DataElemVec>& elem_data,
        ElementVectorBase<evtype, GeoElemVec>& elem_geo,
        ElementVectorBase<evtype, ElemVec>& elem_sol) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    // TODO: Need to make this more robust?
    T max_value = -1e20;

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDEIntegrand
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::FiniteElementSpace coef;
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
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_deriv Output derivative of the functional w.r.t. data
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec, class DataDerivElemVec>
  void add_data_derivative(
      PDEIntegrand& pde, ElementVectorBase<evtype, DataElemVec>& elem_data,
      ElementVectorBase<evtype, GeoElemVec>& elem_geo,
      ElementVectorBase<evtype, ElemVec>& elem_sol,
      ElementVectorBase<evtype, DataDerivElemVec>& elem_deriv) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
      elem_deriv.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Allocate space for the derivative
      QDataSpace deriv;

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDEIntegrand
        double weight = Quadrature::get_weight(j);
        pde.data_derivative(weight * detJ, data.get(j), gref, s, deriv.get(j));
      }

      // Add the derivative of the data back to the data space
      typename DataDerivElemVec::FEDof deriv_dof(i, elem_deriv);
      DataBasis::template add(deriv, deriv_dof);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_deriv.add_element_values(i, deriv_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_deriv.add_values();
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
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_adj Element adjoint vector
   * @param elem_deriv Data derivative
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemAdjVec, class DataDerivElemVec>
  void add_adjoint_residual_data_derivative(
      PDEIntegrand& pde, ElementVectorBase<evtype, DataElemVec>& elem_data,
      ElementVectorBase<evtype, GeoElemVec>& elem_geo,
      ElementVectorBase<evtype, ElemVec>& elem_sol,
      ElementVectorBase<evtype, ElemAdjVec>& elem_adj,
      ElementVectorBase<evtype, DataDerivElemVec>& elem_deriv) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
      elem_adj.get_values();
      elem_deriv.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get data, geo, sol and adj for this element and interpolate them
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      typename ElemAdjVec::FEDof adj_dof(i, elem_adj);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
        elem_adj.get_element_values(i, adj_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;
      QSpace adj;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);
      Basis::template interp(adj_dof, adj);

      // Allocate space for the derivative
      QDataSpace deriv;

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry in the reference domain
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementSpace& aref = adj.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s, a;
        transform.transform(sref, s);
        transform.transform(aref, a);

        // Allocate the Jacobian-adjoint product functor
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::AdjVecProduct ajp(pde, weight * detJ,
                                                 data.get(j), gref, s);

        ajp(a, deriv.get(j));
      }

      // Add the derivative of the data back to the data space
      typename DataDerivElemVec::FEDof deriv_dof(i, elem_deriv);
      DataBasis::template add(deriv, deriv_dof);
      if constexpr (evtype == ElemVecType::Serial) {
        elem_deriv.add_element_values(i, deriv_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_deriv.add_values();
    }
  }

  /**
   * @brief Add the residuals for the finite-element problem
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution
   * @tparam ElemResVec Element vector class for the residual
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_res Element residual vector
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemResVec>
  void add_residual(PDEIntegrand& pde,
                    ElementVectorBase<evtype, DataElemVec>& elem_data,
                    ElementVectorBase<evtype, GeoElemVec>& elem_geo,
                    ElementVectorBase<evtype, ElemVec>& elem_sol,
                    ElementVectorBase<evtype, ElemResVec>& elem_res) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
      elem_res.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      // Evaluate values (and potentially derivatives) for all quadrature
      // points
      // Note: derivatives computed at this point are all w.r.t. computational
      // coordinates!
      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Allocate space for the residual values at each quadrature point
      QSpace res;

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Get the solution/geometry and derivatives in the reference domain
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref =
            geo.get(j);  // this has J

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the coefficients for the weak form of the PDEIntegrand
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::FiniteElementSpace coef;
        pde.weak(weight * detJ, data.get(j), gref, s, coef);

        // Transform the coefficients back to the reference element
        typename PDEIntegrand::FiniteElementSpace& cref = res.get(j);
        transform.rtransform(coef, cref);
      }

      // Add the residual from the quadrature points back to the finite-element
      // mesh
      typename ElemVec::FEDof res_dof(i, elem_res);
      Basis::template add(res, res_dof);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_res.add_element_values(i, res_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_res.add_values();
    }
  }

  /**
   * @brief Compute the matrix-free Jacobian-vector product y = y + J * x
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param pde Instance of the PDEIntegrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_xvec Element solution vector for storing x-components
   * @param elem_yvec Output element solution vector storing y-components
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec>
  void add_jacobian_vector_product(
      PDEIntegrand& pde, ElementVectorBase<evtype, DataElemVec>& elem_data,
      ElementVectorBase<evtype, GeoElemVec>& elem_geo,
      ElementVectorBase<evtype, ElemVec>& elem_sol,
      ElementVectorBase<evtype, ElemVec>& elem_xvec,
      ElementVectorBase<evtype, ElemVec>& elem_yvec) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
      elem_xvec.get_values();
      elem_yvec.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry, solution and input vector for this element and
      // interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      typename ElemVec::FEDof x_dof(i, elem_xvec);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
        elem_xvec.get_element_values(i, x_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;
      QSpace xsol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);
      Basis::template interp(x_dof, xsol);

      // Allocate space for the output vector
      QSpace ysol;

      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementSpace& xref = xsol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace x, s;
        transform.transform(sref, s);
        transform.transform(xref, x);

        // Allocate the Jacobian-vector product functor
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::JacVecProduct jvp(pde, weight * detJ,
                                                 data.get(j), gref, s);

        // Compute the Jacobian-vector product
        typename PDEIntegrand::FiniteElementSpace y;
        jvp(x, y);

        // Transform to back to the reference element
        typename PDEIntegrand::FiniteElementSpace& yref = ysol.get(j);
        transform.rtransform(y, yref);
      }

      // Add the values from the quadrature points back into the finite-element
      // problem
      typename ElemVec::FEDof y_dof(i, elem_yvec);
      Basis::template add(ysol, y_dof);
      if constexpr (evtype == ElemVecType::Serial) {
        elem_yvec.add_element_values(i, y_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_yvec.add_values();
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
   * @param pde The PDEIntegrand instance
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_mat Element matrix output
   *
   * Note: this function uses Jacobian-vector product, which is deprecated and
   * will be removed soon
   */
  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemMat>
  void add_jacobian(PDEIntegrand& pde,
                    ElementVectorBase<evtype, DataElemVec>& elem_data,
                    ElementVectorBase<evtype, GeoElemVec>& elem_geo,
                    ElementVectorBase<evtype, ElemVec>& elem_sol,
                    ElemMat& elem_mat) {
    Timer timer("FiniteElement::add_jacobian()");
    const index_t ncomp = PDEIntegrand::FiniteElementSpace::ncomp;
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Initialize the element matrix
      typename ElemMat::FEMat element_mat(i, elem_mat);

      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // // Allocate the Jacobian-vector product functor
        // double weight = Quadrature::get_weight(j);
        // typename PDEIntegrand::JacVecProduct jvp(pde, weight * detJ,
        // data.get(j), gref, s);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDEIntegrand::QMatType jac;

        // Temporary vectors
        typename PDEIntegrand::FiniteElementSpace pref, p, Jp;

        double weight = Quadrature::get_weight(j);

        for (index_t k = 0; k < ncomp; k++) {
          // Set the value into the matrix
          pref.zero();
          pref[k] = T(1.0);
          transform.transform(pref, p);

          // Allocate the Jacobian-vector product functor
          typename PDEIntegrand::JacVecProduct jvp(pde, weight * detJ,
                                                   data.get(j), gref, s);

          // Compute the Jacobian-vector product
          jvp(p, Jp);

          // Transform to back to the reference element
          transform.rtransform(Jp, pref);

          for (index_t m = 0; m < ncomp; m++) {
            jac(m, k) = pref[m];
          }
        }

        // Add the results of the outer product
        Basis::template add_outer<Quadrature>(j, jac, element_mat);
      }

      elem_mat.add_element_values(i, element_mat);
    }
  }

  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemMat>
  void add_jacobian_new(PDEIntegrand& pde,
                        ElementVectorBase<evtype, DataElemVec>& elem_data,
                        ElementVectorBase<evtype, GeoElemVec>& elem_geo,
                        ElementVectorBase<evtype, ElemVec>& elem_sol,
                        ElemMat& elem_mat) {
    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Initialize the element matrix
      typename ElemMat::FEMat element_mat(i, elem_mat);

      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform from the reference element to the physical space
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the Jacobian for the weak form at the quadrature point
        double weight = Quadrature::get_weight(j);
        typename PDEIntegrand::QMatType jac_ref, jac;
        pde.jacobian(weight * detJ, data.get(j), gref, s, jac);

        // Transform second derivatives from w.r.t. x to w.r.t. xi
        transform.template jtransform<PDEIntegrand>(jac, jac_ref);

        // Add the results of the outer product
        Basis::template add_outer<Quadrature>(j, jac_ref, element_mat);
      }

      elem_mat.add_element_values(i, element_mat);
    }
  }

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void add_geo_derivative(PDEIntegrand& pde, DataElemVec& elem_data,
                          GeoElemVec& elem_geo, ElemVec& elem_sol,
                          GeoElemVec& geo_deriv) {}

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void add_adjoint_residual_geo_derivative(PDEIntegrand& pde,
                                           DataElemVec& elem_data,
                                           GeoElemVec& elem_geo,
                                           ElemVec& elem_sol, ElemVec& elem_adj,
                                           GeoElemVec& geo_deriv) {}
};

template <typename T, class PDEIntegrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class MatrixFree {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename PDEIntegrand::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace =
      QptSpace<Quadrature, typename PDEIntegrand::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace =
      QptSpace<Quadrature, typename PDEIntegrand::FiniteElementSpace>;

  // Quadrature point view of the Jacobian-matrices
  using QMatSpace = QptSpace<Quadrature, typename PDEIntegrand::QMatType>;

  MatrixFree() {}

  template <ElemVecType evtype, class DataElemVec, class GeoElemVec,
            class ElemVec>
  void initialize(PDEIntegrand& pde,
                  ElementVectorBase<evtype, DataElemVec>& elem_data,
                  ElementVectorBase<evtype, GeoElemVec>& elem_geo,
                  ElementVectorBase<evtype, ElemVec>& elem_sol) {
    Timer timer("MatrixFree::initialize()");
    // Re-size the vector as needed
    const index_t num_elements = elem_geo.get_num_elements();
    if (qmat.size() != num_elements) {
      qmat.resize(num_elements);
    }

    // Number of components at the quadrature point
    const index_t ncomp = PDEIntegrand::FiniteElementSpace::ncomp;

    // Get the number of quadrature points
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDEIntegrand::FiniteElementSpace& sref = sol.get(j);
        typename PDEIntegrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename PDEIntegrand::SolutionMapping transform(gref, detJ);

        // Transform the solution the physical element
        typename PDEIntegrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // // Allocate the Jacobian-vector product functor
        // double weight = Quadrature::get_weight(j);
        // typename PDEIntegrand::JacVecProduct jvp(pde, weight * detJ,
        // data.get(j), gref, s);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDEIntegrand::QMatType& jac = qmat[i].get(j);

        // Temporary vectors
        typename PDEIntegrand::FiniteElementSpace pref, p, Jp;

        for (index_t k = 0; k < ncomp; k++) {
          // Set the value into the matrix
          pref.zero();
          pref[k] = T(1.0);
          transform.transform(pref, p);

          // Allocate the Jacobian-vector product functor
          double weight = Quadrature::get_weight(j);
          typename PDEIntegrand::JacVecProduct jvp(pde, weight * detJ,
                                                   data.get(j), gref, s);

          // Compute the Jacobian-vector product
          jvp(p, Jp);

          // Transform to back to the reference element
          transform.rtransform(Jp, pref);

          for (index_t m = 0; m < ncomp; m++) {
            jac(m, k) = pref[m];
          }
        }
      }
    }
  }

  template <ElemVecType evtype, class ElemVec>
  void add_jacobian_vector_product(
      ElementVectorBase<evtype, ElemVec>& elem_xvec,
      ElementVectorBase<evtype, ElemVec>& elem_yvec) {
    const index_t num_elements = qmat.size();
    const index_t num_quadrature_points = Quadrature::get_num_points();
    const index_t ncomp = PDEIntegrand::FiniteElementSpace::ncomp;

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_xvec.get_values();
      elem_yvec.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Kokkos::parallel_for(
      //     num_elements, A2D_LAMBDA(index_t i) {
      // Set up the values for the input vector
      typename ElemVec::FEDof x_dof(i, elem_xvec);
      if constexpr (evtype == ElemVecType::Serial) {
        elem_xvec.get_element_values(i, x_dof);
      }
      QSpace xsol;
      Basis::template interp(x_dof, xsol);

      // Allocate space for the output vector
      QSpace ysol;

      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Transform to the local coordinate system
        typename PDEIntegrand::FiniteElementSpace& yref = ysol.get(j);
        typename PDEIntegrand::FiniteElementSpace& xref = xsol.get(j);

        // The entries of the Jacobian matrix at the quadrature point
        typename PDEIntegrand::QMatType& jac = qmat[i].get(j);

        // Matrix-vector product at the quadrature point
        yref.zero();
        for (index_t ii = 0; ii < ncomp; ii++) {
          for (index_t jj = 0; jj < ncomp; jj++) {
            yref[ii] += jac(ii, jj) * xref[jj];
          }
        }
      }

      // Add to the output-vector for the element
      typename ElemVec::FEDof y_dof(i, elem_yvec);
      Basis::template add(ysol, y_dof);
      if constexpr (evtype == ElemVecType::Serial) {
        elem_yvec.add_element_values(i, y_dof);
      }
    }
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_yvec.add_values();
    }
  }

 private:
  std::vector<QMatSpace> qmat;
};

/**
 * @brief Test the implementation of the PDEIntegrand to check if the
 * derivatives are consistent with the weak form.
 *
 * @tparam T Solution type
 * @tparam PDEIntegrand Type of PDEIntegrand object to test
 * @param pde Instance of the PDEIntegrand object to test
 * @param dh Finite-difference or complex-step step size
 */
template <typename T, class PDEIntegrand>
void TestPDEImplementation(PDEIntegrand& pde, double dh = 1e-7) {
  Timer timer("TestPDEImplementation()");
  typename PDEIntegrand::DataSpace data;
  typename PDEIntegrand::FiniteElementGeometry geo;
  typename PDEIntegrand::FiniteElementSpace s, sref;
  typename PDEIntegrand::FiniteElementSpace p, pref;
  typename PDEIntegrand::FiniteElementSpace coef, cref, cref0;
  typename PDEIntegrand::FiniteElementSpace Jp, Jpref;

  // Generate random data
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> distr(-1.0, 1.0);

  // Set random values for the data
  if constexpr (PDEIntegrand::DataSpace::ncomp > 0) {
    for (index_t i = 0; i < PDEIntegrand::DataSpace::ncomp; i++) {
      data[i] = distr(gen);
    }
  }

  // Set random values for the geometry
  for (index_t i = 0; i < PDEIntegrand::FiniteElementGeometry::ncomp; i++) {
    geo[i] = distr(gen);
  }

  // Set the random values
  for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
    sref[i] = distr(gen);
    pref[i] = distr(gen);
  }

  // Initialize the transform object
  T detJ;
  typename PDEIntegrand::SolutionMapping transform(geo, detJ);

  // Compute the coefficients
  transform.transform(sref, s);
  pde.weak(detJ, data, geo, s, coef);
  transform.rtransform(coef, cref0);

  if constexpr (std::is_same<T, std::complex<double>>::value) {
    for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
      sref[i] = sref[i] + dh * pref[i] * std::complex<double>(0.0, 1.0);
    }
  } else {
    for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
      sref[i] = sref[i] + dh * pref[i];
    }
  }

  // Compute the coefficients
  transform.transform(sref, s);
  pde.weak(detJ, data, geo, s, coef);
  transform.rtransform(coef, cref);

  // Compute the Jacobian-vector product
  typename PDEIntegrand::JacVecProduct jvp(pde, detJ, data, geo, s);

  // Compute the Jacobian-vector product
  transform.transform(pref, p);
  jvp(p, Jp);
  transform.rtransform(Jp, Jpref);

  // Compute the finite-difference value
  typename PDEIntegrand::FiniteElementSpace fd;

  if constexpr (std::is_same<T, std::complex<double>>::value) {
    for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
      fd[i] = std::imag(cref[i]) / dh;
    }
  } else {
    for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
      fd[i] = (cref[i] - cref0[i]) / dh;
    }
  }

  for (index_t i = 0; i < PDEIntegrand::FiniteElementSpace::ncomp; i++) {
    std::cout << "fd[" << std::setw(2) << i << "]: " << std::setw(12)
              << std::real(fd[i]) << " Jpref[" << std::setw(2) << i
              << "]: " << std::setw(12) << std::real(Jpref[i]) << " err["
              << std::setw(2) << i << "]: " << std::setw(12)
              << std::real((fd[i] - Jpref[i]) / fd[i]) << std::endl;
  }
}

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H