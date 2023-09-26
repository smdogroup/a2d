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
#include "multiphysics/fespace.h"
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
template <typename T, class Integrand, class GeoBasis, class Basis>
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
      QptSpace<DOFQuadrature, typename Integrand::FiniteElementGeometry>;

  template <class GeoElemVec, class ElemVec>
  void get_dof_coordinates(GeoElemVec& elem_geo, ElemVec& elem_x,
                           ElemVec& elem_y, ElemVec& elem_z) {
    using same_evtype = have_same_evtype<GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

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
        typename Integrand::FiniteElementGeometry& gref = geo.get(j);
        Vec<T, Integrand::dim>& X = gref.template get<0>().get_value();

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

// A set of element-wise operations that operate on all elements at once
template <typename T, class Integrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class FiniteElement {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename Integrand::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace =
      QptSpace<Quadrature, typename Integrand::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace = QptSpace<Quadrature, typename Integrand::FiniteElementSpace>;

  // Select the Quadrature point object based on the type of object
  template <FEVarType wrt>
  using QSpaceSelect = FEVarSelect<wrt, QDataSpace, QGeoSpace, QSpace>;

  FiniteElement() {}

  /**
   * @brief Compute the value of an integral functional over the finite-element
   * domain
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @param integrand Instance of the Integrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The integral over the element
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec>
  T integrate(const Integrand& integrand, DataElemVec& elem_data,
              GeoElemVec& elem_geo, ElemVec& elem_sol) {
    using same_evtype = have_same_evtype<DataElemVec, GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

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
        double weight = Quadrature::get_weight(j);
        value +=
            integrand.integrand(weight, data.get(j), geo.get(j), sol.get(j));
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
   * @param integrand Instance of the Integrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @return The maximum value over the domain
   */
  template <class DataElemVec, class GeoElemVec, class ElemVec>
  T max(const Integrand& integrand, DataElemVec& elem_data,
        GeoElemVec& elem_geo, ElemVec& elem_sol) {
    using same_evtype = have_same_evtype<DataElemVec, GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

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
        double weight = Quadrature::get_weight(j);
        T value = integrand.max(data.get(j), geo.get(j), sol.get(j));
        if (std::real(value) > std::real(max_value)) {
          max_value = value;
        }
      }
    }

    return max_value;
  }

  /**
   * @brief Add the residuals for the finite-element problem
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution
   * @tparam ElemResVec Element vector class for the residual
   * @param integrand Instance of the Integrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_res Element residual vector
   */
  template <FEVarType wrt, class DataElemVec, class GeoElemVec, class ElemVec,
            class ElemResVec>
  void add_residual(const Integrand& integrand, const T alpha,
                    DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemResVec& elem_res) {
    using same_evtype =
        have_same_evtype<DataElemVec, GeoElemVec, ElemVec, ElemResVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial element vector at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    auto loop_body = KOKKOS_LAMBDA(const index_t i) {
      // Get the data, geometry and solution for this element and
      // interpolate it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);

      // Serialized scatter - only effective for a serialized element vector
      // implementation
      elem_data.get_element_values(i, data_dof);
      elem_geo.get_element_values(i, geo_dof);
      elem_sol.get_element_values(i, sol_dof);

      // Space needed for the data, geo and solution at each quadrature point
      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      // Space for the residual
      QSpaceSelect<wrt> res;

      // Evaluate values (and potentially derivatives) for all quadrature
      // points. Note: derivatives computed at this point are all w.r.t.
      // computational coordinates!
      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      // Compute the weak coefficients at all quadrature points
      for (index_t j = 0; j < num_quadrature_points; j++) {
        T weight = alpha * Quadrature::get_weight(j);
        integrand.template residual<wrt>(weight, data.get(j), geo.get(j),
                                         sol.get(j), res.get(j));
      }

      // Add the residual from the quadrature points back to the
      // finite-element mesh
      typename ElemResVec::FEDof res_dof(i, elem_res);
      if constexpr (wrt == FEVarType::DATA) {
        DataBasis::template add(res, res_dof);
      } else if constexpr (wrt == FEVarType::GEOMETRY) {
        GeoBasis::template add(res, res_dof);
      } else if constexpr (wrt == FEVarType::STATE) {
        Basis::template add(res, res_dof);
      }

      // Serialized gather - only effective for a serialized element vector
      // implementation
      elem_res.add_element_values(i, res_dof);
    };

    // Execution
    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
      elem_res.get_zero_values();

      Kokkos::parallel_for("add_residual", num_elements, loop_body);
      Kokkos::fence();

      elem_res.add_values();
    } else {
      static_assert(evtype == ElemVecType::Serial,
                    "invalid ElemVecType deduced.");
      for (index_t i = 0; i < num_elements; i++) {
        loop_body(i);
      }
    }
  }

  /**
   * @brief Add the residuals for the finite-element problem
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution
   * @tparam ElemResVec Element vector class for the residual
   * @param integrand Instance of the Integrand
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_res Element residual vector
   */
  template <FEVarType of, FEVarType wrt, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemProdVec, class ElemResVec>
  void add_jacobian_product(const Integrand& integrand, const T alpha,
                            DataElemVec& elem_data, GeoElemVec& elem_geo,
                            ElemVec& elem_sol, ElemProdVec& elem_prod,
                            ElemResVec& elem_res) {
    using same_evtype = have_same_evtype<DataElemVec, GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and interpolate
      // it
      typename DataElemVec::FEDof data_dof(i, elem_data);
      typename GeoElemVec::FEDof geo_dof(i, elem_geo);
      typename ElemVec::FEDof sol_dof(i, elem_sol);
      typename ElemProdVec::FEDof prod_dof(i, elem_prod);
      typename ElemResVec::FEDof res_dof(i, elem_res);

      if constexpr (evtype == ElemVecType::Serial) {
        elem_data.get_element_values(i, data_dof);
        elem_geo.get_element_values(i, geo_dof);
        elem_sol.get_element_values(i, sol_dof);
        elem_prod.get_element_values(i, prod_dof);
      }

      QDataSpace data;
      QGeoSpace geo;
      QSpace sol;

      // Space for the residual
      QSpaceSelect<of> res;
      QSpaceSelect<wrt> prod;

      DataBasis::template interp(data_dof, data);
      GeoBasis::template interp(geo_dof, geo);
      Basis::template interp(sol_dof, sol);

      if constexpr (wrt == FEVarType::DATA) {
        DataBasis::template interp(prod_dof, prod);
      } else if constexpr (wrt == FEVarType::GEOMETRY) {
        GeoBasis::template interp(prod_dof, prod);
      } else if constexpr (wrt == FEVarType::STATE) {
        Basis::template interp(prod_dof, prod);
      }

      for (index_t j = 0; j < num_quadrature_points; j++) {
        T weight = alpha * Quadrature::get_weight(j);
        typename Integrand::template FiniteElementJacobian<of, wrt> jac;
        integrand.template jacobian_product<of, wrt>(weight, data.get(j),
                                                     geo.get(j), sol.get(j),
                                                     prod.get(j), res.get(j));
      }

      if constexpr (of == FEVarType::DATA) {
        DataBasis::template add(res, res_dof);
      } else if constexpr (of == FEVarType::GEOMETRY) {
        GeoBasis::template add(res, res_dof);
      } else if constexpr (of == FEVarType::STATE) {
        Basis::template add(res, res_dof);
      }

      // Serialized gather - only effective for a serialized element vector
      // implementation
      elem_res.add_element_values(i, res_dof);
    }
  }

  /**
   * @brief Assemble element Jacobian matrices based on the data, geometry and
   * solution vectors.
   *
   * WARNING: This is intended only for the lowest order elements, ie. p = 1.
   * It scales O(p^9) so it is unsuitable for high-order elements!
   *
   * @tparam DataElemVec Element vector class for the data
   * @tparam GeoElemVec Element vector class for the geometry
   * @tparam ElemVec Element vector class for the solution/residual
   * @tparam ElemMat The element matrix
   * @param integrand The Integrand instance
   * @param elem_data Element vector for the data
   * @param elem_geo Element vector for the geometry
   * @param elem_sol Element solution vector
   * @param elem_mat Element matrix output
   *
   * Note: this function uses Jacobian-vector product, which is deprecated and
   * will be removed soon
   */
  template <FEVarType of, FEVarType wrt, class DataElemVec, class GeoElemVec,
            class ElemVec, class ElemMat>
  void add_jacobian(const Integrand& integrand, const T alpha,
                    DataElemVec& elem_data, GeoElemVec& elem_geo,
                    ElemVec& elem_sol, ElemMat& elem_mat) {
    Timer timer("FiniteElement::add_jacobian()");

    using same_evtype = have_same_evtype<DataElemVec, GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

    const index_t num_elements = elem_geo.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_data.get_values();
      elem_geo.get_values();
      elem_sol.get_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
      // Get the data, geometry and solution for this element and
      // interpolate it
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
        T weight = alpha * Quadrature::get_weight(j);
        typename Integrand::template FiniteElementJacobian<of, wrt> jac;
        integrand.template jacobian<of, wrt>(weight, data.get(j), geo.get(j),
                                             sol.get(j), jac);

        // Add the results of the outer product
        Basis::template add_outer<Quadrature>(j, jac, element_mat);
      }

      elem_mat.add_element_values(i, element_mat);
    }
  }
};

}  // namespace A2D

#endif  // A2D_FE_ELEMENT_H