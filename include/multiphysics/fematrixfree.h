#ifndef FE_MATRIX_FREE_H
#define FE_MATRIX_FREE_H

#include <iomanip>
#include <iostream>
#include <random>
#include <type_traits>

#include "multiphysics/feelementmat.h"
#include "multiphysics/feelementvector.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fesolution.h"
#include "multiphysics/fespace.h"
#include "utils/a2dprofiler.h"

namespace A2D {

template <typename T, FEVarType of, FEVarType wrt, class Integrand,
          class Quadrature, class DataBasis, class GeoBasis, class Basis>
class MatrixFree {
 public:
  // Quadrature point object for the data space
  using QDataSpace = QptSpace<Quadrature, typename Integrand::DataSpace>;

  // Quadrature point object for the geometry
  using QGeoSpace =
      QptSpace<Quadrature, typename Integrand::FiniteElementGeometry>;

  // Quadrature point object for the finite-element space
  using QSpace = QptSpace<Quadrature, typename Integrand::FiniteElementSpace>;

  // Quadrature point view of the Jacobian-matrices
  using QMatSpace =
      QptSpace<Quadrature,
               typename Integrand::template FiniteElementJacobian<of, wrt>>;

  MatrixFree() {}

  template <class DataElemVec, class GeoElemVec, class ElemVec>
  void initialize(const Integrand& integrand, DataElemVec& elem_data,
                  GeoElemVec& elem_geo, ElemVec& elem_sol) {
    using same_evtype = have_same_evtype<DataElemVec, GeoElemVec, ElemVec>;
    static_assert(same_evtype::value,
                  "Cannot mix up different element vector types (e.g. using "
                  "parallel and serial at the same time)");
    constexpr ElemVecType evtype = same_evtype::evtype;

    Timer timer("MatrixFree::initialize()");
    // Re-size the vector as needed
    const index_t num_elements = elem_geo.get_num_elements();
    if (qmat.size() != num_elements) {
      qmat.resize(num_elements);
    }

    // Number of components at the quadrature point
    const index_t ncomp = Integrand::FiniteElementSpace::ncomp;

    // Get the number of quadrature points
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
        typename Integrand::FiniteElementSpace& sref = sol.get(j);
        typename Integrand::FiniteElementGeometry& gref = geo.get(j);

        // Initialize the transform object
        T detJ;
        typename Integrand::SolutionMapping transform(gref, detJ);

        // Transform the solution the physical element
        typename Integrand::FiniteElementSpace s;
        transform.transform(sref, s);

        // Compute the Jacobian for the weak form at the quadrature point
        double weight = Quadrature::get_weight(j);
        typename Integrand::QMatType jac;
        integrand.jacobian(weight * detJ, data.get(j), gref, s, jac);

        // Transform second derivatives from w.r.t. x to w.r.t. xi
        typename Integrand::QMatType& jac_ref = qmat[i].get(j);
        transform.template jtransform<typename Integrand::FiniteElementSpace>(
            jac, jac_ref);
      }
    }
  }

  template <class ElemVec>
  void add_jacobian_vector_product(ElemVec& elem_xvec, ElemVec& elem_yvec) {
    constexpr ElemVecType evtype = ElemVec::evtype;

    const index_t num_elements = qmat.size();
    const index_t num_quadrature_points = Quadrature::get_num_points();
    const index_t ncomp = Integrand::FiniteElementSpace::ncomp;

    if constexpr (evtype == ElemVecType::Parallel) {
      elem_xvec.get_values();
      elem_yvec.get_zero_values();
    }

    for (index_t i = 0; i < num_elements; i++) {
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
        typename Integrand::FiniteElementSpace& yref = ysol.get(j);
        typename Integrand::FiniteElementSpace& xref = xsol.get(j);

        // The entries of the Jacobian matrix at the quadrature point
        typename Integrand::QMatType& jac = qmat[i].get(j);

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

}  // namespace A2D

#endif  // FE_MATRIX_FREE_H