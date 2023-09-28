#ifndef A2D_HELMHOLTZ_H
#define A2D_HELMHOLTZ_H

#include "a2dcore.h"
#include "multiphysics/febase.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femapping.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/fespace.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hypercube_basis.h"

namespace A2D {

template <typename T, index_t D>
class HelmholtzFilter {
 public:
  HelmholtzFilter(T r0) : r0(r0) {}

  // The filter radius
  T r0;

  // Number of dimensions
  static const index_t dim = D;

  // Number of data dimensions
  static const index_t data_dim = 1;

  // Space for the finite-element data
  using DataSpace = FESpace<T, dim, H1Space<T, data_dim, dim>>;

  // Space for the element geometry
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;

  // Finite element space
  using FiniteElementSpace = FESpace<T, dim, H1Space<T, data_dim, dim>>;

  // Define the input or output type based on wrt type
  template <FEVarType wrt>
  using FiniteElementVar =
      FEVarSelect<wrt, DataSpace, FiniteElementGeometry, FiniteElementSpace>;

  // Define the matrix Jacobian type based on the of and wrt types
  template <FEVarType of, FEVarType wrt>
  using FiniteElementJacobian =
      FESymMatSelect<of, wrt, T, DataSpace::ncomp, FiniteElementGeometry::ncomp,
                     FiniteElementSpace::ncomp>;

  /**
   * @brief Find the integral of the compliance over the entire domain
   *
   * @param weight The quadrature weight
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param sref The solution at the quadurature point
   * @return T The integrand contribution
   */
  KOKKOS_FUNCTION T integrand(T weight, const DataSpace& data,
                              const FiniteElementGeometry& geo,
                              const FiniteElementSpace& sref) const {
    T x = data[0];
    FiniteElementSpace s;
    T detJ, dot;
    const T& rho = get_value<0>(s);
    const Vec<T, dim>& rho_grad = get_grad<0>(s);

    RefElementTransform(geo, sref, detJ, s);
    VecDot(rho_grad, rho_grad, dot);
    T output =
        0.5 * weight * detJ * (rho * rho + r0 * r0 * dot - 2.0 * rho * x);

    return output;
  }

  /**
   * @brief Compute the contribution to the residual
   *
   * @tparam wrt Variable type (DATA, GEOMETRY, STATE)
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param res Residual contribution
   */
  template <FEVarType wrt>
  KOKKOS_FUNCTION void residual(T weight, const DataSpace& data,
                                const FiniteElementGeometry& geo_,
                                const FiniteElementSpace& sref_,
                                FiniteElementVar<wrt>& res) const {
    ADObj<T> x(data[0]);
    ADObj<FiniteElementSpace> sref(sref_);
    ADObj<FiniteElementGeometry> geo(geo_);

    ADObj<FiniteElementSpace> s;
    ADObj<T> detJ, dot, output;
    ADObj<T&> rho = get_value<0>(s);
    ADObj<Vec<T, dim>&> rho_grad = get_grad<0>(s);

    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),
        VecDot(rho_grad, rho_grad, dot),
        Eval(0.5 * weight * detJ * (rho * rho + r0 * r0 * dot - 2.0 * rho * x),
             output));

    output.bvalue() = 1.0;
    stack.reverse();

    if constexpr (wrt == FEVarType::DATA) {
      res[0] = x.bvalue();
    } else if constexpr (wrt == FEVarType::GEOMETRY) {
      res.copy(geo.bvalue());
    } else if constexpr (wrt == FEVarType::STATE) {
      res.copy(sref.bvalue());
    }
  }

  /**
   * @brief Compute the Jacobian-vector product
   *
   * @tparam of The residual that we're taking a derivative of
   * @tparam wrt The derivative that we're taking
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param p Direction for Jacobian-vector product
   * @param res Output product
   */
  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian_product(T weight, const DataSpace& data,
                                        const FiniteElementGeometry& geo_,
                                        const FiniteElementSpace& sref_,
                                        const FiniteElementVar<wrt>& p,
                                        FiniteElementVar<of>& res) const {
    A2DObj<T> x(data[0]);
    A2DObj<FiniteElementSpace> sref(sref_);
    A2DObj<FiniteElementGeometry> geo(geo_);

    A2DObj<FiniteElementSpace> s;
    A2DObj<T> detJ, dot, output;
    A2DObj<T&> rho = get_value<0>(s);
    A2DObj<Vec<T, dim>&> rho_grad = get_grad<0>(s);

    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),
        VecDot(rho_grad, rho_grad, dot),
        Eval(0.5 * weight * detJ * (rho * rho + r0 * r0 * dot - 2.0 * rho * x),
             output));

    output.bvalue() = 1.0;
    stack.hproduct();

    if constexpr (of == FEVarType::DATA) {
      res[0] = x.hvalue();
    } else if constexpr (of == FEVarType::GEOMETRY) {
      res.copy(geo.hvalue());
    } else if constexpr (of == FEVarType::STATE) {
      res.copy(sref.hvalue());
    }
  }

  /**
   * @brief Compute the Jacobian at a quadrature point
   *
   * @tparam of The residual that we're taking a derivative of
   * @tparam wrt The derivative that we're taking
   * @param weight Quadrature weight
   * @param data Data at the quadrature point
   * @param geo_ Geometry data at the quadrature point
   * @param sref_ State at the quadrature point
   * @param jac The Jacobian output
   */
  template <FEVarType of, FEVarType wrt>
  KOKKOS_FUNCTION void jacobian(T weight, const DataSpace& data,
                                const FiniteElementGeometry& geo_,
                                const FiniteElementSpace& sref_,
                                FiniteElementJacobian<of, wrt>& jac) const {
    A2DObj<T> x(data[0]);
    A2DObj<FiniteElementSpace> sref(sref_);
    A2DObj<FiniteElementGeometry> geo(geo_);

    A2DObj<FiniteElementSpace> s;
    A2DObj<T> detJ, dot, output;
    A2DObj<T&> rho = get_value<0>(s);
    A2DObj<Vec<T, dim>&> rho_grad = get_grad<0>(s);

    auto stack = MakeStack(
        RefElementTransform(geo, sref, detJ, s),
        VecDot(rho_grad, rho_grad, dot),
        Eval(0.5 * weight * detJ * (rho * rho + r0 * r0 * dot - 2.0 * rho * x),
             output));

    output.bvalue() = 1.0;

    // Extract the Jacobian
    ExtractJacobian<of, wrt>(stack, x, geo, sref, jac);
  }
};

template <class Impl, index_t degree>
class HexHelmholtzFilterElement
    : public ElementIntegrand<
          Impl, HelmholtzFilter<typename Impl::type, 3>,
          HexGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 3, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1HexBasis<typename Impl::type, 1, degree>>> {
 public:
  using T = typename Impl::type;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;

  HexHelmholtzFilterElement(HelmholtzFilter<T, 3> integrand,
                            std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                            std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                            std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const HelmholtzFilter<T, 3>& get_integrand() { return integrand; }

 private:
  HelmholtzFilter<T, 3> integrand;
};

template <class Impl, index_t degree>
class QuadHelmholtzFilterElement
    : public ElementIntegrand<
          Impl, HelmholtzFilter<typename Impl::type, 2>,
          QuadGaussQuadrature<degree + 1>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 1, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 2, degree>>,
          FEBasis<typename Impl::type,
                  LagrangeH1QuadBasis<typename Impl::type, 1, degree>>> {
 public:
  using T = typename Impl::type;
  using Integrand = HelmholtzFilter<T, 2>;
  using DataBasis = FEBasis<T, LagrangeH1QuadBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1QuadBasis<T, 2, degree>>;
  using Basis = FEBasis<T, LagrangeH1QuadBasis<T, 1, degree>>;
  using Vec_t = typename Impl::Vec_t;
  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  QuadHelmholtzFilterElement(HelmholtzFilter<T, 2> integrand,
                             std::shared_ptr<ElementMesh<DataBasis>> data_mesh,
                             std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh,
                             std::shared_ptr<ElementMesh<Basis>> sol_mesh)
      : integrand(integrand) {
    this->set_meshes(data_mesh, geo_mesh, sol_mesh);
  }

  const HelmholtzFilter<T, 2>& get_integrand() { return integrand; }

  void to_vtk(Vec_t& data, Vec_t& geo, Vec_t& sol, const std::string filename) {
    const int num_out = 2;  // Number of outputs to the file
    ElementVector<DataBasis> elem_data(*this->data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*this->geo_mesh, geo);
    ElementVector<Basis> elem_sol(*this->sol_mesh, sol);

    write_quad_to_vtk<num_out, degree, T, DataBasis, GeoBasis, Basis,
                      Integrand>(
        elem_data, elem_geo, elem_sol, filename,
        [](index_t k, typename Integrand::DataSpace& d,
           typename Integrand::FiniteElementGeometry& g,
           typename Integrand::FiniteElementSpace& s) {
          if (k == 0) {
            return get_value<0>(s);
          } else if (k == 1) {
            return get_value<0>(d);
          }
          return 0.0;
        });
  }

 private:
  HelmholtzFilter<T, 2> integrand;
};

}  // namespace A2D

#endif  // A2D_HELMHOLTZ_H