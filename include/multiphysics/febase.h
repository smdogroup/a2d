#ifndef A2D_FE_BASE_H
#define A2D_FE_BASE_H

#include <algorithm>
#include <list>
#include <set>

#include "a2ddefs.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"

namespace A2D {

/**
 * The classes in this file use an Impl object that should contain the following
 * information
 *
 * class Impl {
 *  public: using type = double;
 *
 *  using Vec_t = SolutionVector<type>;
 *
 *  using Mat_t = BSRMat<type, bsize, bsize>;
 *
 *  template <class Basis>
 *  using ElementVector<Basis> = ElementVector_Serial<type, Basis, Vec_t>;
 *
 *  template <class Basis>
 *  using ElementMatrix<Basis> = ElementMat_Serial<type, Basis, Vec_t>;
 * };
 *
 */

/**
 * @brief Element base class.
 *
 * This defines an element that is compatible with the given vector and matrix
 * types.
 */
template <class Impl>
class ElementBase {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;
  using Mat_t = typename Impl::Mat_t;

  virtual ~ElementBase() {}

  virtual const ElementMeshBase& get_mesh(FEVarType wrt) const = 0;

  // Add the residual
  virtual void add_residual(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol,
                            Vec_t& res) = 0;

  // Add the Jacobian
  virtual void add_jacobian(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol,
                            Mat_t& mat) = 0;

  // Add the adjoint-residual product
  virtual void add_adjoint_res_product(FEVarType wrt, T alpha, Vec_t& data,
                                       Vec_t& geo, Vec_t& sol, Vec_t& adjoint,
                                       Vec_t& dfdx) = 0;

  virtual void to_vtk(Vec_t& data, Vec_t& geo, Vec_t& sol,
                      const std::string filename) {}
};

/**
 * @brief Assemble the contributions from a number of elements
 *
 */
template <class Impl>
class ElementAssembler {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;
  using Mat_t = typename Impl::Mat_t;
  using Elem_t = std::shared_ptr<ElementBase<Impl>>;

  ElementAssembler() {}
  virtual ~ElementAssembler() {}

  // Add an element at the specified index
  void add_element(Elem_t element) { elements.push_back(element); }

  // Get the CSR structure for all the elements in the mesh
  void get_bsr_data(const index_t block_size, index_t& nrows,
                    std::vector<index_t>& rowp, std::vector<index_t>& cols) {
    std::set<std::pair<index_t, index_t>> pairs;
    typename std::list<Elem_t>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      const ElementMeshBase& mesh = (*it)->get_mesh(FEVarType::STATE);
      mesh.add_matrix_pairs(block_size, pairs);
    }
    ElementMeshBase::create_block_csr(pairs, nrows, rowp, cols);
  }

  // Add the residual from all the elements
  void add_residual(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol, Vec_t& res) {
    typename std::list<Elem_t>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_residual(alpha, data, geo, sol, res);
    }
  }

  // Add the Jacobian contributions from all the elements
  void add_jacobian(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol, Mat_t& mat) {
    typename std::list<Elem_t>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_jacobian(alpha, data, geo, sol, mat);
    }
  }

  // Add the adjoint-residual product
  virtual void add_adjoint_res_product(FEVarType wrt, T alpha, Vec_t& data,
                                       Vec_t& geo, Vec_t& sol, Vec_t& adjoint,
                                       Vec_t& dfdx) {
    typename std::list<Elem_t>::iterator it;
    for (it = elements.begin(); it != elements.end(); ++it) {
      (*it)->add_adjoint_res_product(wrt, alpha, data, geo, sol, adjoint, dfdx);
    }
  }

  // Write all the vtk files
  virtual void to_vtk(Vec_t& data, Vec_t& geo, Vec_t& sol,
                      const std::string prefix) {
    typename std::list<Elem_t>::iterator it = elements.begin();
    for (int index = 0; it != elements.end(); ++it, index++) {
      std::string filename = prefix + std::string(".vtk");
      (*it)->to_vtk(data, geo, sol, filename);
    }
  }

 private:
  std::list<Elem_t> elements;
};

/*
  Integrand-type finite-element
*/
template <class Impl, class Integrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class ElementIntegrand : public ElementBase<Impl> {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;
  using Mat_t = typename Impl::Mat_t;

  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  template <class Base>
  using ElementMatrix = typename Impl::template ElementMatrix<Base>;

  virtual const Integrand& get_integrand() = 0;

  // Set the meshes
  void set_meshes(std::shared_ptr<ElementMesh<DataBasis>> data,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo,
                  std::shared_ptr<ElementMesh<Basis>> sol) {
    data_mesh = data;
    geo_mesh = geo;
    sol_mesh = sol;
  }

  const ElementMeshBase& get_mesh(FEVarType wrt) const {
    if (wrt == FEVarType::DATA) {
      return *data_mesh;
    } else if (wrt == FEVarType::GEOMETRY) {
      return *geo_mesh;
    } else {
      return *sol_mesh;
    }
  }

  // Add the residual to the residual vector
  void add_residual(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol, Vec_t& res) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<DataBasis> elem_data(*data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*geo_mesh, geo);
    ElementVector<Basis> elem_sol(*sol_mesh, sol);
    ElementVector<Basis> elem_res(*sol_mesh, res);
    element.template add_residual<FEVarType::STATE>(
        integrand, alpha, elem_data, elem_geo, elem_sol, elem_res);
  }

  // Add the Jacobian
  void add_jacobian(T alpha, Vec_t& data, Vec_t& geo, Vec_t& sol, Mat_t& mat) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<DataBasis> elem_data(*data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*geo_mesh, geo);
    ElementVector<Basis> elem_sol(*sol_mesh, sol);
    ElementMatrix<Basis> elem_jac(*sol_mesh, mat);
    element.template add_jacobian<FEVarType::STATE, FEVarType::STATE>(
        integrand, alpha, elem_data, elem_geo, elem_sol, elem_jac);
  }

  // Add the adjoint-residual product
  void add_adjoint_res_product(FEVarType wrt, T alpha, Vec_t& data, Vec_t& geo,
                               Vec_t& sol, Vec_t& adjoint, Vec_t& dfdx) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<DataBasis> elem_data(*data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*geo_mesh, geo);
    ElementVector<Basis> elem_sol(*sol_mesh, sol);
    ElementVector<Basis> elem_adjoint(*sol_mesh, adjoint);

    if (wrt == FEVarType::DATA) {
      ElementVector<DataBasis> elem_dfdx(*data_mesh, dfdx);
      element.template add_jacobian_product<FEVarType::DATA, FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
          elem_dfdx);

    } else if (wrt == FEVarType::GEOMETRY) {
      ElementVector<GeoBasis> elem_dfdx(*geo_mesh, dfdx);
      element
          .template add_jacobian_product<FEVarType::GEOMETRY, FEVarType::STATE>(
              integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
              elem_dfdx);

    } else if (wrt == FEVarType::STATE) {
      ElementVector<Basis> elem_dfdx(*sol_mesh, dfdx);
      element.template add_jacobian_product<FEVarType::STATE, FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_adjoint,
          elem_dfdx);
    }
  }

 protected:
  FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis> element;
  std::shared_ptr<ElementMesh<DataBasis>> data_mesh;
  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh;
  std::shared_ptr<ElementMesh<Basis>> sol_mesh;
};

/**
 * @brief Functional base class
 *
 * This class defines a functional evaluation that is compatible with the
 * given vector and matrix types
 */
template <class Impl>
class FunctionalBase {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;

  virtual ~FunctionalBase() {}

  // Add the residual
  virtual T evaluate(Vec_t& data, Vec_t& geo, Vec_t& sol) = 0;

  // Add to the derivative type
  virtual void add_derivative(FEVarType wrt, T alpha, Vec_t& data, Vec_t& geo,
                              Vec_t& sol, Vec_t& dfdx) = 0;
};

/*
  Integral type functional
*/
template <class Impl, class Integrand, class Quadrature, class DataBasis,
          class GeoBasis, class Basis>
class IntegralFunctional : public FunctionalBase<Impl> {
 public:
  using T = typename Impl::type;
  using Vec_t = typename Impl::Vec_t;
  using Mat_t = typename Impl::Mat_t;

  template <class Base>
  using ElementVector = typename Impl::template ElementVector<Base>;

  template <class Base>
  using ElementMatrix = typename Impl::template ElementMatrix<Base>;

  virtual const Integrand& get_integrand() = 0;

  // Set the meshes
  void set_meshes(std::shared_ptr<ElementMesh<DataBasis>> data,
                  std::shared_ptr<ElementMesh<GeoBasis>> geo,
                  std::shared_ptr<ElementMesh<Basis>> sol) {
    data_mesh = data;
    geo_mesh = geo;
    sol_mesh = sol;
  }

  // Add the residual to the residual vector
  T evaluate(Vec_t& data, Vec_t& geo, Vec_t& sol) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<DataBasis> elem_data(*data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*geo_mesh, geo);
    ElementVector<Basis> elem_sol(*sol_mesh, sol);
    T value = element.integrate(integrand, elem_data, elem_geo, elem_sol);
    return value;
  }

  // Add the adjoint-residual product
  void add_derivative(FEVarType wrt, T alpha, Vec_t& data, Vec_t& geo,
                      Vec_t& sol, Vec_t& dfdx) {
    const Integrand& integrand = this->get_integrand();
    ElementVector<DataBasis> elem_data(*data_mesh, data);
    ElementVector<GeoBasis> elem_geo(*geo_mesh, geo);
    ElementVector<Basis> elem_sol(*sol_mesh, sol);

    if (wrt == FEVarType::DATA) {
      ElementVector<DataBasis> elem_dfdx(*data_mesh, dfdx);
      element.template add_residual<FEVarType::DATA>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);

    } else if (wrt == FEVarType::GEOMETRY) {
      ElementVector<GeoBasis> elem_dfdx(*geo_mesh, dfdx);
      element.template add_residual<FEVarType::GEOMETRY>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);

    } else if (wrt == FEVarType::STATE) {
      ElementVector<Basis> elem_dfdx(*sol_mesh, dfdx);
      element.template add_residual<FEVarType::STATE>(
          integrand, alpha, elem_data, elem_geo, elem_sol, elem_dfdx);
    }
  }

 protected:
  FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis> element;
  std::shared_ptr<ElementMesh<DataBasis>> data_mesh;
  std::shared_ptr<ElementMesh<GeoBasis>> geo_mesh;
  std::shared_ptr<ElementMesh<Basis>> sol_mesh;
};

}  // namespace A2D

#endif  //  A2D_FE_BASE_H