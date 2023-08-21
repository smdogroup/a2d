#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "fem/functional.h"
#include "fem/helmholtz.h"
#include "fem/integrand_elasticity.h"
#include "fem/model.h"
#include "multiarray.h"

namespace py = pybind11;

using namespace A2D;

// Define spatial dim
static const index_t SPATIAL_DIM = 2;

// Define the types that we'll use
typedef A2D::index_t Itype;
typedef double Ttype;

typedef BasisOps<SPATIAL_DIM, QuadBiLinearBasisFunc, Quad4ptQuadrature>
    Basis_CPS4;

typedef ElementBase<Itype, Ttype, ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Elasticity_Element;
typedef NonlinElasticityElement<Itype, Ttype, Basis_CPS4> Elasticity_CPS4;

typedef FEModel<Itype, Ttype, ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Elasticity_Model;
typedef typename ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>::SparseAmg
    Elasticity_Amg;
typedef typename ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>::SparseMat
    Elasticity_Mat;

typedef Functional<Itype, Ttype, ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Elasticity_Functional;
typedef ElementFunctional<Itype, Ttype,
                          ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Elasticity_ElementFunctional;

typedef IntegrandTopoVolume<Itype, Ttype, Basis_CPS4> TopoVolume_CPS4;

typedef IntegrandTopoVonMisesKS<Itype, Ttype, Basis_CPS4>
    TopoVonMisesAggregation_CPS4;

typedef ConstitutiveBase<Itype, Ttype,
                         ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Elasticity_Constitutive;
typedef TopoIsoConstitutive<Itype, Ttype, Basis_CPS4> TopoIsoConstitutive_CPS4;

typedef ElementBase<Itype, Ttype, HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Helmholtz_Element;
typedef HelmholtzElement<Itype, Ttype, Basis_CPS4> Helmholtz_CPS4;

typedef FEModel<Itype, Ttype, HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Helmholtz_Model;
typedef typename HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>::SparseAmg
    Helmholtz_Amg;
typedef typename HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>::SparseMat
    Helmholtz_Mat;

typedef ConstitutiveBase<Itype, Ttype,
                         HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>>
    Helmholtz_Constitutive;
typedef HelmholtzConstitutive<Itype, Ttype, Basis_CPS4>
    HelmholtzConstitutive_CPS4;

template <class multiarray>
void declare_array(py::module& m, const char typestr[]) {
  py::class_<multiarray, std::shared_ptr<multiarray>>(m, typestr,
                                                      py::buffer_protocol())
      .def_buffer([](multiarray& array) -> py::buffer_info {
        std::size_t ndims = array.get_rank();
        std::vector<std::size_t> shape(ndims);
        std::vector<std::size_t> strides(ndims);

        for (std::size_t i = 0; i < ndims; i++) {
          shape[i] = array.extent(i);
          strides[i] = sizeof(typename multiarray::type);
          for (std::size_t j = 0; j < i; j++) {
            strides[j] *= array.extent(i);
          }
        }

        return py::buffer_info(
            array.data, sizeof(typename multiarray::type),
            py::format_descriptor<typename multiarray::type>::format(), ndims,
            shape, strides);
      })
      .def("zero", &multiarray::zero)
      .def("fill", &multiarray::fill)
      .def("scale", &multiarray::scale)
      .def("random", &multiarray::random)
      .def("norm", &multiarray::norm);
}

template <class element, class base, class... args>
void declare_element(py::module& m, const char typestr[]) {
  // Wrap the elasticity elements
  py::class_<element, base, std::shared_ptr<element>>(m, typestr)
      .def(py::init(
          [=](py::array_t<int, py::array::c_style> conn, args... vals) {
            py::buffer_info buf = conn.request();

            if (buf.ndim != 2) {
              throw std::runtime_error("Connectivity dimension must be two");
            }
            if (buf.shape[1] != element::nodes_per_elem) {
              throw std::runtime_error(
                  "Second connectivity dimension must match "
                  "number of nodes per element");
            }

            index_t nelems = buf.shape[0];
            int* ptr = static_cast<int*>(buf.ptr);
            return new element(nelems, ptr, vals...);
          }));
}

template <class model>
void declare_2dmodel(py::module& m, const char typestr[]) {
  // Wrap the model
  py::class_<model, std::shared_ptr<model>>(m, typestr)
      .def(py::init([](py::array_t<Ttype, py::array::c_style> X) {
        py::buffer_info Xbuf = X.request();

        // Check the buffer shapes
        if (Xbuf.ndim != 2) {
          throw std::runtime_error("Model: Node array dimension must be two");
        }
        index_t nnodes = Xbuf.shape[0];
        if (nnodes > 0 && Xbuf.shape[1] != SPATIAL_DIM) {
          throw std::runtime_error(
              "Model: There must be 2 coordinates per node");
        }

        Ttype* Xptr = static_cast<Ttype*>(Xbuf.ptr);
        int* bcs = NULL;
        return new model(nnodes, Xptr, 0, bcs);
      }))
      .def(py::init([](py::array_t<Ttype, py::array::c_style> X,
                       py::array_t<int, py::array::c_style> bcs) {
        py::buffer_info Xbuf = X.request();
        py::buffer_info bcsbuf = bcs.request();

        // Check the buffer shapes
        if (Xbuf.ndim != 2) {
          throw std::runtime_error("Model: Node array dimension must be two");
        }
        index_t nnodes = Xbuf.shape[0];
        if (nnodes > 0 && Xbuf.shape[1] != SPATIAL_DIM) {
          throw std::runtime_error(
              "Model: There must be 2 coordinates per node");
        }
        if (bcsbuf.ndim != 2) {
          throw std::runtime_error(
              "Model: The bcs array must have two dimensions");
        }
        index_t nbcs = bcsbuf.shape[0];
        if (nbcs > 0 && bcsbuf.shape[1] != 2) {
          throw std::runtime_error(
              "Model: The bcs array must have two values index and fixed "
              "dof");
        }

        Ttype* Xptr = static_cast<Ttype*>(Xbuf.ptr);
        int* bcsptr = static_cast<int*>(bcsbuf.ptr);
        return new model(nnodes, Xptr, nbcs, bcsptr);
      }))
      .def("add_constitutive", &model::add_constitutive)
      .def("add_element", &model::add_element)
      .def("new_solution", &model::new_solution)
      .def("init", &model::init)
      .def("set_nodes", &model::set_nodes)
      .def("set_solution", &model::set_solution)
      .def("residual", &model::residual)
      .def("jacobian", &model::jacobian)
      .def("set_design_vars", &model::set_design_vars)
      .def("add_adjoint_dfdx", &model::add_adjoint_dfdx)
      .def("new_matrix", &model::new_matrix)
      .def("new_amg", &model::new_amg, py::arg("num_levels"), py::arg("omega"),
           py::arg("epsilon"), py::arg("mat"), py::arg("print_info") = false);
}

template <class amg>
void declare_amg(py::module& m, const char typestr[]) {
  // Wrap the Amg object
  py::class_<amg, std::shared_ptr<amg>>(m, typestr)
      .def("update", &amg::update)
      .def("mg", &amg::mg, "Multigrid method", py::arg("b"), py::arg("x"),
           py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30)
      .def("cg", &amg::cg, "Conjugate gradient method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30,
           py::arg("iters_per_reset") = 100);
}

PYBIND11_MODULE(pya2d_2d, m) {
  m.doc() = "Wrapping for the a2d model class";

  // Elasticity ----------------------------------------------------------

  // Declare the array types
  declare_array<
      typename ElasticityPDEInfo<SPATIAL_DIM, Itype, Ttype>::SolutionArray>(
      m, "Elasticity_Model::SolutionArray");

  declare_2dmodel<Elasticity_Model>(m, "Elasticity_Model");

  // Virtual base class
  py::class_<Elasticity_Element, std::shared_ptr<Elasticity_Element>>(
      m, "Elasticity_Element");

  // Declare the base class
  declare_element<Elasticity_CPS4, Elasticity_Element>(m, "Elasticity_CPS4");

  // Wrap the Matrix object
  py::class_<Elasticity_Mat, std::shared_ptr<Elasticity_Mat>>(m,
                                                              "Elasticity_Mat");

  // Declare the associated AMG type
  declare_amg<Elasticity_Amg>(m, "Elasticity_Amg");

  // Declare the constitutive classes
  py::class_<Elasticity_Constitutive, std::shared_ptr<Elasticity_Constitutive>>(
      m, "Elasticity_Constitutive");
  py::class_<TopoIsoConstitutive_CPS4, Elasticity_Constitutive,
             std::shared_ptr<TopoIsoConstitutive_CPS4>>(
      m, "TopoIsoConstitutive_CPS4")
      .def(py::init<std::shared_ptr<Elasticity_CPS4>, double, double, double,
                    double, double>());

  py::class_<Elasticity_Functional>(m, "Elasticity_Functional")
      .def(py::init<>())
      .def("add_functional", &Elasticity_Functional::add_functional)
      .def("eval_functional", &Elasticity_Functional::eval_functional)
      .def("eval_dfdu", &Elasticity_Functional::eval_dfdu)
      .def("eval_dfdx", &Elasticity_Functional::eval_dfdx)
      .def("eval_dfdnodes", &Elasticity_Functional::eval_dfdnodes);

  py::class_<Elasticity_ElementFunctional,
             std::shared_ptr<Elasticity_ElementFunctional>>(
      m, "Elasticity_ElementFunctional");

  py::class_<TopoVolume_CPS4, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVolume_CPS4>>(m, "TopoVolume_CPS4")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_CPS4>>());

  py::class_<TopoVonMisesAggregation_CPS4, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVonMisesAggregation_CPS4>>(
      m, "TopoVonMisesAggregation_CPS4")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_CPS4>, double>());

  // Helmholtz ----------------------------------------------------------

  // Declare the array types
  declare_array<
      typename HelmholtzPDEInfo<SPATIAL_DIM, Itype, Ttype>::SolutionArray>(
      m, "Helmholtz_Model::SolutionArray");

  declare_2dmodel<Helmholtz_Model>(m, "Helmholtz_Model");

  // Virtual base class
  py::class_<Helmholtz_Element, std::shared_ptr<Helmholtz_Element>>(
      m, "Helmholtz_Element");

  // Declare the base class
  declare_element<Helmholtz_CPS4, Helmholtz_Element, double>(m,
                                                             "Helmholtz_CPS4");

  // Declare the constitutive classes
  py::class_<Helmholtz_Constitutive, std::shared_ptr<Helmholtz_Constitutive>>(
      m, "Helmholtz_Constitutive");
  py::class_<HelmholtzConstitutive_CPS4, Helmholtz_Constitutive,
             std::shared_ptr<HelmholtzConstitutive_CPS4>>(
      m, "HelmholtzConstitutive_CPS4")
      .def(py::init<std::shared_ptr<Helmholtz_CPS4>>());

  // Wrap the Matrix object
  py::class_<Helmholtz_Mat, std::shared_ptr<Helmholtz_Mat>>(m, "Helmholtz_Mat");

  // Declare the associated AMG type
  declare_amg<Helmholtz_Amg>(m, "Helmholtz_Amg");
}
