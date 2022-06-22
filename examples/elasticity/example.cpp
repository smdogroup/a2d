#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "elasticity3d.h"
#include "functional.h"
#include "helmholtz3d.h"
#include "model.h"
#include "multiarray.h"

namespace py = pybind11;

using namespace A2D;

// Define the types that we'll use
typedef A2D::index_t Itype;
typedef double Ttype;

typedef Basis3D<HexTriLinear, Hex8ptQuadrature> Basis_C3D8;
typedef Basis3D<TetraQuadraticBasis, Tetra5ptQuadrature> Basis_C3D10;

typedef Element<Itype, Ttype, ElasticityPDE<Itype, Ttype>> Elasticity_Element;
typedef LinElasticityElement3D<Itype, Ttype, Basis_C3D8> Elasticity_C3D8;
typedef LinElasticityElement3D<Itype, Ttype, Basis_C3D10> Elasticity_C3D10;

typedef FEModel<Itype, Ttype, ElasticityPDE<Itype, Ttype>> Elasticity_Model;
typedef typename ElasticityPDE<Itype, Ttype>::SparseAmg Elasticity_Amg;
typedef typename ElasticityPDE<Itype, Ttype>::SparseMat Elasticity_Mat;

typedef Functional<Itype, Ttype, ElasticityPDE<Itype, Ttype>>
    Elasticity_Functional;
typedef ElementFunctional<Itype, Ttype, ElasticityPDE<Itype, Ttype>>
    Elasticity_ElementFunctional;

typedef TopoVolume<Itype, Ttype, Basis_C3D8> TopoVolume_C3D8;
typedef TopoVolume<Itype, Ttype, Basis_C3D10> TopoVolume_C3D10;

typedef TopoVonMisesAggregation<Itype, Ttype, Basis_C3D8>
    TopoVonMisesAggregation_C3D8;
typedef TopoVonMisesAggregation<Itype, Ttype, Basis_C3D10>
    TopoVonMisesAggregation_C3D10;

typedef Constitutive<Itype, Ttype, ElasticityPDE<Itype, Ttype>>
    Elasticity_Constitutive;
typedef TopoIsoConstitutive<Itype, Ttype, Basis_C3D8> TopoIsoConstitutive_C3D8;
typedef TopoIsoConstitutive<Itype, Ttype, Basis_C3D10>
    TopoIsoConstitutive_C3D10;

typedef Element<Itype, Ttype, HelmholtzPDE<Itype, Ttype>> Helmholtz_Element;
typedef HelmholtzElement3D<Itype, Ttype, Basis_C3D8> Helmholtz_C3D8;
typedef HelmholtzElement3D<Itype, Ttype, Basis_C3D10> Helmholtz_C3D10;

typedef FEModel<Itype, Ttype, HelmholtzPDE<Itype, Ttype>> Helmholtz_Model;
typedef typename HelmholtzPDE<Itype, Ttype>::SparseAmg Helmholtz_Amg;
typedef typename HelmholtzPDE<Itype, Ttype>::SparseMat Helmholtz_Mat;

typedef Constitutive<Itype, Ttype, HelmholtzPDE<Itype, Ttype>>
    Helmholtz_Constitutive;
typedef HelmholtzConstitutive<Itype, Ttype, Basis_C3D8>
    HelmholtzConstitutive_C3D8;
typedef HelmholtzConstitutive<Itype, Ttype, Basis_C3D10>
    HelmholtzConstitutive_C3D10;

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
void declare_3dmodel(py::module& m, const char typestr[]) {
  // Wrap the model
  py::class_<model, std::shared_ptr<model>>(m, typestr)
      .def(py::init([](py::array_t<Ttype, py::array::c_style> X) {
        py::buffer_info Xbuf = X.request();

        // Check the buffer shapes
        if (Xbuf.ndim != 2) {
          throw std::runtime_error("Model: Node array dimension must be two");
        }
        index_t nnodes = Xbuf.shape[0];
        if (nnodes > 0 && Xbuf.shape[1] != 3) {
          throw std::runtime_error(
              "Model: There must be 3 coordinates per node");
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
        if (nnodes > 0 && Xbuf.shape[1] != 3) {
          throw std::runtime_error(
              "Model: There must be 3 coordinates per node");
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
           py::arg("mat"), py::arg("print_info") = false);
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

PYBIND11_MODULE(example, m) {
  m.doc() = "Wrapping for the a2d model class";

  // Elasticity ----------------------------------------------------------

  // Declare the array types
  declare_array<typename ElasticityPDE<Itype, Ttype>::SolutionArray>(
      m, "Elasticity_Model::SolutionArray");

  declare_3dmodel<Elasticity_Model>(m, "Elasticity_Model");

  // Virtual base class
  py::class_<Elasticity_Element, std::shared_ptr<Elasticity_Element>>(
      m, "Elasticity_Element");

  // Declare the base class
  declare_element<Elasticity_C3D8, Elasticity_Element>(m, "Elasticity_C3D8");
  declare_element<Elasticity_C3D10, Elasticity_Element>(m, "Elasticity_C3D10");

  // Wrap the Matrix object
  py::class_<Elasticity_Mat, std::shared_ptr<Elasticity_Mat>>(m,
                                                              "Elasticity_Mat");

  // Declare the associated AMG type
  declare_amg<Elasticity_Amg>(m, "Elasticity_Amg");

  // Declare the constitutive classes
  py::class_<Elasticity_Constitutive, std::shared_ptr<Elasticity_Constitutive>>(
      m, "Elasticity_Constitutive");
  py::class_<TopoIsoConstitutive_C3D8, Elasticity_Constitutive,
             std::shared_ptr<TopoIsoConstitutive_C3D8>>(
      m, "TopoIsoConstitutive_C3D8")
      .def(py::init<std::shared_ptr<Elasticity_C3D8>, double, double, double,
                    double, double>());
  py::class_<TopoIsoConstitutive_C3D10, Elasticity_Constitutive,
             std::shared_ptr<TopoIsoConstitutive_C3D10>>(
      m, "TopoIsoConstitutive_C3D10")
      .def(py::init<std::shared_ptr<Elasticity_C3D10>, double, double, double,
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

  py::class_<TopoVolume_C3D8, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVolume_C3D8>>(m, "TopoVolume_C3D8")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_C3D8>>());
  py::class_<TopoVolume_C3D10, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVolume_C3D10>>(m, "TopoVolume_C3D10")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_C3D10>>());

  py::class_<TopoVonMisesAggregation_C3D8, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVonMisesAggregation_C3D8>>(
      m, "TopoVonMisesAggregation_C3D8")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_C3D8>, double>());
  py::class_<TopoVonMisesAggregation_C3D10, Elasticity_ElementFunctional,
             std::shared_ptr<TopoVonMisesAggregation_C3D10>>(
      m, "TopoVonMisesAggregation_C3D10")
      .def(py::init<std::shared_ptr<TopoIsoConstitutive_C3D10>, double>());

  // Helmholtz ----------------------------------------------------------

  // Declare the array types
  declare_array<typename HelmholtzPDE<Itype, Ttype>::SolutionArray>(
      m, "Helmholtz_Model::SolutionArray");

  declare_3dmodel<Helmholtz_Model>(m, "Helmholtz_Model");

  // Virtual base class
  py::class_<Helmholtz_Element, std::shared_ptr<Helmholtz_Element>>(
      m, "Helmholtz_Element");

  // Declare the base class
  declare_element<Helmholtz_C3D8, Helmholtz_Element, double>(m,
                                                             "Helmholtz_C3D8");
  declare_element<Helmholtz_C3D10, Helmholtz_Element, double>(
      m, "Helmholtz_C3D10");

  // Declare the constitutive classes
  py::class_<Helmholtz_Constitutive, std::shared_ptr<Helmholtz_Constitutive>>(
      m, "Helmholtz_Constitutive");
  py::class_<HelmholtzConstitutive_C3D8, Helmholtz_Constitutive,
             std::shared_ptr<HelmholtzConstitutive_C3D8>>(
      m, "HelmholtzConstitutive_C3D8")
      .def(py::init<std::shared_ptr<Helmholtz_C3D8>>());
  py::class_<HelmholtzConstitutive_C3D10, Helmholtz_Constitutive,
             std::shared_ptr<HelmholtzConstitutive_C3D10>>(
      m, "HelmholtzConstitutive_C3D10")
      .def(py::init<std::shared_ptr<Helmholtz_C3D10>>());

  // Wrap the Matrix object
  py::class_<Helmholtz_Mat, std::shared_ptr<Helmholtz_Mat>>(m, "Helmholtz_Mat");

  // Declare the associated AMG type
  declare_amg<Helmholtz_Amg>(m, "Helmholtz_Amg");
}
