#include <pybind11/pybind11.h>

#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "model.h"
#include "multiarray.h"

namespace py = pybind11;

using namespace A2D;

// Define the types that we'll use
typedef A2D::index_t Itype;
typedef double Ttype;

typedef Basis3D<HexTriLinear, Hex8ptQuadrature> HexBasis;
typedef Basis3D<TetraQuadraticBasis, Tetra5ptQuadrature> TetBasis;

typedef FEModel<Itype, Ttype, HelmholtzPDE<Itype, Ttype>> HelmholtzModel;
typedef FEModel<Itype, Ttype, ElasticityPDE<Itype, Ttype>> ElasticityModel;

typedef Element<Itype, Ttype, ElasticityPDE<Itype, Ttype>>
    ElasticityElementBase;

typedef LinElasticityElement3D<Itype, Ttype, HexBasis> ElasticityHexElement;
typedef LinElasticityElement3D<Itype, Ttype, TetBasis> ElasticityTetElement;

typedef typename ElasticityPDE<Itype, Ttype>::SparseAmg ElasticityAmg;
typedef typename ElasticityPDE<Itype, Ttype>::SparseMat ElasticityMat;

// typedef typename HelmholtzPDE<Itype, Ttype>::SparseAmg HelmholtzAmg;
// typedef typename HelmholtzPDE<Itype, Ttype>::SparseMat HelmholtzMat;

template <class multiarray>
void declare_array(py::module& m, const char typestr[]) {
  py::class_<multiarray>(m, typestr, py::buffer_protocol())
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
      .def("fill", &multiarray::fill);
}

PYBIND11_MODULE(example, m) {
  m.doc() = "Wrapping for the a2d model class";

  // Declare the array types
  declare_array<typename ElasticityPDE<Itype, Ttype>::SolutionArray>(
      m, "ElasticityModel::SolutionArray");
  declare_array<typename ElasticityPDE<Itype, Ttype>::BCsArray>(
      m, "ElasticityModel::BCsArray");

  declare_array<typename ElasticityHexElement::ConnArray>(
      m, "ElasticityHexElement::ConnArray");
  declare_array<typename ElasticityHexElement::QuadDataArray>(
      m, "ElasticityHexElement::QuadDataArray");

  declare_array<typename ElasticityTetElement::ConnArray>(
      m, "ElasticityTetElement::ConnArray");
  declare_array<typename ElasticityTetElement::QuadDataArray>(
      m, "ElasticityTetElement::QuadDataArray");

  // Virtual base class
  py::class_<ElasticityElementBase>(m, "ElasticityElementBase");

  // Wrap the elasticity elements
  py::class_<ElasticityHexElement, ElasticityElementBase>(
      m, "ElasticityHexElement")
      .def(py::init<const index_t>())
      .def("get_conn", &ElasticityHexElement::get_conn)
      .def("get_quad_data", &ElasticityHexElement::get_quad_data);

  py::class_<ElasticityTetElement, ElasticityElementBase>(
      m, "ElasticityTetElement")
      .def(py::init<const index_t>())
      .def("get_conn", &ElasticityTetElement::get_conn)
      .def("get_quad_data", &ElasticityTetElement::get_quad_data);

  // Wrap the model
  py::class_<ElasticityModel>(m, "ElasticityModel")
      .def(py::init<const index_t, const index_t>())
      .def("add_element", &ElasticityModel::add_element)
      .def("new_solution", &ElasticityModel::new_solution)
      .def("get_nodes", &ElasticityModel::get_nodes,
           py::return_value_policy::reference)
      .def("get_bcs", &ElasticityModel::get_bcs,
           py::return_value_policy::reference)
      .def("set_nodes", &ElasticityModel::set_nodes)
      .def("set_solution", &ElasticityModel::set_solution)
      .def("residual", &ElasticityModel::residual)
      .def("jacobian", &ElasticityModel::jacobian)
      .def("new_matrix", &ElasticityModel::new_matrix)
      .def("new_amg", &ElasticityModel::new_amg, py::arg("num_levels"),
           py::arg("omega"), py::arg("mat"), py::arg("print_info") = false);

  // Wrap the Matrix object
  py::class_<ElasticityMat>(m, "ElasticityMat");

  // Wrap the Amg object
  py::class_<ElasticityAmg>(m, "ElasticityAmg")
      .def("mg", &ElasticityAmg::mg, "Multigrid method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30)
      .def("cg", &ElasticityAmg::cg, "Conjugate gradient method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30,
           py::arg("iters_per_reset") = 100);

  //   // Declare the array types unique for the Helmholtz problem
  //   declare_array<typename Helmholtz::SolutionArray>(m,
  //                                                    "Helmholtz::SolutionArray");
  //   declare_array<typename Helmholtz::QuadDataArray>(m,
  //                                                    "Helmholtz::QuadDataArray");

  //   // Wrap the model function
  //   py::class_<Helmholtz>(m, "Helmholtz")
  //       .def(py::init<const int, const int, const int>())
  //       .def("get_conn", &Helmholtz::get_conn,
  //       py::return_value_policy::reference) .def("get_bcs",
  //       &Helmholtz::get_bcs, py::return_value_policy::reference)
  //       .def("get_nodes", &Helmholtz::get_nodes,
  //            py::return_value_policy::reference)
  //       .def("reset_nodes", &Helmholtz::reset_nodes)
  //       .def("get_solution", &Helmholtz::get_solution,
  //            py::return_value_policy::reference)
  //       .def("reset_solution", &Helmholtz::reset_solution)
  //       .def("get_quad_data", &Helmholtz::get_quad_data,
  //            py::return_value_policy::reference)
  //       .def("add_residual", &Helmholtz::add_residual)
  //       .def("add_jacobian", &Helmholtz::add_jacobian)
  //       .def("new_solution", &Helmholtz::new_solution)
  //       .def("new_matrix", &Helmholtz::new_matrix,
  //            py::return_value_policy::reference)
  //       .def("new_amg", &Helmholtz::new_amg,
  //       py::return_value_policy::reference);

  //   py::class_<HelmholtzMat>(m, "HelmholtzMat").def("zero",
  //   &HelmholtzMat::zero);

  //   py::class_<HelmholtzAmg>(m, "HelmholtzAmg")
  //       .def("mg", &HelmholtzAmg::mg, "Multigrid method", py::arg("b"),
  //            py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") =
  //            500, py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30)
  //       .def("cg", &HelmholtzAmg::cg, "Conjugate gradient method",
  //       py::arg("b"),
  //            py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") =
  //            500, py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30,
  //            py::arg("iters_per_reset") = 100);
}
