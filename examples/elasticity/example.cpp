#include <pybind11/pybind11.h>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "multiarray.h"

typedef A2D::index_t IndexType;
typedef double ScalarType;

typedef HexBasis<HexQuadrature> Basis;
typedef A2D::LinearElasticity3D<IndexType, ScalarType, Basis> Elasticity;
typedef typename Elasticity::base::SparseAmg ElasticityAmg;
typedef typename Elasticity::base::SparseMat ElasticityMat;

typedef A2D::HelmholtzPDE<IndexType, ScalarType, Basis> Helmholtz;
typedef typename Helmholtz::base::SparseMat HelmholtzMat;
typedef typename Helmholtz::base::SparseAmg HelmholtzAmg;

namespace py = pybind11;

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
  declare_array<typename Elasticity::ConnArray>(m, "Elasticity::ConnArray");
  declare_array<typename Elasticity::BCsArray>(m, "Elasticity::BCsArray");
  declare_array<typename Elasticity::SolutionArray>(
      m, "Elasticity::SolutionArray");
  if (!std::is_same<typename Elasticity::SolutionArray,
                    typename Elasticity::NodeArray>::value) {
    declare_array<typename Elasticity::NodeArray>(m, "Elasticity::NodeArray");
  }
  declare_array<typename Elasticity::QuadDataArray>(
      m, "Elasticity::QuadDataArray");

  // Wrap the model function
  py::class_<Elasticity>(m, "Elasticity")
      .def(py::init<const int, const int, const int>())
      .def("get_conn", &Elasticity::get_conn,
           py::return_value_policy::reference)
      .def("get_bcs", &Elasticity::get_bcs, py::return_value_policy::reference)
      .def("get_nodes", &Elasticity::get_nodes,
           py::return_value_policy::reference)
      .def("reset_nodes", &Elasticity::reset_nodes)
      .def("get_solution", &Elasticity::get_solution,
           py::return_value_policy::reference)
      .def("reset_solution", &Elasticity::reset_solution)
      .def("get_quad_data", &Elasticity::get_quad_data,
           py::return_value_policy::reference)
      .def("add_residual", &Elasticity::add_residual)
      .def("add_jacobian", &Elasticity::add_jacobian)
      .def("new_solution", &Elasticity::new_solution)
      .def("new_matrix", &Elasticity::new_matrix,
           py::return_value_policy::reference)
      .def("new_amg", &Elasticity::new_amg, py::return_value_policy::reference);

  py::class_<ElasticityMat>(m, "ElasticityMat")
      .def("zero", &ElasticityMat::zero);

  py::class_<ElasticityAmg>(m, "ElasticityAmg")
      .def("mg", &ElasticityAmg::mg, "Multigrid method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30)
      .def("cg", &ElasticityAmg::cg, "Conjugate gradient method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30,
           py::arg("iters_per_reset") = 100);

  // Declare the array types unique for the Helmholtz problem
  declare_array<typename Helmholtz::SolutionArray>(m,
                                                   "Helmholtz::SolutionArray");
  declare_array<typename Helmholtz::QuadDataArray>(m,
                                                   "Helmholtz::QuadDataArray");

  // Wrap the model function
  py::class_<Helmholtz>(m, "Helmholtz")
      .def(py::init<const int, const int, const int>())
      .def("get_conn", &Helmholtz::get_conn, py::return_value_policy::reference)
      .def("get_bcs", &Helmholtz::get_bcs, py::return_value_policy::reference)
      .def("get_nodes", &Helmholtz::get_nodes,
           py::return_value_policy::reference)
      .def("reset_nodes", &Helmholtz::reset_nodes)
      .def("get_solution", &Helmholtz::get_solution,
           py::return_value_policy::reference)
      .def("reset_solution", &Helmholtz::reset_solution)
      .def("get_quad_data", &Helmholtz::get_quad_data,
           py::return_value_policy::reference)
      .def("add_residual", &Helmholtz::add_residual)
      .def("add_jacobian", &Helmholtz::add_jacobian)
      .def("new_solution", &Helmholtz::new_solution)
      .def("new_matrix", &Helmholtz::new_matrix,
           py::return_value_policy::reference)
      .def("new_amg", &Helmholtz::new_amg, py::return_value_policy::reference);

  py::class_<HelmholtzMat>(m, "HelmholtzMat").def("zero", &HelmholtzMat::zero);

  py::class_<HelmholtzAmg>(m, "HelmholtzAmg")
      .def("mg", &HelmholtzAmg::mg, "Multigrid method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30)
      .def("cg", &HelmholtzAmg::cg, "Conjugate gradient method", py::arg("b"),
           py::arg("x"), py::arg("monitor") = 0, py::arg("max_iters") = 500,
           py::arg("rtol") = 1e-8, py::arg("atol") = 1e-30,
           py::arg("iters_per_reset") = 100);
}
