#include <pybind11/pybind11.h>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "multiarray.h"

typedef int IndexType;
typedef double ScalarType;
typedef HexQuadrature Quadrature;
typedef HexBasis<HexQuadrature> Basis;
typedef NonlinearElasticity3D<IndexType, ScalarType, Basis> Model;

namespace py = pybind11;

template <class multiarray>
void declare_multiarray(py::module& m, const char typestr[]) {
  py::class_<multiarray>(m, typestr, py::buffer_protocol())
      .def_buffer([](multiarray& array) -> py::buffer_info {
        std::size_t ndims = array.rank();
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
      });
}

PYBIND11_MODULE(example, m) {
  m.doc() = "Wrapping for the a2d model class";

  // Declare the array types
  declare_multiarray<typename Model::ConnArray>(m, "ConnArray");
  declare_multiarray<typename Model::SolutionArray>(m, "SolutionArray");
  if (!std::is_same<typename Model::SolutionArray,
                    typename Model::NodeArray>::value) {
    declare_multiarray<typename Model::NodeArray>(m, "NodeArray");
  }
  declare_multiarray<typename Model::QuadDataArray>(m, "QuadDataArray");
  declare_multiarray<typename Model::ElemJacArray>(m, "ElemJacArray");

  // Wrap the model function
  py::class_<Model>(m, "Model")
      .def(py::init<const int, const int>())
      .def("get_conn", &Model::get_conn, py::return_value_policy::reference)
      .def("get_nodes", &Model::get_nodes, py::return_value_policy::reference)
      .def("reset_nodes", &Model::reset_nodes)
      .def("get_solution", &Model::get_solution,
           py::return_value_policy::reference)
      .def("reset_solution", &Model::reset_solution)
      .def("get_quad_data", &Model::get_quad_data,
           py::return_value_policy::reference)
      .def("add_residuals", &Model::add_residuals)
      .def("add_jacobians", &Model::add_jacobians,
           py::return_value_policy::reference)
      .def("get_elem_jac", &Model::get_quad_data,
           py::return_value_policy::reference);
}
