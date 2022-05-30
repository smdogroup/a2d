#include <pybind11/pybind11.h>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "multiarray.h"

typedef int IndexType;
typedef double ScalarType;

typedef HexBasis<HexQuadrature> Basis;
typedef NonlinearElasticity3D<IndexType, ScalarType, Basis> Elasticity;
typedef HelmholtzPDE<IndexType, ScalarType, Basis> Helmholtz;

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
      })
      .def("zero", &multiarray::zero);
}

PYBIND11_MODULE(example, m) {
  m.doc() = "Wrapping for the a2d model class";

  // Declare the array types
  declare_multiarray<typename Elasticity::ConnArray>(m,
                                                     "Elasticity::ConnArray");
  declare_multiarray<typename Elasticity::SolutionArray>(
      m, "Elasticity::SolutionArray");
  if (!std::is_same<typename Elasticity::SolutionArray,
                    typename Elasticity::NodeArray>::value) {
    declare_multiarray<typename Elasticity::NodeArray>(m,
                                                       "Elasticity::NodeArray");
  }
  declare_multiarray<typename Elasticity::QuadDataArray>(
      m, "Elasticity::QuadDataArray");
  declare_multiarray<typename Elasticity::ElemJacArray>(
      m, "Elasticity::ElemJacArray");

  // Wrap the model function
  py::class_<Elasticity>(m, "Elasticity")
      .def(py::init<const int, const int>())
      .def("get_conn", &Elasticity::get_conn,
           py::return_value_policy::reference)
      .def("get_nodes", &Elasticity::get_nodes,
           py::return_value_policy::reference)
      .def("reset_nodes", &Elasticity::reset_nodes)
      .def("get_solution", &Elasticity::get_solution,
           py::return_value_policy::reference)
      .def("reset_solution", &Elasticity::reset_solution)
      .def("get_quad_data", &Elasticity::get_quad_data,
           py::return_value_policy::reference)
      .def("add_residuals", &Elasticity::add_residuals)
      .def("add_jacobians", &Elasticity::add_jacobians,
           py::return_value_policy::reference)
      .def("get_elem_jac", &Elasticity::get_elem_jac,
           py::return_value_policy::reference);

  // Declare the array types unique for the Helmholtz problem
  declare_multiarray<typename Helmholtz::SolutionArray>(
      m, "Helmholtz::SolutionArray");
  declare_multiarray<typename Helmholtz::QuadDataArray>(
      m, "Helmholtz::QuadDataArray");
  declare_multiarray<typename Helmholtz::ElemJacArray>(
      m, "Helmholtz::ElemJacArray");

  // Wrap the model function
  py::class_<Helmholtz>(m, "Helmholtz")
      .def(py::init<const int, const int>())
      .def("get_conn", &Helmholtz::get_conn, py::return_value_policy::reference)
      .def("get_nodes", &Helmholtz::get_nodes,
           py::return_value_policy::reference)
      .def("reset_nodes", &Helmholtz::reset_nodes)
      .def("get_solution", &Helmholtz::get_solution,
           py::return_value_policy::reference)
      .def("reset_solution", &Helmholtz::reset_solution)
      .def("get_quad_data", &Helmholtz::get_quad_data,
           py::return_value_policy::reference)
      .def("add_residuals", &Helmholtz::add_residuals)
      .def("add_jacobians", &Helmholtz::add_jacobians,
           py::return_value_policy::reference)
      .def("get_elem_jac", &Helmholtz::get_elem_jac,
           py::return_value_policy::reference);
}
