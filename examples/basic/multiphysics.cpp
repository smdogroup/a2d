#include <iostream>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"

using namespace A2D;

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  const index_t degree = 2;
  using T = double;

  /*
  using PDE = NonlinearElasticity<T, 3>;
  using Quadrature = HexQuadrature<degree + 1>;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 2, 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  */

  using PDE = MixedPoisson<T, 3>;
  using Quadrature = HexQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;
  using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
                        LagrangeL2HexBasis<T, 1, degree - 1>>;

  constexpr bool use_parallel_elemvec = false;
  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis,
                           use_parallel_elemvec>;

  // Number of elements in each dimension
  const int nx = 25, ny = 25, nz = 25;
  auto node_num = [](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  // Number of edges
  const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = nx * ny * nz;

  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8 * nhex];

  for (int k = 0, e = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, e++) {
        hex[8 * e + 0] = node_num(i, j, k);
        hex[8 * e + 1] = node_num(i + 1, j, k);
        hex[8 * e + 2] = node_num(i + 1, j + 1, k);
        hex[8 * e + 3] = node_num(i, j + 1, k);
        hex[8 * e + 4] = node_num(i, j, k + 1);
        hex[8 * e + 5] = node_num(i + 1, j, k + 1);
        hex[8 * e + 6] = node_num(i + 1, j + 1, k + 1);
        hex[8 * e + 7] = node_num(i, j + 1, k + 1);
      }
    }
  }

  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  index_t nelems = nx * ny * nz;
  index_t nfaces = nx * ny * (nz + 1) + nx * (ny + 1) * nz + (nx + 1) * ny * nz;
  index_t nedges = nx * (ny + 1) * (nz + 1) + (nx + 1) * ny * (nz + 1) +
                   (nx + 1) * (ny + 1) * nz;

  std::cout << "Number of elements: " << nelems << " "
            << conn.get_num_elements() << std::endl;
  std::cout << "Number of faces:    " << nfaces << " " << conn.get_num_faces()
            << std::endl;
  std::cout << "Number of edges:    " << nedges << " " << conn.get_num_edges()
            << std::endl;

  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);
  ElementMesh<DataBasis> datamesh(conn);

  index_t ndof = mesh.get_num_dof();
  SolutionVector<T> global_U(mesh.get_num_dof());
  SolutionVector<T> global_res(mesh.get_num_dof());
  SolutionVector<T> global_geo(geomesh.get_num_dof());
  SolutionVector<T> global_data(datamesh.get_num_dof());

  FE::DataElemVec elem_data(datamesh, global_data);
  FE::GeoElemVec elem_geo(geomesh, global_geo);
  FE::ElemVec elem_sol(mesh, global_U);
  FE::ElemVec elem_res(mesh, global_res);

  SolutionVector<T> global_x(mesh.get_num_dof());
  SolutionVector<T> global_y(mesh.get_num_dof());
  FE::ElemVec elem_x(mesh, global_x);
  FE::ElemVec elem_y(mesh, global_y);

  FE fe;

  // Add the residual
  elem_res.init_zero_values();
  fe.add_residual(elem_data, elem_geo, elem_sol, elem_res);
  elem_res.add_values();

  fe.add_jacobian_vector_product(elem_data, elem_geo, elem_sol, elem_x, elem_y);

  std::cout << "create_block_matrix" << std::endl;
  BSRMat<index_t, T, 3, 3>* mat = mesh.create_block_matrix<T, 3>();
  BSRMat<index_t, T, 3, 3> mat_ref = *mat;

  std::cout << "initialize element matrix" << std::endl;
  ElementMat_Serial<T, Basis, BSRMat<index_t, T, 3, 3>> elem_mat(mesh, mat_ref);

  std::cout << "add_jacobian" << std::endl;
  fe.add_jacobian(elem_data, elem_geo, elem_sol, elem_mat);

  Kokkos::finalize();

  return (0);
}