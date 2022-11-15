#include <iostream>
#include <memory>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"

using namespace A2D;

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const index_t degree = 2;
  using T = double;
  using ET = ElementTypes;

  /*
  using PDE = NonlinearElasticity<T, 3>;
  using Quadrature = HexQuadrature<degree + 1>;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 2, 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, 1>>;
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
  const int nx = 10, ny = 10, nz = 10;
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
        for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  double Xloc[3 * nverts];
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (1.0 * i) / nx;
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j) / ny;
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k) / nz;
      }
    }
  }

  int boundary_verts[(ny + 1) * (nz + 1)];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(0, j, k);
    }
  }

  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);
  ElementMesh<DataBasis> datamesh(conn);

  // Set the boundary conditions on the first field
  index_t basis_select[2] = {1, 0};

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  BoundaryCondition<Basis> bcs(conn, mesh, basis_select, (ny + 1) * (nz + 1),
                               boundary_verts);

  std::cout << "Number of elements:            " << conn.get_num_elements()
            << std::endl;
  std::cout << "Number of degrees of freedom:  " << mesh.get_num_dof()
            << std::endl;

  const index_t *bcs_index;
  std::cout << "Number of boundary conditions: " << bcs.get_bcs(&bcs_index)
            << std::endl;

  index_t ndof = mesh.get_num_dof();
  SolutionVector<T> global_U(mesh.get_num_dof());
  SolutionVector<T> global_res(mesh.get_num_dof());
  SolutionVector<T> global_geo(geomesh.get_num_dof());
  SolutionVector<T> global_data(datamesh.get_num_dof());

  FE::DataElemVec elem_data(datamesh, global_data);
  FE::GeoElemVec elem_geo(geomesh, global_geo);
  FE::ElemVec elem_sol(mesh, global_U);
  FE::ElemVec elem_res(mesh, global_res);

  // Set the global geo values
  for (int k = 0, e = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, e++) {
        // Get the geometry values
        typename FE::GeoElemVec::FEDof geo_dof(e, elem_geo);
        int sign[GeoBasis::ndof];

        for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
          index_t node = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                  j + ET::HEX_VERTS_CART[ii][1],
                                  k + ET::HEX_VERTS_CART[ii][2]);

          // Set the entity DOF
          index_t basis = 0;
          index_t orient = 0;
          GeoBasis::set_entity_dof(basis, ET::VERTEX, ii, orient,
                                   &Xloc[3 * node], geo_dof, sign);
        }

        elem_geo.set_element_values(e, geo_dof);
      }
    }
  }

  SolutionVector<T> global_x(mesh.get_num_dof());
  SolutionVector<T> global_y(mesh.get_num_dof());
  FE::ElemVec elem_x(mesh, global_x);
  FE::ElemVec elem_y(mesh, global_y);

  FE fe;

  // Add the residual
  elem_res.init_zero_values();
  fe.add_residual(elem_data, elem_geo, elem_sol, elem_res);
  elem_res.add_values();

  // fe.add_jacobian_vector_product(elem_data, elem_geo, elem_sol, elem_x,
  // elem_y);

  std::cout << "create_block_matrix" << std::endl;
  BSRMat<index_t, T, 3, 3> *mat = mesh.create_block_matrix<T, 3>();
  BSRMat<index_t, T, 3, 3> mat_ref = *mat;

  std::cout << "initialize element matrix" << std::endl;
  ElementMat_Serial<T, Basis, BSRMat<index_t, T, 3, 3>> elem_mat(mesh, *mat);

  std::cout << "add_jacobian" << std::endl;
  fe.add_jacobian(elem_data, elem_geo, elem_sol, elem_mat);

  // SolutionVector<T> *B1[6];

  // for (index_t i = 0; i < 6; i++ ){
  //   fe.get_near_nullspace(0, elem_data, elem_geo, elem_res);
  //   elem_res.

  // }
  // bcs.set_bcs

  int num_levels = 3;
  T omega = 2.0 / 3.0;
  T epsilon = 0.1;

  // amg = new BSRMatAmg(num_levels, omega, epsilon, bsr, B);

  Kokkos::finalize();

  return (0);
}