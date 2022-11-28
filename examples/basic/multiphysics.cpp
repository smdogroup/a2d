#include <iostream>
#include <memory>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"

using namespace A2D;

//   static void compute_null_space(NodeArray& X, NullSpaceArray& B) {
//     A2D::BLAS::zero(B);
//     if (SPATIAL_DIM == 3) {
//       for (I i = 0; i < B.extent(0); i++) {
//         B(i, 0, 0) = 1.0;
//         B(i, 1, 1) = 1.0;
//         B(i, 2, 2) = 1.0;

//         // Rotation about the x-axis
//         B(i, 1, 3) = X(i, 2);
//         B(i, 2, 3) = -X(i, 1);

//         // Rotation about the y-axis
//         B(i, 0, 4) = X(i, 2);
//         B(i, 2, 4) = -X(i, 0);

//         // Rotation about the z-axis
//         B(i, 0, 5) = X(i, 1);
//         B(i, 1, 5) = -X(i, 0);
//       }
//     } else {
//       for (I i = 0; i < B.extent(0); i++) {
//         B(i, 0, 0) = 1.0;
//         B(i, 1, 1) = 1.0;

//         // Rotation about the z-axis
//         B(i, 0, 2) = X(i, 1);
//         B(i, 1, 2) = -X(i, 0);
//       }
//     }
//   }

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const index_t dim = 3;
  using T = double;
  using ET = ElementTypes;
  // using PDE = MixedPoisson<T, dim>;
  using PDE = Poisson<T, dim>;

  const index_t degree = 8;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, degree>>;
  // using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
  //                       LagrangeL2HexBasis<T, 1, degree - 1>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  using DataElemVec = ElementVector_Serial<T, DataBasis>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis>;
  using ElemVec = ElementVector_Serial<T, Basis>;
  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Matrix-free operator for the problem
  using MatFree = MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  const index_t low_degree = 1;
  using LOrderQuadrature = HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = FEBasis<T>;
  using LOrderGeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, low_degree>>;
  // using LOrderBasis = FEBasis<T, QHdivHexBasis<T, low_degree>,
  //                             LagrangeL2HexBasis<T, 1, low_degree - 1>>;
  using LOrderBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, low_degree>>;
  using LOrderDataElemVec = ElementVector_Serial<T, LOrderDataBasis>;
  using LOrderGeoElemVec = ElementVector_Serial<T, LOrderGeoBasis>;
  using LOrderElemVec = ElementVector_Serial<T, LOrderBasis>;
  using LOrderFE = FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                 LOrderGeoBasis, LOrderBasis>;

  std::cout << "Poisson\n";
  Poisson<std::complex<T>, dim> poisson;
  TestPDEImplementation<std::complex<T>>(poisson);

  std::cout << "Mixed Poisson\n";
  MixedPoisson<std::complex<T>, dim> mixed_poisson;
  TestPDEImplementation<std::complex<T>>(mixed_poisson);

  std::cout << "Nonlinear elasticity\n";
  NonlinearElasticity<std::complex<T>, dim> elasticity;
  TestPDEImplementation<std::complex<T>>(elasticity);

  std::cout << "Heat conduction\n";
  HeatConduction<std::complex<T>, dim> heat_conduction;
  TestPDEImplementation<std::complex<T>>(heat_conduction);

  std::cout << "Mixed heat conduction\n";
  MixedHeatConduction<std::complex<T>, dim> mixed_heat_conduction;
  TestPDEImplementation<std::complex<T>>(mixed_heat_conduction);

  // Number of elements in each dimension
  const int nx = 4, ny = 4, nz = 4;
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

  int boundary1_verts[(ny + 1) * (nz + 1)];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary1_verts[index] = node_num(0, j, k);
    }
  }

  int boundary2_verts[(ny + 1) * (nz + 1)];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary2_verts[index] = node_num(nx, j, k);
    }
  }

  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);
  ElementMesh<DataBasis> datamesh(conn);

  HexProjection<degree, Basis, LOrderBasis> basis_proj;
  HexProjection<degree, GeoBasis, LOrderGeoBasis> geo_proj;
  HexProjection<degree, DataBasis, LOrderDataBasis> data_proj;

  ElementMesh<LOrderBasis> lorder_mesh(mesh, basis_proj);
  ElementMesh<LOrderGeoBasis> lorder_geomesh(geomesh, geo_proj);
  ElementMesh<LOrderDataBasis> lorder_datamesh(datamesh, data_proj);

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select1[2] = {0};
  BoundaryCondition<Basis> bcs1(conn, mesh, basis_select1, (ny + 1) * (nz + 1),
                                boundary1_verts);

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select2[2] = {0};
  BoundaryCondition<Basis> bcs2(conn, mesh, basis_select2, (ny + 1) * (nz + 1),
                                boundary2_verts);

  std::cout << "Number of elements:            " << conn.get_num_elements()
            << std::endl;
  std::cout << "Number of degrees of freedom:  " << mesh.get_num_dof()
            << std::endl;

  // const index_t *bcs_index;
  // std::cout << "Number of boundary conditions: " << bcs1.get_bcs(&bcs_index)
  //           << std::endl;
  // std::cout << "Number of boundary conditions: " << bcs2.get_bcs(&bcs_index)
  //           << std::endl;

  PDE pde;

  index_t ndof = mesh.get_num_dof();
  SolutionVector<T> global_U(mesh.get_num_dof());
  SolutionVector<T> global_res(mesh.get_num_dof());
  SolutionVector<T> global_geo(geomesh.get_num_dof());
  SolutionVector<T> global_data(datamesh.get_num_dof());

  DataElemVec elem_data(datamesh, global_data);
  GeoElemVec elem_geo(geomesh, global_geo);
  ElemVec elem_sol(mesh, global_U);
  ElemVec elem_res(mesh, global_res);

  // Set the geometry from the node locations
  set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  SolutionVector<T> x(mesh.get_num_dof());
  SolutionVector<T> y(mesh.get_num_dof());
  SolutionVector<T> z(mesh.get_num_dof());

  ElemVec elem_x(mesh, x);
  ElemVec elem_y(mesh, y);
  ElemVec elem_z(mesh, z);

  DOFCoordinates<T, PDE, GeoBasis, Basis> coords;
  coords.get_dof_coordinates(elem_geo, elem_x, elem_y, elem_z);

  SolutionVector<T> global_xvec(mesh.get_num_dof());
  SolutionVector<T> global_yvec(mesh.get_num_dof());
  ElemVec elem_xvec(mesh, global_xvec);
  ElemVec elem_yvec(mesh, global_yvec);

  // Create the finite-element model
  FE fe;
  MatFree matfree;

  // Add the residual
  elem_res.init_zero_values();
  fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);
  elem_res.add_values();

  fe.add_jacobian_vector_product(pde, elem_data, elem_geo, elem_sol, elem_xvec,
                                 elem_yvec);

  // Initialize the matrix-free data
  matfree.initialize(pde, elem_data, elem_geo, elem_sol);

  elem_yvec.init_zero_values();
  matfree.add_jacobian_vector_product(elem_xvec, elem_yvec);
  elem_yvec.add_values();

  LOrderFE lorder_fe;

  std::cout << "create_block_matrix" << std::endl;
  index_t nrows;
  std::vector<index_t> rowp, cols;
  lorder_mesh.create_block_csr<1>(nrows, rowp, cols);
  auto mat = std::make_shared<BSRMat<index_t, T, 1, 1>>(
      nrows, nrows, cols.size(), rowp, cols);

  std::cout << "initialize element matrix" << std::endl;
  ElementMat_Serial<T, LOrderBasis, BSRMat<index_t, T, 1, 1>> elem_mat(
      lorder_mesh, *mat);

  // Set the low order element vectors
  LOrderDataElemVec lorder_elem_data(lorder_datamesh, global_data);
  LOrderGeoElemVec lorder_elem_geo(lorder_geomesh, global_geo);
  LOrderElemVec lorder_elem_sol(lorder_mesh, global_U);

  std::cout << "add_jacobian" << std::endl;
  lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                         lorder_elem_sol, elem_mat);

  const index_t *bc_dofs1;
  index_t nbcs1 = bcs1.get_bcs(&bc_dofs1);
  mat->zero_rows(nbcs1, bc_dofs1);

  MultiArrayNew<T *[1][1]> B("B", global_U.get_num_dof());
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    B(i, 0, 0) = 1.0;
  }

  for (index_t i = 0; i < nbcs1; i++) {
    B(bc_dofs1[i], 0, 0) = 0.0;
  }

  index_t num_levels = 3;
  double omega = 3.0 / 4.0;
  double epsilon = 0.0;
  bool print_info = true;
  BSRMatAmg<index_t, T, 1, 1> amg(num_levels, omega, epsilon, mat, B,
                                  print_info);

  MultiArrayNew<T *[1]> xvec("x", global_U.get_num_dof());
  MultiArrayNew<T *[1]> rhs("rhs", global_U.get_num_dof());

  for (index_t i = 0; i < nbcs1; i++) {
    rhs(bc_dofs1[i], 0) = 1.0;
  }

  amg.cg(rhs, xvec, 5);

  // Copy the solution to the solution vector
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    global_U[i] = xvec(i, 0);
  }

  write_hex_to_vtk<1, degree, T, DataBasis, GeoBasis, Basis>(
      pde, elem_data, elem_geo, elem_sol,
      [](index_t k, typename PDE::DataSpace &data,
         typename PDE::FiniteElementGeometry &geo,
         typename PDE::FiniteElementSpace &sol) {
        return sol.template get<0>().get_value();
      });

  // Kokkos::finalize();

  return (0);
}