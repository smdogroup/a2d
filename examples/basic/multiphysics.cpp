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
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, 1>>;
  using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
                        LagrangeL2HexBasis<T, 1, degree - 1>>;

  constexpr bool use_parallel_elemvec = false;
  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis,
                           use_parallel_elemvec>;

  // Number of elements in each dimension
  const int nx = 3, ny = 3, nz = 3;
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

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select1[2] = {1, 0};
  BoundaryCondition<Basis> bcs1(conn, mesh, basis_select1, (ny + 1) * (nz + 1),
                                boundary1_verts);

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select2[2] = {0, 1};
  BoundaryCondition<Basis> bcs2(conn, mesh, basis_select2, (ny + 1) * (nz + 1),
                                boundary2_verts);

  std::cout << "Number of elements:            " << conn.get_num_elements()
            << std::endl;
  std::cout << "Number of degrees of freedom:  " << mesh.get_num_dof()
            << std::endl;

  const index_t *bcs_index;
  std::cout << "Number of boundary conditions: " << bcs1.get_bcs(&bcs_index)
            << std::endl;
  std::cout << "Number of boundary conditions: " << bcs2.get_bcs(&bcs_index)
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

  // Create the finite-element model
  FE fe;

  // Add the residual
  elem_res.init_zero_values();
  fe.add_residual(elem_data, elem_geo, elem_sol, elem_res);
  elem_res.add_values();

  // fe.add_jacobian_vector_product(elem_data, elem_geo, elem_sol, elem_x,
  // elem_y);

  std::cout << "create_block_matrix" << std::endl;
  auto mat = mesh.create_block_matrix<T, 1>();
  BSRMat<index_t, T, 1, 1> &mat_ref = *mat;

  std::cout << "initialize element matrix" << std::endl;
  ElementMat_Serial<T, Basis, BSRMat<index_t, T, 1, 1>> elem_mat(mesh, mat_ref);

  std::cout << "add_jacobian" << std::endl;
  fe.add_jacobian(elem_data, elem_geo, elem_sol, elem_mat);

  const index_t *bc_dofs1;
  index_t nbcs1 = bcs1.get_bcs(&bc_dofs1);
  mat->zero_rows(nbcs1, bc_dofs1);

  // const index_t *bc_dofs2;
  // index_t nbcs2 = bcs1.get_bcs(&bc_dofs2);
  index_t nbcs2 = 1;
  index_t bc_dofs2[] = {33199};

  mat->zero_rows(nbcs2, bc_dofs2);

  MultiArrayNew<T *[1][1]> B("B", global_U.get_num_dof());
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    B(i, 0, 0) = 1.0;
  }

  for (index_t i = 0; i < nbcs1; i++) {
    B(bc_dofs1[i], 0, 0) = 0.0;
  }

  for (index_t i = 0; i < nbcs2; i++) {
    B(bc_dofs2[i], 0, 0) = 0.0;
  }

  // index_t num_levels = 0;
  // double omega = 3.0 / 4.0;
  // double epsilon = 0.1;
  // bool print_info = true;
  // BSRMatAmg<index_t, T, 1, 1> amg(num_levels, omega, epsilon, mat, B,
  //                                 print_info);

  BSRMat<index_t, T, 1, 1> *factor = BSRMatFactorSymbolic(*mat);
  BSRMatCopy(*mat, *factor);
  BSRMatFactor(*factor);

  std::cout << "nnz = " << mat->nnz << std::endl;

  MultiArrayNew<T *[1]> x("x", global_U.get_num_dof());
  MultiArrayNew<T *[1]> rhs("rhs", global_U.get_num_dof());

  for (index_t i = 0; i < nbcs1; i++) {
    rhs(bc_dofs1[i], 0) = 1.0;
  }

  // for (index_t i = 0; i < mat->nbrows; i++) {
  //   for (index_t jp = mat->rowp(i); jp < mat->rowp(i + 1); jp++) {
  //     std::cout << mat->Avals(jp, 0, 0) << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  BSRMatApplyFactor(*factor, rhs, x);

  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    std::cout << "rhs[" << i << "]: " << rhs(i, 0) << std::endl;
  }

  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    std::cout << "u[" << i << "]: " << x(i, 0) << std::endl;
  }

  // amg.cg(rhs, x, 5);

  // Kokkos::finalize();

  return (0);
}