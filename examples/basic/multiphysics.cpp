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
#include "utils/a2dvtk.h"

using namespace A2D;

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const index_t dim = 3;
  const index_t degree = 2;
  const index_t geo_degree = 1;
  using T = double;
  using ET = ElementTypes;

  /*
  using PDE = NonlinearElasticity<T, 3>;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 2, 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, geo_degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, dim, degree>>;
  //*/

  //*
  using PDE = MixedPoisson<T, dim>;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, geo_degree>>;
  using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
                        LagrangeL2HexBasis<T, 1, degree - 1>>;
  // */

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

        for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
          index_t node = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                  j + ET::HEX_VERTS_CART[ii][1],
                                  k + ET::HEX_VERTS_CART[ii][2]);

          // Set the entity DOF
          index_t basis = 0;
          index_t orient = 0;
          GeoBasis::set_entity_dof(basis, ET::VERTEX, ii, orient,
                                   &Xloc[3 * node], geo_dof);
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
  index_t bc_dofs2[] = {global_x.get_num_dof() - 1};

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

  BSRMatApplyFactor(*factor, rhs, x);

  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    std::cout << "rhs[" << i << "]: " << rhs(i, 0) << std::endl;
  }

  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    std::cout << "u[" << i << "]: " << x(i, 0) << std::endl;
  }

  // Copy the solution to the solution vector
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    global_U[i] = x(i, 0);
  }

  const index_t nex = 3;
  index_t nvtk_elems = nhex * nex * nex * nex;
  index_t nvtk_nodes = nhex * (nex + 1) * (nex + 1) * (nex + 1);

  auto vtk_node_num = [](index_t i, index_t j, index_t k) {
    return i + j * (nex + 1) + k * (nex + 1) * (nex + 1);
  };

  MultiArrayNew<int *[8]> vtk_conn("vtk_elems", nvtk_elems);
  MultiArrayNew<double *[3]> vtk_nodes("vtk_nodes", nvtk_nodes);
  MultiArrayNew<double *> vtk_solt("vtk_nodes", nvtk_nodes);
  MultiArrayNew<double *> vtk_solqx("vtk_nodes", nvtk_nodes);
  MultiArrayNew<double *> vtk_solqy("vtk_nodes", nvtk_nodes);
  MultiArrayNew<double *> vtk_solqz("vtk_nodes", nvtk_nodes);

  // for (index_t i = 0; i < global_U.get_num_dof(); i++) {
  //   global_U[i] = 0.0;
  // }

  // T entity_vals[degree * degree];
  // for (index_t k = 0; k < degree * degree; k++) {
  //   entity_vals[k] = 0.0;
  // }
  // entity_vals[2] = 1.0;

  // typename FE::ElemVec::FEDof sol_dof(0, elem_sol);
  // elem_sol.get_element_values(0, sol_dof);

  // // Set the entity DOF
  // index_t basis = 0;
  // index_t orient = 0;
  // Basis::set_entity_dof(basis, ET::FACE, 1, orient, entity_vals, sol_dof);

  // elem_sol.set_element_values(0, sol_dof);

  for (index_t n = 0, counter = 0; n < nhex; n++) {
    typename FE::ElemVec::FEDof sol_dof(n, elem_sol);
    elem_sol.get_element_values(n, sol_dof);

    typename FE::GeoElemVec::FEDof geo_dof(n, elem_geo);
    elem_geo.get_element_values(n, geo_dof);

    for (index_t k = 0; k < nex + 1; k++) {
      for (index_t j = 0; j < nex + 1; j++) {
        for (index_t i = 0; i < nex + 1; i++) {
          index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

          index_t index = vtk_node_num(i, j, k);
          typename PDE::FiniteElementSpace sol;
          typename PDE::FiniteElementGeometry geo;

          Basis::template interp<HexGaussLobattoQuadrature<nex + 1>>(
              index, sol_dof, sol);
          GeoBasis::template interp<HexGaussLobattoQuadrature<nex + 1>>(
              index, geo_dof, geo);

          auto X = geo.template get<0>().get_value();
          auto sigma = sol.template get<0>().get_value();
          auto u = sol.template get<1>().get_value();

          index_t node = off + vtk_node_num(i, j, k);
          vtk_nodes(node, 0) = X(0);
          vtk_nodes(node, 1) = X(1);
          vtk_nodes(node, 2) = X(2);

          vtk_solt(node) = u;
          vtk_solqx(node) = sigma(0);
          vtk_solqy(node) = sigma(1);
          vtk_solqz(node) = sigma(2);
        }
      }
    }

    for (index_t k = 0; k < nex; k++) {
      for (index_t j = 0; j < nex; j++) {
        for (index_t i = 0; i < nex; i++, counter++) {
          index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

          for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
            vtk_conn(counter, ii) =
                off + vtk_node_num(i + ET::HEX_VERTS_CART[ii][0],
                                   j + ET::HEX_VERTS_CART[ii][1],
                                   k + ET::HEX_VERTS_CART[ii][2]);
          }
        }
      }
    }
  }

  ToVTK vtk(vtk_conn, vtk_nodes);
  vtk.write_mesh();
  vtk.write_sol("t", vtk_solt);
  vtk.write_sol("qx", vtk_solqx);
  vtk.write_sol("qy", vtk_solqy);
  vtk.write_sol("qz", vtk_solqz);

  // amg.cg(rhs, x, 5);
  // * /
  // Kokkos::finalize();

  return (0);
}