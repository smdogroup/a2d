#ifndef TOPOLOGY_HEAT_H
#define TOPOLOGY_HEAT_H

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dprofiler.h"

void test_heat_analysis(int argc, char *argv[]) {
  // Magic numbers
  const int degree = 1;            // polynomial degree, 1 = linear
  const int order = 2;             // order = degree + 1
  const int spatial_dim = 3;       // spatial dimension
  const int var_dim = 1;           // solution variable dimension
  const int data_dim = 1;          // design variable dimension
  const int block_size = var_dim;  // block size for BSR matrix

  /* Types */

  // Basic type
  using T = double;
  using I = A2D::index_t;

  // Quadrature and basis
  using Quadrature = A2D::HexGaussQuadrature<order>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using DataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, degree - 1>>;

  // Block sparse compressed row matrix
  using BSRMatType = A2D::BSRMat<I, T, block_size, block_size>;

  // Element-centric views
  using GlobalVecType = A2D::SolutionVector<T>;
  using ElemVec = A2D::ElementVector_Serial<T, Basis, GlobalVecType>;
  using GeoElemVec = A2D::ElementVector_Serial<T, GeoBasis, GlobalVecType>;
  using DataElemVec = A2D::ElementVector_Serial<T, DataBasis, GlobalVecType>;
  using ElemMat = A2D::ElementMat_Serial<T, Basis, BSRMatType>;

  // Physics and functional
  using PDE = A2D::HeatConduction<T, spatial_dim>;
  using FE = A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  /* Load mesh and boundary vertices from vtk */

  // Load vtk
  std::string vtk_name = "3d_hex_small.vtk";
  if (argc > 1) {
    vtk_name = argv[1];
  }
  A2D::ReadVTK3D<I, T> readvtk(vtk_name);

  // Get connectivity for each element type
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;

  // Construct connectivity
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Extract boundary vertices
  std::vector<int> ids{100, 101};
  std::vector<I> bc_verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Construct A2D mesh objects
  A2D::ElementMesh<Basis> mesh(conn);
  A2D::ElementMesh<GeoBasis> geomesh(conn);
  A2D::ElementMesh<DataBasis> datamesh(conn);

  /* Construct solution, data and X vectors and their element views */

  // Allocate vectors
  GlobalVecType sol(mesh.get_num_dof());
  GlobalVecType res(mesh.get_num_dof());
  GlobalVecType data(datamesh.get_num_dof());
  GlobalVecType geo(geomesh.get_num_dof());

  // Allocate elem views
  ElemVec elem_sol(mesh, sol);
  ElemVec elem_res(mesh, res);
  DataElemVec elem_data(datamesh, data);
  GeoElemVec elem_geo(geomesh, geo);

  // Populate data and geo
  data.fill(1.0);
  A2D::set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  /* Construct the Jacobian matrix */

  // Construct the nonzero pattern
  I nrows;
  std::vector<I> rowp, cols;
  mesh.create_block_csr<block_size>(nrows, rowp, cols);
  std::shared_ptr<BSRMatType> mat =
      std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);

  // Populate the stiffness matrix via element view
  ElemMat elem_mat(mesh, *mat);
  FE fe;
  T kappa = 1.0, q = 0.0, heat_source = 1.0, bc_temp = 0.0;
  PDE pde(kappa, q, heat_source);
  fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

  mat->write_mtx("heat_jacobian_nobc.mtx");

  // Apply boundary condition
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data()));
  A2D::DirichletBCs<Basis> bcs(conn, mesh, bcinfo);
  const I *bc_dofs;
  I nbcs = bcs.get_bcs(&bc_dofs);
  mat->zero_rows(nbcs, bc_dofs);

  // Write the Jacobian matrix to mtx for visualization
  mat->write_mtx("heat_jacobian_withbc.mtx");

  /* Construct the residual and rhs */

  fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

  I size = sol.get_num_dof() / block_size;
  A2D::MultiArrayNew<T *[block_size]> sol_vec("sol_vec", size);
  A2D::MultiArrayNew<T *[block_size]> rhs_vec("rhs_vec", size);

  for (I i = 0; i < sol.get_num_dof(); i++) {
    rhs_vec(i / block_size, i % block_size) = -res[i];
  }

  // Set the temperature boundary conditions
  for (I i = 0; i < nbcs; i++) {
    I dof = bc_dofs[i];
    rhs_vec(dof / block_size, dof % block_size) = bc_temp;
  }

  /* Construct the AMG solver and solve the problem */

  // Create amg
  A2D::index_t num_levels = 3;
  double omega = 4.0 / 3.0;
  double epsilon = 0.0;
  bool print_info = false;
  const int null_size = 1;
  A2D::MultiArrayNew<T *[block_size][null_size]> B(
      "B", sol.get_num_dof() / block_size);
  A2D::BLAS::fill(B, 1.0);
  A2D::BSRMatAmg<I, T, block_size, null_size> amg(num_levels, omega, epsilon,
                                                  mat, B, print_info);

  // Solve
  auto mat_vec = [&](A2D::MultiArrayNew<T *[block_size]> &in,
                     A2D::MultiArrayNew<T *[block_size]> &out) -> void {
    A2D::BSRMatVecMult<I, T, block_size, block_size>(*mat, in, out);
  };
  amg.cg(mat_vec, rhs_vec, sol_vec, 5, 100);

  // Record the solution
  for (I i = 0; i < sol.get_num_dof(); i++) {
    sol[i] = sol_vec(i / block_size, i % block_size);
  }

  // Write result to vtk
  A2D::write_hex_to_vtk<2, degree, T, DataBasis, GeoBasis, Basis>(
      pde, elem_data, elem_geo, elem_sol, "heat_analysis.vtk",
      [](A2D::index_t k, typename PDE::DataSpace &d,
         typename PDE::FiniteElementGeometry &g,
         typename PDE::FiniteElementSpace &s) {
        if (k == 0) {
          return (s.template get<0>()).get_value();  // state
        } else {
          return (d.template get<0>()).get_value();  // design
        }
      });
}

#endif