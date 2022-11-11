#include <iostream>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/qhdiv_hex_basis.h"

// #include "sparse/sparse_matrix.h"

#include "multiphysics/femesh.h"

using namespace A2D;

int main(int argc, char* argv[]) {
  const index_t degree = 2;
  using T = double;
  using PDE = NonlinearElasticity<T, 3>;
  using Quadrature = TriQuadrature3;
  using DataBasis = FEBasis<T, LagrangeH1HexBasis<T, 2, 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, 3, 1>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>;

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

  // index_t nelems = nx * ny * nz;
  // index_t nfaces = nx * ny * (nz + 1) + nx * (ny + 1) * nz + (nx + 1) * ny *
  // nz; index_t nedges = nx * (ny + 1) * (nz + 1) + (nx + 1) * ny * (nz + 1) +
  //                  (nx + 1) * (ny + 1) * nz;

  // std::cout << "Number of elements: " << nelems << " "
  //           << conn.get_num_elements() << std::endl;
  // std::cout << "Number of faces:    " << nfaces << " " <<
  // conn.get_num_faces()
  //           << std::endl;
  // std::cout << "Number of edges:    " << nedges << " " <<
  // conn.get_num_edges()
  //           << std::endl;

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

  BSRMat<index_t, T, 3, 3>* mat =
      mesh.create_

          ElementMat_Serial<T, Basis, BSRMat<index_t, T, 3, 3>>
              elem_mat(mesh);

  fe.add_jacobian(elem_data, elem_geo, elem_sol, elem_mat);

  // // Allocate global data, X, U, residual vectors
  // SolutionVector<T> global_data(2 * nnodes);
  // SolutionVector<T> global_X(2 * nnodes);
  // SolutionVector<T> global_U(2 * nnodes);
  // SolutionVector<T> global_res(2 * nnodes);

  // // Create element vector views
  // FE::DataElemVec elem_data(datamesh, global_data);
  // FE::GeoElemVec elem_geo(mesh, global_X);
  // FE::ElemVec elem_sol(mesh, global_U);
  // FE::ElemVec elem_vec(mesh, global_res);

  // // Fabricate data, X, and U
  // // TODO

  // // Populate element views using global vectors
  // elem_data.init_values();
  // elem_geo.init_values();
  // elem_sol.init_values();

  // // Create finite element instance
  // FE fe(elem_data, elem_geo, elem_sol, elem_vec);

  // fe.add_residual(global_U);

  // // fd.add_jacobian_vector_product(pert, res);
  // fe.add_jacobian();

  // // Create the node vector
  // SolutionVector<T> nodes(2 * nnodes);

  // for (index_t j = 0; j < ny + 1; j++) {
  //   for (index_t i = 0; i < nx + 1; i++) {
  //     index_t node = i + j * (nx + 1);
  //     nodes[2 * node] = 1.0 * i;
  //     nodes[2 * node + 1] = 1.0 * j;
  //   }
  // }

  // // Set the connectivity
  // index_t* conn = new index_t[3 * nelems];

  // for (index_t j = 0; j < ny; j++) {
  //   for (index_t i = 0; i < nx; i++) {
  //     index_t elem = 2 * (i + j * nx);
  //     conn[3 * elem] = i + j * (nx + 1);
  //     conn[3 * elem + 1] = i + 1 + j * (nx + 1);
  //     conn[3 * elem + 2] = i + 1 + (j + 1) * (nx + 1);

  //     elem += 1;
  //     conn[3 * elem] = i + j * (nx + 1);
  //     conn[3 * elem + 1] = i + 1 + (j + 1) * (nx + 1);
  //     conn[3 * elem + 2] = i + (j + 1) * (nx + 1);
  //   }
  // }

  // ElementConnectivity connect(nnodes, nelems, conn);
  // delete[] conn;

  // // Create the mesh for the geometry
  // SpaceType geo_space[] = {H1};
  // index_t dims[] = {2};
  // ElementMesh<GeoBasis> geomesh(connect, geo_space, dims);

  // SpaceType data_space[] = {H1};
  // ElementMesh<DataBasis> datamesh(connect, data_space, dims);

  // SpaceType sol_space[] = {H1};
  // ElementMesh<Basis> mesh(connect, sol_space, dims);

  // // Get the total number of degrees of freedom
  // index_t ndof = mesh.get_num_dof();

  // // Allocate global data, X, U, residual vectors
  // SolutionVector<T> global_data(2 * nnodes);
  // SolutionVector<T> global_X(2 * nnodes);
  // SolutionVector<T> global_U(2 * nnodes);
  // SolutionVector<T> global_res(2 * nnodes);

  // // Create element vector views
  // FE::DataElemVec elem_data(datamesh, global_data);
  // FE::GeoElemVec elem_geo(mesh, global_X);
  // FE::ElemVec elem_sol(mesh, global_U);
  // FE::ElemVec elem_vec(mesh, global_res);

  // // Fabricate data, X, and U
  // // TODO

  // // Populate element views using global vectors
  // elem_data.init_values();
  // elem_geo.init_values();
  // elem_sol.init_values();

  // // Create finite element instance
  // FE fe(elem_data, elem_geo, elem_sol, elem_vec);

  // fe.add_residual(global_U);

  // // fd.add_jacobian_vector_product(pert, res);
  // fe.add_jacobian();

  // typedef double T;
  // typedef NonlinearElasticity<T> PDE;

  // typename PDE::FiniteElementGeometry geo;
  // typename PDE::DataSpace data;
  // typename PDE::FiniteElementSpace s, coef;

  // T wdetJ = 0.5842;
  // //   PDE::eval_weak_coef(wdetJ, data, s, coef);

  // typename PDE::JacVecProduct jvp(wdetJ, data, s);

  // typename PDE::FiniteElementSpace p, Jp;
  // jvp(p, Jp);

  // type PoissonMixed2D PDE;
  // typedef FEBasis<T, LagrangeTri1<T, 2>> GeoBasis;
  // typedef FEBasis<T, LagrangeTri0<T, 1>, RT2DTri1<T>> Basis;

  // A2D::Mat<T, Basis::ncomp, Basis::ncomp> mat;
  // A2D::Mat<T, Basis::ndof, Basis::ndof> jac;
  // Basis::add_outer<Quadrature>(0, mat, jac);

  // T dof[Basis::ndof];
  // for (int i = 0; i < Basis::ndof; i++) {
  //   dof[i] = -1.32 + 0.31 * i;
  // }

  // FiniteElementSpace s1, s2;
  // Basis::interp<Quadrature>(0, dof, s1);
  // Basis::interp_basis<Quadrature>(0, dof, s2);

  // for (int i = 0; i < FiniteElementSpace::ncomp; i++) {
  //   std::cout << s1.get_value(i) << "  " << s2.get_value(i) << std::endl;
  // }

  // std::fill(dof, dof + Basis::ndof, T(0.0));
  // Basis::add<Quadrature>(0, s1, dof);

  // for (int i = 0; i < Basis::ndof; i++) {
  //   std::cout << dof[i] << std::endl;
  // }

  // std::fill(dof, dof + Basis::ndof, T(0.0));
  // Basis::add<Quadrature>(0, s2, dof);

  // for (int i = 0; i < Basis::ndof; i++) {
  //   std::cout << dof[i] << std::endl;
  // }

  return (0);
}