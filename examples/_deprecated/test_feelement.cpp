#include <iostream>

#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/integrand_elasticity.h"

using namespace A2D;

void main_body() {
  using T = double;
  using PDEIntegrand = NonlinearElasticity<T, 2>;
  using Quadrature = TriQuadrature3;
  using DataBasis = FEBasis<T, LagrangeTri0<T, 2>>;
  using GeoBasis = FEBasis<T, LagrangeTri1<T, 2>>;
  using Basis = FEBasis<T, LagrangeTri1<T, 2>>;

  // constexpr bool use_parallel_elemvec = false;
  constexpr bool use_parallel_elemvec = true;
  using FE = ElementOps<T, PDEIntegrand, Quadrature, DataBasis, GeoBasis,
                           Basis, use_parallel_elemvec>;

  // Set the node locations
  index_t nx = 10, ny = 10;
  index_t nnodes = (nx + 1) * (ny + 1);
  index_t nelems = 2 * nx * ny;

  // Create the node vector
  SolutionVector<T> nodes(2 * nnodes);

  for (index_t j = 0; j < ny + 1; j++) {
    for (index_t i = 0; i < nx + 1; i++) {
      index_t node = i + j * (nx + 1);
      nodes[2 * node] = 1.0 * i;
      nodes[2 * node + 1] = 1.0 * j;
    }
  }

  // Set the connectivity
  index_t* conn = new index_t[3 * nelems];

  for (index_t j = 0; j < ny; j++) {
    for (index_t i = 0; i < nx; i++) {
      index_t elem = 2 * (i + j * nx);
      conn[3 * elem] = i + j * (nx + 1);
      conn[3 * elem + 1] = i + 1 + j * (nx + 1);
      conn[3 * elem + 2] = i + 1 + (j + 1) * (nx + 1);

      elem += 1;
      conn[3 * elem] = i + j * (nx + 1);
      conn[3 * elem + 1] = i + 1 + (j + 1) * (nx + 1);
      conn[3 * elem + 2] = i + (j + 1) * (nx + 1);
    }
  }

  ElementConnectivity connect(nnodes, nelems, conn);
  delete[] conn;

  // Create the meshes
  SpaceType geo_space[] = {H1};
  index_t dims[] = {2};
  ElementMesh<GeoBasis> geomesh(connect, geo_space, dims);

  SpaceType data_space[] = {H1};
  ElementMesh<DataBasis> datamesh(connect, data_space, dims);

  SpaceType sol_space[] = {L2, EDGE};
  ElementMesh<Basis> mesh(connect, sol_space);

  // Allocate global data, X, U, residual vectors
  SolutionVector<T> global_data(2 * nnodes);
  SolutionVector<T> global_X(2 * nnodes);
  SolutionVector<T> global_U(2 * nnodes);
  SolutionVector<T> global_res(2 * nnodes);

  // Create element vector views
  FE::DataElemVec elem_data(datamesh, global_data);
  FE::GeoElemVec elem_geo(mesh, global_X);
  FE::ElemVec elem_sol(mesh, global_U);
  FE::ElemVec elem_vec(mesh, global_res);

  // Fabricate data, X, and U
  // TODO

  // Populate element views using global vectors
  elem_data.init_values();
  elem_geo.init_values();
  elem_sol.init_values();

  // Create finite element instance
  FE fe(elem_data, elem_geo, elem_sol, elem_vec);

  fe.add_residual(global_res);
}

int main() {
  Kokkos::initialize();
  { main_body(); }
  Kokkos::finalize();
}