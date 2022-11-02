#include <iostream>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/fespace.h"
#include "sparse/sparse_matrix.h"

using namespace A2D;

int main(int argc, char* argv[]) {
  typedef double T;
  typedef NonlinearElasticity2D<T> PDE;
  typedef TriQuadrature3 Quadrature;

  typedef FEBasis<T, LagrangeTri1<T, 2>> GeoBasis;
  typedef FEBasis<T, LagrangeTri1<T, 2>> Basis;

  // type PoissonMixed2D PDE;
  // typedef FEBasis<T, LagrangeTri1<T, 2>> GeoBasis;
  // typedef FEBasis<T, LagrangeTri0<T>, RT2DTri1<T>> Basis;

  typedef FiniteElement<T, Quadrature, PDE, GeoBasis, Basis> FE;

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

  for (index_t j = 0; j < ny + 1; j++) {
    for (index_t i = 0; i < nx + 1; i++) {
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

  // Create the mesh for the geometry
  SpaceType geo_space[] = {H1};
  index_t dims[] = {2};
  ElementMesh<GeoBasis> geomesh(connect, geo_space, dims);

  SpaceType sol_space[] = {L2, EDGE};
  ElementMesh<Basis> mesh(connect, sol_space);

  // Get the total number of degrees of freedom
  index_t ndof = mesh.get_num_dof();

  SolutionVector<T> sol(ndof), res(ndof), pert(ndof);
  FE fe(geomesh, nodes, mesh, sol, res);

  for (index_t i = 0; i < ndof; i++) {
    pert[i] = 1.0;
  }

  fd.add_residual();
  // fd.add_jacobian_vector_product(pert, res);
  // poisson.add_jacobian();

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

  return (0);
}