#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "mpi.h"
#include "multiarray.h"
#include "sparse_amg.h"
#include "sparse_matrix.h"
#include "sparse_numeric.h"
#include "sparse_symbolic.h"

// #define USE_COMPLEX

using namespace A2D;

/*
  Finite-element computations using templating/auto-diff and multi-dimensional
  arrays
*/
int main(int argc, char* argv[]) {
  typedef A2D::index_t IndexType;
  typedef double ScalarType;
  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef HelmholtzPDE<IndexType, ScalarType, Basis> Model;
  // typedef NonlinearElasticity3D<IndexType, ScalarType, Basis> Model;
  // typedef LinearElasticity3D<IndexType, ScalarType, Basis> Model;

  const index_t nx = 64;
  const index_t ny = 64;
  const index_t nz = 64;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);
  const index_t vars_per_node = Model::NUM_VARS;
  const index_t nodes_per_elem = Basis::NUM_NODES;

  Model model(nelems, nnodes, nbcs);

  typename Model::base::ConnArray& conn = model.get_conn();
  typename Model::base::BCsArray& bcs = model.get_bcs();
  typename Model::base::NodeArray& X = model.get_nodes();
  typename Model::base::SolutionArray& U = model.get_solution();

  // Create a new solution array
  typename Model::base::SolutionArray& P = *model.new_solution();
  typename Model::base::ElemSolnArray& Pe = *model.new_elem_solution();

  // Residual vector
  typename Model::base::SolutionArray& residual = *model.new_solution();
  typename Model::base::SolutionArray& solution = *model.new_solution();
  typename Model::base::SolutionArray& res = *model.new_solution();

#ifdef USE_COMPLEX
  double dh = 1e-30;
  ScalarType perturb = ScalarType(0.0, dh);
#else
  double dh = 1e-6;
  ScalarType perturb = dh;
#endif  // USE_COMPLEX

  // Set the boundary conditions
  index_t index = 0;
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      int i = 0;
      int node = i + (nx + 1) * (j + (ny + 1) * k);

      // Set the boundary conditions
      bcs(index, 0) = node;
      for (int ii = 0; ii < 3; ii++) {
        bcs(index, 1) |= 1U << ii;
      }
      index++;
    }
  }

  // Set the node locations
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        X(node, 0) = 1.0 * i / nx;
        X(node, 1) = 1.0 * j / ny;
        X(node, 2) = 1.0 * k / nz;
      }
    }
  }

  // Set the element connectivity
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int elem = i + nx * (j + ny * k);

        int conn_coord[8];
        for (int kk = 0, index = 0; kk < 2; kk++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int ii = 0; ii < 2; ii++, index++) {
              conn_coord[index] =
                  (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
            }
          }
        }

        // Convert to the correct connectivity
        conn(elem, 0) = conn_coord[0];
        conn(elem, 1) = conn_coord[1];
        conn(elem, 2) = conn_coord[3];
        conn(elem, 3) = conn_coord[2];

        conn(elem, 4) = conn_coord[4];
        conn(elem, 5) = conn_coord[5];
        conn(elem, 6) = conn_coord[7];
        conn(elem, 7) = conn_coord[6];
      }
    }
  }

  typename Model::base::QuadDataArray& data = model.get_quad_data();

  double pt_data[] = {1.23, 2.45};
  for (A2D::index_t i = 0; i < data.extent(0); i++) {
    for (A2D::index_t j = 0; j < data.extent(1); j++) {
      for (A2D::index_t k = 0; k < data.extent(2); k++) {
        data(i, j, k) = pt_data[k];
      }
    }
  }

  // Now reset the nodes
  model.reset_nodes();

  // Compute the Jacobian matrix
  double t0 = MPI_Wtime();
  typename Model::base::SparseMat* J = model.new_matrix();
  J->zero();
  model.add_jacobian(*J);
  t0 = MPI_Wtime() - t0;

  std::cout << "Jacobian time: " << t0 << std::endl;

  double t1 = MPI_Wtime();
  int num_levels = 4;
  double omega = 1.333;
  bool print_info = true;

  typename Model::base::SparseAmg* amg =
      model.new_amg(num_levels, omega, J, print_info);

  t1 = MPI_Wtime() - t1;
  std::cout << "Set up time for AMG: " << t1 << std::endl;

  // Set the residuals and apply the boundary conditions
  solution.fill(1.0);
  residual.zero();
  VecZeroBCRows(bcs, solution);
  BSRMatVecMult(*J, solution, residual);

  // Set the solution back to zero
  solution.zero();

  IndexType monitor = 10;
  IndexType max_iters = 80;
  double t2 = MPI_Wtime();
  amg->cg(residual, solution, monitor, max_iters);
  t2 = MPI_Wtime() - t2;
  std::cout << "Conjugate gradient solution time: " << t2 << std::endl;

  double t3 = MPI_Wtime();
  amg->mg(residual, solution, monitor, max_iters);
  t3 = MPI_Wtime() - t3;
  std::cout << "Multigrid solution time: " << t3 << std::endl;

  return (0);
}
