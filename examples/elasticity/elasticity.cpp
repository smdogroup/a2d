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
#ifdef USE_COMPLEX
  typedef std::complex<double> ScalarType;
#else
  typedef double ScalarType;
#endif  // USE_COMPLEX

  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  // typedef HelmholtzPDE<IndexType, ScalarType, Basis> Model;
  // typedef NonlinearElasticity3D<IndexType, ScalarType, Basis> Model;
  typedef LinearElasticity3D<IndexType, ScalarType, Basis> Model;

  const int nx = 64;
  const int ny = 64;
  const int nz = 64;
  const int nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const int nelems = nx * ny * nz;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;

  Model model(nelems, nnodes);

  typename Model::base::ConnArray& conn = model.get_conn();
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

  // Set the values of the solution
  for (int i = 0; i < nnodes; i++) {
    for (int j = 0; j < vars_per_node; j++) {
      U(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      P(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      U(i, j) = U(i, j) + perturb * P(i, j);
    }
  }

  // Reset the soultion
  model.reset_solution();

  element_scatter(conn, P, Pe);

  typename Model::base::ElemResArray& elem_res = model.get_elem_res();
  typename Model::base::ElemJacArray& elem_jac = model.get_elem_jac();

  residual.zero();
  elem_res.zero();
  elem_jac.zero();

  model.add_residuals(residual);

  double t0 = MPI_Wtime();
  model.add_jacobians();
  t0 = MPI_Wtime() - t0;

  std::cout << "Jacobian time: " << t0 << std::endl;

#ifdef USE_COMPLEX

  for (int i = 0; i < 3 && i < elem_jac.extent(0); i++) {
    for (int jy = 0; jy < nodes_per_elem; jy++) {
      for (int iy = 0; iy < vars_per_node; iy++) {
        double fd = elem_res(i, jy, iy).imag() / dh;

        ScalarType result = 0.0;
        for (int jx = 0; jx < nodes_per_elem; jx++) {
          for (int ix = 0; ix < vars_per_node; ix++) {
            result += elem_jac(i, jy, jx, iy, ix) * Pe(i, jx, ix);
          }
        }

        std::cout << "fd: " << std::setw(20) << fd
                  << " result: " << std::setw(20) << result.real()
                  << " error: " << std::setw(20) << (fd - result.real()) / fd
                  << std::endl;
      }
    }
  }

#endif  // USE_COMPLEX

  CLayout<2> bcs_layout((nz + 1) * (ny + 1));
  MultiArray<index_t, CLayout<2>> bcs(bcs_layout);
  bcs.zero();

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

  typedef A2D::BSRMat<IndexType, ScalarType, vars_per_node, vars_per_node>
      SparseMat;

  // Create the sparse matrix representing the Jacobian matrix
  SparseMat& J =
      *BSRMatFromConnectivity2<IndexType, ScalarType, vars_per_node>(conn);

  // Form the Jacobian matrix
  J.zero();
  BSRMatAddElementMatrices(conn, model.get_elem_jac(), J);

  // // Zero rows associated with the boundary conditions
  BSRMatZeroBCRows(bcs, J);

  // The near null-space basis for 3D elasticity
  const index_t null_space_basis = 6;
  CLayout<vars_per_node, null_space_basis> near_nullspace_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node, null_space_basis>> B(
      near_nullspace_layout);

  // Form the near null - space basis
  B.zero();
  for (IndexType i = 0; i < nnodes; i++) {
    B(i, 0, 0) = 1.0;
    B(i, 1, 1) = 1.0;
    B(i, 2, 2) = 1.0;

    if (null_space_basis == 6) {
      // Rotation about the x-axis
      B(i, 1, 3) = X(i, 2);
      B(i, 2, 3) = -X(i, 1);

      // Rotation about the y-axis
      B(i, 0, 4) = X(i, 2);
      B(i, 2, 4) = -X(i, 0);

      // Rotation about the z-axis
      B(i, 0, 5) = X(i, 1);
      B(i, 1, 5) = -X(i, 0);
    }
  }

  // const index_t null_space_basis = 1;
  // CLayout<vars_per_node, null_space_basis> near_nullspace_layout(nnodes);
  // MultiArray<ScalarType, CLayout<vars_per_node, null_space_basis>> B(
  //     near_nullspace_layout);

  // // Form the near null - space basis
  // B.zero();
  // for (IndexType i = 0; i < nnodes; i++) {
  //   B(i, 0, 0) = 1.0;
  // }

  // Apply the boundary conditions to the null space
  VecZeroBCRows(bcs, B);

  double t1 = MPI_Wtime();

  int num_levels = 4;
  double omega = 1.333;
  bool print_info = true;
  BSRMatAmg<IndexType, ScalarType, vars_per_node, null_space_basis>* amg =
      new BSRMatAmg<IndexType, ScalarType, vars_per_node, null_space_basis>(
          num_levels, omega, &J, &B, print_info);

  t1 = MPI_Wtime() - t1;
  std::cout << "Set up time for AMG: " << t1 << std::endl;

  // Set the residuals and apply the boundary conditions
  for (IndexType i = 0; i < solution.extent(0); i++) {
    for (IndexType j = 0; j < solution.extent(1); j++) {
      solution(i, j) = 1.0;
    }
  }
  residual.zero();
  VecZeroBCRows(bcs, solution);
  BSRMatVecMult(J, solution, residual);

  // Set the solution back to zero
  solution.zero();

  IndexType monitor = 5;
  IndexType max_iters = 80;
  amg->cg(residual, solution, monitor, max_iters);
  amg->mg(residual, solution, monitor, max_iters);

  return (0);
}
