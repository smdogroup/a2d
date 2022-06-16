#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "model.h"
#include "mpi.h"

using namespace A2D;

/*
  Finite-element computations using templating/auto-diff and multi-dimensional
  arrays
*/
int main(int argc, char* argv[]) {
  typedef index_t I;
  typedef double T;
  typedef Basis3D<HexTriLinear, Hex8ptQuadrature> Basis;
  typedef ElasticityPDE<I, T> PDE;

  const index_t nx = 64;
  const index_t ny = 64;
  const index_t nz = 64;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);

  FEModel<I, T, PDE> model(nnodes, nbcs);
  LinElasticityElement3D<I, T, Basis> element(nelems);

  model.add_element(&element);

  // Set the boundary conditions
  auto bcs = model.get_bcs();
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

  // Set the connectivity
  auto conn = element.get_conn();
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

  // Set the node locations
  auto X = model.get_nodes();
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

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model.set_nodes(X);

  // Set the element data
  auto data = element.get_quad_data();
  double pt_data[] = {1.23, 2.45};
  for (A2D::index_t i = 0; i < data.extent(0); i++) {
    for (A2D::index_t j = 0; j < data.extent(1); j++) {
      for (A2D::index_t k = 0; k < data.extent(2); k++) {
        data(i, j, k) = pt_data[k];
      }
    }
  }

  // Compute the Jacobian matrix
  double t0 = MPI_Wtime();
  auto J = model.new_matrix();
  t0 = MPI_Wtime() - t0;
  std::cout << "Jacobian initialization time: " << t0 << std::endl;

  double t1 = MPI_Wtime();
  model.jacobian(*J);
  t1 = MPI_Wtime() - t1;
  std::cout << "Jacobian computational time: " << t1 << std::endl;

  double t2 = MPI_Wtime();
  int num_levels = 4;
  double omega = 1.333;
  bool print_info = true;
  auto amg = model.new_amg(num_levels, omega, J, print_info);
  t2 = MPI_Wtime() - t2;
  std::cout << "Set up time for AMG: " << t2 << std::endl;

  // Set the residuals and apply the boundary conditions
  auto solution = model.new_solution();
  auto residual = model.new_solution();
  solution->fill(1.0);
  residual->zero();
  BSRMatVecMult(*J, *solution, *residual);

  // Set the solution back to zero
  solution->zero();

  index_t monitor = 10;
  index_t max_iters = 80;
  double t3 = MPI_Wtime();
  amg->cg(*residual, *solution, monitor, max_iters);
  t3 = MPI_Wtime() - t3;
  std::cout << "Conjugate gradient solution time: " << t3 << std::endl;

  double t4 = MPI_Wtime();
  amg->mg(*residual, *solution, monitor, max_iters);
  t4 = MPI_Wtime() - t4;
  std::cout << "Multigrid solution time: " << t4 << std::endl;

  model.set_solution(*solution);
  T energy = model.energy();
  std::cout << "Model energy: " << energy << std::endl;

  // Allocate the stress functional
  StressIntegral3D<I, T, Basis> stress_functional(element, 1.0);

  Functional<I, T, PDE> functional;
  functional.add_functional(&stress_functional);

  std::cout << "Stress integral " << functional.eval_functional() << std::endl;

  return (0);
}
