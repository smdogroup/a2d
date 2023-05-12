#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cassert>

#include "a2d.h"
#include "amgx_solver.h"
#include "toolkit.h"

using namespace A2D;

void block_element_random_update(int block_dimx, int block_dimy, int i, void* v,
                                 AMGX_Mode mode) {
  int j1, j2, block_size;
  block_size = block_dimy * block_dimx;

  for (j1 = 0; j1 < block_dimy; j1++) {
    for (j2 = 0; j2 < block_dimx; j2++) {
      if ((AMGX_GET_MODE_VAL(AMGX_MatPrecision, mode) == AMGX_matFloat)) {
        float* t = (float*)v;
        t[i * block_size + j1 * block_dimx + j2] *=
            (1.0f + (rand() % 10) * (1e-6f));
      } else {
        double* t = (double*)v;
        t[i * block_size + j1 * block_dimx + j2] *=
            (1.0 + (rand() % 10) * (1e-12));
      }
    }
  }
}

typedef index_t I;
typedef double T;
static const int SPATIAL_DIM = 3;
typedef BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature> Basis;
typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> ElasticityPDE;

template <class ConnArray, class XArray>
void populate_X_conn_3d(const int nx, const int ny, const int nz, XArray& X,
                        ConnArray& conn) {
  // Set X
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

  // Set connectivity
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
}

template <class BcsArray>
void populate_bcs_3d(const int nx, const int ny, const int nz, BcsArray& bcs) {
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
}

template <class DvArray>
void set_dv_3d(const int nx, const int ny, const int nz, DvArray& x) {
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);
        x(node, 0) = 1.0;
        // if (i < nx / 2 and j < ny / 2 and k < nz / 2) {
        //   x(node, 0) = 1e-3;
        // } else {
        //   x(node, 0) = 1.0;
        // }
      }
    }
  }
}

template <class Model, class RhsArray>
void set_force_3d(const int nx, const int ny, const int nz, Model model,
                  RhsArray& residual) {
  residual->zero();
  for (int k = nz / 4; k < 3 * nz / 4; k++) {
    int node = nx + (nx + 1) * (0 + (ny + 1) * k);
    (*residual)(node, 1) = -1e2;
  }
  model->zero_bcs(residual);
}

int main(int argc, char** argv) {
  const index_t nx = 16;
  const index_t ny = 64;
  const index_t nz = 128;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);

  std::shared_ptr<FEModel<I, T, ElasticityPDE>> model =
      std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  populate_bcs_3d(nx, ny, nz, bcs);

  // Set the connectivity
  auto conn = element->get_conn();

  // Set the node locations
  auto X = model->get_nodes();

  populate_X_conn_3d(nx, ny, nz, X, conn);

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();

  // Set the element
  T q = 5.0, E = 70e3, nu = 0.3;
  T density = 1.0, design_stress = 1e3;
  auto constitutive = std::make_shared<TopoIsoConstitutive<I, T, Basis>>(
      element, q, E, nu, density, design_stress);
  model->add_constitutive(constitutive);

  // Create the design vector
  A2D::CLayout<1> design_layout(model->nnodes);
  auto x = std::make_shared<A2D::MultiArray<T, A2D::CLayout<1>>>(design_layout);

  // Set the design variable values
  // set_dv_3d(nx, ny, nz, *x);
  x->fill(1.0);
  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, ElasticityPDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // ================================================================
  // Set up input parameters for the solver
  // Compute the Jacobian matrix
  auto J = model->new_matrix();
  model->jacobian(J);
  // Set the residuals and apply the boundary conditions
  auto residual = model->new_solution();
  // set_force_3d(nx, ny, nz, model, residual);
  // Initialize the solution vector
  auto solution = model->new_solution();
  solution->fill(1.0);
  model->zero_bcs(solution);
  residual->zero();
  BSRMatVecMult(*J, *solution, *residual);
  model->zero_bcs(residual);
  solution->zero();

  // ================================================================
  // Apply the preconditioned conjugate gradient method.
  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.0;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);
  index_t monitor = 50;
  index_t max_iters = 1000;
  // Compute the solution
  // TIMER("MG");
  // amg->mg(*residual, *solution, monitor, max_iters);
  // TIMER("MG");

  // ================================================================
  // Apply multigrid method.
  solution->zero();
  TIMER("PCG");
  monitor = 10;
  max_iters = 100;
  amg->cg(*residual, *solution, monitor, max_iters);
  TIMER("PCG");

  // ================================================================
  // Compute the solution using AMGX
  solution->zero();
  amgxSolver(J, residual, solution, argc, argv);

  // ================================================================
  // print the solution
  // for (int i = 0; i < 8; i++) {
  //   if (((double*)solution->data)[i] < 1e-6) {
  //     ((double*)solution->data)[i] = 0.0;
  //   }
  //   std::cout << " " << ((double*)solution->data)[i];
  // }
  // std::cout << std::endl;

  // ================================================================
  // ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())>
  // vtk(conn,
  //                                                                        X);
  // vtk.write_mesh();
  // vtk.write_sol("x", *x, 0);
  // vtk.write_sol("ux", *solution, 0);
  // vtk.write_sol("uy", *solution, 1);
  // vtk.write_sol("uz", *solution, 2);

  return 0;
}
