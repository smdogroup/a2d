#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp2d.h"
#include "a2dtmp3d.h"
#include "a2dvtk.h"
#include "elasticity.h"
#include "helmholtz.h"
#include "model.h"
#include "mpi.h"

using namespace A2D;

template <class ConnArray, class XArray>
void populate_X_conn_2d(const int nx, const int ny, XArray& X,
                        ConnArray& conn) {
  // Set X
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      int node = i + (nx + 1) * j;

      X(node, 0) = 1.0 * i / nx;
      X(node, 1) = 1.0 * j / ny;
    }
  }

  // Set connectivity
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int elem = i + nx * j;

      int conn_coord[4];
      for (int jj = 0, index = 0; jj < 2; jj++) {
        for (int ii = 0; ii < 2; ii++, index++) {
          conn_coord[index] = (i + ii) + (nx + 1) * (j + jj);
        }
      }

      // Convert to the correct connectivity
      conn(elem, 0) = conn_coord[0];
      conn(elem, 1) = conn_coord[1];
      conn(elem, 2) = conn_coord[3];
      conn(elem, 3) = conn_coord[2];
    }
  }
}

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
void populate_bcs_2d(const int nx, const int ny, BcsArray& bcs) {
  index_t index = 0;
  for (int j = 0; j < ny + 1; j++) {
    int i = 0;
    int node = i + (nx + 1) * j;

    // Set the boundary conditions
    bcs(index, 0) = node;
    for (int ii = 0; ii < 2; ii++) {
      bcs(index, 1) |= 1U << ii;
    }
    index++;
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
void set_dv_2d(const int nx, const int ny, DvArray& x) {
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      int node = i + (nx + 1) * j;
      if (i < nx / 2 and j < ny / 2) {
        x(node, 0) = 1e-3;
      } else {
        x(node, 0) = 1.0;
      }
    }
  }
}

template <class DvArray>
void set_dv_3d(const int nx, const int ny, const int nz, DvArray& x) {
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);
        if (i < nx / 2 and j < ny / 2 and k < nz / 2) {
          x(node, 0) = 1e-3;
        } else {
          x(node, 0) = 1.0;
        }
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
    residual(node, 1) = -1e2;
  }
  model->zero_bcs(residual);
}

template <class Model, class RhsArray>
void set_force_2d(const int nx, const int ny, Model& model,
                  RhsArray& residual) {
  residual->zero();
  (*residual)(nx, 1) = -1e2;
  model->zero_bcs(residual);
}

int main(int argc, char* argv[]) {
  static const int SPATIAL_DIM = 2;
  typedef index_t I;
  typedef double T;
  // typedef BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature>
  // Basis;
  typedef BasisOps<SPATIAL_DIM, QuadBiLinearBasisFunc, Quad4ptQuadrature> Basis;
  typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> ElasticityPDE;
  typedef HelmholtzPDEInfo<SPATIAL_DIM, I, T> HelmholtzPDE;

  const index_t nx = 32;
  const index_t ny = 32;
  const index_t nnodes = (nx + 1) * (ny + 1);
  const index_t nelems = nx * ny;
  const index_t nbcs = (ny + 1);

  auto model = std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  double r0 = 0.1;
  auto model_helm = std::make_shared<FEModel<I, T, HelmholtzPDE>>(nnodes, 0);
  auto element_helm =
      std::make_shared<HelmholtzElement<I, T, Basis>>(nelems, r0);
  model_helm->add_element(element_helm);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  populate_bcs_2d(nx, ny, bcs);

  // Set the connectivity
  auto conn = element->get_conn();
  auto conn_helm = element_helm->get_conn();

  // Set the node locations
  auto X = model->get_nodes();
  auto X_helm = model_helm->get_nodes();

  populate_X_conn_2d(nx, ny, X, conn);
  populate_X_conn_2d(nx, ny, X_helm, conn_helm);

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();
  model_helm->init();

  // Set the element
  T q = 5.0, E = 70e3, nu = 0.3;
  T density = 1.0, design_stress = 1e3;
  auto constitutive = std::make_shared<TopoIsoConstitutive<I, T, Basis>>(
      element, q, E, nu, density, design_stress);
  model->add_constitutive(constitutive);

  auto constitutive_helm =
      std::make_shared<HelmholtzConstitutive<I, T, Basis>>(element_helm);
  model_helm->add_constitutive(constitutive_helm);

  // Create the design vector
  A2D::CLayout<1> design_layout(model->nnodes);
  auto x = std::make_shared<A2D::MultiArray<T, A2D::CLayout<1>>>(design_layout);

  // Set the design variable values
  set_dv_2d(nx, ny, *x);

  model->set_design_vars(x);
  model_helm->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, ElasticityPDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // Compute the Jacobian matrix
  auto J = model->new_matrix();
  model->jacobian(J);

  auto Jh = model_helm->new_matrix();
  model_helm->jacobian(Jh);

  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.0;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);

  auto amg_helm =
      model_helm->new_amg(num_levels, omega, epsilon, Jh, print_info);

  // Set the residuals and apply the boundary conditions
  auto solution = model->new_solution();
  auto residual = model->new_solution();

  set_force_2d(nx, ny, model, residual);

  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 80;
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);

  auto rho = model_helm->new_solution();
  auto residual_helm = model_helm->new_solution();
  model_helm->residual(residual_helm);
  rho->zero();
  amg_helm->cg(*residual_helm, *rho, monitor, max_iters);

  // Set the solution
  model->set_solution(solution);

  agg_functional->compute_offset();
  functional->eval_functional();

#if 0
  // Compute the adjoint right-hand-side
  auto dfdu = model->new_solution();
  functional->eval_dfdu(dfdu);
  dfdu->scale(-1.0);
  model->zero_bcs(dfdu);

  // Compute the adjoint variables
  auto adjoint = model->new_solution();
  double t4 = MPI_Wtime();
  amg->mg(*dfdu, *adjoint, monitor, max_iters);
  t4 = MPI_Wtime() - t4;
  std::cout << "Adjoint solution time: " << t3 << std::endl;

  // Complete the adjoint derivative
  auto dfdx =
      std::shared_ptr<A2D::MultiArray<T, A2D::CLayout<1>>>(x->duplicate());
  dfdx->zero();
  functional->eval_dfdx(dfdx);
  model->add_adjoint_dfdx(adjoint, dfdx);

  // Compute a projected derivative and test against complex step
  auto px =
      std::shared_ptr<A2D::MultiArray<T, A2D::CLayout<1>>>(x->duplicate());
  px->random();
  T ans = dfdx->dot(*px);

  // Set the new design variable values
  double dh = 1e-30;
  x->fill(1.0);
  x->axpy(T(0.0, dh), *px);
  model->set_design_vars(x);
  model->jacobian(J);
  amg->update();

  double t5 = MPI_Wtime();
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);
  t5 = MPI_Wtime() - t5;
  std::cout << "Conjugate gradient solution time: " << t5 << std::endl;

  // Set the solution
  model->set_solution(solution);

  // Compute the complex-step result
  T fd = functional->eval_functional().imag() / dh;

  std::cout << "Complex-step result: " << std::setw(20) << std::setprecision(16)
            << fd.real() << std::endl;
  std::cout << "Adjoint-data result: " << std::setw(20) << std::setprecision(16)
            << ans.real() << std::endl;
  std::cout << "Relative error:      " << std::setw(20) << std::setprecision(16)
            << (ans.real() - fd.real()) / ans.real() << std::endl;
#endif

  ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(conn,
                                                                         X);
  vtk.write_mesh();
  vtk.write_sol("x", *x, 0);
  vtk.write_sol("rho", *rho, 0);
  vtk.write_sol("ux", *solution, 0);
  vtk.write_sol("uy", *solution, 1);

  return (0);
}
