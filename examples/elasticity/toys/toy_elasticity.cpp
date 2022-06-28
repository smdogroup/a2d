#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dprofiler.h"
#include "a2dtmp3d.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "model.h"

using namespace A2D;

int main(int argc, char* argv[]) {
  static const index_t SPATIAL_DIM = 2;
  typedef index_t I;
  typedef std::complex<double> T;
  // typedef std::complex<double> T;
  typedef BasisOps<SPATIAL_DIM, QuadBiLinearBasisFunc, Quad4ptQuadrature> Basis;
  typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> PDE;

  const index_t nx = 32;
  const index_t ny = 32;
  const index_t nnodes = (nx + 1) * (ny + 1);
  const index_t nelems = nx * ny;
  const index_t nbcs = (ny + 1);

  auto model = std::make_shared<FEModel<I, T, PDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  index_t index = 0;
  for (int j = 0; j < ny + 1; j++) {
    int i = 0;
    int node = i + (nx + 1) * j;

    // Set the boundary conditions
    bcs(index, 0) = node;
    for (int ii = 0; ii < SPATIAL_DIM; ii++) {
      bcs(index, 1) |= 1U << ii;
    }
    index++;
  }

  // Set the connectivity
  auto conn = element->get_conn();
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

  // Set the node locations
  auto X = model->get_nodes();
  for (int j = 0; j < ny + 1; j++) {
    for (int i = 0; i < nx + 1; i++) {
      int node = i + (nx + 1) * j;

      X(node, 0) = 1.0 * i / nx;
      X(node, 1) = 1.0 * j / ny;
    }
  }

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
  x->fill(1.0);
  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, PDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // Compute the Jacobian matrix
  Timer* t;
  t = new Timer("Jacobian initialization");
  auto J = model->new_matrix();
  delete t;

  t = new Timer("Jacobian computation");
  model->jacobian(J);
  delete t;

  t = new Timer("AMG Set up");
  int num_levels = 3;
  double omega = 1.333;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, J, print_info);
  delete t;

  // Set the residuals and apply the boundary conditions
  auto solution = model->new_solution();
  auto residual = model->new_solution();
  solution->fill(1.0);
  model->zero_bcs(solution);

  residual->zero();
  BSRMatVecMult(*J, *solution, *residual);
  model->zero_bcs(residual);

  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 80;
  t = new Timer("CG solution");
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);
  delete t;

  // Set the solution
  model->set_solution(solution);

  agg_functional->compute_offset();
  functional->eval_functional();

  // Compute the adjoint right-hand-side
  auto dfdu = model->new_solution();
  functional->eval_dfdu(dfdu);
  dfdu->scale(-1.0);
  model->zero_bcs(dfdu);

  // Compute the adjoint variables
  auto adjoint = model->new_solution();
  t = new Timer("Adjoint solution");
  amg->mg(*dfdu, *adjoint, monitor, max_iters);
  delete t;

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

  // Set the new desigh variable values
  double dh = 1e-30;
  x->fill(1.0);
  x->axpy(T(0.0, dh), *px);
  model->set_design_vars(x);
  model->jacobian(J);
  amg->update();

  t = new Timer("CG solution");
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);
  delete t;

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

  return (0);
}
