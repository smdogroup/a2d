#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2d.h"

using namespace A2D;

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

  double lx = 1.0, ly = 1.0;
  MesherRect2D mesher(nx, ny, lx, ly);

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
  mesher.set_bcs(bcs);

  // Set the connectivity
  auto conn = element->get_conn();
  auto conn_helm = element_helm->get_conn();

  // Set the node locations
  auto X = model->get_nodes();
  auto X_helm = model_helm->get_nodes();

  mesher.set_X_conn(X, conn);
  mesher.set_X_conn(X_helm, conn_helm);

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
  mesher.set_dv(*x);

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

  mesher.set_force(model, residual);

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
  amg->mg(*dfdu, *adjoint, monitor, max_iters);

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

  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);

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
