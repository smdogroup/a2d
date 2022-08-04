#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2d.h"

using namespace A2D;

typedef index_t I;
typedef double T;

void main_body(int argc, char* argv[]) {
  // Define problem dimension
  static const int SPATIAL_DIM = 3;
  const index_t nx = 160;
  const index_t ny = 16;
  const index_t nz = 16;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const index_t nelems = nx * ny * nz;
  const index_t nbcs = (ny + 1) * (nz + 1);
  typedef MesherBrick3D<nx, ny, nz> Mesher;
  typedef BasisOps<SPATIAL_DIM, HexTriLinearBasisFunc, Hex8ptQuadrature> Basis;
  typedef ElasticityPDEInfo<SPATIAL_DIM, I, T> ElasticityPDE;

  // Set PDE model and element
  auto model = std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the boundary conditions
  auto bcs = model->get_bcs();
  Mesher::set_bcs(bcs);

  // Set the connectivity
  auto conn = element->get_conn();

  // Set the node locations
  auto X = model->get_nodes();
  Mesher::set_X_conn(X, conn, 10.0, 1.0, 1.0);

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
  Mesher::set_dv(*x);
  model->set_design_vars(x);

  // Set up the stress functional
  auto functional = std::make_shared<Functional<I, T, ElasticityPDE>>();
  auto agg_functional =
      std::make_shared<TopoVonMisesAggregation<I, T, Basis>>(constitutive);
  functional->add_functional(agg_functional);

  // Compute the Jacobian matrix
  auto J = model->new_matrix();
  model->jacobian(J);

  int num_levels = 3;
  double omega = 0.6667;
  double epsilon = 0.01;
  bool print_info = true;
  auto amg = model->new_amg(num_levels, omega, epsilon, J, print_info);

  // Set the residuals and apply the boundary conditions
  auto solution = model->new_solution();
  auto residual = model->new_solution();

  Mesher::set_force(model, residual);

  // Compute the solution
  index_t monitor = 10;
  index_t max_iters = 3000;
  solution->zero();
  amg->cg(*residual, *solution, monitor, max_iters);

  ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(conn,
                                                                         X);
  vtk.write_mesh();
  vtk.write_sol("x", *x, 0);
  vtk.write_sol("ux", *solution, 0);
  vtk.write_sol("uy", *solution, 1);
  vtk.write_sol("uz", *solution, 2);
}

int main(int argc, char* argv[]) {
#ifdef A2D_USE_KOKKOS
  Kokkos::initialize();
#endif

  Timer main_timer("main");
  main_body(argc, argv);

#ifdef A2D_USE_KOKKOS
  Kokkos::finalize();
#endif
  return (0);
}
