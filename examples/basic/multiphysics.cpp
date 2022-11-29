#include <iostream>
#include <memory>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"

void test_febasis() {
  using T = double;
  const A2D::index_t degree = 2;
  const A2D::index_t dim = 3;

  using Quadrature = A2D::HexGaussQuadrature<degree + 1>;
  using Space =
      A2D::FESpace<T, dim, A2D::HdivSpace<T, dim>, A2D::H1Space<T, dim, dim>>;
  using Basis = A2D::FEBasis<T, A2D::QHdivHexBasis<T, degree>,
                             A2D::LagrangeH1HexBasis<T, dim, degree>>;
  const A2D::index_t ncomp = Space::ncomp;
  using MatType = A2D::Mat<T, Basis::ndof, Basis::ndof>;
  using QMatType = A2D::Mat<T, ncomp, ncomp>;

  // Generate random data
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> distr(-1.0, 1.0);

  // Set the degrees of freedom and basis
  T dof[Basis::ndof], res[Basis::ndof], result[Basis::ndof];
  for (A2D::index_t i = 0; i < Basis::ndof; i++) {
    dof[i] = distr(gen);
    res[i] = 0.0;
    result[i] = 0.0;
  }

  Space s, p;
  MatType mat;
  QMatType qmat;

  for (A2D::index_t i = 0; i < ncomp; i++) {
    for (A2D::index_t j = 0; j < ncomp; j++) {
      qmat(i, j) = distr(gen);
    }
  }

  A2D::index_t pt = degree + 4;
  Basis::add_outer<Quadrature>(pt, qmat, mat);

  Basis::interp_basis<Quadrature>(pt, dof, s);
  for (A2D::index_t i = 0; i < ncomp; i++) {
    p[i] = 0.0;
    for (A2D::index_t j = 0; j < ncomp; j++) {
      p[i] += qmat(i, j) * s[j];
    }
  }
  Basis::add_basis<Quadrature>(pt, p, res);

  for (A2D::index_t i = 0; i < Basis::ndof; i++) {
    for (A2D::index_t j = 0; j < Basis::ndof; j++) {
      result[i] += mat(i, j) * dof[j];
    }
  }

  std::cout << std::setw(15) << "add_outer " << std::setw(15) << "basis"
            << std::setw(15) << "rel_err" << std::endl;
  for (A2D::index_t i = 0; i < Basis::ndof; i++) {
    std::cout << std::setw(15) << result[i] << std::setw(15) << res[i]
              << std::setw(15) << (result[i] - res[i]) / result[i] << std::endl;
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  test_febasis();

  using T = double;
  const A2D::index_t dim = 3;

  std::cout << "Poisson\n";
  A2D::Poisson<std::complex<T>, dim> poisson;
  A2D::TestPDEImplementation<std::complex<T>>(poisson);

  std::cout << "Mixed Poisson\n";
  A2D::MixedPoisson<std::complex<T>, dim> mixed_poisson;
  A2D::TestPDEImplementation<std::complex<T>>(mixed_poisson);

  std::cout << "Topology linear elasticity\n";
  T E = 70e3, nu = 0.3, q = 5.0;
  A2D::TopoLinearElasticity<std::complex<T>, dim> elasticity(E, nu, q);
  A2D::TestPDEImplementation<std::complex<T>>(elasticity);

  std::cout << "Heat conduction\n";
  A2D::HeatConduction<std::complex<T>, dim> heat_conduction;
  A2D::TestPDEImplementation<std::complex<T>>(heat_conduction);

  std::cout << "Mixed heat conduction\n";
  A2D::MixedHeatConduction<std::complex<T>, dim> mixed_heat_conduction;
  A2D::TestPDEImplementation<std::complex<T>>(mixed_heat_conduction);

  return (0);
}