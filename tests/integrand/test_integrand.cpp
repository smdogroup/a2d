#include <complex>
#include <iostream>
#include <vector>

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/integrand_elasticity.h"
// #include "multiphysics/integrand_heat_conduction.h"
// #include "multiphysics/integrand_poisson.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_test.h"

using namespace A2D;

void box(index_t b_nx = 5, index_t b_ny = 5, index_t b_nz = 5,
         double b_lx = 1.0, double b_ly = 1.0, double b_lz = 1.0) {
  using T = std::complex<double>;

  // helper functor
  auto node_num = [&](int i, int j, int k) {
    return i + j * (b_nx + 1) + k * (b_nx + 1) * (b_ny + 1);
  };

  // Compute number of vertices and hex elements
  index_t nverts = (b_nx + 1) * (b_ny + 1) * (b_nz + 1);
  index_t nhex = b_nx * b_ny * b_nz;

  // Allocate temporary arrays
  index_t *hex = new index_t[8 * nhex];
  T *Xloc = new T[3 * nverts];

  // Populate hex
  using ET = ElementTypes;
  for (int k = 0, e = 0; k < b_nz; k++) {
    for (int j = 0; j < b_ny; j++) {
      for (int i = 0; i < b_nx; i++, e++) {
        for (int ii = 0; ii < ET::HEX_NVERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  // Populate Xloc
  for (int k = 0; k < b_nz + 1; k++) {
    for (int j = 0; j < b_ny + 1; j++) {
      for (int i = 0; i < b_nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (b_lx * i) / b_nx;
        Xloc[3 * node_num(i, j, k) + 1] = (b_ly * j) / b_ny;
        Xloc[3 * node_num(i, j, k) + 2] = (b_lz * k) / b_nz;
      }
    }
  }

  // Create A2D's connectivity object
  index_t ntets = 0, nwedge = 0, npyrmd = 0;
  index_t *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  const index_t block_size = 3;
  using MatType = BSRMat<T, block_size, block_size>;
  using VecType = SolutionVector<T>;
  const int dim = 3;
  const int degree = 2;
  constexpr GreenStrainType etype = GreenStrainType::LINEAR;
  using Elem = HexTopoElement<T, dim, etype, degree, VecType, MatType>;
  auto data_mesh = std::make_shared<ElementMesh<Elem::DataBasis>>(conn);
  auto geo_mesh = std::make_shared<ElementMesh<Elem::GeoBasis>>(conn);
  auto sol_mesh = std::make_shared<ElementMesh<Elem::Basis>>(conn);

  T E = 70.0, nu = 0.3, q = 5.0;

  TopoElasticityIntegrand<T, dim, etype> integrand(E, nu, q);
  Elem element(integrand, data_mesh, geo_mesh, sol_mesh);

  // Create the matrix
  index_t nrows;
  std::vector<index_t> rowp, cols;
  sol_mesh->template create_block_csr<block_size>(nrows, rowp, cols);

  // Create the shared pointer
  SolutionVector<T> data(data_mesh->get_num_dof());
  SolutionVector<T> geo(geo_mesh->get_num_dof());
  SolutionVector<T> sol(sol_mesh->get_num_dof());
  SolutionVector<T> res(sol_mesh->get_num_dof());

  // Set the geometry
  ElementVector_Serial<T, Elem::GeoBasis, VecType> elem_geo(*geo_mesh, geo);
  set_geo_from_hex_nodes<Elem::GeoBasis>(nhex, hex, Xloc, elem_geo);

  auto mat = std::make_shared<MatType>(nrows, nrows, cols.size(), rowp, cols);
  element.add_residual(data, geo, sol, res);
  element.add_jacobian(data, geo, sol, *mat);
}

template <typename T>
using TopoIntegrand = TopoElasticityIntegrand<T, 3, GreenStrainType::LINEAR>;

template <typename T>
using NonlinearTopoIntegrand =
    TopoElasticityIntegrand<T, 3, GreenStrainType::NONLINEAR>;

// IntegrandTopoLinearElasticity<T, 3>;

// template <typename T>
// using PoissonIntegrand = Poisson<T, 3>;

// template <typename T>
// using MixedPoissonIntegrand = MixedPoisson<T, 3>;

// template <typename T>
// using HeatIntegrand = HeatConduction<T, 3>;

bool TestIntegrands(bool component, bool write_output) {
  bool passed = true;

  using Tc = std::complex<double>;

  TopoIntegrand<Tc> integrand1(Tc(0.7), Tc(0.3), Tc(5.0));
  A2D::Test::A2DIntegrandAnalysisTest<TopoIntegrand, Tc> test1(integrand1);
  passed = passed && A2D::Test::Run(test1, component, write_output);

  NonlinearTopoIntegrand<Tc> integrand2(Tc(0.7), Tc(0.3), Tc(5.0));
  A2D::Test::A2DIntegrandAnalysisTest<NonlinearTopoIntegrand, Tc> test2(
      integrand2);
  passed = passed && A2D::Test::Run(test2, component, write_output);

  // PoissonIntegrand<Tc> integrand2;
  // A2D::Test::A2DIntegrandAnalysisTest<PoissonIntegrand, Tc>
  // test2(integrand2); passed = passed && A2D::Test::Run(test2, component,
  // write_output);

  // MixedPoissonIntegrand<Tc> integrand3;
  // A2D::Test::A2DIntegrandAnalysisTest<MixedPoissonIntegrand, Tc> test3(
  //     integrand3);
  // passed = passed && A2D::Test::Run(test3, component, write_output);

  // HeatIntegrand<Tc> integrand4(Tc(0.7), Tc(0.3), Tc(5.0));
  // A2D::Test::A2DIntegrandAnalysisTest<HeatIntegrand, Tc> test4(integrand4);
  // passed = passed && A2D::Test::Run(test4, component, write_output);

  return passed;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  bool component = false;     // Default to a projection test
  bool write_output = false;  // Don't write output;

  // Check for the write_output flag
  for (int i = 0; i < argc; i++) {
    std::string str(argv[i]);
    if (str.compare("--write_output") == 0) {
      write_output = true;
    }
    if (str.compare("--component") == 0) {
      component = true;
    }
  }

  box();

  typedef std::function<bool(bool, bool)> TestFunc;
  std::vector<TestFunc> tests;

  tests.push_back(TestIntegrands);

  bool passed = true;
  for (int i = 0; i < tests.size(); i++) {
    bool test_passed = tests[i](component, write_output);

    // If the test fails, perform a component test and print the results
    // to screen
    if (!test_passed) {
      bool comp = true;
      bool write = true;
      tests[i](comp, write);
    }

    passed = test_passed && passed;
  }

  // Convert to a fail flag
  int fail = 0;
  if (!passed) {
    fail = 1;
  }

  Kokkos::finalize();
  return fail;
}
