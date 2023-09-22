#include <complex>
#include <iostream>
#include <vector>

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/integrand_elasticity.h"
// #include "multiphysics/integrand_heat_conduction.h"
// #include "multiphysics/integrand_poisson.h"
#include "multiphysics/feanalysis.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_test.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"

using namespace A2D;

// void solve() {
//   // Compute the residual and zero the boundary conditions
//   res.zero();
//   assembler->add_residual(data, geo, sol, res);
//   zero_bcs(res);

//   // Solve for the update
//   factor();
//   chol->solve(res.data());

//   // Update the solution
//   for (index_t i = 0; i < sol.get_num_dof(); i++) {
//     sol[i] -= res[i];
//   }
// }

// void tovtk(const std::string filename) {
//   if constexpr (spatial_dim == 2) {
//     write_quad_to_vtk<3, degree, T, DataBasisElas, GeoBasisElas, BasisElas,
//                       IntegrandElas>(
//         elem_data_elas, elem_geo_elas, elem_sol_elas, filename,
//         [](index_t k, typename IntegrandElas::DataSpace &d,
//            typename IntegrandElas::FiniteElementGeometry &g,
//            typename IntegrandElas::FiniteElementSpace &s) {
//           if (k == 2) {  // write data
//             return (d.template get<0>()).get_value();
//           } else {  // write solution components
//             auto u = (s.template get<0>()).get_value();
//             return u(k);
//           }
//         });
//   } else {  // spatial_dim == 3
//     write_hex_to_vtk<4, degree, T, DataBasisElas, GeoBasisElas, BasisElas,
//                      IntegrandElas>(
//         elem_data_elas, elem_geo_elas, elem_sol_elas, filename,
//         [](index_t k, typename IntegrandElas::DataSpace &d,
//            typename IntegrandElas::FiniteElementGeometry &g,
//            typename IntegrandElas::FiniteElementSpace &s) {
//           if (k == 3) {  // write data
//             return (d.template get<0>()).get_value();
//           } else {  // write solution components
//             auto u = (s.template get<0>()).get_value();
//             return u(k);
//           }
//         });
//   }
// }
// };

void box(index_t b_nx = 5, index_t b_ny = 5, index_t b_nz = 5,
         double b_lx = 1.0, double b_ly = 1.0, double b_lz = 1.0) {
  using T = double;  // std::complex<double>;

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
  using BSRMat_t = BSRMat<T, block_size, block_size>;
  using Vec_t = SolutionVector<T>;
  const int dim = 3;
  const int degree = 2;
  constexpr GreenStrainType etype = GreenStrainType::LINEAR;
  using HexElem = HexTopoElement<T, dim, etype, degree, Vec_t, BSRMat_t>;

  // Create the meshes for the Hex elements
  auto data_mesh = std::make_shared<ElementMesh<HexElem::DataBasis>>(conn);
  auto geo_mesh = std::make_shared<ElementMesh<HexElem::GeoBasis>>(conn);
  auto sol_mesh = std::make_shared<ElementMesh<HexElem::Basis>>(conn);

  // Get the sizes of the different meshes
  index_t ndata = data_mesh->get_num_dof();
  index_t ngeo = data_mesh->get_num_dof();
  index_t ndof = data_mesh->get_num_dof();

  // Create the element
  T E = 70.0, nu = 0.3, q = 5.0;
  T design_stress = 1.0, ks_param = 10.0;
  TopoElasticityIntegrand<T, dim, etype> elem_integrand(E, nu, q);

  auto assembler = std::make_shared<ElementAssembler<Vec_t, BSRMat_t>>();
  assembler->add_element(
      std::make_shared<HexElem>(elem_integrand, data_mesh, geo_mesh, sol_mesh));

  // Create the assembler object
  DirectCholeskyAnalysis<T, block_size> analysis(ndata, ngeo, ndof, assembler);

  // Set the geometry
  ElementVector_Serial<T, HexElem::GeoBasis, Vec_t> elem_geo(
      *geo_mesh, analysis.get_geo());
  set_geo_from_hex_nodes<HexElem::GeoBasis>(nhex, hex, Xloc, elem_geo);

  analysis.factor();

  // using HexFunc = HexTopoVonMises<T, dim, etype, degree, Vec_t>;
  // Analysis<T> analysis(assembler, bsr_mat);

  // TopoVonMisesKS<T, dim, etype> func_integrand(E, nu, q, design_stress,
  //                                              ks_param);
  // HexFunc functional(func_integrand, data_mesh, geo_mesh, sol_mesh);

  // assembler.add_residual(data, geo, sol, res);
  // assembler.add_jacobian(data, geo, sol, *mat);

  // T value = functional.evaluate(data, geo, sol);
  // functional.add_derivative(FEVarType::STATE, data, geo, sol, res);
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
