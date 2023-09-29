#include <complex>
#include <iostream>
#include <vector>

#include "a2dcore.h"
#include "ad/a2dtest.h"
#include "multiphysics/feanalysis.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"

using namespace A2D;

void box(index_t b_nx = 20, index_t b_ny = 20, index_t b_nz = 20,
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

  index_t num_boundary_verts = (b_ny + 1) * (b_nz + 1);
  index_t *boundary_verts = new index_t[num_boundary_verts];
  for (int k = 0, count = 0; k < b_nz + 1; k++) {
    for (int j = 0; j < b_ny + 1; j++, count++) {
      int i = 0;
      boundary_verts[count] = node_num(i, j, k);
    }
  }

  // Create A2D's connectivity object
  index_t ntets = 0, nwedge = 0, npyrmd = 0;
  index_t *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  // Label the element boundaries on which to apply boundary conditions
  DirichletBCInfo bcinfo;
  index_t bc_label =
      conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);
  bcinfo.add_boundary_condition(bc_label);

  // Set up the type of implementation we're using
  const index_t degree = 1;
  const index_t dim = 3;
  const index_t block_size = dim;
  using Impl_t = typename DirectCholeskyAnalysis<T, block_size>::Impl_t;
  using Vec_t = typename Impl_t::Vec_t;

  // Set the type of element
  constexpr GreenStrainType etype = GreenStrainType::LINEAR;
  using HexElem = HexTopoElement<Impl_t, dim, etype, degree>;

  // Create the meshes for the Hex elements
  auto data_mesh = std::make_shared<ElementMesh<HexElem::DataBasis>>(conn);
  auto geo_mesh = std::make_shared<ElementMesh<HexElem::GeoBasis>>(conn);
  auto sol_mesh = std::make_shared<ElementMesh<HexElem::Basis>>(conn);

  // Get the number of different dof
  index_t ndata = data_mesh->get_num_dof();
  index_t ngeo = geo_mesh->get_num_dof();
  index_t ndof = sol_mesh->get_num_dof();
  auto data = std::make_shared<Vec_t>(ndata);
  auto geo = std::make_shared<Vec_t>(ngeo);
  auto sol = std::make_shared<Vec_t>(ndof);
  auto res = std::make_shared<Vec_t>(ndof);

  // Set the data
  for (index_t i = 0; i < ndata; i++) {
    (*data)[i] = 1.0;
  }

  // Create the element
  T E = 70.0, nu = 0.3, q = 5.0;
  T design_stress = 1.0, ks_param = 10.0;
  TopoElasticityIntegrand<T, dim, etype> elem_integrand(E, nu, q);

  auto assembler = std::make_shared<ElementAssembler<Impl_t>>();
  assembler->add_element(
      std::make_shared<HexElem>(elem_integrand, data_mesh, geo_mesh, sol_mesh));

  // Set up the boundary conditions
  auto bcs = std::make_shared<DirichletBCs<T>>();
  bcs->add_bcs(std::make_shared<DirichletBasis<T, HexElem::Basis>>(
      conn, *sol_mesh, bcinfo, 0.0));

  // Create the assembler object
  DirectCholeskyAnalysis<T, block_size> analysis(data, geo, sol, res, assembler,
                                                 bcs);

  // Set the geometry
  typename Impl_t::ElementVector<HexElem::GeoBasis> elem_geo(
      *geo_mesh, *analysis.get_geo());
  set_geo_from_hex_nodes<HexElem::GeoBasis>(nhex, hex, Xloc, elem_geo);

  analysis.linear_solve();

  using HexFunc = HexTopoVonMises<Impl_t, dim, etype, degree>;
  TopoVonMisesKS<T, dim, etype> func_integrand(E, nu, q, design_stress,
                                               ks_param);
  HexFunc functional(func_integrand, data_mesh, geo_mesh, sol_mesh);

  constexpr FEVarType wrt = FEVarType::GEOMETRY;
  typename Impl_t::Vec_t dfdx(ngeo);
  T value = analysis.evaluate(functional);
  analysis.eval_adjoint_derivative(functional, wrt, dfdx);
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();
  int fail = 0;

  box();

  Kokkos::finalize();
  return fail;
}
