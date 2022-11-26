#include <iostream>
#include <memory>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dvtk.h"

using namespace A2D;

template <class GeoBasis, typename I, typename T, class GeoElemVec>
void set_geo_from_hex_nodes(const index_t nhex, const I hex[], const T Xloc[],
                            GeoElemVec &elem_geo) {
  for (int e = 0; e < nhex; e++) {
    // Get the geometry values
    typename GeoElemVec::FEDof geo_dof(e, elem_geo);

    for (int ii = 0; ii < GeoBasis::ndof; ii++) {
      double pt[3];
      GeoBasis::get_dof_point(ii, pt);

      double N[8];
      N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
      N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

      // Interpolate to find the basis
      if (ii % 3 == 0) {
        T x = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          x += N[kk] * Xloc[3 * hex[8 * e + kk]];
        }
        geo_dof[ii] = x;
      } else if (ii % 3 == 1) {
        T y = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          y += N[kk] * Xloc[3 * hex[8 * e + kk] + 1];
        }
        geo_dof[ii] = y;
      } else if (ii % 3 == 2) {
        T z = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          z += N[kk] * Xloc[3 * hex[8 * e + kk] + 2];
        }
        geo_dof[ii] = z;
      }
    }

    elem_geo.set_element_values(e, geo_dof);
  }
}

template <index_t outputs, index_t degree, typename T, class DataBasis,
          class GeoBasis, class Basis, class PDE, class DataElemVec,
          class GeoElemVec, class ElemVec, class FunctorType>
void write_hex_to_vtk(PDE &pde, DataElemVec &elem_data, GeoElemVec &elem_geo,
                      ElemVec &elem_sol, const FunctorType &func) {
  using ET = ElementTypes;
  const index_t nex = degree;
  using QuadPts = HexGaussLobattoQuadrature<nex + 1>;

  const index_t nhex = elem_sol.get_num_elements();
  index_t nvtk_elems = nhex * nex * nex * nex;
  index_t nvtk_nodes = nhex * (nex + 1) * (nex + 1) * (nex + 1);

  auto vtk_node_num = [](index_t i, index_t j, index_t k) {
    return i + (nex + 1) * (j + (nex + 1) * k);
  };

  MultiArrayNew<int *[8]> vtk_conn("vtk_elems", nvtk_elems);
  MultiArrayNew<T *[3]> vtk_nodes("vtk_nodes", nvtk_nodes);
  MultiArrayNew<T *[outputs]> vtk_outputs("vtk_outputs", nvtk_nodes);

  for (index_t n = 0, counter = 0; n < nhex; n++) {
    // Get the data values
    typename DataElemVec::FEDof data_dof(n, elem_data);
    elem_data.get_element_values(n, data_dof);
    QptSpace<QuadPts, typename PDE::DataSpace> data;
    DataBasis::template interp(data_dof, data);

    // Get the geometry values
    typename GeoElemVec::FEDof geo_dof(n, elem_geo);
    elem_geo.get_element_values(n, geo_dof);
    QptSpace<QuadPts, typename PDE::FiniteElementGeometry> geo;
    GeoBasis::template interp(geo_dof, geo);

    // Get the degrees of freedom for the element
    typename ElemVec::FEDof sol_dof(n, elem_sol);
    elem_sol.get_element_values(n, sol_dof);
    QptSpace<QuadPts, typename PDE::FiniteElementSpace> sol;
    Basis::template interp(sol_dof, sol);

    const index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

    for (index_t k = 0; k < nex + 1; k++) {
      for (index_t j = 0; j < nex + 1; j++) {
        for (index_t i = 0; i < nex + 1; i++) {
          const index_t index = vtk_node_num(i, j, k);
          const index_t node = off + index;

          typename PDE::FiniteElementSpace &sref = sol.get(index);
          typename PDE::FiniteElementGeometry &gref = geo.get(index);

          // Get the Jacobian transformation
          A2D::Mat<T, PDE::dim, PDE::dim> &J =
              gref.template get<0>().get_grad();

          // Compute the inverse of the transformation
          A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
          A2D::MatInverse(J, Jinv);

          // Compute the determinant of the Jacobian matrix
          T detJ;
          A2D::MatDet(J, detJ);

          // Transform the solution the physical element
          typename PDE::FiniteElementSpace x, s;
          sref.transform(detJ, J, Jinv, s);

          auto X = gref.template get<0>().get_value();
          vtk_nodes(node, 0) = X(0);
          vtk_nodes(node, 1) = X(1);
          vtk_nodes(node, 2) = X(2);

          for (index_t kk = 0; kk < outputs; kk++) {
            vtk_outputs(node, kk) = func(kk, data.get(index), gref, s);
          }
        }
      }
    }

    for (index_t k = 0; k < nex; k++) {
      for (index_t j = 0; j < nex; j++) {
        for (index_t i = 0; i < nex; i++, counter++) {
          index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

          for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
            vtk_conn(counter, ii) =
                off + vtk_node_num(i + ET::HEX_VERTS_CART[ii][0],
                                   j + ET::HEX_VERTS_CART[ii][1],
                                   k + ET::HEX_VERTS_CART[ii][2]);
          }
        }
      }
    }
  }

  ToVTK vtk(vtk_conn, vtk_nodes);
  vtk.write_mesh();

  MultiArrayNew<T *> vtk_vec("vtk_vec", nvtk_nodes);
  for (index_t i = 0; i < outputs; i++) {
    for (index_t j = 0; j < nvtk_nodes; j++) {
      vtk_vec(j) = vtk_outputs(j, i);
    }
    char name[256];
    std::snprintf(name, sizeof(name), "solution%d", i + 1);
    vtk.write_sol(name, vtk_vec);
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const index_t dim = 3;
  using T = double;
  using ET = ElementTypes;
  // using PDE = MixedPoisson<T, dim>;
  using PDE = Poisson<T, dim>;

  const index_t degree = 8;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, degree>>;
  // using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
  //                       LagrangeL2HexBasis<T, 1, degree - 1>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  using DataElemVec = ElementVector_Serial<T, DataBasis>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis>;
  using ElemVec = ElementVector_Serial<T, Basis>;
  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Matrix-free operator for the problem
  using MatFree = MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  const index_t low_degree = 1;
  using LOrderQuadrature = HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = FEBasis<T>;
  using LOrderGeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, low_degree>>;
  // using LOrderBasis = FEBasis<T, QHdivHexBasis<T, low_degree>,
  //                             LagrangeL2HexBasis<T, 1, low_degree - 1>>;
  using LOrderBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, low_degree>>;
  using LOrderDataElemVec = ElementVector_Serial<T, LOrderDataBasis>;
  using LOrderGeoElemVec = ElementVector_Serial<T, LOrderGeoBasis>;
  using LOrderElemVec = ElementVector_Serial<T, LOrderBasis>;
  using LOrderFE = FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                 LOrderGeoBasis, LOrderBasis>;

  std::cout << "Poisson\n";
  Poisson<std::complex<T>, dim> poisson;
  TestPDEImplementation<std::complex<T>>(poisson);

  std::cout << "Mixed Poisson\n";
  MixedPoisson<std::complex<T>, dim> mixed_poisson;
  TestPDEImplementation<std::complex<T>>(mixed_poisson);

  std::cout << "Nonlinear elasticity\n";
  NonlinearElasticity<std::complex<T>, dim> elasticity;
  TestPDEImplementation<std::complex<T>>(elasticity);

  std::cout << "Heat conduction\n";
  HeatConduction<std::complex<T>, dim> heat_conduction;
  TestPDEImplementation<std::complex<T>>(heat_conduction);

  std::cout << "Mixed heat conduction\n";
  MixedHeatConduction<std::complex<T>, dim> mixed_heat_conduction;
  TestPDEImplementation<std::complex<T>>(mixed_heat_conduction);

  // Number of elements in each dimension
  const int nx = 4, ny = 4, nz = 4;
  auto node_num = [](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  // Number of edges
  const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = nx * ny * nz;

  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8 * nhex];

  for (int k = 0, e = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, e++) {
        for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  double Xloc[3 * nverts];
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (1.0 * i) / nx;
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j) / ny;
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k) / nz;
      }
    }
  }

  int boundary1_verts[(ny + 1) * (nz + 1)];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary1_verts[index] = node_num(0, j, k);
    }
  }

  int boundary2_verts[(ny + 1) * (nz + 1)];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary2_verts[index] = node_num(nx, j, k);
    }
  }

  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);
  ElementMesh<DataBasis> datamesh(conn);

  HexProjection<degree, Basis, LOrderBasis> basis_proj;
  HexProjection<degree, GeoBasis, LOrderGeoBasis> geo_proj;
  HexProjection<degree, DataBasis, LOrderDataBasis> data_proj;

  ElementMesh<LOrderBasis> lorder_mesh(mesh, basis_proj);
  ElementMesh<LOrderGeoBasis> lorder_geomesh(geomesh, geo_proj);
  ElementMesh<LOrderDataBasis> lorder_datamesh(datamesh, data_proj);

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select1[2] = {0};
  BoundaryCondition<Basis> bcs1(conn, mesh, basis_select1, (ny + 1) * (nz + 1),
                                boundary1_verts);

  // Set boundary conditions based on the vertex indices and finite-element
  // space
  index_t basis_select2[2] = {0};
  BoundaryCondition<Basis> bcs2(conn, mesh, basis_select2, (ny + 1) * (nz + 1),
                                boundary2_verts);

  std::cout << "Number of elements:            " << conn.get_num_elements()
            << std::endl;
  std::cout << "Number of degrees of freedom:  " << mesh.get_num_dof()
            << std::endl;

  // const index_t *bcs_index;
  // std::cout << "Number of boundary conditions: " << bcs1.get_bcs(&bcs_index)
  //           << std::endl;
  // std::cout << "Number of boundary conditions: " << bcs2.get_bcs(&bcs_index)
  //           << std::endl;

  PDE pde;

  index_t ndof = mesh.get_num_dof();
  SolutionVector<T> global_U(mesh.get_num_dof());
  SolutionVector<T> global_res(mesh.get_num_dof());
  SolutionVector<T> global_geo(geomesh.get_num_dof());
  SolutionVector<T> global_data(datamesh.get_num_dof());

  DataElemVec elem_data(datamesh, global_data);
  GeoElemVec elem_geo(geomesh, global_geo);
  ElemVec elem_sol(mesh, global_U);
  ElemVec elem_res(mesh, global_res);

  // Set the geometry from the node locations
  set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  SolutionVector<T> global_x(mesh.get_num_dof());
  SolutionVector<T> global_y(mesh.get_num_dof());
  ElemVec elem_x(mesh, global_x);
  ElemVec elem_y(mesh, global_y);

  // Create the finite-element model
  FE fe;
  MatFree matfree;

  // Add the residual
  elem_res.init_zero_values();
  fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);
  elem_res.add_values();

  fe.add_jacobian_vector_product(pde, elem_data, elem_geo, elem_sol, elem_x,
                                 elem_y);

  // Initialize the matrix-free data
  matfree.initialize(pde, elem_data, elem_geo, elem_sol);

  elem_y.init_zero_values();
  matfree.add_jacobian_vector_product(elem_x, elem_y);
  elem_y.add_values();

  LOrderFE lorder_fe;

  std::cout << "create_block_matrix" << std::endl;
  auto mat = lorder_mesh.create_block_matrix<T, 1>();
  BSRMat<index_t, T, 1, 1> &mat_ref = *mat;

  std::cout << "initialize element matrix" << std::endl;
  ElementMat_Serial<T, LOrderBasis, BSRMat<index_t, T, 1, 1>> elem_mat(
      lorder_mesh, mat_ref);

  // Set the low order element vectors
  LOrderDataElemVec lorder_elem_data(lorder_datamesh, global_data);
  LOrderGeoElemVec lorder_elem_geo(lorder_geomesh, global_geo);
  LOrderElemVec lorder_elem_sol(lorder_mesh, global_U);

  std::cout << "add_jacobian" << std::endl;
  lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                         lorder_elem_sol, elem_mat);

  const index_t *bc_dofs1;
  index_t nbcs1 = bcs1.get_bcs(&bc_dofs1);
  mat->zero_rows(nbcs1, bc_dofs1);

  // const index_t *bc_dofs2;
  // index_t nbcs2 = bcs1.get_bcs(&bc_dofs2);
  index_t nbcs2 = 1;
  index_t bc_dofs2[] = {global_x.get_num_dof() - 1};

  mat->zero_rows(nbcs2, bc_dofs2);

  MultiArrayNew<T *[1][1]> B("B", global_U.get_num_dof());
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    B(i, 0, 0) = 1.0;
  }

  for (index_t i = 0; i < nbcs1; i++) {
    B(bc_dofs1[i], 0, 0) = 0.0;
  }

  for (index_t i = 0; i < nbcs2; i++) {
    B(bc_dofs2[i], 0, 0) = 0.0;
  }

  index_t num_levels = 3;
  double omega = 3.0 / 4.0;
  double epsilon = 0.0;
  bool print_info = true;
  BSRMatAmg<index_t, T, 1, 1> amg(num_levels, omega, epsilon, mat, B,
                                  print_info);

  MultiArrayNew<T *[1]> x("x", global_U.get_num_dof());
  MultiArrayNew<T *[1]> rhs("rhs", global_U.get_num_dof());

  for (index_t i = 0; i < nbcs1; i++) {
    rhs(bc_dofs1[i], 0) = 1.0;
  }

  amg.cg(rhs, x, 5);

  // Copy the solution to the solution vector
  for (index_t i = 0; i < global_U.get_num_dof(); i++) {
    global_U[i] = x(i, 0);
  }

  write_hex_to_vtk<1, degree, T, DataBasis, GeoBasis, Basis>(
      pde, elem_data, elem_geo, elem_sol,
      [](index_t k, typename PDE::DataSpace &data,
         typename PDE::FiniteElementGeometry &geo,
         typename PDE::FiniteElementSpace &sol) {
        return sol.template get<0>().get_value();
      });

  // Kokkos::finalize();

  return (0);
}