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

using namespace A2D;

extern "C" {
void dggev_(const char *jobvl, const char *jobvr, int *N, double *A, int *lda,
            double *B, int *ldb, double *alphar, double *alphai, double *beta,
            double *vl, int *ldvl, double *vr, int *ldvr, double *work,
            int *lwork, int *info);
}

/**
 * @brief Compute the eigenvalues of the problem
 *
 * A * v(j) = lambda(j) * B * v(j)
 *
 * @param n The matrix size
 * @param A The A matrix
 * @param B The B matrix
 * @param eigvals The eigenvalues of the problem
 */
void compute_eigenvalues(index_t n, double *A, double *B, double *eigvals,
                         bool write = false) {
  int N = n;

  double *alphar = new double[n];
  double *alphai = new double[n];
  double *beta = new double[n];
  int lwork = 20 * n;
  double *work = new double[lwork];

  int info;
  dggev_("N", "N", &N, A, &N, B, &N, alphar, alphai, beta, NULL, &N, NULL, &N,
         work, &lwork, &info);

  const double eps = 1e-12;
  for (int i = 0; i < n; i++) {
    if (write) {
      std::cout << std::setw(12) << alphar[i] / beta[i] << " " << std::setw(12)
                << alphar[i] << "  " << std::setw(12) << alphai[i] << " "
                << std::setw(12) << beta[i] << std::endl;
    }

    if (beta[i] > eps || beta[i] < -eps) {
      eigvals[i] = alphar[i] / beta[i];
    } else {
      eigvals[i] = 0.0;
    }
  }

  delete[] alphar;
  delete[] alphai;
  delete[] beta;
  delete[] work;
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const index_t dim = 3;
  using T = double;
  using ET = ElementTypes;

  const index_t degree = 10;
  const index_t low_degree = 1;

  // using PDE = MixedPoisson<T, dim>;
  // using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
  //                       LagrangeL2HexBasis<T, 1, degree - 1>>;
  // using LOrderBasis = FEBasis<T, QHdivHexBasis<T, low_degree>,
  //                             LagrangeL2HexBasis<T, 1, low_degree - 1>>;

  using BasisVecType = A2D::SolutionVector<T>;

  using PDE = Poisson<T, dim>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  using LOrderBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, low_degree>>;

  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, degree>>;
  using DataElemVec = A2D::ElementVector_Serial<T, DataBasis, BasisVecType>;
  using GeoElemVec = A2D::ElementVector_Serial<T, GeoBasis, BasisVecType>;
  using ElemVec = A2D::ElementVector_Serial<T, Basis, BasisVecType>;
  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  using LOrderQuadrature = HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = FEBasis<T>;
  using LOrderGeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, low_degree>>;
  using LOrderDataElemVec =
      A2D::ElementVector_Serial<T, LOrderDataBasis, BasisVecType>;
  using LOrderGeoElemVec =
      A2D::ElementVector_Serial<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = A2D::ElementVector_Serial<T, LOrderBasis, BasisVecType>;

  using LOrderFE = FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                 LOrderGeoBasis, LOrderBasis>;

  // Number of elements in each dimension
  auto node_num = [](int i, int j, int k) { return i + 2 * (j + 2 * k); };

  // Number of edges
  const int nverts = 8;
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = 1;

  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8];

  for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
    hex[ii] = node_num(ET::HEX_VERTS_CART[ii][0], ET::HEX_VERTS_CART[ii][1],
                       ET::HEX_VERTS_CART[ii][2]);
  }

  double Xloc[3 * 8];
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        Xloc[3 * node_num(i, j, k)] = (1.0 * i);
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j);
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k);
      }
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

  SolutionVector<T> sol(mesh.get_num_dof());
  SolutionVector<T> geo(geomesh.get_num_dof());
  SolutionVector<T> data(datamesh.get_num_dof());

  DataElemVec elem_data(datamesh, data);
  GeoElemVec elem_geo(geomesh, geo);
  ElemVec elem_sol(mesh, sol);

  // Set the geometry from the node locations
  set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  LOrderDataElemVec lorder_elem_data(lorder_datamesh, data);
  LOrderGeoElemVec lorder_elem_geo(lorder_geomesh, geo);
  LOrderElemVec lorder_elem_sol(lorder_mesh, sol);

  // Create the finite-element model
  FE fe;
  LOrderFE lorder_fe;
  PDE pde;

  const index_t block_size = 1;
  using BSRMatType = BSRMat<index_t, T, block_size, block_size>;

  index_t nrows = mesh.get_num_dof();
  std::vector<index_t> rowp(nrows + 1), cols(nrows * nrows);

  // mesh.create_block_csr<block_size>(nrows, rowp, cols);
  index_t nnz = nrows * nrows;
  for (index_t i = 0; i < nrows; i++) {
    rowp[i] = i * nrows;
    for (index_t j = 0; j < nrows; j++) {
      cols[rowp[i] + j] = j;
    }
  }
  rowp[nrows] = nrows * nrows;

  BSRMatType horder_mat(nrows, nrows, rowp[nrows], rowp, cols);
  ElementMat_Serial<T, Basis, BSRMatType> elem_mat(mesh, horder_mat);
  fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

  lorder_mesh.create_block_csr<block_size>(nrows, rowp, cols);
  BSRMatType lorder_mat(nrows, nrows, rowp[nrows], rowp, cols);
  ElementMat_Serial<T, LOrderBasis, BSRMatType> lorder_elem_mat(lorder_mesh,
                                                                lorder_mat);
  lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                         lorder_elem_sol, lorder_elem_mat);

  index_t n, m;
  T *Ap, *Ah;
  horder_mat.to_dense(&n, &m, &Ap);
  lorder_mat.to_dense(&n, &m, &Ah);

  double *eigvals = new double[n];
  compute_eigenvalues(n, Ap, Ah, eigvals);

  double eig_min = 1e20;
  double eig_max = -1e20;
  double eig_min_abs = 1e20;

  for (int i = 0; i < n; i++) {
    if (eigvals[i] != 0.0) {
      if (std::fabs(eigvals[i]) < eig_min_abs) {
        eig_min_abs = std::fabs(eigvals[i]);
      }
      if (eigvals[i] < eig_min) {
        eig_min = eigvals[i];
      }
      if (eigvals[i] > eig_max) {
        eig_max = eigvals[i];
      }
    }
  }

  std::cout << "dof:         " << mesh.get_num_dof() << std::endl;
  std::cout << "dof:         " << (degree + 1) * (degree + 1) * (degree + 1)
            << std::endl;
  std::cout << "degree:      " << degree << std::endl;
  std::cout << "eig_min:     " << eig_min << std::endl;
  std::cout << "eig_max:     " << eig_max << std::endl;
  std::cout << "eig_min_abs: " << eig_min_abs << std::endl;

  delete[] Ap;
  delete[] Ah;
  delete[] eigvals;

  // write_hex_to_vtk<1, degree, T, DataBasis, GeoBasis, Basis>(
  //     pde, elem_data, elem_geo, elem_sol,
  //     [](index_t k, typename PDE::DataSpace &data,
  //        typename PDE::FiniteElementGeometry &geo,
  //        typename PDE::FiniteElementSpace &sol) {
  //       return sol.template get<0>().get_value();
  //     });

  return (0);
}