#include <iostream>
#include <memory>
#include <string>

#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/integrand_heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/integrand_poisson.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dparser.h"
#include "utils/a2dvtk.h"

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

enum class PDE_TYPE { POISSON, ELASTICITY };

template <int degree, PDE_TYPE integrand_type>
void main_body(std::string type, int ar = 1, double h = 1.0,
               bool write_vtk = false) {
  using I = index_t;
  using T = double;

  constexpr int spatial_dim = 3;
  constexpr int low_degree = 1;

  // Switch between poisson and elasticity
  using Integrand = typename std::conditional<
      integrand_type == PDE_TYPE::POISSON, Poisson<T, spatial_dim>,
      IntegrandTopoLinearElasticity<T, spatial_dim>>::type;
  using Basis = typename std::conditional<
      integrand_type == PDE_TYPE::POISSON,
      FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>,
      FEBasis<T, LagrangeH1HexBasis<T, 3, degree>>>::type;
  using LOrderBasis = typename std::conditional<
      integrand_type == PDE_TYPE::POISSON,
      FEBasis<T, LagrangeH1HexBasis<T, 1, low_degree>>,
      FEBasis<T, LagrangeH1HexBasis<T, 3, low_degree>>>::type;

  using BasisVecType = A2D::SolutionVector<T>;

  // using Integrand = MixedPoisson<T, dim>;
  // using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
  //                       LagrangeL2HexBasis<T, 1, degree - 1>>;
  // using LOrderBasis = FEBasis<T, QHdivHexBasis<T, low_degree>,
  //                             LagrangeL2HexBasis<T, 1, low_degree - 1>>;

  // using Integrand = Poisson<T, dim>;
  // using Basis = FEBasis<T, LagrangeH1HexBasis<T, 1, degree>>;
  // using LOrderBasis = FEBasis<T, LagrangeH1HexBasis<T, 1, 1>>;

  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using DataElemVec = A2D::ElementVector_Serial<T, DataBasis, BasisVecType>;
  using GeoElemVec = A2D::ElementVector_Serial<T, GeoBasis, BasisVecType>;
  using ElemVec = A2D::ElementVector_Serial<T, Basis, BasisVecType>;
  using FE =
      FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis>;

  using LOrderQuadrature = HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = FEBasis<T>;
  using LOrderGeoBasis =
      FEBasis<T, LagrangeH1HexBasis<T, spatial_dim, low_degree>>;
  using LOrderDataElemVec =
      A2D::ElementVector_Serial<T, LOrderDataBasis, BasisVecType>;
  using LOrderGeoElemVec =
      A2D::ElementVector_Serial<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = A2D::ElementVector_Serial<T, LOrderBasis, BasisVecType>;

  using LOrderFE = FiniteElement<T, Integrand, LOrderQuadrature,
                                 LOrderDataBasis, LOrderGeoBasis, LOrderBasis>;

  /** Create mesh for a single element
   *
   *        7 --------------- 6
   *       / |              / |
   *      /  |             /  |
   *     /   |            /   |
   *    4 -------------- 5    |
   *    |    |           |    |
   *    |    3 ----------|--- 2
   *    |   /            |   /
   *    |  /             |  /
   *    | /              | /
   *    0 -------------- 1
   */

  // - Number of elements in each dimension
  auto node_num = [](I i, I j, I k) { return i + 2 * (j + 2 * k); };

  // - Number of edges
  const I nverts = 8;
  I ntets = 0, nwedge = 0, npyrmd = 0;
  const I nhex = 1;

  I *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  I hex[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  // - Nodal location
  double Xloc[3 * 8];

  if (type == "box") {
    Xloc[0] = 0.0, Xloc[1] = 0.0, Xloc[2] = 0.0;          // pt0
    Xloc[3] = ar * 1.0, Xloc[4] = 0.0, Xloc[5] = 0.0;     // pt1
    Xloc[6] = ar * 1.0, Xloc[7] = 1.0, Xloc[8] = 0.0;     // pt2
    Xloc[9] = 0.0, Xloc[10] = 1.0, Xloc[11] = 0.0;        // pt3
    Xloc[12] = 0.0, Xloc[13] = 0.0, Xloc[14] = 1.0;       // pt4
    Xloc[15] = ar * 1.0, Xloc[16] = 0.0, Xloc[17] = 1.0;  // pt5
    Xloc[18] = ar * 1.0, Xloc[19] = 1.0, Xloc[20] = 1.0;  // pt6
    Xloc[21] = 0.0, Xloc[22] = 1.0, Xloc[23] = 1.0;       // pt7

  } else if (type == "distortion") {
    Xloc[0] = 0.5 - 0.5 / ar, Xloc[1] = 0.0, Xloc[2] = 0.0;    // pt0
    Xloc[3] = 0.5 + 0.5 / ar, Xloc[4] = 0.0, Xloc[5] = 0.0;    // pt1
    Xloc[6] = 0.5 + 0.5 / ar, Xloc[7] = 1.0, Xloc[8] = 0.0;    // pt2
    Xloc[9] = 0.5 - 0.5 / ar, Xloc[10] = 1.0, Xloc[11] = 0.0;  // pt3
    Xloc[12] = 0.0, Xloc[13] = 0.5 - 0.5 / ar, Xloc[14] = h;   // pt4
    Xloc[15] = 1.0, Xloc[16] = 0.5 - 0.5 / ar, Xloc[17] = h;   // pt5
    Xloc[18] = 1.0, Xloc[19] = 0.5 + 0.5 / ar, Xloc[20] = h;   // pt6
    Xloc[21] = 0.0, Xloc[22] = 0.5 + 0.5 / ar, Xloc[23] = h;   // pt7
  } else {
    std::printf("Invalid type %s\n", type.c_str());
    exit(-1);
  }

  // Save element to vtk
  if (write_vtk) {
    A2D::ToVTK3D(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd, pyrmd,
                 Xloc, "spectral_element.vtk");
  }

  // Build element connectivity that A2D needs
  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);

  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);
  ElementMesh<DataBasis> datamesh(conn);

  ElementMesh<LOrderBasis> lorder_mesh(mesh);
  ElementMesh<LOrderGeoBasis> lorder_geomesh(geomesh);
  ElementMesh<LOrderDataBasis> lorder_datamesh(datamesh);

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
  std::shared_ptr<Integrand> integrand;
  if constexpr (integrand_type == PDE_TYPE::POISSON) {
    integrand = std::make_shared<Integrand>();
  } else {
    integrand = std::make_shared<Integrand>(1.0, 0.3, 5.0);
  }

  const index_t block_size = 1;
  using BSRMatType = BSRMat<T, block_size, block_size>;

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
  fe.add_jacobian(*integrand, elem_data, elem_geo, elem_sol, elem_mat);

  lorder_mesh.template create_block_csr<block_size>(nrows, rowp, cols);
  BSRMatType lorder_mat(nrows, nrows, rowp[nrows], rowp, cols);
  ElementMat_Serial<T, LOrderBasis, BSRMatType> lorder_elem_mat(lorder_mesh,
                                                                lorder_mat);
  lorder_fe.add_jacobian(*integrand, lorder_elem_data, lorder_elem_geo,
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

  for (I i = 0; i < n; i++) {
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

  int dof = mesh.get_num_dof();
  double cond = std::fabs(eig_max) / std::fabs(eig_min);

  if (degree == 1 && ar == 1) {
    std::printf("%10s%10s%10s%15s%15s%15s\n", "degree", "AR", "dof", "eig_min",
                "eig_max", "cond");
  }
  std::printf("%10d%10d%10d%15.5e%15.5e%15.5e\n", degree, ar, dof, eig_min,
              eig_max, cond);

  delete[] Ap;
  delete[] Ah;
  delete[] eigvals;

  // write_hex_to_vtk<1, degree, T, DataBasis, GeoBasis, Basis>(
  //     integrand, elem_data, elem_geo, elem_sol,
  //     [](index_t k, typename Integrand::DataSpace &data,
  //        typename Integrand::FiniteElementGeometry &geo,
  //        typename Integrand::FiniteElementSpace &sol) {
  //       return sol.template get<0>().get_value();
  //     });
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();
  {
    // Get cmd arguments
    ArgumentParser parser(argc, argv);
    std::string type = parser.parse_option("--type", std::string("distortion"));
    std::string integrand =
        parser.parse_option("--integrand", std::string("poisson"));
    double h = parser.parse_option("--h", 1.0);
    parser.help_info();

    // Check option validity
    std::vector<std::string> valid_types = {"distortion", "box"};
    assert_option_in(type, valid_types);
    std::vector<std::string> valid_integrands = {"poisson", "elasticity"};
    assert_option_in(integrand, valid_integrands);

    int ar[] = {1, 2, 4, 8, 16, 32};

    if (integrand == "poisson") {
      for (auto a : ar) {
        main_body<1, PDE_TYPE::POISSON>(type, a, h);
        main_body<2, PDE_TYPE::POISSON>(type, a, h);
        main_body<3, PDE_TYPE::POISSON>(type, a, h);
        main_body<4, PDE_TYPE::POISSON>(type, a, h);
        main_body<5, PDE_TYPE::POISSON>(type, a, h);
        main_body<6, PDE_TYPE::POISSON>(type, a, h);
        main_body<7, PDE_TYPE::POISSON>(type, a, h);
        main_body<8, PDE_TYPE::POISSON>(type, a, h);
        main_body<9, PDE_TYPE::POISSON>(type, a, h);
        main_body<10, PDE_TYPE::POISSON>(type, a, h, true);
      }
    } else {
      for (auto a : ar) {
        main_body<1, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<2, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<3, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<4, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<5, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<6, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<7, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<8, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<9, PDE_TYPE::ELASTICITY>(type, a, h);
        main_body<10, PDE_TYPE::ELASTICITY>(type, a, h, true);
      }
    }
  }
  Kokkos::finalize();
}