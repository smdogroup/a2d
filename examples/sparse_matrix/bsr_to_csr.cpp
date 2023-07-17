#include <vector>

#include "a2dobjs.h"
#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

template <typename T>
class FEProb {
 public:
  static constexpr int spatial_dim = 3;
  static constexpr int data_dim = 1;
  static constexpr int var_dim = 3;
  static constexpr int block_size = var_dim;
  static constexpr int degree = 1;

  using Vec_t = SolutionVector<T>;
  using BSRMat_t = BSRMat<T, block_size, block_size>;

  using PDE = TopoLinearElasticity<T, spatial_dim>;
  using Quadrature = HexGaussQuadrature<degree + 1>;
  using DataBasis = FEBasis<T, LagrangeL2HexBasis<T, data_dim, degree - 1>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = FEBasis<T, LagrangeH1HexBasis<T, var_dim, degree>>;

  using FE = FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;
  using DataElemVec = ElementVector_Serial<T, DataBasis, Vec_t>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis, Vec_t>;
  using ElemVec = ElementVector_Serial<T, Basis, Vec_t>;

  FEProb(MeshConnectivity3D &conn, int nhex, const int hex[],
         const double Xloc[], T E = 70e3, T nu = 0.3, T q = 5.0)
      : E(E),
        nu(nu),
        q(q),
        datamesh(conn),
        geomesh(conn),
        mesh(conn),
        data(datamesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        sol(mesh.get_num_dof()),
        elem_data(datamesh, data),
        elem_geo(geomesh, geo),
        elem_sol(mesh, sol) {
    // Set geometry
    set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);
  }

  BSRMat_t create_fea_bsr_matrix(MeshConnectivity3D &conn) {
    // Symbolically create block CSR matrix

    index_t nrows;
    std::vector<index_t> rowp, cols;
    mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    BSRMat_t bsr_mat(nrows, nrows, cols.size(), rowp, cols);

    // Populate the CSR matrix
    FE fe;
    PDE pde(E, nu, q);

    ElementMat_Serial<T, Basis, BSRMat_t> elem_mat(mesh, bsr_mat);
    fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

    return bsr_mat;
  }

 private:
  T E, nu, q;
  ElementMesh<DataBasis> datamesh;
  ElementMesh<GeoBasis> geomesh;
  ElementMesh<Basis> mesh;
  Vec_t data;
  Vec_t geo;
  Vec_t sol;
  DataElemVec elem_data;
  GeoElemVec elem_geo;
  ElemVec elem_sol;
};

int main() {
  Kokkos::initialize();
  {
    using T = double;

    // Number of elements in each dimension
    const int nx = 2, ny = 1, nz = 1;
    const double lx = 1.5, ly = 2.0, lz = 1.0;

    // Set up mesh
    const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
    int ntets = 0, nwedge = 0, npyrmd = 0;
    const int nhex = nx * ny * nz;
    int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
    int hex[8 * nhex];
    double Xloc[3 * nverts];
    A2D::MesherBrick3D mesher(nx, ny, nz, lx, ly, lz);
    mesher.set_X_conn<int, double>(Xloc, hex);
    A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                                 npyrmd, pyrmd);

    FEProb<T> prob(conn, nhex, hex, Xloc);

    // Create bsr mat
    FEProb<T>::BSRMat_t bsr_mat = prob.create_fea_bsr_matrix(conn);

    // Create csr mat
    CSRMat<T> csr_mat = bsr_to_csr(bsr_mat);

    // Create csc mat
    CSCMat<T> csc_mat = bsr_to_csc(bsr_mat);

    // Write to mtx
    bsr_mat.write_mtx("bsr_mat.mtx");
    csr_mat.write_mtx("csr_mat.mtx");
    csc_mat.write_mtx("csc_mat.mtx");
  }

  Kokkos::finalize();
}