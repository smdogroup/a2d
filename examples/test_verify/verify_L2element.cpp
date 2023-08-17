#include "a2dobjs.h"
#include "a2dtypes.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/feelementtypes.h"
#include "multiphysics/feelementvector.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/fespace.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "utils/a2dvtk.h"

using namespace A2D;

template <typename T, index_t D>
class ElementTesterPDE {
 public:
  static constexpr index_t dim = D;
  using FiniteElementSpace = FESpace<T, dim, L2Space<T, 1, dim>>;
  using FiniteElementGeometry = FESpace<T, dim, H1Space<T, dim, dim>>;
};

template <index_t order>
class SamplerQuadrature {
 public:
  static const bool is_tensor_product = true;
  static constexpr index_t num_quad_points = order * order * order;
  static const index_t tensor_dim0 = order;
  static const index_t tensor_dim1 = order;
  static const index_t tensor_dim2 = order;

  static double get_tensor_point(const index_t dim, const index_t pt) {
    double tp = 2 * double(pt) / double(order - 1) - 1.0;
    return tp;
  }

  static double get_tensor_weight(const index_t dim, const index_t pt) {
    return 1.0;
  }

  static index_t get_tensor_index(const index_t q0, const index_t q1,
                                  const index_t q2) {
    return q0 + order * (q1 + order * q2);
  }

  static constexpr index_t get_num_points() { return num_quad_points; }

  static void get_point(const index_t n, double pt[]) {
    index_t i = n % order;
    index_t j = (n % (order * order)) / order;
    index_t k = n / (order * order);

    pt[0] = 2 * double(i) / double(order - 1) - 1.0;
    pt[1] = 2 * double(j) / double(order - 1) - 1.0;
    pt[2] = 2 * double(k) / double(order - 1) - 1.0;
    printf("[%2d]%10.5f %10.5f %10.5f\n", n, pt[0], pt[1], pt[2]);
  }

  static double get_weight(const index_t n) { return 1.0; }
};

void main_body() {
  // Set up magic numbers
  constexpr index_t dim = 3;
  constexpr index_t geo_degree = 1;  // use quadratic element for the mesh
  constexpr index_t degree = 3;      // use quadratic element for "physics"
  constexpr index_t nsample_per_dim = 20;

  // Set up type aliases
  using T = double;
  using ET = ElementTypes;
  using SamplingQuadrature = SamplerQuadrature<nsample_per_dim>;
  using Basis = FEBasis<T, LagrangeL2HexBasis<T, 1, degree>>;
  using GeoBasis = FEBasis<T, LagrangeH1HexBasis<T, dim, geo_degree>>;
  using PDE = ElementTesterPDE<T, dim>;

  const index_t nx = 1, ny = 1, nz = 1;  // number of elements in each dimension
  const index_t nelems = nx * ny * nz;
  const index_t nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  auto node_num = [](index_t i, index_t j, index_t k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };  // compute element-local node numbering given (i, j, k)

  const index_t nverts = (nx + 1) * (ny + 1) * (nz + 1);
  index_t ntets = 0, nwedge = 0, npyrmd = 0;
  const index_t nhex = nx * ny * nz;

  index_t *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  index_t hex[8 * nhex];

  for (index_t k = 0, e = 0; k < nz; k++) {
    for (index_t j = 0; j < ny; j++) {
      for (index_t i = 0; i < nx; i++, e++) {
        for (index_t ii = 0; ii < ET::HEX_NVERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  // Create mesh connectivity
  MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                          pyrmd);
  ElementMesh<Basis> mesh(conn);
  ElementMesh<GeoBasis> geomesh(conn);

  // Create mesh nodal locations
  double Xloc[3 * nverts];
  for (index_t k = 0; k < nz + 1; k++) {
    for (index_t j = 0; j < ny + 1; j++) {
      for (index_t i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (1.0 * i);
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j);
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k);
      }
    }
  }

  // Write mesh to vtk
  ToVTK3D vtk(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd, pyrmd,
              Xloc, "mesh.vtk");

  // Create and populate fake solution variable
  SolutionVector<T> sol(mesh.get_num_dof());
  std::printf("num dof: %d\n", mesh.get_num_dof());
  for (index_t i = 0; i < mesh.get_num_dof(); i++) {
    sol[i] = std::rand() / (double)RAND_MAX + 0.1;
  }

  printf("geo dof: %d\n", geomesh.get_num_dof());
  printf("dim: %d\n", dim);
  printf("nverts: %d\n", nverts);
  assert(geomesh.get_num_dof() == dim * nverts);
  SolutionVector<T> sol_X(geomesh.get_num_dof());
  for (index_t i = 0; i < geomesh.get_num_dof(); i++) {
    sol_X[i] = Xloc[i];
  }

  // Allocate an inplace element view
  using ElemVec = ElementVector_Serial<T, Basis, SolutionVector<T>>;
  using GeoElemVec = ElementVector_Serial<T, GeoBasis, SolutionVector<T>>;
  ElemVec elemvec(mesh, sol);
  GeoElemVec geoelemvec(geomesh, sol_X);

  // Allocate interpolation scalar and vector arrays
  index_t nsamples =
      nsample_per_dim * nsample_per_dim * nsample_per_dim * nelems;
  std::vector<PDE::FiniteElementSpace> spaces;
  std::vector<PDE::FiniteElementGeometry> geo_spaces;
  spaces.reserve(nsamples);
  geo_spaces.reserve(nsamples);

  std::vector<double> vector_field;
  vector_field.reserve(dim * nsamples);

  std::vector<double> vector_field_Xloc;
  vector_field_Xloc.reserve(dim * nsamples);

  // Within each element, evaluate the basis functions at sampling points
  // Loop over elements
  for (index_t k = 0, elem = 0, index = 0; k < nz; k++) {
    for (index_t j = 0; j < ny; j++) {
      for (index_t i = 0; i < nx; i++, elem++) {
        ElemVec::FEDof dof(elem, elemvec);
        GeoElemVec::FEDof dof_geo(elem, geoelemvec);

        elemvec.get_element_values(elem, dof);

        // Get the geometry values
        for (index_t ii = 0; ii < ET::HEX_NVERTS; ii++) {
          index_t node = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                  j + ET::HEX_VERTS_CART[ii][1],
                                  k + ET::HEX_VERTS_CART[ii][2]);

          // Set the entity DOF
          index_t basis = 0;
          index_t orient = 0;
          GeoBasis::set_entity_dof(basis, ET::Vertex, ii, orient,
                                   &Xloc[3 * node], dof_geo);
        }

        geoelemvec.set_element_values(elem, dof_geo);

        using QptSpaceSol =
            QptSpace<SamplingQuadrature, PDE::FiniteElementSpace>;
        QptSpaceSol qptspace;
        using QptSpaceGeo =
            QptSpace<SamplingQuadrature, PDE::FiniteElementGeometry>;
        QptSpaceGeo qptspace_geo;

        Basis::template interp<SamplingQuadrature, ElemVec::FEDof,
                               PDE::FiniteElementSpace>(dof, qptspace);
        GeoBasis::template interp<SamplingQuadrature, GeoElemVec::FEDof,
                                  PDE::FiniteElementGeometry>(dof_geo,
                                                              qptspace_geo);

        // Loop over sample points in each element
        for (index_t kk = 0, pt = 0; kk < nsample_per_dim; kk++) {
          for (index_t jj = 0; jj < nsample_per_dim; jj++) {
            for (index_t ii = 0; ii < nsample_per_dim; ii++, pt++, index++) {
              T u = qptspace.get(pt).get<0>().get_value();
              // T& div = qptspace.get(pt).get<0>().get_div();

              Vec<T, dim>& x = qptspace_geo.get(pt).get<0>().get_value();

              // printf(
              //     "elem: %2d, sample: %3d, x: %10.5f y: %10.5f z: %10.5f u: "
              //     "%10.5f v: %10.5f w:%10.5f "
              //     "div:%10.5f\n",
              //     elem, pt, x(0), x(1), x(2), u(0), u(1), u(2), div);

              vector_field[dim * index] = u;
              for (index_t idim = 0; idim < dim; idim++) {
                vector_field_Xloc[dim * index + idim] = x(idim);
              }
            }
          }
        }
      }
    }
  }

  // Write vector field to vtk
  VectorFieldToVTK vtk_vec(nsamples, vector_field_Xloc.data(),
                           vector_field.data());

  return;
}

int main() {
  Kokkos::initialize();
  { main_body(); }
  Kokkos::finalize();
  return 0;
}
