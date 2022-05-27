#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"
#include "elasticity3d.h"
#include "helmholtz3d.h"
#include "multiarray.h"

#define USE_COMPLEX 1

using namespace A2D;

/*
  Finite-element computations using templating/auto-diff and multi-dimensional
  arrays
*/
int main(int argc, char* argv[]) {
  typedef int32_t IndexType;
#ifdef USE_COMPLEX
  typedef std::complex<double> ScalarType;
#else
  typedef double ScalarType;
#endif  // USE_COMPLEX

  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  // typedef HelmholtzPDE<IndexType, ScalarType, Basis> Model;
  typedef LinearElasticity3D<IndexType, ScalarType, Basis> Model;

  const int nx = 24;
  const int ny = 64;
  const int nz = 64;
  const int nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const int nelems = nx * ny * nz;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;

  Model model(nelems, nnodes);

  typename Model::base::ConnArray& conn = model.get_conn();
  typename Model::base::NodeArray& X = model.get_nodes();
  typename Model::base::SolutionArray& U = model.get_solution();

  // Create a new solution array
  typename Model::base::SolutionArray& P = *model.new_solution();
  typename Model::base::ElemSolnArray& Pe = *model.new_elem_solution();

  // Residual vector
  typename Model::base::SolutionArray& residual = *model.new_solution();

#ifdef USE_COMPLEX
  double dh = 1e-30;
  ScalarType perturb = ScalarType(0.0, dh);
#else
  double dh = 1e-6;
  ScalarType perturb = dh;
#endif  // USE_COMPLEX

  // Set the node locations
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        X(node, 0) = 1.0 * i / nx;
        X(node, 1) = 1.0 * j / ny;
        X(node, 2) = 1.0 * k / nz;
      }
    }
  }

  // Set the element connectivity
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int elem = i + nx * (j + ny * k);

        for (int kk = 0, index = 0; kk < 2; kk++) {
          for (int jj = 0; jj < 2; jj++) {
            for (int ii = 0; ii < 2; ii++, index++) {
              conn(elem, index) =
                  (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
            }
          }
        }
      }
    }
  }

  typename Model::base::QuadDataArray& data = model.get_quad_data();

  double pt_data[] = {1.23, 2.45};
  for (std::size_t i = 0; i < data.extent(0); i++) {
    for (std::size_t j = 0; j < data.extent(1); j++) {
      for (std::size_t k = 0; k < data.extent(2); k++) {
        data(i, j, k) = pt_data[k];
      }
    }
  }

  // Now reset the nodes
  model.reset_nodes();

  // Set the values of the solution
  for (int i = 0; i < nnodes; i++) {
    for (int j = 0; j < vars_per_node; j++) {
      U(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      P(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      U(i, j) = U(i, j) + perturb * P(i, j);
    }
  }

  // Reset the soultion
  model.reset_solution();

  element_scatter(conn, P, Pe);

  typename Model::base::ElemResArray& res = model.get_elem_res();
  typename Model::base::ElemJacArray& jac = model.get_elem_jac();

  residual.zero();
  res.zero();
  jac.zero();

  model.add_residuals(residual);
  model.add_jacobians();

#ifdef USE_COMPLEX

  for (int i = 0; i < 10; i++) {
    for (int jy = 0; jy < nodes_per_elem; jy++) {
      for (int iy = 0; iy < vars_per_node; iy++) {
        double fd = res(i, jy, iy).imag() / dh;

        ScalarType result = 0.0;
        for (int jx = 0; jx < nodes_per_elem; jx++) {
          for (int ix = 0; ix < vars_per_node; ix++) {
            result += jac(i, jy, jx, iy, ix) * Pe(i, jx, ix);
          }
        }

        std::cout << "fd: " << std::setw(20) << fd
                  << " result: " << std::setw(20) << result.real()
                  << " error: " << std::setw(20) << (fd - result.real()) / fd
                  << std::endl;
      }
    }
  }

#endif  // USE_COMPLEX

  return (0);
}
