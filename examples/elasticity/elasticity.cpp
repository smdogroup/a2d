#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "a2dtmp.h"
#include "elasticity3d.h"
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
  typedef NonlinearElasticity3D<Basis> Model;

  const int nx = 24;
  const int ny = 64;
  const int nz = 64;
  const int nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const int nelems = nx * ny * nz;

  const int vars_per_node = 3;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  IndexType* conn_data = new IndexType[nodes_per_elem * nelems];
  CLayout<8> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  ScalarType* X_data = new ScalarType[3 * nnodes];
  CLayout<3> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<3> > X(node_layout, X_data);

  ScalarType* U_data = new ScalarType[3 * nnodes];
  MultiArray<ScalarType, CLayout<3> > U(node_layout, U_data);

  ScalarType* P_data = new ScalarType[3 * nnodes];
  MultiArray<ScalarType, CLayout<3> > P(node_layout, P_data);

  ScalarType* Xe_data = new ScalarType[3 * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, 3> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xe(node_element_layout,
                                                         Xe_data);

  ScalarType* Ue_data = new ScalarType[3 * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Ue(node_element_layout,
                                                         Ue_data);

  ScalarType* Pe_data = new ScalarType[3 * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Pe(node_element_layout,
                                                         Pe_data);

  ScalarType* Xq_data = new ScalarType[3 * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xq(node_element_layout,
                                                         Xq_data);

#ifdef USE_COMPLEX
  double dh = 1e-30;
  ScalarType perturb = ScalarType(0.0, dh);
#else
  double dh = 1e-6;
  ScalarType perturb = dh;
#endif  // USE_COMPLEX

  // Set the values of the solution
  for (int i = 0; i < nnodes; i++) {
    for (int j = 0; j < 3; j++) {
      U(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;
      P(i, j) = -1.0 + 2.0 * rand() / RAND_MAX;

      U(i, j) = U(i, j) + perturb * P(i, j);
    }
  }

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

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);
  element_scatter(conn, P, Pe);

  // Store the data for Jinv
  ScalarType* Jinv_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_grad_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Jinv(element_grad_layout,
                                                            Jinv_data);

  // Store data for detJ
  ScalarType* detJ_data = new ScalarType[num_quad_pts * nelems];
  CLayout<num_quad_pts> element_detJ_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts> > detJ(element_detJ_layout,
                                                      detJ_data);

  // Store data for Uxi
  ScalarType* Uxi_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_Uxi_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Uxi(element_Uxi_layout,
                                                           Uxi_data);

  // Allocate space for the residuals
  ScalarType* res_data = new ScalarType[3 * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, 3> res_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > res(res_layout, res_data);

  // Allocate space for the data
  ScalarType* mat_data = new ScalarType[2 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 2> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 2> > data(element_mat_layout,
                                                         mat_data);

  ScalarType* jac_data =
      new ScalarType[nelems * nodes_per_elem * nodes_per_elem * 3 * 3];
  CLayout<nodes_per_elem, nodes_per_elem, 3, 3> jac_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, nodes_per_elem, 3, 3> > jac(
      jac_layout, jac_data);

  double mu = 0.3;
  double lambda = 1.2;

  for (int i = 0; i < nelems; i++) {
    for (int j = 0; j < num_quad_pts; j++) {
      data(i, j, 0) = mu;
      data(i, j, 1) = lambda;
    }
  }

  Basis::interp<vars_per_node>(Xe, Xq);
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);

  ScalarType energy;
  Model::energy<ScalarType>(data, detJ, Jinv, Uxi, energy);
  Model::residuals<ScalarType>(data, detJ, Jinv, Uxi, res);
  Model::jacobians<ScalarType>(data, detJ, Jinv, Uxi, jac);

#ifdef USE_COMPLEX

  for (int i = 0; i < 10; i++) {
    for (int jy = 0; jy < nodes_per_elem; jy++) {
      for (int iy = 0; iy < 3; iy++) {
        double fd = res(i, jy, iy).imag() / dh;

        ScalarType result = 0.0;
        for (int jx = 0; jx < nodes_per_elem; jx++) {
          for (int ix = 0; ix < 3; ix++) {
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
