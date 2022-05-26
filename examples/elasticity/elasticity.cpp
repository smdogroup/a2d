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
  typedef HelmholtzPDE<Basis> Model;
  // typedef NonlinearElasticity3D<Basis> Model;

  const int data_per_point = Model::NUM_DATA;
  const int spatial_dim = Model::SPATIAL_DIM;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  const int nx = 24;
  const int ny = 64;
  const int nz = 64;
  const int nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const int nelems = nx * ny * nz;

  IndexType* conn_data = new IndexType[nodes_per_elem * nelems];
  CLayout<nodes_per_elem> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  ScalarType* X_data = new ScalarType[spatial_dim * nnodes];
  CLayout<spatial_dim> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<spatial_dim> > X(node_layout, X_data);

  ScalarType* U_data = new ScalarType[vars_per_node * nnodes];
  CLayout<vars_per_node> solution_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node> > U(solution_layout, U_data);

  ScalarType* P_data = new ScalarType[vars_per_node * nnodes];
  MultiArray<ScalarType, CLayout<vars_per_node> > P(solution_layout, P_data);

  ScalarType* Xe_data = new ScalarType[spatial_dim * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, spatial_dim> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, spatial_dim> > Xe(
      node_element_layout, Xe_data);

  ScalarType* Ue_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> solution_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Ue(
      solution_element_layout, Ue_data);

  ScalarType* Pe_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Pe(
      solution_element_layout, Pe_data);

  // Allocate space for the data
  ScalarType* mat_data = new ScalarType[data_per_point * num_quad_pts * nelems];
  CLayout<num_quad_pts, data_per_point> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, data_per_point> > data(
      element_mat_layout, mat_data);

#ifdef USE_COMPLEX
  double dh = 1e-30;
  ScalarType perturb = ScalarType(0.0, dh);
#else
  double dh = 1e-6;
  ScalarType perturb = dh;
#endif  // USE_COMPLEX

  // Set the values of the solution
  for (int i = 0; i < nnodes; i++) {
    for (int j = 0; j < vars_per_node; j++) {
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

  element_scatter(conn, U, Ue);
  element_scatter(conn, X, Xe);
  element_scatter(conn, P, Pe);

  // Store data for Uq
  ScalarType* Uq_data = new ScalarType[vars_per_node * num_quad_pts * nelems];
  CLayout<num_quad_pts, vars_per_node> solution_quadpt_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, vars_per_node> > Uq(
      solution_quadpt_layout, Uq_data);

  Basis::interp<vars_per_node>(Ue, Uq);

  // Store data for detJ
  ScalarType* detJ_data = new ScalarType[num_quad_pts * nelems];
  CLayout<num_quad_pts> element_detJ_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts> > detJ(element_detJ_layout,
                                                      detJ_data);

  // Store the data for Jinv
  ScalarType* Jinv_data =
      new ScalarType[spatial_dim * spatial_dim * num_quad_pts * nelems];
  CLayout<num_quad_pts, spatial_dim, spatial_dim> element_grad_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, spatial_dim, spatial_dim> > Jinv(
      element_grad_layout, Jinv_data);

  // Store data for Uxi
  ScalarType* Uxi_data =
      new ScalarType[vars_per_node * spatial_dim * num_quad_pts * nelems];
  CLayout<num_quad_pts, vars_per_node, spatial_dim> element_Uxi_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, vars_per_node, spatial_dim> >
      Uxi(element_Uxi_layout, Uxi_data);

  // Allocate space for the residuals
  ScalarType* res_data =
      new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> res_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > res(
      res_layout, res_data);

  ScalarType* jac_data =
      new ScalarType[nelems * nodes_per_elem * nodes_per_elem * vars_per_node *
                     vars_per_node];
  CLayout<nodes_per_elem, nodes_per_elem, vars_per_node, vars_per_node>
      jac_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, nodes_per_elem, vars_per_node,
                                 vars_per_node> >
      jac(jac_layout, jac_data);

  double pt_data[] = {0.3, 1.2};

  for (int i = 0; i < nelems; i++) {
    for (int j = 0; j < num_quad_pts; j++) {
      for (int k = 0; k < data_per_point; k++) {
        data(i, j, k) = pt_data[k];
      }
    }
  }

  // Zero the residual
  res.zero();

  // Zero the Jacobian
  jac.zero();

  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);

  // ScalarType energy;
  // Model::energy<ScalarType>(data, detJ, Jinv, Uxi, energy);
  // Model::residuals<ScalarType>(data, detJ, Jinv, Uxi, res);
  // Model::jacobians<ScalarType>(data, detJ, Jinv, Uxi, jac);
  Model::residuals<ScalarType>(data, detJ, Jinv, Uq, Uxi, res);
  Model::jacobians<ScalarType>(data, detJ, Jinv, Uq, Uxi, jac);

#ifdef USE_COMPLEX

  for (int i = 0; i < 10; i++) {
    for (int jy = 0; jy < nodes_per_elem; jy++) {
      for (int iy = 0; iy < vars_per_node; iy++) {
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
