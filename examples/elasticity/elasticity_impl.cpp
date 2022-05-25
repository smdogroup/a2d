#include "elasticity_impl.h"

using namespace A2D;

void compute_residual(int nelems, int nnodes, int* conn_data, double* X_data,
                      double* mat_data, double* U_data, double* res_data) {
  typedef int IndexType;

  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef NonlinearElasticity3D Model;
  typedef double ScalarType;

  const int vars_per_node = 3;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<8> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<3> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<3> > X(node_layout, X_data);
  MultiArray<ScalarType, CLayout<3> > U(node_layout, U_data);

  CLayout<num_quad_pts, 2> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 2> > data(element_mat_layout,
                                                         mat_data);

  CLayout<nodes_per_elem, 3> res_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > res(res_layout, res_data);

  ScalarType* Xe_data = new ScalarType[3 * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, 3> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xe(node_element_layout,
                                                         Xe_data);

  ScalarType* Ue_data = new ScalarType[3 * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Ue(node_element_layout,
                                                         Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

  // Store data for detJ
  ScalarType* detJ_data = new ScalarType[num_quad_pts * nelems];
  CLayout<num_quad_pts> element_detJ_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts> > detJ(element_detJ_layout,
                                                      detJ_data);

  // Store the data for Jinv
  ScalarType* Jinv_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_grad_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Jinv(element_grad_layout,
                                                            Jinv_data);

  // Store data for Uxi
  ScalarType* Uxi_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_Uxi_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Uxi(element_Uxi_layout,
                                                           Uxi_data);

  // Zero the residual array
  for (int i = 0; i < res.extent(0); i++) {
    for (int j = 0; j < res.extent(1); j++) {
      for (int k = 0; k < res.extent(2); k++) {
        res(i, j, k) = 0.0;
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Basis::residuals<ScalarType, Model>(data, detJ, Jinv, Uxi, res);

  // Free data
  delete[] Xe_data;
  delete[] Ue_data;
  delete[] detJ_data;
  delete[] Jinv_data;
  delete[] Uxi_data;
}

void compute_jacobian(int nelems, int nnodes, int* conn_data, double* X_data,
                      double* mat_data, double* U_data, double* jac_data) {
  typedef int IndexType;

  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef NonlinearElasticity3D Model;
  typedef double ScalarType;

  const int vars_per_node = 3;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<8> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<3> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<3> > X(node_layout, X_data);
  MultiArray<ScalarType, CLayout<3> > U(node_layout, U_data);

  CLayout<num_quad_pts, 2> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 2> > data(element_mat_layout,
                                                         mat_data);

  CLayout<nodes_per_elem, nodes_per_elem, 3, 3> jac_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, nodes_per_elem, 3, 3> > jac(
      jac_layout, jac_data);

  ScalarType* Xe_data = new ScalarType[3 * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, 3> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xe(node_element_layout,
                                                         Xe_data);

  ScalarType* Ue_data = new ScalarType[3 * nodes_per_elem * nelems];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Ue(node_element_layout,
                                                         Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

  // Store data for detJ
  ScalarType* detJ_data = new ScalarType[num_quad_pts * nelems];
  CLayout<num_quad_pts> element_detJ_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts> > detJ(element_detJ_layout,
                                                      detJ_data);

  // Store the data for Jinv
  ScalarType* Jinv_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_grad_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Jinv(element_grad_layout,
                                                            Jinv_data);

  // Store data for Uxi
  ScalarType* Uxi_data = new ScalarType[3 * 3 * num_quad_pts * nelems];
  CLayout<num_quad_pts, 3, 3> element_Uxi_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Uxi(element_Uxi_layout,
                                                           Uxi_data);

  // Zero the Jacobian
  for (int i = 0; i < jac.extent(0); i++) {
    for (int j = 0; j < jac.extent(1); j++) {
      for (int k = 0; k < jac.extent(2); k++) {
        for (int l = 0; l < jac.extent(3); l++) {
          for (int m = 0; m < jac.extent(4); m++) {
            jac(i, j, k, l, m) = 0.0;
          }
        }
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Basis::jacobians<ScalarType, Model>(data, detJ, Jinv, Uxi, jac);

  // Free data
  delete[] Xe_data;
  delete[] Ue_data;
  delete[] detJ_data;
  delete[] Jinv_data;
  delete[] Uxi_data;
}