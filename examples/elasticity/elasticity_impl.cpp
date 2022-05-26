#include "elasticity_impl.h"

using namespace A2D;

void compute_residual(int nelems, int nnodes, int* conn_data, double* X_data,
                      double* mat_data, double* U_data, double* res_data) {
  typedef int IndexType;
  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef LinearElasticity3D<Basis> Model;
  typedef double ScalarType;

  const int data_per_point = Model::NUM_DATA;
  const int spatial_dim = Model::SPATIAL_DIM;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<nodes_per_elem> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<spatial_dim> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<spatial_dim> > X(node_layout, X_data);

  CLayout<vars_per_node> solution_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node> > U(solution_layout, U_data);

  CLayout<num_quad_pts, data_per_point> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, data_per_point> > data(
      element_mat_layout, mat_data);

  CLayout<nodes_per_elem, vars_per_node> res_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > res(
      res_layout, res_data);

  ScalarType* Xe_data = new ScalarType[spatial_dim * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, spatial_dim> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, spatial_dim> > Xe(
      node_element_layout, Xe_data);

  ScalarType* Ue_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> solution_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Ue(
      solution_element_layout, Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

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

  // Zero the residual array
  for (std::size_t i = 0; i < res.extent(0); i++) {
    for (std::size_t j = 0; j < res.extent(1); j++) {
      for (std::size_t k = 0; k < res.extent(2); k++) {
        res(i, j, k) = 0.0;
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Model::residuals<ScalarType>(data, detJ, Jinv, Uxi, res);

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
  typedef LinearElasticity3D<Basis> Model;
  typedef double ScalarType;

  const int data_per_point = Model::NUM_DATA;
  const int spatial_dim = Model::SPATIAL_DIM;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<nodes_per_elem> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<spatial_dim> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<spatial_dim> > X(node_layout, X_data);

  CLayout<vars_per_node> solution_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node> > U(solution_layout, U_data);

  CLayout<num_quad_pts, data_per_point> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, data_per_point> > data(
      element_mat_layout, mat_data);

  CLayout<nodes_per_elem, nodes_per_elem, vars_per_node, vars_per_node>
      jac_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, nodes_per_elem, vars_per_node,
                                 vars_per_node> >
      jac(jac_layout, jac_data);

  ScalarType* Xe_data = new ScalarType[spatial_dim * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, spatial_dim> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, spatial_dim> > Xe(
      node_element_layout, Xe_data);

  ScalarType* Ue_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> solution_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Ue(
      solution_element_layout, Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

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

  // Zero the Jacobian
  for (std::size_t i = 0; i < jac.extent(0); i++) {
    for (std::size_t j = 0; j < jac.extent(1); j++) {
      for (std::size_t k = 0; k < jac.extent(2); k++) {
        for (std::size_t l = 0; l < jac.extent(3); l++) {
          for (std::size_t m = 0; m < jac.extent(4); m++) {
            jac(i, j, k, l, m) = 0.0;
          }
        }
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Model::jacobians<ScalarType>(data, detJ, Jinv, Uxi, jac);

  // Free data
  delete[] Xe_data;
  delete[] Ue_data;
  delete[] detJ_data;
  delete[] Jinv_data;
  delete[] Uxi_data;
}

void compute_helmholtz_residual(int nelems, int nnodes, int* conn_data,
                                double* X_data, double* mat_data,
                                double* U_data, double* res_data) {
  typedef int IndexType;
  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef HelmholtzPDE<Basis> Model;
  typedef double ScalarType;

  const int data_per_point = Model::NUM_DATA;
  const int spatial_dim = Model::SPATIAL_DIM;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<nodes_per_elem> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<spatial_dim> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<spatial_dim> > X(node_layout, X_data);

  CLayout<vars_per_node> solution_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node> > U(solution_layout, U_data);

  CLayout<num_quad_pts, data_per_point> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, data_per_point> > data(
      element_mat_layout, mat_data);

  CLayout<nodes_per_elem, vars_per_node> res_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > res(
      res_layout, res_data);

  ScalarType* Xe_data = new ScalarType[spatial_dim * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, spatial_dim> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, spatial_dim> > Xe(
      node_element_layout, Xe_data);

  ScalarType* Ue_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> solution_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Ue(
      solution_element_layout, Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

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

  // Zero the residual array
  for (std::size_t i = 0; i < res.extent(0); i++) {
    for (std::size_t j = 0; j < res.extent(1); j++) {
      for (std::size_t k = 0; k < res.extent(2); k++) {
        res(i, j, k) = 0.0;
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Model::residuals<ScalarType>(data, detJ, Jinv, Uq, Uxi, res);

  // Free data
  delete[] Xe_data;
  delete[] Ue_data;
  delete[] detJ_data;
  delete[] Jinv_data;
  delete[] Uxi_data;
}

void compute_helmholtz_jacobian(int nelems, int nnodes, int* conn_data,
                                double* X_data, double* mat_data,
                                double* U_data, double* jac_data) {
  typedef int IndexType;
  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef HelmholtzPDE<Basis> Model;
  typedef double ScalarType;

  const int data_per_point = Model::NUM_DATA;
  const int spatial_dim = Model::SPATIAL_DIM;
  const int vars_per_node = Model::NUM_VARS;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  CLayout<nodes_per_elem> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  CLayout<spatial_dim> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<spatial_dim> > X(node_layout, X_data);

  CLayout<vars_per_node> solution_layout(nnodes);
  MultiArray<ScalarType, CLayout<vars_per_node> > U(solution_layout, U_data);

  CLayout<num_quad_pts, data_per_point> element_mat_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, data_per_point> > data(
      element_mat_layout, mat_data);

  CLayout<nodes_per_elem, nodes_per_elem, vars_per_node, vars_per_node>
      jac_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, nodes_per_elem, vars_per_node,
                                 vars_per_node> >
      jac(jac_layout, jac_data);

  ScalarType* Xe_data = new ScalarType[spatial_dim * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, spatial_dim> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, spatial_dim> > Xe(
      node_element_layout, Xe_data);

  ScalarType* Ue_data = new ScalarType[vars_per_node * nodes_per_elem * nelems];
  CLayout<nodes_per_elem, vars_per_node> solution_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, vars_per_node> > Ue(
      solution_element_layout, Ue_data);

  element_scatter(conn, X, Xe);
  element_scatter(conn, U, Ue);

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

  // Zero the Jacobian
  for (std::size_t i = 0; i < jac.extent(0); i++) {
    for (std::size_t j = 0; j < jac.extent(1); j++) {
      for (std::size_t k = 0; k < jac.extent(2); k++) {
        for (std::size_t l = 0; l < jac.extent(3); l++) {
          for (std::size_t m = 0; m < jac.extent(4); m++) {
            jac(i, j, k, l, m) = 0.0;
          }
        }
      }
    }
  }

  // Interpolate the material data to the nodes
  Basis::compute_jtrans<ScalarType>(Xe, detJ, Jinv);
  Basis::gradient<ScalarType, vars_per_node>(Ue, Uxi);
  Model::jacobians<ScalarType>(data, detJ, Jinv, Uq, Uxi, jac);

  // Free data
  delete[] Xe_data;
  delete[] Ue_data;
  delete[] detJ_data;
  delete[] Jinv_data;
  delete[] Uxi_data;
}
