#include "multiarray.h"
#include "elasticity3d.h"
#include <cstdint>
#include <iostream>

/*
  Create a slice of a multi-dimensional array
*/
int main( int argc, char *argv[] ){
  typedef int32_t IndexType;
  typedef double ScalarType;
  typedef HexQuadrature Quadrature;
  typedef HexBasis<HexQuadrature> Basis;
  typedef NonlinearElasticity3D Model;

  const int nx = 128;
  const int ny = 64;
  const int nz = 64;
  const int nnodes = (nx + 1) * (ny + 1) * (nz + 1);
  const int nelems = nx * ny * nz;

  const int vars_per_node = 3;
  const int nodes_per_elem = Basis::NUM_NODES;
  const int num_quad_pts = Quadrature::NUM_QUAD_PTS;

  IndexType* conn_data = new IndexType[ nodes_per_elem * nelems ];
  CLayout<8> conn_layout(nelems);
  MultiArray<IndexType, CLayout<nodes_per_elem> > conn(conn_layout, conn_data);

  ScalarType* X_data = new ScalarType[ 3 * nnodes ];
  CLayout<3> node_layout(nnodes);
  MultiArray<ScalarType, CLayout<3> > X(node_layout, X_data);

  ScalarType* Xe_data = new ScalarType[ 3 * nodes_per_elem * nelems ];
  CLayout<nodes_per_elem, 3> node_element_layout(nelems);
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xe(node_element_layout, Xe_data);

  ScalarType* Xq_data = new ScalarType[ 3 * nodes_per_elem * nelems ];
  MultiArray<ScalarType, CLayout<nodes_per_elem, 3> > Xq(node_element_layout, Xq_data);

  // Set the node locations
  for ( int k = 0; k < nz + 1; k++ ){
    for ( int j = 0; j < ny + 1; j++ ){
      for ( int i = 0; i < nx + 1; i++ ){
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        X(node, 0) = 1.0 * i/nx;
        X(node, 1) = 1.0 * j/ny;
        X(node, 2) = 1.0 * k/nz;
      }
    }
  }

  // Set the element connectivity
  for ( int k = 0; k < nz; k++ ){
    for ( int j = 0; j < ny; j++ ){
      for ( int i = 0; i < nx; i++ ){
        int elem = i + nx * (j + ny * k);

        for ( int kk = 0, index = 0; kk < 2; kk++ ){
          for ( int jj = 0; jj < 2; jj++ ){
            for ( int ii = 0; ii < 2; ii++, index++ ){
              conn(elem, index) =
                (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
            }
          }
        }
      }
    }
  }

  element_scatter(conn, X, Xe);
  // element_scatter(conn, U, Ue);

  // Store the data for Jinv
  ScalarType* Jinv_data = new ScalarType[ 3 * 3 * num_quad_pts * nelems ];
  CLayout<num_quad_pts, 3, 3> element_grad_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Jinv(element_grad_layout, Jinv_data);

  // Store data for detJ
  ScalarType* detJ_data = new ScalarType[ num_quad_pts * nelems ];
  CLayout<num_quad_pts> element_detJ_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts> > detJ(element_detJ_layout, detJ_data);

  // Store data for Ux
  ScalarType* Ux_data = new ScalarType[ 3 * 3 * num_quad_pts * nelems ];
  CLayout<num_quad_pts, 3, 3> element_Ux_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 3, 3> > Ux(element_Ux_layout, Ux_data);

  // Allocate space for the data
  ScalarType* Eq_data = new ScalarType[ 2 * num_quad_pts * nelems ];
  CLayout<num_quad_pts, 2> element_Eq_layout(nelems);
  MultiArray<ScalarType, CLayout<num_quad_pts, 2> > Eq(element_Eq_layout, Eq_data);

  Basis::interp<vars_per_node>(Xe, Xq);
  Basis::compute_jtrans<ScalarType>(Xe, Jinv, detJ);
  Basis::gradient<ScalarType, vars_per_node>(Xe, Jinv, Ux);

  ScalarType energy;
  Basis::energy<ScalarType, Model>(Ux, Eq, Jinv, detJ, energy);

  return (0);
}
