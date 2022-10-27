#ifndef A2D_FE_MESH_H
#define A2D_FE_MESH_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"

namespace A2D {

const int tri_edge_nodes[3][2] = {{2, 1}, {2, 0}, {0, 1}};

/*
  Compute a node to triangle or node to quad data structure
*/
void compute_nodes_to_elem(A2D::index_t nnodes, A2D::index_t nelems,
                           A2D::index_t num_elem_nodes,
                           const A2D::index_t conn[], A2D::index_t** _ptr,
                           A2D::index_t** _node_to_elems) {
  // Set the pointer
  A2D::index_t* ptr = new A2D::index_t[nnodes + 1];
  std::fill(ptr, ptr + nnodes + 1, 0);

  // Count up the references
  const A2D::index_t conn_size = nelems * num_elem_nodes;
  for (A2D::index_t i = 0; i < conn_size; i++) {
    ptr[conn[i] + 1]++;
  }

  // Set the pointer into the array
  for (A2D::index_t i = 0; i < nnodes; i++) {
    ptr[i + 1] += ptr[i];
  }

  // Compute the node to quads
  A2D::index_t* node_to_elems = new A2D::index_t[ptr[nnodes]];
  const A2D::index_t* conn_ptr = conn;
  for (A2D::index_t i = 0; i < nelems; i++) {
    for (A2D::index_t j = 0; j < num_elem_nodes; j++) {
      A2D::index_t node = conn_ptr[0];
      if (node >= 0) {
        node_to_elems[ptr[node]] = i;
        ptr[node]++;
        conn_ptr++;
      }
    }
  }

  // Reset the pointer array
  for (A2D::index_t i = nnodes; i > 0; i--) {
    ptr[i] = ptr[i - 1];
  }
  ptr[0] = 0;

  // Set the output points
  *_ptr = ptr;
  *_node_to_elems = node_to_elems;
}

/*
  Compute all of the edges within the triangular mesh
*/
A2D::index_t compute_planar_edges(A2D::index_t nnodes, A2D::index_t ntris,
                                  const A2D::index_t tris[],
                                  A2D::index_t* tri_edge_nums,
                                  A2D::index_t* tri_orient) {
  // Compute the edges in the triangular mesh
  A2D::index_t* ptr;
  A2D::index_t* node_to_tris;
  compute_nodes_to_elem(nnodes, ntris, 3, tris, &ptr, &node_to_tris);

  // Set the no-label
  const A2D::index_t no_label = std::numeric_limits<A2D::index_t>::max();

  // Now compute the neighbors for each triangle
  for (A2D::index_t i = 0; i < 3 * ntris; i++) {
    tri_edge_nums[i] = no_label;
  }

  A2D::index_t ne = 0;
  for (A2D::index_t i = 0; i < ntris; i++) {
    // Search through each edge of the each triangle
    for (A2D::index_t j = 0; j < 3; j++) {
      if (tri_edge_nums[3 * i + j] == no_label) {
        tri_edge_nums[3 * i + j] = ne;
        tri_orient[3 * i + j] = 1;

        A2D::index_t e0[2];
        e0[0] = tris[3 * i + tri_edge_nodes[j][0]];
        e0[1] = tris[3 * i + tri_edge_nodes[j][1]];

        // Search for the neighboring that shares this edge
        A2D::index_t kp = ptr[e0[0]];
        A2D::index_t kpend = ptr[e0[0] + 1];
        for (; kp < kpend; kp++) {
          // Find the potential neighbor
          A2D::index_t n = node_to_tris[kp];

          // Don't count the same edge twice
          if (n == i) {
            continue;
          }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is
          // planar in parameter space)
          A2D::index_t quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for (A2D::index_t e = 0; e < 3; e++) {
            A2D::index_t e1[2];
            e1[0] = tris[3 * n + tri_edge_nodes[e][0]];
            e1[1] = tris[3 * n + tri_edge_nodes[e][1]];

            // Check if the adjacent edge matches in either direction
            if ((e0[0] == e1[0] && e0[1] == e1[1]) ||
                (e0[0] == e1[1] && e0[1] == e1[0])) {
              // Label the other edge that shares this same node
              tri_edge_nums[3 * n + e] = ne;

              // Opposite orientation
              tri_orient[3 * i + j] = 0;

              quit = 1;
            }
          }
          if (quit) {
            break;
          }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the data that is no longer required
  delete[] ptr;
  delete[] node_to_tris;

  return ne;
}

/*
  The element connectivity class

  This class will need a substantial overhaul
  1. The connectivity is fixed at this point
  2.
*/
class ElementConnectivity {
 public:
  ElementConnectivity(A2D::index_t nnodes, A2D::index_t nelems,
                      A2D::index_t* conn_)
      : nnodes(nnodes), nelems(nelems) {
    conn = new A2D::index_t[3 * nelems];
    for (A2D::index_t i = 0; i < 3 * nelems; i++) {
      conn[i] = conn_[i];
    }

    edges = new A2D::index_t[3 * nelems];
    orient = new A2D::index_t[3 * nelems];
    nedges = compute_planar_edges(nnodes, nelems, conn, edges, orient);
  }

  A2D::index_t get_face_dof(A2D::index_t elem, A2D::index_t index) {
    return elem;
  }

  A2D::index_t get_edge_dof(A2D::index_t elem, A2D::index_t index, int& ort) {
    if (orient[3 * elem + index]) {
      ort = 1;
    } else {
      ort = -1;
    }

    return edges[3 * elem + index];
  }

  A2D::index_t get_node_dof(A2D::index_t elem, A2D::index_t index) {
    return conn[3 * elem + index];
  }

  A2D::index_t get_num_elements() { return nelems; }
  A2D::index_t get_num_edges() { return nedges; }
  A2D::index_t get_num_nodes() { return nnodes; }

 private:
  A2D::index_t nnodes;
  A2D::index_t nedges;
  A2D::index_t nelems;
  A2D::index_t* conn;
  A2D::index_t* edges;
  A2D::index_t* orient;  // edge orientation 1 for +ve, 0 for -ve
};

/*
  ElementMesh - Map from an element to the global to element local degrees
  of freedom
*/
enum SpaceType { L2, H1, EDGE };

template <class Basis>
class ElementMesh {
 public:
  ElementMesh(A2D::ElementConnectivity& conn, SpaceType spaces_[],
              A2D::index_t dim_[] = NULL)
      : conn(conn) {
    for (A2D::index_t i = 0; i < nspaces; i++) {
      spaces[i] = spaces_[i];
    }
    if (dim_) {
      for (A2D::index_t i = 0; i < nspaces; i++) {
        dim[i] = dim_[i];
      }
    } else {
      for (A2D::index_t i = 0; i < nspaces; i++) {
        dim[i] = 1;
      }
    }

    offset[0] = 0;

    for (A2D::index_t index = 0; index < nspaces; index++) {
      if (spaces[index] == L2) {
        counts[index] = dim[index] * conn.get_num_elements();
      } else if (spaces[index] == H1) {
        counts[index] = dim[index] * conn.get_num_nodes();
      } else if (spaces[index] == EDGE) {
        counts[index] = dim[index] * conn.get_num_edges();
      }

      offset[index + 1] = offset[index] + counts[index];
    }
  }

  static const A2D::index_t nspaces = Basis::nbasis;

  A2D::index_t get_num_elements() { return conn.get_num_elements(); }
  A2D::index_t get_num_dof() { return offset[nspaces]; }

  int get_global_dof_sign(A2D::index_t elem, A2D::index_t space,
                          A2D::index_t index) {
    if (spaces[space] == EDGE) {
      int orient;
      conn.get_edge_dof(elem, index, orient);
      return orient;
    }
    return 1;
  }
  A2D::index_t get_global_dof(A2D::index_t elem, A2D::index_t space,
                              A2D::index_t index) {
    if (spaces[space] == L2) {
      return offset[space] +
             dim[space] * conn.get_face_dof(elem, index / dim[space]) +
             index % dim[space];
    } else if (spaces[space] == H1) {
      return offset[space] +
             dim[space] * conn.get_node_dof(elem, index / dim[space]) +
             index % dim[space];
    } else if (spaces[space] == EDGE) {
      int orient;
      return offset[space] +
             dim[space] * conn.get_edge_dof(elem, index / dim[space], orient) +
             index % dim[space];
    }
    return 0;
  }

 private:
  A2D::ElementConnectivity conn;
  A2D::index_t counts[nspaces];
  A2D::index_t offset[nspaces + 1];
  SpaceType spaces[nspaces];
  A2D::index_t dim[nspaces];
};

}  // namespace A2D

#endif  // A2D_FE_MESH_H