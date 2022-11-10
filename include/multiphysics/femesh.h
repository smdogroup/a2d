#ifndef A2D_FE_MESH_H
#define A2D_FE_MESH_H

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"
#include "multiphysics/feelementtypes.h"

namespace A2D {

/**
 * @brief Mesh connecivity class for 3D meshes composed of tetrahedral,
 * hexahedral, wedge and pyramid elements
 *
 */
class MeshConnectivity3D {
 public:
  static constexpr index_t NO_LABEL = std::numeric_limits<index_t>::max();

  // Use the definitions from the element types
  using ET = ElementTypes;

  template <typename I>
  MeshConnectivity3D(I nverts, I ntets, I* tets, I nhex, I* hex, I nwedge,
                     I* wedge, I npyrmd, I* pyrmd)
      : nverts(nverts),
        ntets(ntets),
        nhex(nhex),
        nwedge(nwedge),
        npyrmd(npyrmd) {
    // Allocate space for the element -> vert connectivity
    tet_verts = new index_t[ET::TET_VERTS * ntets];
    hex_verts = new index_t[ET::HEX_VERTS * nhex];
    wedge_verts = new index_t[ET::WEDGE_VERTS * nwedge];
    pyrmd_verts = new index_t[ET::PYRMD_VERTS * npyrmd];

    // vert -> element connectivity
    vert_element_ptr = NULL;
    vert_elements = NULL;

    // element -> face connectivity
    tet_faces = NULL;
    hex_faces = NULL;
    wedge_faces = NULL;
    pyrmd_faces = NULL;
    tri_face_elements = NULL;
    quad_face_elements = NULL;

    // face -> edge connectivity
    tri_edges = NULL;
    quad_edges = NULL;

    // Count the total number of elements
    nelems = ntets + nhex + nwedge + npyrmd;

    // Set the connectivity: element -> verts
    for (index_t i = 0; i < ET::TET_VERTS * ntets; i++) {
      tet_verts[i] = tets[i];
    }
    for (index_t i = 0; i < ET::HEX_VERTS * nhex; i++) {
      hex_verts[i] = hex[i];
    }
    for (index_t i = 0; i < ET::WEDGE_VERTS * nwedge; i++) {
      wedge_verts[i] = wedge[i];
    }
    for (index_t i = 0; i < ET::PYRMD_VERTS * npyrmd; i++) {
      pyrmd_verts[i] = pyrmd[i];
    }

    init_vert_element_data();
    init_face_data();
    init_edge_data();
  }
  ~MeshConnectivity3D() {
    if (tet_verts) {
      delete[] tet_verts;
    }
    if (hex_verts) {
      delete[] hex_verts;
    }
    if (wedge_verts) {
      delete[] wedge_verts;
    }
    if (pyrmd_verts) {
      delete[] pyrmd_verts;
    }

    if (vert_element_ptr) {
      delete[] vert_element_ptr;
    }
    if (vert_elements) {
      delete[] vert_elements;
    }

    if (tet_faces) {
      delete[] tet_faces;
    }
    if (hex_faces) {
      delete[] hex_faces;
    }
    if (wedge_faces) {
      delete[] wedge_faces;
    }
    if (pyrmd_faces) {
      delete[] pyrmd_faces;
    }
    if (tri_face_elements) {
      delete[] tri_face_elements;
    }
    if (quad_face_elements) {
      delete[] quad_face_elements;
    }

    if (tri_edges) {
      delete[] tri_edges;
    }
    if (quad_edges) {
      delete[] quad_edges;
    }
  }

  // Get global element counts
  index_t get_num_elements() { return nelems; }
  index_t get_num_faces() { return nfaces; }
  index_t get_num_edges() { return nedges; }
  index_t get_num_verts() { return nverts; }

  /**
   * @brief Get the element faces array
   *
   * @param elem The element index
   * @param faces The face indicies associated with the element
   * @return The number of faces
   */
  index_t get_element_faces(index_t elem, const index_t* faces[]) {
    if (elem < ntets) {
      *faces = &tet_faces[ET::TET_FACES * elem];
      return ET::TET_FACES;
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      *faces = &hex_faces[ET::HEX_FACES * elem];
      return ET::HEX_FACES;
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      *faces = &wedge_faces[ET::WEDGE_FACES * elem];
      return ET::WEDGE_FACES;
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      *faces = &pyrmd_faces[ET::PYRMD_FACES * elem];
      return ET::PYRMD_FACES;
    }
    *faces = NULL;
    return 0;
  }

  /**
   * @brief Get the element verts
   *
   * @param elem The element index
   * @param verts The vertex indices
   * @return The number of vertices
   */
  index_t get_element_verts(index_t elem, const index_t* verts[]) {
    if (elem < ntets) {
      *verts = &tet_verts[ET::TET_VERTS * elem];
      return ET::TET_VERTS;
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      *verts = &hex_verts[ET::HEX_VERTS * elem];
      return ET::HEX_VERTS;
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      *verts = &wedge_verts[ET::WEDGE_VERTS * elem];
      return ET::WEDGE_VERTS;
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      *verts = &pyrmd_verts[ET::PYRMD_VERTS * elem];
      return ET::PYRMD_VERTS;
    }
    *verts = NULL;
    return 0;
  }

  /**
   * @brief Get the adjacent elements from the vert index
   *
   * @param vert The index of the vertex
   * @param elems Array of the adjacent element indices
   * @return The number of elements
   */
  index_t get_adjacent_elements_from_vert(const index_t vert,
                                          const index_t* elems[]) {
    if (vert < nverts) {
      *elems = &vert_elements[vert_element_ptr[vert]];
      return vert_element_ptr[vert + 1] - vert_element_ptr[vert];
    }
    elems = NULL;
    return 0;
  }

  /**
   * @brief Get the indices of the face in its local order
   *
   * @param elem The element index
   * @param f The local face index
   * @param verts The vertices
   * @return The number of vertices
   */
  index_t get_element_face_verts(index_t elem, index_t f, index_t verts[]) {
    if (elem < ntets) {
      verts[1] = tet_verts[ET::TET_VERTS * elem + ET::TET_FACE_VERTS[f][1]];
      verts[2] = tet_verts[ET::TET_VERTS * elem + ET::TET_FACE_VERTS[f][2]];
      verts[0] = tet_verts[ET::TET_VERTS * elem + ET::TET_FACE_VERTS[f][0]];
      return 3;
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      verts[0] = hex_verts[ET::HEX_VERTS * elem + ET::HEX_FACE_VERTS[f][0]];
      verts[1] = hex_verts[ET::HEX_VERTS * elem + ET::HEX_FACE_VERTS[f][1]];
      verts[2] = hex_verts[ET::HEX_VERTS * elem + ET::HEX_FACE_VERTS[f][2]];
      verts[3] = hex_verts[ET::HEX_VERTS * elem + ET::HEX_FACE_VERTS[f][3]];
      return 4;
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;

      if (f < ET::WEDGE_TRI_FACES) {
        verts[0] = wedge_verts[ET::WEDGE_VERTS * elem +
                               ET::WEDGE_TRI_FACE_VERTS[f][0]];
        verts[1] = wedge_verts[ET::WEDGE_VERTS * elem +
                               ET::WEDGE_TRI_FACE_VERTS[f][1]];
        verts[2] = wedge_verts[ET::WEDGE_VERTS * elem +
                               ET::WEDGE_TRI_FACE_VERTS[f][2]];
        return 3;
      } else {
        verts[0] =
            wedge_verts[ET::WEDGE_VERTS * elem +
                        ET::WEDGE_QUAD_FACE_VERTS[f - ET::WEDGE_TRI_FACES][0]];
        verts[1] =
            wedge_verts[ET::WEDGE_VERTS * elem +
                        ET::WEDGE_QUAD_FACE_VERTS[f - ET::WEDGE_TRI_FACES][1]];
        verts[2] =
            wedge_verts[ET::WEDGE_VERTS * elem +
                        ET::WEDGE_QUAD_FACE_VERTS[f - ET::WEDGE_TRI_FACES][2]];
        verts[3] =
            wedge_verts[ET::WEDGE_VERTS * elem +
                        ET::WEDGE_QUAD_FACE_VERTS[f - ET::WEDGE_TRI_FACES][3]];
        return 4;
      }
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;

      if (f < ET::PYRMD_TRI_FACES) {
        verts[0] = pyrmd_verts[ET::PYRMD_VERTS * elem +
                               ET::PYRMD_TRI_FACE_VERTS[f][0]];
        verts[1] = pyrmd_verts[ET::PYRMD_VERTS * elem +
                               ET::PYRMD_TRI_FACE_VERTS[f][1]];
        verts[2] = pyrmd_verts[ET::PYRMD_VERTS * elem +
                               ET::PYRMD_TRI_FACE_VERTS[f][2]];
        return 3;
      } else {
        verts[0] =
            pyrmd_verts[ET::PYRMD_VERTS * elem + ET::PYRMD_QUAD_FACE_VERTS[0]];
        verts[1] =
            pyrmd_verts[ET::PYRMD_VERTS * elem + ET::PYRMD_QUAD_FACE_VERTS[1]];
        verts[2] =
            pyrmd_verts[ET::PYRMD_VERTS * elem + ET::PYRMD_QUAD_FACE_VERTS[2]];
        verts[3] =
            pyrmd_verts[ET::PYRMD_VERTS * elem + ET::PYRMD_QUAD_FACE_VERTS[3]];
        return 4;
      }
    }
    return 0;
  }

  /**
   * @brief Get the face element vertices in a sorted order
   *
   * This call will return the vertices in the same order regardless of the
   * orientation of the face. This is used to compare faces.
   *
   * @param elem The element index
   * @param f The local face index
   * @param verts The vertices
   * @return The number of vertices
   */
  index_t get_element_global_face_verts(index_t elem, index_t local_face,
                                        index_t verts[]) {
    index_t t[4];
    index_t n = get_element_face_verts(elem, local_face, t);

    if (n == 3) {
      // Find the smallest vertex index
      if (t[0] < t[1] && t[0] < t[2]) {  // t[0] is the smallest
        if (t[1] < t[2]) {
          verts[0] = t[0];
          verts[1] = t[1];
          verts[2] = t[2];
        } else {  // Reverse the order
          verts[0] = t[0];
          verts[1] = t[2];
          verts[2] = t[1];
        }
      } else if (t[1] < t[0] && t[1] < t[2]) {  // t[1] is the smallest
        if (t[2] < t[0]) {
          verts[0] = t[1];
          verts[1] = t[2];
          verts[2] = t[0];
        } else {  // Reverse the order
          verts[0] = t[1];
          verts[1] = t[0];
          verts[2] = t[2];
        }
      } else {  // t[2] is the smallest
        if (t[0] < t[1]) {
          verts[0] = t[2];
          verts[1] = t[0];
          verts[2] = t[1];
        } else {  // Reverse the order
          verts[0] = t[2];
          verts[1] = t[1];
          verts[2] = t[0];
        }
      }

      return 3;
    } else {                                            // n == 4
      if (t[0] < t[1] && t[0] < t[2] && t[0] < t[3]) {  // t[0] is smallest
        if (t[1] < t[3]) {
          verts[0] = t[0];
          verts[1] = t[1];
          verts[2] = t[2];
          verts[3] = t[3];
        } else {
          verts[0] = t[0];
          verts[1] = t[3];
          verts[2] = t[2];
          verts[3] = t[1];
        }
      } else if (t[1] < t[0] && t[1] < t[2] &&
                 t[1] < t[3]) {  // t[1] is smallest
        if (t[2] < t[0]) {
          verts[0] = t[1];
          verts[1] = t[2];
          verts[2] = t[3];
          verts[3] = t[0];
        } else {
          verts[0] = t[1];
          verts[1] = t[0];
          verts[2] = t[3];
          verts[3] = t[2];
        }
      } else if (t[2] < t[0] && t[2] < t[1] &&
                 t[2] < t[3]) {  // t[2] is smallest
        if (t[3] < t[1]) {
          verts[0] = t[2];
          verts[1] = t[3];
          verts[2] = t[0];
          verts[3] = t[1];
        } else {
          verts[0] = t[2];
          verts[1] = t[1];
          verts[2] = t[0];
          verts[3] = t[3];
        }
      } else {  // t[3] is smallest
        if (t[0] < t[2]) {
          verts[0] = t[3];
          verts[1] = t[0];
          verts[2] = t[1];
          verts[3] = t[2];
        } else {
          verts[0] = t[3];
          verts[1] = t[2];
          verts[2] = t[1];
          verts[3] = t[0];
        }
      }

      return 4;
    }
  }

  /**
   * @brief Check for equality between two faces
   *
   * @return True if the faces match, false otherwise
   */
  bool global_face_equality(index_t na, const index_t a[], index_t nb,
                            const index_t b[]) {
    if (na == nb) {
      if (na == 3) {
        if (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]) {
          return true;
        }
      } else {
        if (a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3]) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * @brief Get the face elements associated with the given face
   *
   * @param face Global face index
   * @param e1 Returned value of the first element
   * @param e2 Returned value of the second element (if it exists)
   * @return Boolean if this face is on the boundary
   */
  bool get_face_elements(index_t face, index_t* e1, index_t* e2) {
    if (face < ntri_faces) {
      *e1 = tri_face_elements[2 * face];
      *e2 = tri_face_elements[2 * face + 1];
      return (tri_face_elements[2 * face + 1] == NO_LABEL);
    } else {
      face = face - ntri_faces;
      *e1 = quad_face_elements[2 * face];
      *e2 = quad_face_elements[2 * face + 1];
      return (quad_face_elements[2 * face + 1] == NO_LABEL);
    }
  }

  /**
   * @brief Get the edges associated with the global face index
   *
   * @param face Input face index
   * @param edges The edges associated with the face index
   * @return The number of edges
   */
  index_t get_face_edges(index_t face, const index_t* edges[]) {
    if (face < ntri_faces) {
      *edges = &tri_edges[ET::TRI_EDGES * face];
      return ET::TRI_EDGES;
    } else {
      face = face - ntri_faces;
      *edges = &quad_edges[ET::QUAD_EDGES * face];
      return ET::QUAD_EDGES;
    }
  }

  /**
   * @brief Get the number of elements associated with an edge
   *
   */
  index_t get_num_element_edges(index_t elem) {
    if (elem < ntets) {
      return ET::TET_EDGES;
    } else if (elem < ntets + nhex) {
      return ET::HEX_EDGES;
    } else if (elem < ntets + nhex + nwedge) {
      return ET::WEDGE_EDGES;
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      return ET::PYRMD_EDGES;
    }
    return 0;
  }

  /**
   * @brief Get the element edges
   *
   * @param elem The element index
   * @param edges The edge indices
   * @return The number of edges for this element
   */
  index_t get_element_edges(index_t elem, index_t edges[]) {
    const index_t* faces;
    index_t nf = get_element_faces(elem, &faces);
    if (elem < ntets) {
      for (index_t i = 0; i < ET::TET_EDGES; i++) {
        index_t f0 = faces[ET::TET_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::TET_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);
        edges[i] = f0_edges[e0];
      }
      return ET::TET_EDGES;
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      for (index_t i = 0; i < ET::HEX_EDGES; i++) {
        index_t f0 = faces[ET::HEX_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::HEX_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);
        edges[i] = f0_edges[e0];
      }
      return ET::HEX_EDGES;
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      for (index_t i = 0; i < ET::WEDGE_EDGES; i++) {
        index_t f0 = faces[ET::WEDGE_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::WEDGE_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);
        edges[i] = f0_edges[e0];
      }
      return ET::WEDGE_EDGES;
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      for (index_t i = 0; i < ET::PYRMD_EDGES; i++) {
        index_t f0 = faces[ET::PYRMD_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::PYRMD_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);
        edges[i] = f0_edges[e0];
      }
      return ET::PYRMD_EDGES;
    }
    return 0;
  }

  /**
   * @brief Get the vertices of the associated element edge
   */
  void get_element_edge_verts(index_t elem, index_t edge, index_t verts[]) {
    if (elem < ntets) {
      verts[0] = tet_verts[ET::TET_EDGES * elem + ET::TET_EDGE_VERTS[edge][0]];
      verts[1] = tet_verts[ET::TET_EDGES * elem + ET::TET_EDGE_VERTS[edge][1]];
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      verts[0] = hex_verts[ET::HEX_EDGES * elem + ET::HEX_EDGE_VERTS[edge][0]];
      verts[1] = hex_verts[ET::HEX_EDGES * elem + ET::HEX_EDGE_VERTS[edge][1]];
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      verts[0] =
          wedge_verts[ET::WEDGE_EDGES * elem + ET::WEDGE_EDGE_VERTS[edge][0]];
      verts[1] =
          wedge_verts[ET::WEDGE_EDGES * elem + ET::WEDGE_EDGE_VERTS[edge][1]];
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      verts[0] =
          pyrmd_verts[ET::PYRMD_EDGES * elem + ET::PYRMD_EDGE_VERTS[edge][0]];
      verts[1] =
          pyrmd_verts[ET::PYRMD_EDGES * elem + ET::PYRMD_EDGE_VERTS[edge][1]];
    }
  }

  /**
   * @brief Get the edges associated with the given element and local face index
   *
   * @param elem The element index
   * @param f The local face number
   * @param edges The edges associated with the local face
   * @return The number of edges for the given face
   */
  index_t get_element_face_edges(index_t elem, index_t f,
                                 const index_t* edges[]) {
    const index_t* faces;
    get_element_faces(elem, &faces);
    index_t face = faces[f];

    if (face < ntri_faces) {
      *edges = &tri_edges[ET::TRI_EDGES * face];
      return ET::TRI_EDGES;
    } else {
      face = face - ntri_faces;
      *edges = &quad_edges[ET::QUAD_EDGES * face];
      return ET::QUAD_EDGES;
    }
  }

 private:
  /**
   * @brief Initialize the data for the connection between vertices and
   * the elements
   *
   * Given a vertex index, this enables retrieval of all elements that reference
   * that vertex.
   */
  void init_vert_element_data() {
    // Construct vert -> element data
    vert_element_ptr = new index_t[nverts + 1];
    std::fill(vert_element_ptr, vert_element_ptr + nverts + 1, 0);

    for (index_t elem = 0; elem < nelems; elem++) {
      const index_t* verts;
      index_t nv = get_element_verts(elem, &verts);
      for (index_t i = 0; i < nv; i++) {
        vert_element_ptr[verts[i] + 1]++;
      }
    }

    for (index_t i = 0; i < nverts; i++) {
      vert_element_ptr[i + 1] += vert_element_ptr[i];
    }

    vert_elements = new index_t[vert_element_ptr[nverts]];

    for (index_t elem = 0; elem < nelems; elem++) {
      const index_t* verts;
      index_t nv = get_element_verts(elem, &verts);
      for (index_t i = 0; i < nv; i++) {
        vert_elements[vert_element_ptr[verts[i]]] = elem;
        vert_element_ptr[verts[i]]++;
      }
    }

    for (index_t i = nverts; i > 0; i--) {
      vert_element_ptr[i] = vert_element_ptr[i - 1];
    }
    vert_element_ptr[0] = 0;
  }

  /**
   * @brief Get a non-constant element faces array
   *
   * @param elem The element index
   * @param faces The face indicies associated with the element
   * @return The number of faces
   */
  index_t get_element_faces(index_t elem, index_t* faces[]) {
    if (elem < ntets) {
      *faces = &tet_faces[ET::TET_FACES * elem];
      return ET::TET_FACES;
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      *faces = &hex_faces[ET::HEX_FACES * elem];
      return ET::HEX_FACES;
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      *faces = &wedge_faces[ET::WEDGE_FACES * elem];
      return ET::WEDGE_FACES;
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      *faces = &pyrmd_faces[ET::PYRMD_FACES * elem];
      return ET::PYRMD_FACES;
    }
    *faces = NULL;
    return 0;
  }

  /**
   * @brief Get the edges associated with the global face index
   *
   * @param face Input face index
   * @param edges The edges associated with the face index
   * @return The number of edges
   */
  index_t get_face_edges(index_t face, index_t* edges[]) {
    if (face < ntri_faces) {
      *edges = &tri_edges[ET::TRI_EDGES * face];
      return ET::TRI_EDGES;
    } else {
      face = face - ntri_faces;
      *edges = &quad_edges[ET::QUAD_EDGES * face];
      return ET::QUAD_EDGES;
    }
  }

  /**
   * @brief Initialize data associated with the face information
   *
   * This code uniquely orders the faces assocaited with each element
   */
  void init_face_data() {
    if (tet_faces) {
      delete[] tet_faces;
    }
    if (hex_faces) {
      delete[] hex_faces;
    }
    if (wedge_faces) {
      delete[] wedge_faces;
    }
    if (pyrmd_faces) {
      delete[] pyrmd_faces;
    }

    // Allocate the face data
    tet_faces = new index_t[ET::TET_FACES * ntets];
    hex_faces = new index_t[ET::HEX_FACES * nhex];
    wedge_faces = new index_t[ET::WEDGE_FACES * nwedge];
    pyrmd_faces = new index_t[ET::PYRMD_FACES * npyrmd];

    // Fill the face arrays with NO_LABEL
    std::fill(tet_faces, tet_faces + ET::TET_FACES * ntets, NO_LABEL);
    std::fill(hex_faces, hex_faces + ET::HEX_FACES * nhex, NO_LABEL);
    std::fill(wedge_faces, wedge_faces + ET::WEDGE_FACES * nwedge, NO_LABEL);
    std::fill(pyrmd_faces, pyrmd_faces + ET::PYRMD_FACES * npyrmd, NO_LABEL);

    // Prepare to count and number the number of faces. This keeps track of
    // separate triangle and quadrilateral face counts
    ntri_faces = 0;
    nquad_faces = 0;
    nfaces = 0;

    for (index_t elem = 0; elem < nelems; elem++) {
      // Loop over the elements of this face
      index_t* faces;
      const index_t nf = get_element_faces(elem, &faces);

      for (index_t face = 0; face < nf; face++) {
        if (faces[face] == NO_LABEL) {
          // Get the unique set of verts corresponding to this face
          index_t face_verts[4];
          index_t nface_verts =
              get_element_global_face_verts(elem, face, face_verts);

          bool face_located = false;

          // For all adjacent elements check if there is a match
          for (index_t m = 0; m < nface_verts; m++) {
            const index_t* adj_elems;
            const index_t nadj_elems =
                get_adjacent_elements_from_vert(face_verts[m], &adj_elems);

            // Loop over all adjacent elements that are not the current element
            for (index_t k = 0; k < nadj_elems; k++) {
              if (adj_elems[k] != elem) {
                // Loop over the elements of this face
                index_t* adj_faces;
                const index_t adj_nf =
                    get_element_faces(adj_elems[k], &adj_faces);

                for (index_t adj_face = 0; adj_face < adj_nf; adj_face++) {
                  // Get the unique set of vertices corresponding to this face
                  index_t adj_face_verts[4];
                  index_t adj_nface_verts = get_element_global_face_verts(
                      adj_elems[k], adj_face, adj_face_verts);

                  if (adj_faces[adj_face] == NO_LABEL &&
                      global_face_equality(nface_verts, face_verts,
                                           adj_nface_verts, adj_face_verts)) {
                    if (nface_verts == 3) {
                      adj_faces[adj_face] = ntri_faces;
                      faces[face] = ntri_faces;
                      ntri_faces++;
                    } else {
                      adj_faces[adj_face] = nquad_faces;
                      faces[face] = nquad_faces;
                      nquad_faces++;
                    }

                    face_located = true;
                    break;
                  }
                }

                if (face_located) {
                  break;
                }
              }
            }

            if (face_located) {
              break;
            }
          }

          // No adjacent face was found. This is a boundary face
          if (!face_located) {
            if (nface_verts == 3) {
              faces[face] = ntri_faces;
              ntri_faces++;
            } else {
              faces[face] = nquad_faces;
              nquad_faces++;
            }

            // Add this face to the boundary face list???
          }
        }
      }
    }

    // Sum up the total number of faces
    nfaces = ntri_faces + nquad_faces;

    // Set the face indices associated with the two adjacent elements. At this
    // point the face indices stored are with respect to separate lists of
    // triangular and quadrilateral faces. We order the triangular faces first,
    // so we have to add the total number of triangular faces to each
    // quadrilateral face index to get the global index.
    tri_face_elements = new index_t[2 * ntri_faces];
    std::fill(tri_face_elements, tri_face_elements + 2 * ntri_faces, NO_LABEL);

    quad_face_elements = new index_t[2 * nquad_faces];
    std::fill(quad_face_elements, quad_face_elements + 2 * nquad_faces,
              NO_LABEL);

    for (index_t elem = 0; elem < nelems; elem++) {
      // Loop over the elements of this face
      index_t* faces;
      const index_t nf = get_element_faces(elem, &faces);

      for (index_t face = 0; face < nf; face++) {
        index_t face_verts[4];
        index_t nface_verts =
            get_element_global_face_verts(elem, face, face_verts);

        if (nface_verts == 3) {
          if (tri_face_elements[2 * faces[face]] == NO_LABEL) {
            tri_face_elements[2 * faces[face]] = elem;
          } else if (tri_face_elements[2 * faces[face] + 1] == NO_LABEL) {
            tri_face_elements[2 * faces[face] + 1] = elem;
          }
        } else {
          if (quad_face_elements[2 * faces[face]] == NO_LABEL) {
            quad_face_elements[2 * faces[face]] = elem;
          } else if (quad_face_elements[2 * faces[face] + 1] == NO_LABEL) {
            quad_face_elements[2 * faces[face] + 1] = elem;
          }

          // Reset the face index into the global face index - add the number of
          // triangle faces. Now the faces will index from a global face number.
          faces[face] += ntri_faces;
        }
      }
    }
  }

  /**
   * @brief Make the edge numbering consistent on an element by copying any
   * values not set to NO_LABEL to the other common element edges. If no common
   * edges have NO_LABEL then do nothing.
   *
   * @param elem The element index
   */
  void make_element_edge_consistent(index_t elem) {
    const index_t* faces;
    index_t nf = get_element_faces(elem, &faces);
    if (elem < ntets) {
      for (index_t i = 0; i < ET::TET_EDGES; i++) {
        index_t f0 = faces[ET::TET_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::TET_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);

        index_t f1 = faces[ET::TET_EDGE_TO_ADJ_FACES[i][1]];
        index_t e1 = ET::TET_EDGE_TO_ADJ_FACE_EDGE[i][1];
        index_t* f1_edges;
        get_face_edges(f1, &f1_edges);

        if (f0_edges[e0] != f1_edges[e1]) {
          if (f0_edges[e0] == NO_LABEL) {
            f0_edges[e0] = f1_edges[e1];
          } else {
            f1_edges[e1] = f0_edges[e0];
          }
        }
      }
    } else if (elem < ntets + nhex) {
      elem = elem - ntets;
      for (index_t i = 0; i < ET::HEX_EDGES; i++) {
        index_t f0 = faces[ET::HEX_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::HEX_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);

        index_t f1 = faces[ET::HEX_EDGE_TO_ADJ_FACES[i][1]];
        index_t e1 = ET::HEX_EDGE_TO_ADJ_FACE_EDGE[i][1];
        index_t* f1_edges;
        get_face_edges(f1, &f1_edges);

        if (f0_edges[e0] != f1_edges[e1]) {
          if (f0_edges[e0] == NO_LABEL) {
            f0_edges[e0] = f1_edges[e1];
          } else {
            f1_edges[e1] = f0_edges[e0];
          }
        }
      }
    } else if (elem < ntets + nhex + nwedge) {
      elem = elem - ntets - nhex;
      for (index_t i = 0; i < ET::WEDGE_EDGES; i++) {
        index_t f0 = faces[ET::WEDGE_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::WEDGE_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);

        index_t f1 = faces[ET::WEDGE_EDGE_TO_ADJ_FACES[i][1]];
        index_t e1 = ET::WEDGE_EDGE_TO_ADJ_FACE_EDGE[i][1];
        index_t* f1_edges;
        get_face_edges(f1, &f1_edges);

        if (f0_edges[e0] != f1_edges[e1]) {
          if (f0_edges[e0] == NO_LABEL) {
            f0_edges[e0] = f1_edges[e1];
          } else {
            f1_edges[e1] = f0_edges[e0];
          }
        }
      }
    } else if (elem < ntets + nhex + nwedge + npyrmd) {
      elem = elem - ntets - nhex - nwedge;
      for (index_t i = 0; i < ET::PYRMD_EDGES; i++) {
        index_t f0 = faces[ET::PYRMD_EDGE_TO_ADJ_FACES[i][0]];
        index_t e0 = ET::PYRMD_EDGE_TO_ADJ_FACE_EDGE[i][0];
        index_t* f0_edges;
        get_face_edges(f0, &f0_edges);

        index_t f1 = faces[ET::PYRMD_EDGE_TO_ADJ_FACES[i][1]];
        index_t e1 = ET::PYRMD_EDGE_TO_ADJ_FACE_EDGE[i][1];
        index_t* f1_edges;
        get_face_edges(f1, &f1_edges);

        if (f0_edges[e0] != f1_edges[e1]) {
          if (f0_edges[e0] == NO_LABEL) {
            f0_edges[e0] = f1_edges[e1];
          } else {
            f1_edges[e1] = f0_edges[e0];
          }
        }
      }
    }
  }

  /**
   * @brief Initialize and order the edge information.
   *
   * This relies on the face connectivity data - so that must be initialized
   * first.
   */
  void init_edge_data() {
    nedges = 0;

    if (tri_edges) {
      delete[] tri_edges;
    }
    if (quad_edges) {
      delete[] quad_edges;
    }

    tri_edges = new index_t[ET::TRI_EDGES * ntri_faces];
    quad_edges = new index_t[ET::QUAD_EDGES * nquad_faces];

    // Fill the face arrays with NO_LABEL
    std::fill(tri_edges, tri_edges + ET::TRI_EDGES * ntri_faces, NO_LABEL);
    std::fill(quad_edges, quad_edges + ET::QUAD_EDGES * nquad_faces, NO_LABEL);

    for (index_t elem = 0; elem < nelems; elem++) {
      // Ensure that the faces for this element are ordered consistently.
      make_element_edge_consistent(elem);

      index_t edges[ET::MAX_ELEMENT_EDGES];
      index_t ne = get_element_edges(elem, edges);

      for (index_t e0 = 0; e0 < ne; e0++) {
        if (edges[e0] == NO_LABEL) {
          // Try and find adjacent elements with this edge
          index_t v0[2];
          get_element_edge_verts(elem, e0, v0);

          // Search for all adjacent edges
          const index_t* adj_elems;
          const index_t nadj_elems =
              get_adjacent_elements_from_vert(v0[0], &adj_elems);

          for (index_t i = 0; i < nadj_elems; i++) {
            index_t adj_ne = get_num_element_edges(adj_elems[i]);

            // Check for a matching vertex
            for (index_t e1 = 0; e1 < adj_ne; e1++) {
              index_t v1[2];
              get_element_edge_verts(elem, e1, v1);

              // We have a match
              if ((v0[0] == v1[0] && v0[1] == v1[1]) ||
                  (v0[0] == v1[1] && v0[1] == v1[0])) {
              }
            }
          }
        }
      }
    }
  }

  // Input counts of the verts, tet, hex, wedge and pyramid elements
  index_t nverts;
  index_t ntets, nhex, nwedge, npyrmd;

  // Derived count information
  index_t nelems;  // Total number of elements
  index_t ntri_faces, nquad_faces;
  index_t nfaces;
  index_t nedges;

  // Element -> vert connectivity
  index_t* tet_verts;
  index_t* hex_verts;
  index_t* wedge_verts;
  index_t* pyrmd_verts;

  // Global data for pointing from verts to elements
  index_t* vert_element_ptr;
  index_t* vert_elements;

  // Element -> face connectivity
  index_t* tet_faces;
  index_t* hex_faces;
  index_t* wedge_faces;
  index_t* pyrmd_faces;

  // Face -> element connectivity
  index_t* tri_face_elements;
  index_t* quad_face_elements;

  // Face -> edge connectivity
  index_t* tri_edges;
  index_t* quad_edges;
};

// const int tri_edge_nodes[3][2] = {{2, 1}, {2, 0}, {0, 1}};

// /*
//   Compute a node to triangle or node to quad data structure
// */
// void compute_nodes_to_elem(A2D::index_t nnodes, A2D::index_t nelems,
//                            A2D::index_t num_elem_nodes,
//                            const A2D::index_t conn[], A2D::index_t**
//                            _ptr, A2D::index_t** _node_to_elems) {
//   // Set the pointer
//   A2D::index_t* ptr = new A2D::index_t[nnodes + 1];
//   std::fill(ptr, ptr + nnodes + 1, 0);

//   // Count up the references
//   const A2D::index_t conn_size = nelems * num_elem_nodes;
//   for (A2D::index_t i = 0; i < conn_size; i++) {
//     ptr[conn[i] + 1]++;
//   }

//   // Set the pointer into the array
//   for (A2D::index_t i = 0; i < nnodes; i++) {
//     ptr[i + 1] += ptr[i];
//   }

//   // Compute the node to quads
//   A2D::index_t* node_to_elems = new A2D::index_t[ptr[nnodes]];
//   const A2D::index_t* conn_ptr = conn;
//   for (A2D::index_t i = 0; i < nelems; i++) {
//     for (A2D::index_t j = 0; j < num_elem_nodes; j++) {
//       A2D::index_t node = conn_ptr[0];
//       if (node >= 0) {
//         node_to_elems[ptr[node]] = i;
//         ptr[node]++;
//         conn_ptr++;
//       }
//     }
//   }

//   // Reset the pointer array
//   for (A2D::index_t i = nnodes; i > 0; i--) {
//     ptr[i] = ptr[i - 1];
//   }
//   ptr[0] = 0;

//   // Set the output points
//   *_ptr = ptr;
//   *_node_to_elems = node_to_elems;
// }

// /*
//   Compute all of the edges within the triangular mesh
// */
// A2D::index_t compute_planar_edges(A2D::index_t nnodes, A2D::index_t
// ntris,
//                                   const A2D::index_t tris[],
//                                   A2D::index_t* tri_edge_nums,
//                                   A2D::index_t* tri_orient) {
//   // Compute the edges in the triangular mesh
//   A2D::index_t* ptr;
//   A2D::index_t* node_to_tris;
//   compute_nodes_to_elem(nnodes, ntris, 3, tris, &ptr, &node_to_tris);

//   // Set the no-label
//   const A2D::index_t no_label = std::numeric_limits<A2D::index_t>::max();

//   // Now compute the neighbors for each triangle
//   for (A2D::index_t i = 0; i < 3 * ntris; i++) {
//     tri_edge_nums[i] = no_label;
//   }

//   A2D::index_t ne = 0;
//   for (A2D::index_t i = 0; i < ntris; i++) {
//     // Search through each edge of the each triangle
//     for (A2D::index_t j = 0; j < 3; j++) {
//       if (tri_edge_nums[3 * i + j] == no_label) {
//         tri_edge_nums[3 * i + j] = ne;
//         tri_orient[3 * i + j] = 1;

//         A2D::index_t e0[2];
//         e0[0] = tris[3 * i + tri_edge_nodes[j][0]];
//         e0[1] = tris[3 * i + tri_edge_nodes[j][1]];

//         // Search for the neighboring that shares this edge
//         A2D::index_t kp = ptr[e0[0]];
//         A2D::index_t kpend = ptr[e0[0] + 1];
//         for (; kp < kpend; kp++) {
//           // Find the potential neighbor
//           A2D::index_t n = node_to_tris[kp];

//           // Don't count the same edge twice
//           if (n == i) {
//             continue;
//           }

//           // Flag to indicate that we have found the other edge (there
//           // will only be at most one other match since this is
//           // planar in parameter space)
//           A2D::index_t quit = 0;

//           // Search over all the edges on this quad, and see
//           // if they match
//           for (A2D::index_t e = 0; e < 3; e++) {
//             A2D::index_t e1[2];
//             e1[0] = tris[3 * n + tri_edge_nodes[e][0]];
//             e1[1] = tris[3 * n + tri_edge_nodes[e][1]];

//             // Check if the adjacent edge matches in either direction
//             if ((e0[0] == e1[0] && e0[1] == e1[1]) ||
//                 (e0[0] == e1[1] && e0[1] == e1[0])) {
//               // Label the other edge that shares this same node
//               tri_edge_nums[3 * n + e] = ne;

//               // Opposite orientation
//               tri_orient[3 * i + j] = 0;

//               quit = 1;
//             }
//           }
//           if (quit) {
//             break;
//           }
//         }

//         // Increment the edge number
//         ne++;
//       }
//     }
//   }

//   // Free the data that is no longer required
//   delete[] ptr;
//   delete[] node_to_tris;

//   return ne;
// }

/*
  The element connectivity class

  This class will need a substantial overhaul
  1. The connectivity is fixed at this point
  2.
*/
// class ElementConnectivity {
//  public:
//   ElementConnectivity(A2D::index_t nnodes, A2D::index_t nelems,
//                       A2D::index_t* conn_)
//       : nnodes(nnodes), nelems(nelems) {
//     conn = new A2D::index_t[3 * nelems];
//     for (A2D::index_t i = 0; i < 3 * nelems; i++) {
//       conn[i] = conn_[i];
//     }

//     edges = new A2D::index_t[3 * nelems];
//     orient = new A2D::index_t[3 * nelems];
//     nedges = compute_planar_edges(nnodes, nelems, conn, edges, orient);
//   }

//   A2D::index_t get_face_dof(A2D::index_t elem, A2D::index_t index) {
//     return elem;
//   }

//   A2D::index_t get_edge_dof(A2D::index_t elem, A2D::index_t index, int&
//   ort)
//   {
//     if (orient[3 * elem + index]) {
//       ort = 1;
//     } else {
//       ort = -1;
//     }

//     return edges[3 * elem + index];
//   }

//   A2D::index_t get_node_dof(A2D::index_t elem, A2D::index_t index) {
//     return conn[3 * elem + index];
//   }

//   A2D::index_t get_num_elements() { return nelems; }
//   A2D::index_t get_num_edges() { return nedges; }
//   A2D::index_t get_num_nodes() { return nnodes; }

//  private:
//   A2D::index_t nnodes;
//   A2D::index_t nedges;
//   A2D::index_t nelems;
//   A2D::index_t* conn;
//   A2D::index_t* edges;
//   A2D::index_t* orient;  // edge orientation 1 for +ve, 0 for -ve
// };

/*
  ElementMesh - Map from an element to the global to element local degrees
  of freedom
*/
// enum SpaceType { L2, H1, EDGE };

// template <class Basis>
// class ElementMesh {
//  public:
//   ElementMesh(A2D::ElementConnectivity& conn, SpaceType spaces_[],
//               A2D::index_t dim_[] = NULL)
//       : conn(conn) {
//     for (A2D::index_t i = 0; i < nspaces; i++) {
//       spaces[i] = spaces_[i];
//     }
//     if (dim_) {
//       for (A2D::index_t i = 0; i < nspaces; i++) {
//         dim[i] = dim_[i];
//       }
//     } else {
//       for (A2D::index_t i = 0; i < nspaces; i++) {
//         dim[i] = 1;
//       }
//     }

//     offset[0] = 0;

//     for (A2D::index_t index = 0; index < nspaces; index++) {
//       if (spaces[index] == L2) {
//         counts[index] = dim[index] * conn.get_num_elements();
//       } else if (spaces[index] == H1) {
//         counts[index] = dim[index] * conn.get_num_nodes();
//       } else if (spaces[index] == EDGE) {
//         counts[index] = dim[index] * conn.get_num_edges();
//       }

//       offset[index + 1] = offset[index] + counts[index];
//     }
//   }

//   static const A2D::index_t nspaces = Basis::nbasis;

//   A2D::index_t get_num_elements() { return conn.get_num_elements(); }
//   A2D::index_t get_num_dof() { return offset[nspaces]; }

//   int get_global_dof_sign(A2D::index_t elem, A2D::index_t space,
//                           A2D::index_t index) {
//     if (spaces[space] == EDGE) {
//       int orient;
//       conn.get_edge_dof(elem, index, orient);
//       return orient;
//     }
//     return 1;
//   }
//   A2D::index_t get_global_dof(A2D::index_t elem, A2D::index_t space,
//                               A2D::index_t index) {
//     if (spaces[space] == L2) {
//       return offset[space] +
//              dim[space] * conn.get_face_dof(elem, index / dim[space]) +
//              index % dim[space];
//     } else if (spaces[space] == H1) {
//       return offset[space] +
//              dim[space] * conn.get_node_dof(elem, index / dim[space]) +
//              index % dim[space];
//     } else if (spaces[space] == EDGE) {
//       int orient;
//       return offset[space] +
//              dim[space] * conn.get_edge_dof(elem, index / dim[space],
//              orient)
//              + index % dim[space];
//     }
//     return 0;
//   }

//  private:
//   A2D::ElementConnectivity conn;
//   A2D::index_t counts[nspaces];
//   A2D::index_t offset[nspaces + 1];
//   SpaceType spaces[nspaces];
//   A2D::index_t dim[nspaces];
// };

}  // namespace A2D

#endif  // A2D_FE_MESH_H
