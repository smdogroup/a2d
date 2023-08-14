#ifndef A2D_FE_MESH_INL_H
#define A2D_FE_MESH_INL_H

#include "multiphysics/femesh.h"

namespace A2D {

/**
 * @brief Initialize the data for the connection between vertices and
 * the elements
 *
 * Given a vertex index, this enables retrieval of all elements that
 * reference that vertex.
 */
inline void MeshConnectivityBase::init_vert_element_data() {
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

inline MeshConnectivityBase::MeshConnectivityBase(index_t nverts,
                                                  index_t nelems)
    : nverts(nverts), nelems(nelems) {
  nfaces = 0;
  nedges = 0;

  nline_faces = 0;
  ntri_faces = 0;
  nquad_faces = 0;

  // vert -> element connectivity
  vert_element_ptr = NULL;
  vert_elements = NULL;

  // Face -> element connectivity
  line_face_elements = NULL;
  tri_face_elements = NULL;
  quad_face_elements = NULL;

  // Set to NULL all the boundary info
  num_boundary_labels = 1;
  num_boundary_faces = 0;
  boundary_labels = NULL;
  boundary_faces = NULL;
}

inline MeshConnectivityBase::~MeshConnectivityBase() {
  if (vert_element_ptr) {
    DELETE_ARRAY(vert_element_ptr);
  }
  if (vert_elements) {
    DELETE_ARRAY(vert_elements);
  }

  if (line_face_elements) {
    DELETE_ARRAY(line_face_elements);
  }
  if (tri_face_elements) {
    DELETE_ARRAY(tri_face_elements);
  }
  if (quad_face_elements) {
    DELETE_ARRAY(quad_face_elements);
  }

  if (boundary_faces) {
    DELETE_ARRAY(boundary_faces);
  }
  if (boundary_labels) {
    DELETE_ARRAY(boundary_labels);
  }
}

/**
 * @brief Initialize all derived connectivity data and set up boundary
 *
 * Note: allocation must be performed explicitly before calling this function
 */
inline void MeshConnectivityBase::initialize() {
  init_vert_element_data();
  init_face_data();
  init_edge_data();

  // Count up all the faces that are on the boundary
  num_boundary_faces = 0;
  for (index_t face = 0; face < nfaces; face++) {
    index_t e1, e2;
    get_face_elements(face, &e1, &e2);
    if (e1 == NO_LABEL || e2 == NO_LABEL) {
      num_boundary_faces++;
    }
  }

  // Set the boundary labels - 0 for anything on the boundary
  boundary_faces = new index_t[num_boundary_faces];
  boundary_labels = new index_t[num_boundary_faces];
  for (index_t face = 0, count = 0; face < nfaces; face++) {
    index_t e1, e2;
    get_face_elements(face, &e1, &e2);
    if (e1 == NO_LABEL || e2 == NO_LABEL) {
      boundary_faces[count] = face;
      boundary_labels[count] = 0;
      count++;
    }
  }
}

/**
 * @brief Get the boundary faces and labels
 *
 * @param boundary_faces_ Array of the boundary faces
 * @param boundary_labels_ Array of the boundary labels
 * @return The number of boundary faces
 */
inline index_t MeshConnectivityBase::get_boundary_faces(
    const index_t* boundary_faces_[], const index_t* boundary_labels_[]) {
  if (boundary_faces_) {
    *boundary_faces_ = boundary_faces;
  }
  if (boundary_labels_) {
    *boundary_labels_ = boundary_labels;
  }
  return num_boundary_faces;
}

/**
 * @brief Count up the number of boundary faces with the given label
 *
 * @param label The label index
 * @return The number
 */
inline index_t MeshConnectivityBase::get_num_boundary_faces_with_label(
    const index_t label) {
  index_t count = 0;
  for (index_t i = 0; i < num_boundary_faces; i++) {
    if (boundary_labels[i] == label) {
      count++;
    }
  }

  return count;
}

/**
 * @brief Add a boundary label from the vertices
 *
 * Any boundary face with all nodes that touch the given set of vertices is
 * labeled
 *
 * @return index_t
 */
template <typename IdxType>
inline index_t MeshConnectivityBase::add_boundary_label_from_verts(
    index_t nv, const IdxType vert_list[]) {
  index_t* vert_labels = new index_t[nverts];
  std::fill(vert_labels, vert_labels + nverts, 0);

  index_t count = 0;

  // Apply a specific label to all the verts that touch the given set of nodes
  const index_t VERT_LABEL = 1;
  for (index_t i = 0; i < nv; i++) {
    if (vert_list[i] >= 0 && vert_list[i] < nverts) {
      vert_labels[vert_list[i]] = VERT_LABEL;
    }
  }

  for (index_t i = 0; i < num_boundary_faces; i++) {
    index_t face = boundary_faces[i];

    // Get the element corresponding to the face
    index_t elem, e2;
    get_face_elements(face, &elem, &e2);

    // Find the face index for the the boundary face;
    index_t face_index = 0;
    const index_t* faces;
    index_t nf = get_element_faces(elem, &faces);
    for (index_t j = 0; j < nf; j++) {
      if (faces[j] == face) {
        face_index = j;
        break;
      }
    }

    // Check if all the vertices are labeled or not
    index_t v[4];
    index_t nv = get_element_face_verts(elem, face_index, v);
    bool labeled = true;
    for (index_t k = 0; k < nv; k++) {
      if (vert_labels[v[k]] != VERT_LABEL) {
        labeled = false;
      }
    }

    if (labeled) {
      boundary_labels[i] = num_boundary_labels;
      count++;
    }
  }

  DELETE_ARRAY(vert_labels);

  // Return the label that was applied to the faces
  if (count > 0) {
    num_boundary_labels++;
    return num_boundary_labels - 1;
  }

  // No faces were found
  return NO_LABEL;
}

/**
 * @brief Get the element faces array
 *
 * @param elem The element index
 * @param faces The face indicies associated with the element
 * @return The number of faces
 */
inline index_t MeshConnectivityBase::get_element_faces(index_t elem,
                                                       const index_t* faces[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  *faces = nullptr;
  if (*meta.faces) {
    *faces = &(*meta.faces)[meta.NFACES * elem];
  }
  return meta.NFACES;
}

/**
 * @brief Get the element verts
 *
 * @param elem The element index
 * @param verts The vertex indices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_verts(index_t elem,
                                                       const index_t* verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  *verts = nullptr;
  if (*meta.verts) {
    *verts = &(*meta.verts)[meta.NVERTS * elem];
  }
  return meta.NVERTS;
}

/**
 * @brief Get the adjacent elements from the vert index
 *
 * @param vert The index of the vertex
 * @param elems Array of the adjacent element indices
 * @return The number of elements
 */
inline index_t MeshConnectivityBase::get_adjacent_elements_from_vert(
    const index_t vert, const index_t* elems[]) {
  if (vert < nverts) {
    *elems = &vert_elements[vert_element_ptr[vert]];
    return vert_element_ptr[vert + 1] - vert_element_ptr[vert];
  }
  *elems = NULL;
  return 0;
}

/**
 * @brief Get the global indices of the face in its local order
 *
 * @param elem The element index
 * @param f The local face index
 * @param verts The vertices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_face_verts(index_t elem,
                                                            index_t f,
                                                            index_t verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nverts = meta.FACE_NVERTS[f];
  for (index_t i = 0; i < nverts; i++) {
    verts[i] = (*meta.verts)[meta.NVERTS * elem + meta.FACE_VERTS[f][i]];
  }
  return nverts;
}

/**
 * @brief Get the local indices of the face in its local order
 *
 * @param elem The element index
 * @param f The local face index
 * @param verts The vertices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_face_vert_indices(
    index_t elem, index_t f, index_t verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nverts = meta.FACE_NVERTS[f];
  for (index_t i = 0; i < nverts; i++) {
    verts[i] = meta.FACE_VERTS[f][i];
  }
  return nverts;
}

/**
 * @brief Get the edges of the face in its local order
 *
 * @param elem The element index
 * @param f The local face index
 * @param edges The edges
 * @return The number of edges
 */
inline index_t MeshConnectivityBase::get_element_face_edges(index_t elem,
                                                            index_t f,
                                                            index_t edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nedges = meta.FACE_NEDGES[f];

  for (index_t i = 0; i < nedges; i++) {
    edges[i] = (*meta.edges)[meta.NEDGES * elem + meta.FACE_EDGES[f][i]];
  }

  return nedges;
}

/**
 * @brief Get the indices of the face in its local order
 *
 * @param elem The element index
 * @param f The local face index
 * @param edges The vertices
 * @return The number of edges
 */
inline index_t MeshConnectivityBase::get_element_face_edge_indices(
    index_t elem, index_t f, index_t edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nedges = meta.FACE_NEDGES[f];

  for (index_t i = 0; i < nedges; i++) {
    edges[i] = meta.FACE_EDGES[f][i];
  }

  return nedges;
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
inline index_t MeshConnectivityBase::get_element_global_face_verts(
    index_t elem, index_t local_face, index_t verts[]) {
  index_t t[ET::MAX_FACE_VERTS];
  index_t n = get_element_face_verts(elem, local_face, t);

  if (n == 2) {
    // Find the smallest vertex index
    if (t[0] < t[1]) {  // t[0] is the smallest
      verts[0] = t[0];
      verts[1] = t[1];
    } else {  // t[1] is the smallest
      verts[0] = t[1];
      verts[1] = t[0];
    }
  } else if (n == 3) {
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
  } else if (n == 4) {                                // n == 4
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
    } else if (t[1] < t[0] && t[1] < t[2] && t[1] < t[3]) {  // t[1] is smallest
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
    } else if (t[2] < t[0] && t[2] < t[1] && t[2] < t[3]) {  // t[2] is smallest
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
  } else {
    char msg[256];
    snprintf(msg, 256, "%d is an invalid number of vertices for a face", n);
    throw std::runtime_error(msg);
  }
  return n;
}

/**
 * @brief Check for equality between two faces
 *
 * @return True if the faces match, false otherwise
 */
inline bool MeshConnectivityBase::global_face_equality(index_t na,
                                                       const index_t a[],
                                                       index_t nb,
                                                       const index_t b[]) {
  if (na == nb) {
    if (na == 2) {
      if (a[0] == b[0] && a[1] == b[1]) {
        return true;
      }
    } else if (na == 3) {
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
 * @brief Check for equality between two edges
 *
 * @return True if the edges match, false otherwise
 */
inline bool MeshConnectivityBase::global_edge_equality(const index_t a[],
                                                       const index_t b[]) {
  if ((a[0] == b[0] && a[1] == b[1]) || (a[0] == b[1] && a[1] == b[0])) {
    return true;
  }
  return false;
}

/**
 * @brief Get the face elements associated with the given face
 *
 * @param face Global face index
 * @param e1 Returned value of the first element
 * @param e2 Returned value of the second element, NO_INDEX indicates second
 * element doesn't exist, and this face is a boundary face
 * @return Boolean if this face is on the boundary
 */
inline bool MeshConnectivityBase::get_face_elements(index_t face, index_t* e1,
                                                    index_t* e2) {
  if (face < nline_faces) {
    *e1 = line_face_elements[2 * face];
    *e2 = line_face_elements[2 * face + 1];
    return (line_face_elements[2 * face + 1] == NO_LABEL);
  } else if (face < ntri_faces) {
    face = face - nline_faces;
    *e1 = tri_face_elements[2 * face];
    *e2 = tri_face_elements[2 * face + 1];
    return (tri_face_elements[2 * face + 1] == NO_LABEL);
  } else if (face < nquad_faces) {
    face = face - nline_faces - ntri_faces;
    *e1 = quad_face_elements[2 * face];
    *e2 = quad_face_elements[2 * face + 1];
    return (quad_face_elements[2 * face + 1] == NO_LABEL);
  } else {
    char msg[256];
    std::snprintf(msg, 256, "%d is not a valid global face index", face);
    throw std::runtime_error(msg);
  }
}

/**
 * @brief Get the element edges
 *
 * @param elem The element index
 * @param edges The edge indices
 * @return The number of edges for this element
 */
inline index_t MeshConnectivityBase::get_element_edges(index_t elem,
                                                       const index_t* edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  *edges = nullptr;
  if (*meta.edges) {
    *edges = &(*meta.edges)[meta.NEDGES * elem];
  }
  return meta.NEDGES;
}

/**
 * @brief Get the vertices of the associated element edge
 *
 * This returns the edge oriented by the local direction
 *
 * @param elem The element index
 * @param edge The local element edge index
 * @param verts The global vertices, ordered in
 */
inline void MeshConnectivityBase::get_element_edge_verts(index_t elem,
                                                         index_t edge,
                                                         index_t verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);

  if (meta.is_valid_element) {
    verts[0] = (*meta.verts)[meta.NVERTS * elem + meta.EDGE_VERTS[edge][0]];
    verts[1] = (*meta.verts)[meta.NVERTS * elem + meta.EDGE_VERTS[edge][1]];
  }
}

template <typename I>
inline MeshConnectivity2D::MeshConnectivity2D(I nverts, I ntri, I* tri, I nquad,
                                              I* quad)
    : MeshConnectivityBase(nverts, ntri + nquad),
      ntri(ntri),
      nquad(nquad),
      meta_tri{true,
               ET::TRI_VERTS,
               ET::TRI_EDGES,
               ET::TRI_NFACES,
               ET::TRI_FACE_NVERTS,
               ET::TRI_FACE_VERTS,
               ET::TRI_FACE_NEDGES,
               ET::TRI_FACE_EDGES,
               ET::TRI_EDGE_VERTS,
               &tri_verts,
               &tri_edges,
               &tri_faces},
      meta_quad{true,
                ET::QUAD_VERTS,
                ET::QUAD_EDGES,
                ET::QUAD_FACES,
                ET::QUAD_FACE_NVERTS,
                ET::QUAD_FACE_VERTS,
                ET::QUAD_FACE_NEDGES,
                ET::QUAD_FACE_EDGES,
                ET::QUAD_EDGE_VERTS,
                &quad_verts,
                &quad_edges,
                &quad_faces},
      meta_none{false,
                0,
                0,
                0,
                ET::NONE_FACE_NQUANTS,
                nullptr,
                ET::NONE_FACE_NQUANTS,
                nullptr,
                nullptr,
                nullptr,
                nullptr,
                nullptr} {
  Timer timer("MeshConnectivity2D");
  // Allocate space for the element -> vert connectivity
  tri_verts = new index_t[ET::TRI_VERTS * ntri];
  quad_verts = new index_t[ET::QUAD_VERTS * nquad];

  // Set the connectivity: element -> verts
  for (index_t i = 0; i < ET::TRI_VERTS * ntri; i++) {
    tri_verts[i] = tri[i];
  }
  for (index_t i = 0; i < ET::QUAD_VERTS * nquad; i++) {
    quad_verts[i] = quad[i];
  }

  // Allocate the face data
  tri_faces = new index_t[ET::TRI_FACES * ntri];
  quad_faces = new index_t[ET::QUAD_FACES * nquad];

  // Fill the face arrays with NO_LABEL
  std::fill(tri_faces, tri_faces + ET::TRI_FACES * ntri, NO_LABEL);
  std::fill(quad_faces, quad_faces + ET::QUAD_FACES * nquad, NO_LABEL);

  // element -> edge connectivity
  tri_edges = new index_t[ET::TRI_EDGES * ntri];
  quad_edges = new index_t[ET::QUAD_EDGES * nquad];

  // Fill the face arrays with NO_LABEL
  std::fill(tri_edges, tri_edges + ET::TRI_EDGES * ntri, NO_LABEL);
  std::fill(quad_edges, quad_edges + ET::QUAD_EDGES * nquad, NO_LABEL);

  initialize();
}

inline MeshConnectivity2D::~MeshConnectivity2D() {
  if (tri_verts) {
    DELETE_ARRAY(tri_verts);
  }
  if (quad_verts) {
    DELETE_ARRAY(quad_verts);
  }

  if (tri_faces) {
    DELETE_ARRAY(tri_faces);
  }
  if (quad_faces) {
    DELETE_ARRAY(quad_faces);
  }

  if (tri_edges) {
    DELETE_ARRAY(tri_edges);
  }
  if (quad_edges) {
    DELETE_ARRAY(quad_edges);
  }
}

/**
 * @brief Shift elem to get local index within its element type (tri or quad)
 * and return the meta data for this element type
 *
 * @param elem global element index, on exit, elem is modified such that it
 * becomes the index within the element type (tri or quad)
 * @return a const reference to the meta data struct of the element type
 */
const inline ElemConnMetaData& MeshConnectivity2D::get_local_elem_and_meta(
    index_t& elem) {
  if (elem < ntri) {
    return meta_tri;
  } else if (elem < ntri + nquad) {
    elem -= ntri;
    return meta_quad;
  } else {
    return meta_none;
  }
}

template <typename I>
inline MeshConnectivity3D::MeshConnectivity3D(I nverts, I ntets, I* tets,
                                              I nhex, I* hex, I nwedge,
                                              I* wedge, I npyrmd, I* pyrmd)
    : MeshConnectivityBase(nverts, ntets + nhex + nwedge + npyrmd),
      ntets(ntets),
      nhex(nhex),
      nwedge(nwedge),
      npyrmd(npyrmd),
      meta_tet{true,
               ET::TET_VERTS,
               ET::TET_EDGES,
               ET::TET_FACES,
               ET::TET_FACE_NVERTS,
               ET::TET_FACE_VERTS,
               ET::TET_FACE_NEDGES,
               ET::TET_FACE_EDGES,
               ET::TET_EDGE_VERTS,
               &tet_verts,
               &tet_edges,
               &tet_faces},
      meta_hex{true,
               ET::HEX_VERTS,
               ET::HEX_EDGES,
               ET::HEX_FACES,
               ET::HEX_FACE_NVERTS,
               ET::HEX_FACE_VERTS,
               ET::HEX_FACE_NEDGES,
               ET::HEX_FACE_EDGES,
               ET::HEX_EDGE_VERTS,
               &hex_verts,
               &hex_edges,
               &hex_faces},
      meta_wedge{true,
                 ET::WEDGE_VERTS,
                 ET::WEDGE_EDGES,
                 ET::WEDGE_FACES,
                 ET::WEDGE_FACE_NVERTS,
                 ET::WEDGE_FACE_VERTS,
                 ET::WEDGE_FACE_NEDGES,
                 ET::WEDGE_FACE_EDGES,
                 ET::WEDGE_EDGE_VERTS,
                 &wedge_verts,
                 &wedge_edges,
                 &wedge_faces},
      meta_pyrmd{true,
                 ET::PYRMD_VERTS,
                 ET::PYRMD_EDGES,
                 ET::PYRMD_FACES,
                 ET::PYRMD_FACE_NVERTS,
                 ET::PYRMD_FACE_VERTS,
                 ET::PYRMD_FACE_NEDGES,
                 ET::PYRMD_FACE_EDGES,
                 ET::PYRMD_EDGE_VERTS,
                 &pyrmd_verts,
                 &pyrmd_edges,
                 &pyrmd_faces},
      meta_none{false,
                0,
                0,
                0,
                ET::NONE_FACE_NQUANTS,
                nullptr,
                ET::NONE_FACE_NQUANTS,
                nullptr,
                nullptr,
                nullptr,
                nullptr,
                nullptr} {
  Timer timer("MeshConnectivity3D");
  // Allocate space for the element -> vert connectivity
  tet_verts = new index_t[ET::TET_VERTS * ntets];
  hex_verts = new index_t[ET::HEX_VERTS * nhex];
  wedge_verts = new index_t[ET::WEDGE_VERTS * nwedge];
  pyrmd_verts = new index_t[ET::PYRMD_VERTS * npyrmd];

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

  // element -> edge connectivity
  tet_edges = new index_t[ET::TET_EDGES * ntets];
  hex_edges = new index_t[ET::HEX_EDGES * nhex];
  wedge_edges = new index_t[ET::WEDGE_EDGES * nwedge];
  pyrmd_edges = new index_t[ET::PYRMD_EDGES * npyrmd];

  // Fill the face arrays with NO_LABEL
  std::fill(tet_edges, tet_edges + ET::TET_EDGES * ntets, NO_LABEL);
  std::fill(hex_edges, hex_edges + ET::HEX_EDGES * nhex, NO_LABEL);
  std::fill(wedge_edges, wedge_edges + ET::WEDGE_EDGES * nwedge, NO_LABEL);
  std::fill(pyrmd_edges, pyrmd_edges + ET::PYRMD_EDGES * npyrmd, NO_LABEL);

  initialize();
}

inline MeshConnectivity3D::~MeshConnectivity3D() {
  if (tet_verts) {
    DELETE_ARRAY(tet_verts);
  }
  if (hex_verts) {
    DELETE_ARRAY(hex_verts);
  }
  if (wedge_verts) {
    DELETE_ARRAY(wedge_verts);
  }
  if (pyrmd_verts) {
    DELETE_ARRAY(pyrmd_verts);
  }

  if (tet_faces) {
    DELETE_ARRAY(tet_faces);
  }
  if (hex_faces) {
    DELETE_ARRAY(hex_faces);
  }
  if (wedge_faces) {
    DELETE_ARRAY(wedge_faces);
  }
  if (pyrmd_faces) {
    DELETE_ARRAY(pyrmd_faces);
  }

  if (tet_edges) {
    DELETE_ARRAY(tet_edges);
  }
  if (hex_edges) {
    DELETE_ARRAY(hex_edges);
  }
  if (wedge_edges) {
    DELETE_ARRAY(wedge_edges);
  }
  if (pyrmd_edges) {
    DELETE_ARRAY(pyrmd_edges);
  }
}

/**
 * @brief Shift elem to get local index within its element type (tet, hex, etc.)
 * and return the meta data for this element type
 *
 * @param elem global element index, on exit, elem is modified such that it
 * becomes the index within the element type (tet, hex, etc.)
 * @return a const reference to the meta data struct of the element type
 */
const inline ElemConnMetaData& MeshConnectivity3D::get_local_elem_and_meta(
    index_t& elem) {
  if (elem < ntets) {
    return meta_tet;
  } else if (elem < ntets + nhex) {
    elem -= ntets;
    return meta_hex;
  } else if (elem < ntets + nhex + nwedge) {
    elem -= ntets + nhex;
    return meta_wedge;
  } else if (elem < ntets + nhex + nwedge + npyrmd) {
    elem -= ntets + nhex + nwedge;
    return meta_pyrmd;
  } else {
    return meta_none;
  }
}

/**
 * @brief Label the verts, edges and faces that touch the list of vertices
 *
 * @param nv The number of input vertices
 * @param verts The vertex numbers
 * @param vert_labels An array of length nverts
 * @param edge_labels An array of length nedges
 * @param face_labels An array of length nfaces
 */
template <typename IdxType>
inline void MeshConnectivityBase::get_labels_from_verts(const index_t nv,
                                                        const IdxType verts[],
                                                        index_t vert_labels[],
                                                        index_t edge_labels[],
                                                        index_t face_labels[]) {
  std::fill(vert_labels, vert_labels + nverts, NO_LABEL);
  std::fill(edge_labels, edge_labels + nedges, NO_LABEL);
  std::fill(face_labels, face_labels + nfaces, NO_LABEL);

  // Label the vertices
  for (index_t i = 0; i < nv; i++) {
    if (verts[i] < nverts) {
      vert_labels[verts[i]] = i;
    }
  }

  // Loop over elements and label the edges and faces
  for (index_t elem = 0; elem < nelems; elem++) {
    // Check whether the edges are labeled
    const index_t* elem_edges;
    index_t ne = get_element_edges(elem, &elem_edges);
    for (index_t e = 0; e < ne; e++) {
      index_t v[2];
      get_element_edge_verts(elem, e, v);

      // Both verts of this element edge are labeled - label the edge
      if (vert_labels[v[0]] != NO_LABEL && vert_labels[v[1]] != NO_LABEL) {
        edge_labels[elem_edges[e]] = elem_edges[e];
      }
    }

    const index_t* elem_faces;
    index_t nf = get_element_faces(elem, &elem_faces);
    for (index_t f = 0; f < nf; f++) {
      index_t fv[ET::MAX_FACE_VERTS];  // Face vertices
      index_t nfv = get_element_face_verts(elem, f, fv);

      // Check if all the face vertices are labeled
      if (nfv == 2) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL) {
          face_labels[elem_faces[f]] = elem_faces[f];
        }
      } else if (nfv == 3) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL &&
            vert_labels[fv[2]] != NO_LABEL) {
          face_labels[elem_faces[f]] = elem_faces[f];
        }
      } else if (nfv == 4) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL &&
            vert_labels[fv[2]] != NO_LABEL && vert_labels[fv[3]] != NO_LABEL) {
          face_labels[elem_faces[f]] = elem_faces[f];
        }
      }
    }
  }
}

/**
 * @brief Get a non-constant element faces array
 *
 * @param elem The element index
 * @param faces The face indicies associated with the element
 * @return The number of faces
 */
inline index_t MeshConnectivityBase::get_element_faces(index_t elem,
                                                       index_t* faces[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  *faces = nullptr;
  if (*meta.faces) {
    *faces = &(*meta.faces)[meta.NFACES * elem];
  }
  return meta.NFACES;
}

/**
 * @brief Get the element edges
 *
 * @param elem The element index
 * @param edges The edge indices
 * @return The number of edges for this element
 */
inline index_t MeshConnectivityBase::get_element_edges(index_t elem,
                                                       index_t* edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  *edges = nullptr;
  if (*meta.edges) {
    *edges = &(*meta.edges)[meta.NEDGES * elem];
  }
  return meta.NEDGES;
}

/**
 * @brief Initialize data associated with the face information
 *
 * This code uniquely orders the faces assocaited with each element
 */
inline void MeshConnectivityBase::init_face_data() {
  // Prepare to count and number the number of faces. This keeps track of
  // counts of separate face types (line, triangle, quadrilateral)
  for (index_t elem = 0; elem < nelems; elem++) {
    // Loop over the elements of this face
    index_t* faces;
    const index_t nf = get_element_faces(elem, &faces);

    for (index_t face = 0; face < nf; face++) {
      if (faces[face] == NO_LABEL) {
        // Get the unique set of verts corresponding to this face
        index_t face_verts[ET::MAX_FACE_VERTS];
        index_t nface_verts =
            get_element_global_face_verts(elem, face, face_verts);

        bool face_located = false;

        // For all adjacent elements check if there is a match
        for (index_t m = 0; m < nface_verts; m++) {
          const index_t* adj_elems;
          const index_t nadj_elems =
              get_adjacent_elements_from_vert(face_verts[m], &adj_elems);

          // Loop over all adjacent elements that are not the current
          // element
          for (index_t k = 0; k < nadj_elems; k++) {
            if (adj_elems[k] != elem) {
              // Loop over the elements of this face
              index_t* adj_faces;
              const index_t adj_nf =
                  get_element_faces(adj_elems[k], &adj_faces);

              for (index_t adj_face = 0; adj_face < adj_nf; adj_face++) {
                // Get the unique set of vertices corresponding to this
                // face
                index_t adj_face_verts[4];
                index_t adj_nface_verts = get_element_global_face_verts(
                    adj_elems[k], adj_face, adj_face_verts);

                if (adj_faces[adj_face] == NO_LABEL &&
                    global_face_equality(nface_verts, face_verts,
                                         adj_nface_verts, adj_face_verts)) {
                  if (nface_verts == 2) {
                    adj_faces[adj_face] = nline_faces;
                    faces[face] = nline_faces;
                    nline_faces++;
                  } else if (nface_verts == 3) {
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
          if (nface_verts == 2) {
            faces[face] = nline_faces;
            nline_faces++;
          } else if (nface_verts == 3) {
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
  nfaces = nline_faces + ntri_faces + nquad_faces;

  // Set the face indices associated with the two adjacent elements. At
  // this point the face indices stored are with respect to separate lists
  // of triangular and quadrilateral faces. We order the triangular faces
  // first, so we have to add the total number of triangular faces to each
  // quadrilateral face index to get the global index.
  line_face_elements = new index_t[2 * nline_faces];
  std::fill(line_face_elements, line_face_elements + 2 * nline_faces, NO_LABEL);

  tri_face_elements = new index_t[2 * ntri_faces];
  std::fill(tri_face_elements, tri_face_elements + 2 * ntri_faces, NO_LABEL);

  quad_face_elements = new index_t[2 * nquad_faces];
  std::fill(quad_face_elements, quad_face_elements + 2 * nquad_faces, NO_LABEL);

  for (index_t elem = 0; elem < nelems; elem++) {
    // Loop over the elements of this face
    index_t* faces;
    const index_t nf = get_element_faces(elem, &faces);

    for (index_t face = 0; face < nf; face++) {
      index_t face_verts[ET::MAX_FACE_VERTS];
      index_t nface_verts =
          get_element_global_face_verts(elem, face, face_verts);

      if (nface_verts == 2) {
        if (line_face_elements[2 * faces[face]] == NO_LABEL) {
          line_face_elements[2 * faces[face]] = elem;
        } else if (line_face_elements[2 * faces[face] + 1] == NO_LABEL) {
          line_face_elements[2 * faces[face] + 1] = elem;
        }
      } else if (nface_verts == 3) {
        if (tri_face_elements[2 * faces[face]] == NO_LABEL) {
          tri_face_elements[2 * faces[face]] = elem;
        } else if (tri_face_elements[2 * faces[face] + 1] == NO_LABEL) {
          tri_face_elements[2 * faces[face] + 1] = elem;
        }

        // Reset the face index into the global face index - add the
        // number of line faces. Now the faces will index from a
        // global face number.
        faces[face] += nline_faces;
      } else {
        if (quad_face_elements[2 * faces[face]] == NO_LABEL) {
          quad_face_elements[2 * faces[face]] = elem;
        } else if (quad_face_elements[2 * faces[face] + 1] == NO_LABEL) {
          quad_face_elements[2 * faces[face] + 1] = elem;
        }

        // Reset the face index into the global face index - add the
        // number of line and triangle faces. Now the faces will index from a
        // global face number.
        faces[face] += nline_faces + ntri_faces;
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
inline void MeshConnectivityBase::init_edge_data() {
  for (index_t elem = 0; elem < nelems; elem++) {
    // Get the number of element edges
    index_t* elem_edges;
    index_t ne = get_element_edges(elem, &elem_edges);

    for (index_t e0 = 0; e0 < ne; e0++) {
      if (elem_edges[e0] == NO_LABEL) {
        index_t edge_num = nedges;
        nedges++;

        elem_edges[e0] = edge_num;

        // Find adjacent elements with this edge
        index_t v0[2];
        get_element_edge_verts(elem, e0, v0);

        // Loop over all adjacent elements to find touching edges
        const index_t* adj_elems;
        const index_t nadj_elems =
            get_adjacent_elements_from_vert(v0[0], &adj_elems);

        // Loop back over the elements and edges and set the global
        //   numbering
        for (index_t i = 0; i < nadj_elems; i++) {
          if (adj_elems[i] != elem) {
            index_t* adj_elem_edges;
            index_t adj_ne = get_element_edges(adj_elems[i], &adj_elem_edges);

            // Check for a matching vertex
            for (index_t e1 = 0; e1 < adj_ne; e1++) {
              index_t v1[2];
              get_element_edge_verts(adj_elems[i], e1, v1);

              // Check if we have a match
              if (global_edge_equality(v0, v1)) {
                adj_elem_edges[e1] = edge_num;
              }
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Add a boundary condition
 *
 * @param label The label for the boundary condition in the connectivity
 * @param basis Index of the basis functions in the FEBasis object
 */
inline void DirichletBCInfo::add_boundary_condition(index_t label,
                                                    index_t basis) {
  index_t info = 0;
  for (index_t i = 0; i < MAX_FIXED_FIELDS; i++) {
    info |= 1 << i;
  }
  data.push_back(std::make_tuple(label, basis, info));
}

/**
 * @brief Add a boundary condition that fixes only the specified subset of
 * dof
 *
 * @param label The label for the boundary condition in the connectivity
 * @param basis Index of the basis functions in the FEBasis object
 */
inline void DirichletBCInfo::add_boundary_condition(index_t label,
                                                    index_t basis,
                                                    index_t nfixed,
                                                    const index_t fixed[]) {
  index_t info = 0;
  for (index_t i = 0; i < nfixed; i++) {
    info |= 1 << fixed[i];
  }
  data.push_back(std::make_tuple(label, basis, info));
}

/**
 * @brief Search if a boundary condition is active for the specified label
 *
 * @param label The label index
 * @return True if the label is used, false otherwise
 */
inline bool DirichletBCInfo::active_for_label(index_t label) {
  for (index_t i = 0; i < data.size(); i++) {
    if (label == std::get<0>(data[i])) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Search if a boundary condition is active for the specified label
 * and basis pair
 *
 * @param label The label index
 * @param basis The basis index
 * @return True if the label is used, false otherwise
 */
inline bool DirichletBCInfo::active_for_basis(index_t label, index_t basis) {
  for (index_t i = 0; i < data.size(); i++) {
    if (label == std::get<0>(data[i]) && basis == std::get<1>(data[i])) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Search if a boundary condition is active for the specified label and
 * basis pair
 *
 * @param label The label index
 * @param basis The basis index
 * @param dof The dof index
 * @return True if the label is used, false otherwise
 */
inline bool DirichletBCInfo::active_for_dof(index_t label, index_t basis,
                                            index_t dof) {
  for (index_t i = 0; i < data.size(); i++) {
    if (label == std::get<0>(data[i]) && basis == std::get<1>(data[i]) &&
        ((1 << dof) & std::get<2>(data[i]))) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Construct a new Element Mesh object
 *
 * Order the degrees of freedom associated with the finite-element basis
 * across the entire mesh
 *
 * @param conn The mesh connectivity
 */
template <class Basis>
ElementMesh<Basis>::ElementMesh(MeshConnectivityBase& conn)
    : nelems(conn.get_num_elements()), num_dof(0) {
  // Count up the number of degrees of freedom
  element_dof = new index_t[nelems * ndof_per_element];
  element_sign = new int[nelems * ndof_per_element];
  std::fill(element_dof, element_dof + nelems * ndof_per_element, NO_INDEX);

  // Perform a sweep of the elements
  std::vector<index_t> ids(nelems, NO_INDEX), stack(nelems);

  index_t start = 0, end = 1, level = 0;
  stack[0] = 0;
  ids[0] = level;

  while (start < end) {
    index_t next = end;
    level++;

    for (index_t i = start; i < end; i++) {
      // Loop over
      const index_t* faces;
      index_t nf = conn.get_element_faces(stack[i], &faces);

      for (index_t j = 0; j < nf; j++) {
        // Look at the adjacent elements
        index_t e1, e2;
        conn.get_face_elements(faces[j], &e1, &e2);

        if (e1 < nelems && ids[e1] == NO_INDEX) {
          ids[e1] = level;
          stack[next] = e1;
          next++;
        }
        if (e2 < nelems && ids[e2] == NO_INDEX) {
          ids[e2] = level;
          stack[next] = e2;
          next++;
        }
      }
    }

    start = end;
    end = next;
  }

  std::vector<index_t> face_owners(conn.get_num_faces(), NO_INDEX);
  std::vector<index_t> edge_owners(conn.get_num_edges(), NO_INDEX);
  std::vector<index_t> vert_owners(conn.get_num_verts(), NO_INDEX);

  index_t dof_counter = 0;
  index_t dof[ndof_per_element];

  for (index_t basis = 0; basis < Basis::nbasis; basis++) {
    for (index_t counter = nelems; counter > 0; counter--) {
      index_t elem = stack[counter - 1];

      index_t* elem_dof = &element_dof[elem * ndof_per_element];
      int* elem_sign = &element_sign[elem * ndof_per_element];

      // The volume DOF are always owned by the element - no need to check
      // for the element that owns them
      index_t ndof;
      if constexpr (dim == 3) {
        ndof = Basis::get_entity_ndof(basis, ET::VOLUME, 0);
      } else if constexpr (dim == 2) {
        ndof = Basis::get_entity_ndof(basis, ET::FACE, 0);
      } else if constexpr (dim == 1) {
        ndof = Basis::get_entity_ndof(basis, ET::EDGE, 0);
      }

      for (index_t i = 0; i < ndof; i++, dof_counter++) {
        dof[i] = dof_counter;
      }
      if constexpr (dim == 3) {
        Basis::set_entity_dof(basis, ET::VOLUME, 0, 0, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::VOLUME, 0, 0, elem_sign);
      } else if constexpr (dim == 2) {
        Basis::set_entity_dof(basis, ET::FACE, 0, 0, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::FACE, 0, 0, elem_sign);
      } else if constexpr (dim == 1) {
        Basis::set_entity_dof(basis, ET::EDGE, 0, 0, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::EDGE, 0, 0, elem_sign);
      }

      // Order the faces
      const index_t* faces;
      index_t nf = conn.get_element_faces(elem, &faces);
      for (index_t index = 0; index < nf; index++) {
        index_t face = faces[index];
        index_t orient = 0;
        if (face_owners[face] == NO_INDEX || face_owners[face] == elem) {
          face_owners[face] = elem;

          if constexpr (dim == 3) {
            ndof = Basis::get_entity_ndof(basis, ET::FACE, index);
          } else if constexpr (dim == 2) {
            ndof = Basis::get_entity_ndof(basis, ET::EDGE, index);
          } else {
            ndof = 0;
          }
          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = face_owners[face];
          const index_t* owner_faces;
          index_t nf_owner = conn.get_element_faces(owner_elem, &owner_faces);

          for (index_t i = 0; i < nf_owner; i++) {
            if (owner_faces[i] == face) {
              index_t ref[4], verts[4];
              index_t nverts = conn.get_element_face_verts(owner_elem, i, ref);
              conn.get_element_face_verts(elem, index, verts);

              if constexpr (dim == 3) {
                Basis::get_entity_dof(
                    basis, ET::FACE, i,
                    &element_dof[owner_elem * ndof_per_element], dof);
              } else if constexpr (dim == 2) {
                Basis::get_entity_dof(
                    basis, ET::EDGE, i,
                    &element_dof[owner_elem * ndof_per_element], dof);
              }

              if (nverts == 4) {
                orient = ET::get_quad_face_orientation(ref, verts);
              }
              break;
            }
          }
        }

        if constexpr (dim == 3) {
          Basis::set_entity_dof(basis, ET::FACE, index, orient, dof, elem_dof);
          Basis::set_entity_signs(basis, ET::FACE, index, orient, elem_sign);
        } else if constexpr (dim == 2) {
          Basis::set_entity_dof(basis, ET::EDGE, index, orient, dof, elem_dof);
          Basis::set_entity_signs(basis, ET::EDGE, index, orient, elem_sign);
        }
      }

      // Order the edges
      const index_t* edges;
      index_t ne = conn.get_element_edges(elem, &edges);
      for (index_t index = 0; index < ne; index++) {
        index_t edge = edges[index];
        index_t orient = 0;
        if (edge_owners[edge] == NO_INDEX || edge_owners[edge] == elem) {
          edge_owners[edge] = elem;

          if constexpr (dim == 3) {
            ndof = Basis::get_entity_ndof(basis, ET::EDGE, index);
          } else {
            ndof = 0;
          }
          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = edge_owners[edge];
          const index_t* owner_edges;
          index_t ne_owner = conn.get_element_edges(owner_elem, &owner_edges);

          for (index_t i = 0; i < ne_owner; i++) {
            if (owner_edges[i] == edge) {
              index_t ref[2], verts[2];
              conn.get_element_edge_verts(owner_elem, i, ref);
              conn.get_element_edge_verts(elem, index, verts);

              if constexpr (dim == 3) {
                Basis::get_entity_dof(
                    basis, ET::EDGE, i,
                    &element_dof[owner_elem * ndof_per_element], dof);
              }

              if (ref[0] == verts[1] && ref[1] == verts[0]) {
                orient = 1;
              }
              break;
            }
          }
        }

        if constexpr (dim == 3) {
          Basis::set_entity_dof(basis, ET::EDGE, index, orient, dof, elem_dof);
          Basis::set_entity_signs(basis, ET::EDGE, index, orient, elem_sign);
        }
      }

      // Order the vertices
      const index_t* verts;
      index_t nv = conn.get_element_verts(elem, &verts);
      for (index_t index = 0; index < nv; index++) {
        index_t vert = verts[index];
        index_t orient = 0;
        if (vert_owners[vert] == NO_INDEX || vert_owners[vert] == elem) {
          vert_owners[vert] = elem;

          ndof = Basis::get_entity_ndof(basis, ET::VERTEX, index);
          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = vert_owners[vert];
          const index_t* owner_verts;
          index_t nv_owner = conn.get_element_verts(owner_elem, &owner_verts);

          for (index_t i = 0; i < nv_owner; i++) {
            if (owner_verts[i] == vert) {
              Basis::get_entity_dof(basis, ET::VERTEX, i,
                                    &element_dof[owner_elem * ndof_per_element],
                                    dof);
              break;
            }
          }
        }

        Basis::set_entity_dof(basis, ET::VERTEX, index, orient, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::VERTEX, index, orient, elem_sign);
      }
    }

    num_dof_offset[basis] = dof_counter;
  }

  // Set the number of degrees of freedom
  num_dof = dof_counter;
}

/**
 * @brief Construct a new ElementMesh object
 *
 * Take the degrees of freedom for only thoes surface elements with the
 * specified label
 *
 * @tparam InteriorBasis FEBasis type used in the interior
 * @param label Surface label index
 * @param conn The mesh connectivity information associated with the mesh
 * @param mesh Mesh object generated with the InteriorBasis
 */
template <class Basis>
template <class InteriorBasis>
ElementMesh<Basis>::ElementMesh(const index_t label, MeshConnectivityBase& conn,
                                ElementMesh<InteriorBasis>& mesh)
    : nelems(conn.get_num_boundary_faces_with_label(label)) {
  element_dof = new index_t[nelems * ndof_per_element];
  element_sign = new int[nelems * ndof_per_element];

  // Get the number of boundary faces
  const index_t* boundary_faces;
  const index_t* boundary_labels;
  index_t num_boundary_faces =
      conn.get_boundary_faces(&boundary_faces, &boundary_labels);

  for (index_t i = 0, elem_count = 0; i < num_boundary_faces; i++) {
    if (boundary_labels[i] == label) {
      index_t face = boundary_faces[i];

      // Get the element adjacent to the boundary face (e2 can be ignored)
      index_t elem, e2;
      conn.get_face_elements(face, &elem, &e2);

      // Find the face index for the element
      index_t face_index = 0;
      const index_t* faces;
      index_t nf = conn.get_element_faces(elem, &faces);
      for (index_t index = 0; index < nf; index++) {
        if (face == faces[index]) {
          face_index = index;
        }
      }

      // Get the signs associated wtih the original mesh
      const index_t* dof;
      const int* signs;
      mesh.get_element_dof(elem, &dof);
      mesh.get_element_signs(elem, &signs);

      // Set pointers for the entity dof
      index_t* surf_dof = &element_dof[elem_count * ndof_per_element];
      int* surf_signs = &element_sign[elem_count * ndof_per_element];

      // The degree of freedom indices for each entity
      index_t entity_dof[InteriorBasis::ndof];

      // Set the vertices associated with the face
      index_t v[4];
      index_t nv = conn.get_element_face_vert_indices(elem, face_index, v);
      for (index_t j = 0; j < nv; j++) {
        for (index_t basis = 0; basis < Basis::nbasis; basis++) {
          index_t vert_index = v[j];

          // Get the vertex dof from the interior element
          InteriorBasis::get_entity_dof(basis, ET::VERTEX, vert_index, dof,
                                        entity_dof);

          index_t surf_vert_index = j;
          index_t orient = 0;

          // Set the same vertex index on the corresponding surface
          Basis::set_entity_dof(basis, ET::VERTEX, surf_vert_index, orient,
                                entity_dof, surf_dof);
          Basis::set_entity_signs(basis, ET::VERTEX, surf_vert_index, orient,
                                  surf_signs);
        }
      }

      // Set the edges associated with the face
      index_t e[4];
      index_t ne = conn.get_element_face_edge_indices(elem, face_index, e);
      for (index_t j = 0; j < ne; j++) {
        for (index_t basis = 0; basis < Basis::nbasis; basis++) {
          index_t edge_index = e[j];

          // Get the vertex dof from the interior element
          if constexpr (dim == 3) {
            InteriorBasis::get_entity_dof(basis, ET::EDGE, edge_index, dof,
                                          entity_dof);
          }

          // Get the edge orientation relative to the face
          index_t orient = 0;

          // Get the index of the edge on the face
          index_t surf_edge_index = j;

          // Set the same edge on the corresponding surface
          if constexpr (dim == 3) {
            Basis::set_entity_dof(basis, ET::EDGE, surf_edge_index, orient,
                                  entity_dof, surf_dof);
            Basis::set_entity_signs(basis, ET::EDGE, surf_edge_index, orient,
                                    surf_signs);
          }
        }
      }

      // Loop over all the boundary faces
      for (index_t basis = 0; basis < InteriorBasis::nbasis; basis++) {
        // Get the degrees of freeom from the element
        index_t orient = 0;
        if constexpr (dim == 3) {
          InteriorBasis::get_entity_dof(basis, ET::FACE, face_index, dof,
                                        entity_dof);
        } else if constexpr (dim == 2) {
          InteriorBasis::get_entity_dof(basis, ET::EDGE, face_index, dof,
                                        entity_dof);
        }

        // Set the degrees of freedom - the face element has the same
        // orientation as its interior owner
        if constexpr (dim == 3) {
          Basis::set_entity_dof(basis, ET::FACE, 0, orient, entity_dof,
                                surf_dof);
          Basis::set_entity_signs(basis, ET::FACE, 0, orient, surf_signs);
        } else if constexpr (dim == 2) {
          Basis::set_entity_dof(basis, ET::EDGE, 0, orient, entity_dof,
                                surf_dof);
          Basis::set_entity_signs(basis, ET::EDGE, 0, orient, surf_signs);
        }
      }

      elem_count++;
    }
  }
}

/**
 * @brief Construct a new Element Mesh object
 *
 * This element mesh is automatically created from a high-order element
 * discretization
 *
 * @tparam HOrderBasis The high-order finite-element basis
 * @param mesh The mesh class created for the high-order basis
 */
template <class Basis>
template <class HOrderBasis>
ElementMesh<Basis>::ElementMesh(ElementMesh<HOrderBasis>& mesh)
    : nelems(HOrderBasis::get_num_lorder_elements() * mesh.get_num_elements()),
      num_dof(mesh.get_num_dof()) {
  element_dof = new index_t[nelems * ndof_per_element];
  element_sign = new int[nelems * ndof_per_element];

  for (index_t i = 0; i < Basis::nbasis; i++) {
    num_dof_offset[i] = mesh.get_num_cumulative_dof(i);
  }

  // Loop over all the high-order elements and generate the low-order
  // representation
  for (index_t j = 0; j < mesh.get_num_elements(); j++) {
    // Get the signs associated wtih the high-order mesh
    const index_t* horder_dof;
    const int* horder_signs;
    mesh.get_element_dof(j, &horder_dof);
    mesh.get_element_signs(j, &horder_signs);

    for (index_t i = 0; i < HOrderBasis::get_num_lorder_elements(); i++) {
      index_t elem = i + j * HOrderBasis::get_num_lorder_elements();

      // Index into the high-order element
      index_t* lorder_dof = &element_dof[ndof_per_element * elem];
      HOrderBasis::get_lorder_dof(i, horder_dof, lorder_dof);

      // Signs indicating any orientation flip
      int* lorder_signs = &element_sign[ndof_per_element * elem];
      HOrderBasis::get_lorder_signs(i, horder_signs, lorder_signs);
    }
  }
}

template <class Basis>
template <index_t basis>
int ElementMesh<Basis>::get_global_dof_sign(index_t elem, index_t index) {
  return element_sign[ndof_per_element * elem +
                      Basis::template get_dof_offset<basis>() + index];
}

template <class Basis>
template <index_t basis>
index_t ElementMesh<Basis>::get_global_dof(index_t elem, index_t index) {
  return element_dof[ndof_per_element * elem +
                     Basis::template get_dof_offset<basis>() + index];
}

template <class Basis>
template <index_t M, index_t basis_offset>
void ElementMesh<Basis>::create_block_csr(index_t& nrows,
                                          std::vector<index_t>& rowp,
                                          std::vector<index_t>& cols) {
  std::set<std::pair<index_t, index_t>> node_set;

  const constexpr index_t ndof_per_elem =
      Basis::template get_dof_offset<basis_offset>();

  for (index_t i = 0; i < nelems; i++) {
    index_t dof_reduced[Basis::ndof];
    index_t n = 0;
    for (index_t j1 = 0; j1 < ndof_per_elem; j1++, n++) {
      index_t row = element_dof[i * Basis::ndof + j1] / M;
      while (j1 + 1 < ndof_per_elem &&
             row == (element_dof[i * Basis::ndof + j1 + 1] / M)) {
        j1++;
      }
      dof_reduced[n] = row;
    }

    for (index_t j1 = 0; j1 < n; j1++) {
      for (index_t j2 = 0; j2 < n; j2++) {
        node_set.insert(std::make_pair(dof_reduced[j1], dof_reduced[j2]));
      }
    }
  }

  nrows = num_dof_offset[basis_offset - 1] / M;
  if (num_dof_offset[basis_offset - 1] % M > 0) {
    nrows += 1;
  }

  // Find the number of nodes referenced by other nodes
  rowp.resize(nrows + 1);
  std::fill(rowp.begin(), rowp.end(), 0);

  typename std::set<std::pair<index_t, index_t>>::iterator it;
  for (it = node_set.begin(); it != node_set.end(); it++) {
    rowp[it->first + 1] += 1;
  }

  // Set the pointer into the rows
  rowp[0] = 0;
  for (index_t i = 0; i < nrows; i++) {
    rowp[i + 1] += rowp[i];
  }

  index_t nnz = rowp[nrows];
  cols.resize(nnz);

  for (it = node_set.begin(); it != node_set.end(); it++) {
    cols[rowp[it->first]] = it->second;
    rowp[it->first]++;
  }

  // Reset the pointer into the nodes
  for (index_t i = nrows; i > 0; i--) {
    rowp[i] = rowp[i - 1];
  }
  rowp[0] = 0;

  // Sort the cols array
  SortCSRData(nrows, rowp, cols);
}

/**
 * @brief Construct a new Boundary Condition object by constraining the
 * degrees of freedom from the vertices, edges and faces that touch the
 * specified degrees of freedom
 *
 * @param conn The connectivity
 * @param mesh The mesh object with the ordered degrees of freedom
 */
template <class Basis>
DirichletBCs<Basis>::DirichletBCs(MeshConnectivityBase& conn,
                                  ElementMesh<Basis>& mesh,
                                  DirichletBCInfo& bcinfo) {
  const index_t* boundary_faces;
  const index_t* boundary_labels;
  index_t num_boundary_faces =
      conn.get_boundary_faces(&boundary_faces, &boundary_labels);

  index_t num_dof = mesh.get_num_dof();
  index_t boundary_dof_counter = 0;
  index_t* dof_flags = new index_t[num_dof];
  index_t* boundary_dof = new index_t[num_dof];

  const index_t NO_LABEL = 1;
  std::fill(dof_flags, dof_flags + num_dof, NO_LABEL);

  // Loop over all the boundaries
  for (index_t boundary_face = 0; boundary_face < num_boundary_faces;
       boundary_face++) {
    index_t face = boundary_faces[boundary_face];
    index_t label = boundary_labels[boundary_face];

    // Find the edges and
    if (bcinfo.active_for_label(label)) {
      // Find the element for this face
      index_t elem, e2;
      conn.get_face_elements(face, &elem, &e2);

      // Get the degrees of freedom for this element
      const index_t* elem_dof;
      mesh.get_element_dof(elem, &elem_dof);

      // Array of dof that may be extracted
      index_t dof[Basis::ndof];

      // Find the face index for the active face
      index_t face_index = 0;
      const index_t* faces;
      index_t nf = conn.get_element_faces(elem, &faces);
      for (index_t i = 0; i < nf; i++) {
        if (faces[i] == face) {
          face_index = i;
          break;
        }
      }

      for (index_t basis = 0; basis < Basis::nbasis; basis++) {
        if (bcinfo.active_for_basis(label, basis)) {
          index_t nface_dof;
          if constexpr (dim == 3) {
            nface_dof = Basis::get_entity_ndof(basis, ET::FACE, face_index);
            Basis::get_entity_dof(basis, ET::FACE, face_index, elem_dof, dof);
          } else if constexpr (dim == 2) {
            nface_dof = Basis::get_entity_ndof(basis, ET::EDGE, face_index);
            Basis::get_entity_dof(basis, ET::EDGE, face_index, elem_dof, dof);
          } else {
            nface_dof = 0;
          }

          for (index_t k = 0; k < nface_dof; k++) {
            bool active = bcinfo.active_for_dof(label, basis,
                                                k % Basis::get_stride(basis));

            if (active && dof_flags[dof[k]] == NO_LABEL) {
              dof_flags[dof[k]] = 0;

              boundary_dof[boundary_dof_counter] = dof[k];
              boundary_dof_counter++;
            }
          }

          // Set the boundary conditions on the element face edges and verts
          index_t e[4];
          int ne = conn.get_element_face_edge_indices(elem, face_index, e);

          for (index_t j = 0; j < ne; j++) {
            index_t edge_index = e[j];

            index_t nedge_dof;
            if constexpr (dim == 3) {
              nedge_dof = Basis::get_entity_ndof(basis, ET::EDGE, edge_index);
              Basis::get_entity_dof(basis, ET::EDGE, edge_index, elem_dof, dof);
            } else {
              nedge_dof = 0;
            }

            for (index_t k = 0; k < nedge_dof; k++) {
              bool active = bcinfo.active_for_dof(label, basis,
                                                  k % Basis::get_stride(basis));

              if (active && dof_flags[dof[k]] == NO_LABEL) {
                dof_flags[dof[k]] = 0;

                boundary_dof[boundary_dof_counter] = dof[k];
                boundary_dof_counter++;
              }
            }
          }

          // Get the vertices associated with the face
          index_t v[4];
          index_t nv = conn.get_element_face_vert_indices(elem, face_index, v);

          for (index_t j = 0; j < nv; j++) {
            index_t vert_index = v[j];

            index_t nvert_dof =
                Basis::get_entity_ndof(basis, ET::VERTEX, vert_index);
            Basis::get_entity_dof(basis, ET::VERTEX, vert_index, elem_dof, dof);

            for (index_t k = 0; k < nvert_dof; k++) {
              bool active = bcinfo.active_for_dof(label, basis,
                                                  k % Basis::get_stride(basis));
              if (active && dof_flags[dof[k]] == NO_LABEL) {
                dof_flags[dof[k]] = 0;

                boundary_dof[boundary_dof_counter] = dof[k];
                boundary_dof_counter++;
              }
            }
          }
        }
      }
    }
  }

  DELETE_ARRAY(dof_flags);

  ndof = boundary_dof_counter;
  dof = new index_t[ndof];
  for (index_t i = 0; i < ndof; i++) {
    dof[i] = boundary_dof[i];
  }

  // Free the boundary_dof array
  DELETE_ARRAY(boundary_dof);
}

}  // namespace A2D

#endif  // A2D_FE_MESH_INL_H