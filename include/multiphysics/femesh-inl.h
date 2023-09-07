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
    : nelems(nelems), nbounds(0), nedges(0), nverts(nverts) {
  // Number of bounds of all types
  nline_bounds = 0;
  ntri_bounds = 0;
  nquad_bounds = 0;

  // vert -> element connectivity
  vert_element_ptr = nullptr;
  vert_elements = nullptr;

  // Bound -> element connectivity
  line_bound_elements = nullptr;
  tri_bound_elements = nullptr;
  quad_bound_elements = nullptr;

  // Set to NULL all the boundary info
  num_boundary_labels = 1;
  num_boundary_bounds = 0;
  boundary_labels = nullptr;
  boundary_bounds = nullptr;
}

inline MeshConnectivityBase::~MeshConnectivityBase() {
  if (vert_element_ptr) {
    DELETE_ARRAY(vert_element_ptr);
  }
  if (vert_elements) {
    DELETE_ARRAY(vert_elements);
  }

  if (line_bound_elements) {
    DELETE_ARRAY(line_bound_elements);
  }
  if (tri_bound_elements) {
    DELETE_ARRAY(tri_bound_elements);
  }
  if (quad_bound_elements) {
    DELETE_ARRAY(quad_bound_elements);
  }

  if (boundary_bounds) {
    DELETE_ARRAY(boundary_bounds);
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
  init_bound_data();
  init_edge_data();

  // Count up all the bounds that are on the boundary
  num_boundary_bounds = 0;
  for (index_t bound = 0; bound < nbounds; bound++) {
    index_t e1, e2;
    get_bound_elements(bound, &e1, &e2);
    if (e1 == NO_LABEL || e2 == NO_LABEL) {
      num_boundary_bounds++;
    }
  }

  // Set the boundary labels - 0 for anything on the external boundary of the
  // mesh
  boundary_bounds = new index_t[num_boundary_bounds];
  boundary_labels = new index_t[num_boundary_bounds];
  for (index_t bound = 0, count = 0; bound < nbounds; bound++) {
    index_t e1, e2;
    get_bound_elements(bound, &e1, &e2);
    if (e1 == NO_LABEL || e2 == NO_LABEL) {
      boundary_bounds[count] = bound;
      boundary_labels[count] = 0;
      count++;
    }
  }
}

/**
 * @brief Get the boundary bounds and labels
 *
 * @param boundary_bounds_ Array of the boundary bounds
 * @param boundary_labels_ Array of the boundary labels
 * @return The number of boundary bounds
 */
inline index_t MeshConnectivityBase::get_boundary_bounds(
    const index_t* boundary_bounds_[], const index_t* boundary_labels_[]) {
  if (boundary_bounds_) {
    *boundary_bounds_ = boundary_bounds;
  }
  if (boundary_labels_) {
    *boundary_labels_ = boundary_labels;
  }
  return num_boundary_bounds;
}

/**
 * @brief Count up the number of boundary bounds with the given label
 *
 * @param label The label index
 * @return The number
 */
inline index_t MeshConnectivityBase::get_num_boundary_bounds_with_label(
    const index_t label) {
  index_t count = 0;
  for (index_t i = 0; i < num_boundary_bounds; i++) {
    if (boundary_labels[i] == label) {
      count++;
    }
  }

  return count;
}

/**
 * @brief Add a boundary label from the vertices
 *
 * Any boundary bound with all nodes that touch the given set of vertices is
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

  for (index_t i = 0; i < num_boundary_bounds; i++) {
    index_t bound = boundary_bounds[i];

    // Get the element corresponding to the bound
    index_t elem, e2;
    get_bound_elements(bound, &elem, &e2);

    // Find the bound index for the the boundary bound;
    index_t bound_index = 0;
    const index_t* bounds;
    index_t nf = get_element_bounds(elem, &bounds);
    for (index_t j = 0; j < nf; j++) {
      if (bounds[j] == bound) {
        bound_index = j;
        break;
      }
    }

    // Check if all the vertices are labeled or not
    index_t v[4];
    index_t nv = get_element_bound_verts(elem, bound_index, v);
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

  // Return the label that was applied to the bounds
  if (count > 0) {
    num_boundary_labels++;
    return num_boundary_labels - 1;
  }

  // No bounds were found
  return NO_LABEL;
}

/**
 * @brief Get the element bounds array
 *
 * @param elem The element index
 * @param bounds The bounds indicies associated with the element
 * @return The number of bounds
 */
inline index_t MeshConnectivityBase::get_element_bounds(
    index_t elem, const index_t* bounds[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nbounds = meta.get_nbounds();
  *bounds = &(*meta.bounds)[nbounds * elem];
  return nbounds;
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
  *verts = &(*meta.verts)[meta.get_nverts() * elem];
  return meta.get_nverts();
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
 * @brief Get the global indices of the bound in its local order
 *
 * @param elem The element index
 * @param b The local bound index
 * @param verts The vertices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_bound_verts(index_t elem,
                                                             index_t b,
                                                             index_t verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nverts = meta.get_bound_nverts(b);
  for (index_t i = 0; i < nverts; i++) {
    verts[i] =
        (*meta.verts)[meta.get_nverts() * elem + meta.get_bound_vert(b, i)];
  }
  return nverts;
}

/**
 * @brief Get the local indices of the bound in its local order
 *
 * @param elem The element index
 * @param b The local bound index
 * @param verts The vertices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_bound_vert_indices(
    index_t elem, index_t b, index_t verts[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nverts = meta.get_bound_nverts(b);
  for (index_t i = 0; i < nverts; i++) {
    verts[i] = meta.get_bound_vert(b, i);
  }
  return nverts;
}

/**
 * @brief Get the edges of the bound in its local order
 *
 * @param elem The element index
 * @param b The local bound index
 * @param edges The edges
 * @return The number of edges
 */
inline index_t MeshConnectivityBase::get_element_bound_edges(index_t elem,
                                                             index_t b,
                                                             index_t edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nedges = meta.get_bound_nedges(b);

  if (nedges == 0) {
    return nedges;
  }

  for (index_t i = 0; i < nedges; i++) {
    edges[i] =
        (*meta.edges)[meta.get_nedges() * elem + meta.get_bound_edge(b, i)];
  }

  return nedges;
}

/**
 * @brief Get the indices of the bound in its local order
 *
 * @param elem The element index
 * @param b The local bound index
 * @param edges The vertices
 * @return The number of edges
 */
inline index_t MeshConnectivityBase::get_element_bound_edge_indices(
    index_t elem, index_t b, index_t edges[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nedges = meta.get_bound_nedges(b);

  if (nedges == 0) {
    return nedges;
  }

  for (index_t i = 0; i < nedges; i++) {
    edges[i] = meta.get_bound_edge(b, i);
  }

  return nedges;
}

/**
 * @brief Get the bound element vertices in a sorted order
 *
 * This call will return the vertices in the same order regardless of the
 * orientation of the bound. This is used to compare bounds.
 *
 * @param elem The element index
 * @param local_bound The local bound index
 * @param verts The vertices
 * @return The number of vertices
 */
inline index_t MeshConnectivityBase::get_element_global_bound_verts(
    index_t elem, index_t local_bound, index_t verts[]) {
  index_t t[ET::MAX_BOUND_VERTS];
  index_t n = get_element_bound_verts(elem, local_bound, t);

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
    snprintf(msg, 256, "%d is an invalid number of vertices for a bound", n);
    throw std::runtime_error(msg);
  }
  return n;
}

/**
 * @brief Check for equality between two bounds
 *
 * @return True if the bounds match, false otherwise
 */
inline bool MeshConnectivityBase::global_bound_equality(index_t na,
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
 * @brief Get the bound elements associated with the given bound
 *
 * @param bound Global bound index
 * @param e1 Returned value of the first element
 * @param e2 Returned value of the second element, NO_INDEX indicates second
 * element doesn't exist, and this bound is an external boundary bound
 * @return Boolean if this bound is on the boundary
 */
inline bool MeshConnectivityBase::get_bound_elements(index_t bound, index_t* e1,
                                                     index_t* e2) {
  if (bound < nline_bounds) {
    *e1 = line_bound_elements[2 * bound];
    *e2 = line_bound_elements[2 * bound + 1];
    return (line_bound_elements[2 * bound + 1] == NO_LABEL);
  } else if (bound < ntri_bounds) {
    bound = bound - nline_bounds;
    *e1 = tri_bound_elements[2 * bound];
    *e2 = tri_bound_elements[2 * bound + 1];
    return (tri_bound_elements[2 * bound + 1] == NO_LABEL);
  } else if (bound < nquad_bounds) {
    bound = bound - nline_bounds - ntri_bounds;
    *e1 = quad_bound_elements[2 * bound];
    *e2 = quad_bound_elements[2 * bound + 1];
    return (quad_bound_elements[2 * bound + 1] == NO_LABEL);
  } else {
    char msg[256];
    std::snprintf(msg, 256, "%d is not a valid global bound index", bound);
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
  index_t nedges = meta.get_nedges();
  if (nedges == 0) {
    *edges = nullptr;
  } else {
    *edges = &(*meta.edges)[nedges * elem];
  }
  return nedges;
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
  verts[0] =
      (*meta.verts)[meta.get_nverts() * elem + meta.get_edge_vert(edge, 0)];
  verts[1] =
      (*meta.verts)[meta.get_nverts() * elem + meta.get_edge_vert(edge, 1)];
}

template <typename I>
inline MeshConnectivity2D::MeshConnectivity2D(I nverts, I ntri, I* tri, I nquad,
                                              I* quad)
    : MeshConnectivityBase(nverts, ntri + nquad),
      ntri(ntri),
      nquad(nquad),
      meta_tri(MetaDataFactory::create_2d_meta(ET::Element::Tri, &tri_bounds,
                                               &tri_verts)),
      meta_quad(MetaDataFactory::create_2d_meta(ET::Element::Quad, &quad_bounds,
                                                &quad_verts)) {
  Timer timer("MeshConnectivity2D");
  // Allocate space for the element -> vert connectivity
  tri_verts = new index_t[ET::TRI_NVERTS * ntri];
  quad_verts = new index_t[ET::QUAD_NVERTS * nquad];

  // Set the connectivity: element -> verts
  for (index_t i = 0; i < ET::TRI_NVERTS * ntri; i++) {
    tri_verts[i] = tri[i];
  }
  for (index_t i = 0; i < ET::QUAD_NVERTS * nquad; i++) {
    quad_verts[i] = quad[i];
  }

  // Allocate the bound data
  tri_bounds = new index_t[ET::TRI_NBOUNDS * ntri];
  quad_bounds = new index_t[ET::QUAD_NBOUNDS * nquad];

  // Fill the bound arrays with NO_LABEL
  std::fill(tri_bounds, tri_bounds + ET::TRI_NBOUNDS * ntri, NO_LABEL);
  std::fill(quad_bounds, quad_bounds + ET::QUAD_NBOUNDS * nquad, NO_LABEL);

  initialize();
}

inline MeshConnectivity2D::~MeshConnectivity2D() {
  if (tri_verts) {
    DELETE_ARRAY(tri_verts);
  }
  if (quad_verts) {
    DELETE_ARRAY(quad_verts);
  }

  if (tri_bounds) {
    DELETE_ARRAY(tri_bounds);
  }
  if (quad_bounds) {
    DELETE_ARRAY(quad_bounds);
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
    char msg[256];
    std::snprintf(msg, 256, "element index %d is beyond the range [0, %d)",
                  elem, get_num_elements());
    throw std::runtime_error(msg);
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
      meta_tet(MetaDataFactory::create_3d_meta(ET::Element::Tet, &tet_bounds,
                                               &tet_edges, &tet_verts)),
      meta_hex(MetaDataFactory::create_3d_meta(ET::Element::Hex, &hex_bounds,
                                               &hex_edges, &hex_verts)),
      meta_wedge(MetaDataFactory::create_3d_meta(
          ET::Element::Wedge, &wedge_bounds, &wedge_edges, &wedge_verts)),
      meta_pyrmd(MetaDataFactory::create_3d_meta(
          ET::Element::Pyrmd, &pyrmd_bounds, &pyrmd_edges, &pyrmd_verts)) {
  Timer timer("MeshConnectivity3D");
  // Allocate space for the element -> vert connectivity
  tet_verts = new index_t[ET::TET_NVERTS * ntets];
  hex_verts = new index_t[ET::HEX_NVERTS * nhex];
  wedge_verts = new index_t[ET::WEDGE_NVERTS * nwedge];
  pyrmd_verts = new index_t[ET::PYRMD_NVERTS * npyrmd];

  // Set the connectivity: element -> verts
  for (index_t i = 0; i < ET::TET_NVERTS * ntets; i++) {
    tet_verts[i] = tets[i];
  }
  for (index_t i = 0; i < ET::HEX_NVERTS * nhex; i++) {
    hex_verts[i] = hex[i];
  }
  for (index_t i = 0; i < ET::WEDGE_NVERTS * nwedge; i++) {
    wedge_verts[i] = wedge[i];
  }
  for (index_t i = 0; i < ET::PYRMD_NVERTS * npyrmd; i++) {
    pyrmd_verts[i] = pyrmd[i];
  }

  // Allocate the bound data
  tet_bounds = new index_t[ET::TET_NBOUNDS * ntets];
  hex_bounds = new index_t[ET::HEX_NBOUNDS * nhex];
  wedge_bounds = new index_t[ET::WEDGE_NBOUNDS * nwedge];
  pyrmd_bounds = new index_t[ET::PYRMD_NBOUNDS * npyrmd];

  // Fill the bound arrays with NO_LABEL
  std::fill(tet_bounds, tet_bounds + ET::TET_NBOUNDS * ntets, NO_LABEL);
  std::fill(hex_bounds, hex_bounds + ET::HEX_NBOUNDS * nhex, NO_LABEL);
  std::fill(wedge_bounds, wedge_bounds + ET::WEDGE_NBOUNDS * nwedge, NO_LABEL);
  std::fill(pyrmd_bounds, pyrmd_bounds + ET::PYRMD_NBOUNDS * npyrmd, NO_LABEL);

  // element -> edge connectivity
  tet_edges = new index_t[ET::TET_NEDGES * ntets];
  hex_edges = new index_t[ET::HEX_NEDGES * nhex];
  wedge_edges = new index_t[ET::WEDGE_NEDGES * nwedge];
  pyrmd_edges = new index_t[ET::PYRMD_NEDGES * npyrmd];

  // Fill the bound arrays with NO_LABEL
  std::fill(tet_edges, tet_edges + ET::TET_NEDGES * ntets, NO_LABEL);
  std::fill(hex_edges, hex_edges + ET::HEX_NEDGES * nhex, NO_LABEL);
  std::fill(wedge_edges, wedge_edges + ET::WEDGE_NEDGES * nwedge, NO_LABEL);
  std::fill(pyrmd_edges, pyrmd_edges + ET::PYRMD_NEDGES * npyrmd, NO_LABEL);

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

  if (tet_bounds) {
    DELETE_ARRAY(tet_bounds);
  }
  if (hex_bounds) {
    DELETE_ARRAY(hex_bounds);
  }
  if (wedge_bounds) {
    DELETE_ARRAY(wedge_bounds);
  }
  if (pyrmd_bounds) {
    DELETE_ARRAY(pyrmd_bounds);
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
    char msg[256];
    std::snprintf(msg, 256, "element index %d is beyond the range [0, %d)",
                  elem, get_num_elements());
    throw std::runtime_error(msg);
  }
}

/**
 * @brief Label the verts, edges and bounds that touch the list of vertices
 *
 * @param nv The number of input vertices
 * @param verts The vertex numbers
 * @param vert_labels An array of length nverts
 * @param edge_labels An array of length nedges
 * @param bound_labels An array of length nbounds
 */
template <typename IdxType>
inline void MeshConnectivityBase::get_labels_from_verts(
    const index_t nv, const IdxType verts[], index_t vert_labels[],
    index_t edge_labels[], index_t bound_labels[]) {
  std::fill(vert_labels, vert_labels + nverts, NO_LABEL);
  std::fill(edge_labels, edge_labels + nedges, NO_LABEL);
  std::fill(bound_labels, bound_labels + nbounds, NO_LABEL);

  // Label the vertices
  for (index_t i = 0; i < nv; i++) {
    if (verts[i] < nverts) {
      vert_labels[verts[i]] = i;
    }
  }

  // Loop over elements and label the edges and bounds
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

    const index_t* elem_bounds;
    index_t nf = get_element_bounds(elem, &elem_bounds);
    for (index_t f = 0; f < nf; f++) {
      index_t fv[ET::MAX_BOUND_VERTS];  // bound vertices
      index_t nfv = get_element_bound_verts(elem, f, fv);

      // Check if all the bound vertices are labeled
      if (nfv == 2) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL) {
          bound_labels[elem_bounds[f]] = elem_bounds[f];
        }
      } else if (nfv == 3) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL &&
            vert_labels[fv[2]] != NO_LABEL) {
          bound_labels[elem_bounds[f]] = elem_bounds[f];
        }
      } else if (nfv == 4) {
        if (vert_labels[fv[0]] != NO_LABEL && vert_labels[fv[1]] != NO_LABEL &&
            vert_labels[fv[2]] != NO_LABEL && vert_labels[fv[3]] != NO_LABEL) {
          bound_labels[elem_bounds[f]] = elem_bounds[f];
        }
      }
    }
  }
}

/**
 * @brief Get a non-constant element bounds array
 *
 * @param elem The element index
 * @param bounds The bound indicies associated with the element
 * @return The number of bounds
 */
inline index_t MeshConnectivityBase::get_element_bounds(index_t elem,
                                                        index_t* bounds[]) {
  const ElemConnMetaData& meta = get_local_elem_and_meta(elem);
  index_t nbounds = meta.get_nbounds();
  *bounds = &(*meta.bounds)[nbounds * elem];
  return nbounds;
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
  index_t nedges = meta.get_nedges();
  if (nedges == 0) {
    *edges = nullptr;
  } else {
    *edges = &(*meta.edges)[nedges * elem];
  }
  return nedges;
}

/**
 * @brief Initialize data associated with the bound information
 *
 * This code uniquely orders the bounds associated with each element
 */
inline void MeshConnectivityBase::init_bound_data() {
  // Prepare to count and number the number of bounds. This keeps track of
  // counts of separate bound types (line, triangle, quadrilateral)
  for (index_t elem = 0; elem < nelems; elem++) {
    // Loop over the elements of this bound
    index_t* bounds;
    const index_t nf = get_element_bounds(elem, &bounds);

    for (index_t bound = 0; bound < nf; bound++) {
      if (bounds[bound] == NO_LABEL) {
        // Get the unique set of verts corresponding to this bound
        index_t bound_verts[ET::MAX_BOUND_VERTS];
        index_t nbound_verts =
            get_element_global_bound_verts(elem, bound, bound_verts);

        bool bound_located = false;

        // For all adjacent elements check if there is a match
        for (index_t m = 0; m < nbound_verts; m++) {
          const index_t* adj_elems;
          const index_t nadj_elems =
              get_adjacent_elements_from_vert(bound_verts[m], &adj_elems);

          // Loop over all adjacent elements that are not the current
          // element
          for (index_t k = 0; k < nadj_elems; k++) {
            if (adj_elems[k] != elem) {
              // Loop over the elements of this bound
              index_t* adj_bounds;
              const index_t adj_nf =
                  get_element_bounds(adj_elems[k], &adj_bounds);

              for (index_t adj_bound = 0; adj_bound < adj_nf; adj_bound++) {
                // Get the unique set of vertices corresponding to this
                // bound
                index_t adj_bound_verts[4];
                index_t adj_nbound_verts = get_element_global_bound_verts(
                    adj_elems[k], adj_bound, adj_bound_verts);

                if (adj_bounds[adj_bound] == NO_LABEL &&
                    global_bound_equality(nbound_verts, bound_verts,
                                          adj_nbound_verts, adj_bound_verts)) {
                  if (nbound_verts == 2) {
                    adj_bounds[adj_bound] = nline_bounds;
                    bounds[bound] = nline_bounds;
                    nline_bounds++;
                  } else if (nbound_verts == 3) {
                    adj_bounds[adj_bound] = ntri_bounds;
                    bounds[bound] = ntri_bounds;
                    ntri_bounds++;
                  } else {
                    adj_bounds[adj_bound] = nquad_bounds;
                    bounds[bound] = nquad_bounds;
                    nquad_bounds++;
                  }

                  bound_located = true;
                  break;
                }
              }

              if (bound_located) {
                break;
              }
            }
          }

          if (bound_located) {
            break;
          }
        }

        // No adjacent bound was found. This is a boundary bound
        if (!bound_located) {
          if (nbound_verts == 2) {
            bounds[bound] = nline_bounds;
            nline_bounds++;
          } else if (nbound_verts == 3) {
            bounds[bound] = ntri_bounds;
            ntri_bounds++;
          } else {
            bounds[bound] = nquad_bounds;
            nquad_bounds++;
          }

          // Add this bound to the boundary bound list???
        }
      }
    }
  }

  // Sum up the total number of bounds
  nbounds = nline_bounds + ntri_bounds + nquad_bounds;

  // Set the bound indices associated with the two adjacent elements. At
  // this point the bound indices stored are with respect to separate lists
  // of triangular and quadrilateral bounds. We order the triangular bounds
  // first, so we have to add the total number of triangular bounds to each
  // quadrilateral bound index to get the global index.
  line_bound_elements = new index_t[2 * nline_bounds];
  std::fill(line_bound_elements, line_bound_elements + 2 * nline_bounds,
            NO_LABEL);

  tri_bound_elements = new index_t[2 * ntri_bounds];
  std::fill(tri_bound_elements, tri_bound_elements + 2 * ntri_bounds, NO_LABEL);

  quad_bound_elements = new index_t[2 * nquad_bounds];
  std::fill(quad_bound_elements, quad_bound_elements + 2 * nquad_bounds,
            NO_LABEL);

  for (index_t elem = 0; elem < nelems; elem++) {
    // Loop over the elements of this bound
    index_t* bounds;
    const index_t nf = get_element_bounds(elem, &bounds);

    for (index_t bound = 0; bound < nf; bound++) {
      index_t bound_verts[ET::MAX_BOUND_VERTS];
      index_t nbound_verts =
          get_element_global_bound_verts(elem, bound, bound_verts);

      if (nbound_verts == 2) {
        if (line_bound_elements[2 * bounds[bound]] == NO_LABEL) {
          line_bound_elements[2 * bounds[bound]] = elem;
        } else if (line_bound_elements[2 * bounds[bound] + 1] == NO_LABEL) {
          line_bound_elements[2 * bounds[bound] + 1] = elem;
        }
      } else if (nbound_verts == 3) {
        if (tri_bound_elements[2 * bounds[bound]] == NO_LABEL) {
          tri_bound_elements[2 * bounds[bound]] = elem;
        } else if (tri_bound_elements[2 * bounds[bound] + 1] == NO_LABEL) {
          tri_bound_elements[2 * bounds[bound] + 1] = elem;
        }

        // Reset the bound index into the global bound index - add the
        // number of line bounds. Now the bounds will index from a
        // global bound number.
        bounds[bound] += nline_bounds;
      } else {
        if (quad_bound_elements[2 * bounds[bound]] == NO_LABEL) {
          quad_bound_elements[2 * bounds[bound]] = elem;
        } else if (quad_bound_elements[2 * bounds[bound] + 1] == NO_LABEL) {
          quad_bound_elements[2 * bounds[bound] + 1] = elem;
        }

        // Reset the bound index into the global bound index - add the
        // number of line and triangle bounds. Now the bounds will index from a
        // global bound number.
        bounds[bound] += nline_bounds + ntri_bounds;
      }
    }
  }
}

/**
 * @brief Initialize and order the edge information.
 *
 * This relies on the bound connectivity data - so that must be initialized
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
  num_dof_offset = DofOffsetArray("num_dof_offset");
  element_dof = ElementDofArray("element_dof", nelems);
  element_sign = ElementSignArray("element_sign", nelems);
  BLAS::fill(element_dof, NO_INDEX);

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
      const index_t* bounds;
      index_t nf = conn.get_element_bounds(stack[i], &bounds);

      for (index_t j = 0; j < nf; j++) {
        // Look at the adjacent elements
        index_t e1, e2;
        conn.get_bound_elements(bounds[j], &e1, &e2);

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

  std::vector<index_t> bound_owners(conn.get_num_bounds(), NO_INDEX);
  std::vector<index_t> edge_owners(conn.get_num_edges(), NO_INDEX);
  std::vector<index_t> vert_owners(conn.get_num_verts(), NO_INDEX);

  index_t dof_counter = 0;
  index_t dof[ndof_per_element];

  for (index_t basis = 0; basis < Basis::nbasis; basis++) {
    for (index_t counter = nelems; counter > 0; counter--) {
      index_t elem = stack[counter - 1];

      auto elem_dof = Kokkos::subview(element_dof, elem, Kokkos::ALL);
      auto elem_sign = Kokkos::subview(element_sign, elem, Kokkos::ALL);

      // The domain DOF are always owned by the element - no need to check
      // for the element that owns them
      index_t ndof = Basis::get_entity_ndof(basis, ET::Domain, 0);

      for (index_t i = 0; i < ndof; i++, dof_counter++) {
        dof[i] = dof_counter;
      }
      Basis::set_entity_dof(basis, ET::Domain, 0, 0, dof, elem_dof);
      Basis::set_entity_signs(basis, ET::Domain, 0, 0, elem_sign);

      // Order the bounds
      const index_t* bounds;
      index_t nf = conn.get_element_bounds(elem, &bounds);
      for (index_t index = 0; index < nf; index++) {
        index_t bound = bounds[index];
        index_t orient = 0;
        if (bound_owners[bound] == NO_INDEX || bound_owners[bound] == elem) {
          bound_owners[bound] = elem;

          ndof = Basis::get_entity_ndof(basis, ET::Bound, index);

          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = bound_owners[bound];
          const index_t* owner_bounds;
          index_t nf_owner = conn.get_element_bounds(owner_elem, &owner_bounds);

          auto owner_elem_dof =
              Kokkos::subview(element_dof, owner_elem, Kokkos::ALL);

          for (index_t i = 0; i < nf_owner; i++) {
            if (owner_bounds[i] == bound) {
              index_t ref[ET::MAX_BOUND_VERTS], verts[ET::MAX_BOUND_VERTS];
              index_t nverts = conn.get_element_bound_verts(owner_elem, i, ref);
              conn.get_element_bound_verts(elem, index, verts);

              Basis::get_entity_dof(basis, ET::Bound, i, owner_elem_dof, dof);

              if (nverts == 4) {
                orient = ET::get_quad_domain_orientation(ref, verts);
              } else if (nverts == 2) {
                if (ref[0] == verts[1] && ref[1] == verts[0]) {
                  orient = 1;
                }
              }
              break;
            }
          }
        }

        Basis::set_entity_dof(basis, ET::Bound, index, orient, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::Bound, index, orient, elem_sign);
      }

      // Order the edges, only effective for dim == 3
      const index_t* edges;
      index_t ne = conn.get_element_edges(elem, &edges);
      for (index_t index = 0; index < ne; index++) {
        index_t edge = edges[index];
        index_t orient = 0;
        if (edge_owners[edge] == NO_INDEX || edge_owners[edge] == elem) {
          edge_owners[edge] = elem;

          ndof = Basis::get_entity_ndof(basis, ET::Edge, index);
          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = edge_owners[edge];
          const index_t* owner_edges;
          index_t ne_owner = conn.get_element_edges(owner_elem, &owner_edges);

          auto owner_elem_dof =
              Kokkos::subview(element_dof, owner_elem, Kokkos::ALL);

          for (index_t i = 0; i < ne_owner; i++) {
            if (owner_edges[i] == edge) {
              index_t ref[2], verts[2];
              conn.get_element_edge_verts(owner_elem, i, ref);
              conn.get_element_edge_verts(elem, index, verts);

              Basis::get_entity_dof(basis, ET::Edge, i, owner_elem_dof, dof);

              if (ref[0] == verts[1] && ref[1] == verts[0]) {
                orient = 1;
              }
              break;
            }
          }
        }

        Basis::set_entity_dof(basis, ET::Edge, index, orient, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::Edge, index, orient, elem_sign);
      }

      // Order the vertices (only effective for dim == 2 or 3)
      const index_t* verts;
      index_t nv = conn.get_element_verts(elem, &verts);
      for (index_t index = 0; index < nv; index++) {
        index_t vert = verts[index];
        index_t orient = 0;
        if (vert_owners[vert] == NO_INDEX || vert_owners[vert] == elem) {
          vert_owners[vert] = elem;

          ndof = Basis::get_entity_ndof(basis, ET::Vertex, index);
          for (index_t i = 0; i < ndof; i++, dof_counter++) {
            dof[i] = dof_counter;
          }
        } else {
          index_t owner_elem = vert_owners[vert];
          const index_t* owner_verts;
          index_t nv_owner = conn.get_element_verts(owner_elem, &owner_verts);

          auto owner_elem_dof =
              Kokkos::subview(element_dof, owner_elem, Kokkos::ALL);

          for (index_t i = 0; i < nv_owner; i++) {
            if (owner_verts[i] == vert) {
              Basis::get_entity_dof(basis, ET::Vertex, i, owner_elem_dof, dof);
              break;
            }
          }
        }

        Basis::set_entity_dof(basis, ET::Vertex, index, orient, dof, elem_dof);
        Basis::set_entity_signs(basis, ET::Vertex, index, orient, elem_sign);
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
 * Note: typically Basis in this case is in a lower dimension (dim=2, e.g.) and
 * InteroirBasis is in a higher dimension (dim=3, e.g.)
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
    : nelems(conn.get_num_boundary_bounds_with_label(label)) {
  num_dof_offset = DofOffsetArray("num_dof_offset");
  element_dof = ElementDofArray("element_dof", nelems);
  element_sign = ElementSignArray("element_sign", nelems);

  // Get the number of boundary bounds
  const index_t* boundary_bounds;
  const index_t* boundary_labels;
  index_t num_boundary_bounds =
      conn.get_boundary_bounds(&boundary_bounds, &boundary_labels);

  for (index_t i = 0, elem_count = 0; i < num_boundary_bounds; i++) {
    if (boundary_labels[i] == label) {
      index_t bound = boundary_bounds[i];

      // Get the element adjacent to the boundary bound (e2 can be ignored)
      index_t elem, e2;
      conn.get_bound_elements(bound, &elem, &e2);

      // Find the bound index for the element
      index_t bound_index = 0;
      const index_t* bounds;
      index_t nf = conn.get_element_bounds(elem, &bounds);
      for (index_t index = 0; index < nf; index++) {
        if (bound == bounds[index]) {
          bound_index = index;
        }
      }

      // Get the signs associated wtih the original mesh
      auto dof = mesh.get_element_dof(elem);

      // Set pointers for the entity dof
      auto surf_dof = Kokkos::subview(element_dof, elem_count, Kokkos::ALL);
      auto surf_signs = Kokkos::subview(element_sign, elem_count, Kokkos::ALL);

      // The degree of freedom indices for each entity
      index_t entity_dof[InteriorBasis::ndof];

      // If this dim == 2:
      //   - set vertices from interior vertices
      // If this dim == 1:
      //   - set bounds from interior vertices
      index_t v[4];
      index_t nv = conn.get_element_bound_vert_indices(elem, bound_index, v);
      for (index_t j = 0; j < nv; j++) {
        for (index_t basis = 0; basis < Basis::nbasis; basis++) {
          index_t vert_index = v[j];

          // Get the vertex dof from the interior element
          InteriorBasis::get_entity_dof(basis, ET::Vertex, vert_index, dof,
                                        entity_dof);

          index_t surf_vert_index = j;
          index_t orient = 0;

          // Set the same vertex index on the corresponding surface
          ET::ElementEntity dest_entity;
          if constexpr (dim == 2) {
            dest_entity = ET::Vertex;
          } else if constexpr (dim == 1) {
            dest_entity = ET::Bound;
          }
          Basis::set_entity_dof(basis, dest_entity, surf_vert_index, orient,
                                entity_dof, surf_dof);
          Basis::set_entity_signs(basis, dest_entity, surf_vert_index, orient,
                                  surf_signs);
        }
      }

      // Set bounds from interior edges
      // Note: this is only effective when dim of this basis is 2
      index_t e[4];
      index_t ne = conn.get_element_bound_edge_indices(elem, bound_index, e);
      for (index_t j = 0; j < ne; j++) {
        for (index_t basis = 0; basis < Basis::nbasis; basis++) {
          index_t edge_index = e[j];

          // Get domain dof for this basis from the bonud dof of the interior
          // element
          InteriorBasis::get_entity_dof(basis, ET::Edge, edge_index, dof,
                                        entity_dof);

          // Get the edge orientation relative to the bound
          index_t orient = 0;

          // Get the index of the edge on the bound
          index_t surf_edge_index = j;

          // Set the same edge on the corresponding surface
          Basis::set_entity_dof(basis, ET::Bound, surf_edge_index, orient,
                                entity_dof, surf_dof);
          Basis::set_entity_signs(basis, ET::Bound, surf_edge_index, orient,
                                  surf_signs);
        }
      }

      // Set domain from interior bound
      for (index_t basis = 0; basis < InteriorBasis::nbasis; basis++) {
        // Get the degrees of freeom from the element
        index_t orient = 0;
        InteriorBasis::get_entity_dof(basis, ET::Bound, bound_index, dof,
                                      entity_dof);

        // Set the degrees of freedom - the bound element has the same
        // orientation as its interior owner
        Basis::set_entity_dof(basis, ET::Domain, 0, orient, entity_dof,
                              surf_dof);
        Basis::set_entity_signs(basis, ET::Domain, 0, orient, surf_signs);
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
ElementMesh<Basis>::ElementMesh(const ElementMesh<HOrderBasis>& mesh)
    : nelems(HOrderBasis::get_num_lorder_elements() * mesh.get_num_elements()),
      num_dof(mesh.get_num_dof()) {
  num_dof_offset = DofOffsetArray("num_dof_offset");
  element_dof = ElementDofArray("element_dof", nelems);
  element_sign = ElementSignArray("element_sign", nelems);

  for (index_t i = 0; i < Basis::nbasis; i++) {
    num_dof_offset[i] = mesh.get_num_cumulative_dof(i);
  }

  // Loop over all the high-order elements and generate the low-order
  // representation
  for (index_t j = 0; j < mesh.get_num_elements(); j++) {
    // Get the signs associated wtih the high-order mesh
    auto horder_dof = mesh.get_element_dof(j);
    auto horder_signs = mesh.get_element_signs(j);

    for (index_t i = 0; i < HOrderBasis::get_num_lorder_elements(); i++) {
      index_t elem = i + j * HOrderBasis::get_num_lorder_elements();

      // Index into the high-order element
      auto lorder_dof = Kokkos::subview(element_dof, elem, Kokkos::ALL);
      auto lorder_signs = Kokkos::subview(element_sign, elem, Kokkos::ALL);
      HOrderBasis::get_lorder_dof(i, horder_dof, lorder_dof);

      // Signs indicating any orientation flip
      HOrderBasis::get_lorder_signs(i, horder_signs, lorder_signs);
    }
  }
}

template <class Basis>
template <index_t basis>
int ElementMesh<Basis>::get_global_dof_sign(index_t elem, index_t index) const {
  return element_sign(elem, Basis::template get_dof_offset<basis>() + index);
}

template <class Basis>
template <index_t basis>
index_t ElementMesh<Basis>::get_global_dof(index_t elem, index_t index) const {
  return element_dof(elem, Basis::template get_dof_offset<basis>() + index);
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
      index_t row = element_dof(i, j1) / M;
      while (j1 + 1 < ndof_per_elem && row == (element_dof(i, j1 + 1) / M)) {
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
 * degrees of freedom from the vertices, edges and bounds that touch the
 * specified degrees of freedom
 *
 * @param conn The connectivity
 * @param mesh The mesh object with the ordered degrees of freedom
 */
template <class Basis>
DirichletBCs<Basis>::DirichletBCs(MeshConnectivityBase& conn,
                                  ElementMesh<Basis>& mesh,
                                  DirichletBCInfo& bcinfo) {
  const index_t* boundary_bounds;
  const index_t* boundary_labels;
  index_t num_boundary_bounds =
      conn.get_boundary_bounds(&boundary_bounds, &boundary_labels);

  index_t num_dof = mesh.get_num_dof();
  index_t boundary_dof_counter = 0;
  index_t* dof_flags = new index_t[num_dof];
  index_t* boundary_dof = new index_t[num_dof];

  const index_t NO_LABEL = 1;
  std::fill(dof_flags, dof_flags + num_dof, NO_LABEL);

  // Loop over all the boundaries
  for (index_t boundary_bound = 0; boundary_bound < num_boundary_bounds;
       boundary_bound++) {
    index_t bound = boundary_bounds[boundary_bound];
    index_t label = boundary_labels[boundary_bound];

    // Find the edges and
    if (bcinfo.active_for_label(label)) {
      // Find the element for this bound
      index_t elem, e2;
      conn.get_bound_elements(bound, &elem, &e2);

      // Get the degrees of freedom for this element
      auto elem_dof = mesh.get_element_dof(elem);

      // Array of dof that may be extracted
      index_t dof[Basis::ndof];

      // Find the bound index for the active bound
      index_t bound_index = 0;
      const index_t* bounds;
      index_t nf = conn.get_element_bounds(elem, &bounds);
      for (index_t i = 0; i < nf; i++) {
        if (bounds[i] == bound) {
          bound_index = i;
          break;
        }
      }

      for (index_t basis = 0; basis < Basis::nbasis; basis++) {
        if (bcinfo.active_for_basis(label, basis)) {
          // Set the boundary conditions on the bounds
          index_t nbound_dof =
              Basis::get_entity_ndof(basis, ET::Bound, bound_index);
          Basis::get_entity_dof(basis, ET::Bound, bound_index, elem_dof, dof);

          for (index_t k = 0; k < nbound_dof; k++) {
            bool active = bcinfo.active_for_dof(label, basis,
                                                k % Basis::get_stride(basis));

            if (active && dof_flags[dof[k]] == NO_LABEL) {
              dof_flags[dof[k]] = 0;

              boundary_dof[boundary_dof_counter] = dof[k];
              boundary_dof_counter++;
            }
          }

          // Set the boundary conditions on the edges, only effective for dim =
          // 3
          index_t e[4];
          int ne = conn.get_element_bound_edge_indices(elem, bound_index, e);
          for (index_t j = 0; j < ne; j++) {
            index_t edge_index = e[j];

            index_t nedge_dof;
            nedge_dof = Basis::get_entity_ndof(basis, ET::Edge, edge_index);
            Basis::get_entity_dof(basis, ET::Edge, edge_index, elem_dof, dof);

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

          // Set the boundary conditions on the vertices, only effective for dim
          // = 2 or 3
          index_t v[4];
          index_t nv =
              conn.get_element_bound_vert_indices(elem, bound_index, v);

          for (index_t j = 0; j < nv; j++) {
            index_t vert_index = v[j];

            index_t nvert_dof =
                Basis::get_entity_ndof(basis, ET::Vertex, vert_index);
            Basis::get_entity_dof(basis, ET::Vertex, vert_index, elem_dof, dof);

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