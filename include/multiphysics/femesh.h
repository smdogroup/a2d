#ifndef A2D_FE_MESH_H
#define A2D_FE_MESH_H

#include <algorithm>
#include <set>

#include "a2ddefs.h"
#include "multiphysics/feelementtypes.h"
#include "sparse/sparse_matrix.h"
#include "sparse/sparse_symbolic.h"
#include "utils/a2dprofiler.h"

namespace A2D {

// Meta data for a specific type of elements
// Note: const qualifier is used as many as possible
struct ElemConnMetaData {
 private:
  using ET = ElementTypes;

  ElemConnMetaData(index_t nbounds, index_t nedges, index_t nverts,
                   const index_t* bound_nverts,
                   const index_t (*bound_verts)[ET::MAX_BOUND_VERTS],
                   const index_t* bound_nedges,
                   const index_t (*bound_edges)[ET::MAX_BOUND_EDGES],
                   const index_t (*edge_verts)[2], index_t* const* bounds,
                   index_t* const* edges, const index_t* const* verts)
      : NBOUNDS(nbounds),
        NEDGES(nedges),
        NVERTS(nverts),
        BOUND_NVERTS(bound_nverts),
        BOUND_VERTS(bound_verts),
        BOUND_NEDGES(bound_nedges),
        BOUND_EDGES(bound_edges),
        EDGE_VERTS(edge_verts),
        bounds(bounds),
        edges(edges),
        verts(verts) {}

  // Meta data objects can only be created by the factory class
  friend class MetaDataFactory;

 public:
  index_t get_nbounds() const { return NBOUNDS; };
  index_t get_nedges() const { return NEDGES; };
  index_t get_nverts() const { return NVERTS; };

  // Get number of vertices of a given bound
  index_t get_bound_nverts(index_t b) const { return BOUND_NVERTS[b]; }

  // Get index of vertex v in bound b
  index_t get_bound_vert(index_t b, index_t v) const {
    return BOUND_VERTS[b][v];
  }

  // Get number of edges of a given bound, return 0 if bound has no edge at all
  index_t get_bound_nedges(index_t b) const {
    return BOUND_NEDGES ? BOUND_NEDGES[b] : 0;
  }

  // Get index of edge e in bound b
  index_t get_bound_edge(index_t b, index_t e) const {
    return BOUND_EDGES ? BOUND_EDGES[b][e] : NO_INDEX;
  }

  // Get index of vert v in edge e
  index_t get_edge_vert(index_t e, index_t v) const {
    return EDGE_VERTS ? EDGE_VERTS[e][v] : NO_INDEX;
  }

 private:
  // Number of quantities for a single element of this type
  const index_t NBOUNDS;
  const index_t NEDGES;  // 3D element only
  const index_t NVERTS;

  // Bound -> vertices
  const index_t* const BOUND_NVERTS;
  const index_t (*const BOUND_VERTS)[ET::MAX_BOUND_VERTS];

  // Bound -> edges (3D element only)
  const index_t* const BOUND_NEDGES;
  const index_t (*const BOUND_EDGES)[ET::MAX_BOUND_EDGES];

  // Edge -> vertices (3D element only)
  const index_t (*const EDGE_VERTS)[2];

 public:
  // Connectivity: Elem -> ...
  index_t* const* const bounds;       // element -> bounds
  index_t* const* const edges;        // element -> edge (3D only)
  const index_t* const* const verts;  // element -> vertex
};

// A factory class to safely cerate meta data
class MetaDataFactory {
 private:
  using ET = ElementTypes;

 public:
  static ElemConnMetaData create_2d_meta(ET::Element element,
                                         index_t* const* bounds,
                                         const index_t* const* verts) {
    switch (element) {
      case ET::Element::Tri:
        return create_tri(bounds, verts);
      case ET::Element::Quad:
        return create_quad(bounds, verts);
      default:
        char msg[256];
        std::snprintf(msg, 256,
                      "element enumerator %d does not represent a 2d element",
                      element);
        throw std::runtime_error(msg);
    }
  }

  static ElemConnMetaData create_3d_meta(ET::Element element, index_t** bounds,
                                         index_t** edges, index_t** verts) {
    switch (element) {
      case ET::Element::Tet:
        return create_tet(bounds, edges, verts);
      case ET::Element::Hex:
        return create_hex(bounds, edges, verts);
      case ET::Element::Wedge:
        return create_wedge(bounds, edges, verts);
      case ET::Element::Pyrmd:
        return create_pyrmd(bounds, edges, verts);
      default:
        char msg[256];
        std::snprintf(msg, 256,
                      "element enumerator %d does not represent a 3d element",
                      element);
        throw std::runtime_error(msg);
    }
  }

 private:
  static ElemConnMetaData create_tri(index_t* const* tri_bounds,
                                     const index_t* const* tri_verts) {
    return ElemConnMetaData(ET::TRI_NBOUNDS,       // NBOUNDS
                            0,                     // NEDGES
                            ET::TRI_NVERTS,        // NVERTS
                            ET::TRI_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::TRI_BOUND_VERTS,   // BOUND_VERTS
                            nullptr,               // BOUND_NEDGES
                            nullptr,               // BOUND_EDGES
                            nullptr,               // EDGE_VERTS
                            tri_bounds,            // bounds
                            nullptr,               // edges
                            tri_verts              // verts
    );
  }

  static ElemConnMetaData create_quad(index_t* const* quad_bounds,
                                      const index_t* const* quad_verts) {
    return ElemConnMetaData(ET::QUAD_NBOUNDS,       // NBOUNDS
                            0,                      // NEDGES
                            ET::QUAD_NVERTS,        // NVERTS
                            ET::QUAD_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::QUAD_BOUND_VERTS,   // BOUND_VERTS
                            nullptr,                // BOUND_NEDGES
                            nullptr,                // BOUND_EDGES
                            nullptr,                // EDGE_VERTS
                            quad_bounds,            // bounds
                            nullptr,                // edges
                            quad_verts              // verts
    );
  }

  static ElemConnMetaData create_tet(index_t* const* tet_bounds,
                                     index_t* const* tet_edges,
                                     const index_t* const* tet_verts) {
    return ElemConnMetaData(ET::TET_NBOUNDS,       // NBOUNDS
                            ET::TET_NEDGES,        // NEDGES
                            ET::TET_NVERTS,        // NVERTS
                            ET::TET_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::TET_BOUND_VERTS,   // BOUND_VERTS
                            ET::TET_BOUND_NEDGES,  // BOUND_NEDGES
                            ET::TET_BOUND_EDGES,   // BOUND_EDGES
                            ET::TET_EDGE_VERTS,    // EDGE_VERTS
                            tet_bounds,            // bounds
                            tet_edges,             // edges
                            tet_verts              // verts
    );
  }
  static ElemConnMetaData create_hex(index_t* const* hex_bounds,
                                     index_t* const* hex_edges,
                                     const index_t* const* hex_verts) {
    return ElemConnMetaData(ET::HEX_NBOUNDS,       // NBOUNDS
                            ET::HEX_NEDGES,        // NEDGES
                            ET::HEX_NVERTS,        // NVERTS
                            ET::HEX_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::HEX_BOUND_VERTS,   // BOUND_VERTS
                            ET::HEX_BOUND_NEDGES,  // BOUND_NEDGES
                            ET::HEX_BOUND_EDGES,   // BOUND_EDGES
                            ET::HEX_EDGE_VERTS,    // EDGE_VERTS
                            hex_bounds,            // bounds
                            hex_edges,             // edges
                            hex_verts              // verts
    );
  }
  static ElemConnMetaData create_wedge(index_t* const* wedge_bounds,
                                       index_t* const* wedge_edges,
                                       const index_t* const* wedge_verts) {
    return ElemConnMetaData(ET::WEDGE_NBOUNDS,       // NBOUNDS
                            ET::WEDGE_NEDGES,        // NEDGES
                            ET::WEDGE_NVERTS,        // NVERTS
                            ET::WEDGE_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::WEDGE_BOUND_VERTS,   // BOUND_VERTS
                            ET::WEDGE_BOUND_NEDGES,  // BOUND_NEDGES
                            ET::WEDGE_BOUND_EDGES,   // BOUND_EDGES
                            ET::WEDGE_EDGE_VERTS,    // EDGE_VERTS
                            wedge_bounds,            // bounds
                            wedge_edges,             // edges
                            wedge_verts              // verts
    );
  }
  static ElemConnMetaData create_pyrmd(index_t* const* pyrmd_bounds,
                                       index_t* const* pyrmd_edges,
                                       const index_t* const* pyrmd_verts) {
    return ElemConnMetaData(ET::PYRMD_NBOUNDS,       // NBOUNDS
                            ET::PYRMD_NEDGES,        // NEDGES
                            ET::PYRMD_NVERTS,        // NVERTS
                            ET::PYRMD_BOUND_NVERTS,  // BOUND_NVERTS
                            ET::PYRMD_BOUND_VERTS,   // BOUND_VERTS
                            ET::PYRMD_BOUND_NEDGES,  // BOUND_NEDGES
                            ET::PYRMD_BOUND_EDGES,   // BOUND_EDGES
                            ET::PYRMD_EDGE_VERTS,    // EDGE_VERTS
                            pyrmd_bounds,            // bounds
                            pyrmd_edges,             // edges
                            pyrmd_verts              // verts
    );
  }
};

class MeshConnectivityBase {
 public:
  using ET = ElementTypes;
  static constexpr index_t NO_LABEL = MAX_INDEX;

  // Get the right meta data based on element index
  virtual inline const ElemConnMetaData& get_local_elem_and_meta(
      index_t& elem) = 0;

  // Get global element counts
  index_t get_num_elements() { return nelems; }
  index_t get_num_bounds() { return nbounds; }
  index_t get_num_edges() { return nedges; }
  index_t get_num_verts() { return nverts; }

  // Get quantities associated to an element
  index_t get_element_verts(index_t elem, const index_t* verts[]);
  index_t get_element_bounds(index_t elem, const index_t* bounds[]);
  index_t get_element_bound_verts(index_t elem, index_t b, index_t verts[]);
  index_t get_element_global_bound_verts(index_t elem, index_t local_bound,
                                         index_t verts[]);
  index_t get_element_bound_vert_indices(index_t elem, index_t b,
                                         index_t verts[]);

  // Only 3D elements have edges
  index_t get_element_bound_edges(index_t elem, index_t b, index_t edges[]);
  index_t get_element_bound_edge_indices(index_t elem, index_t b,
                                         index_t edges[]);
  index_t get_element_edges(index_t elem, const index_t* edges[]);
  void get_element_edge_verts(index_t elem, index_t edge, index_t verts[]);

  // Get the adjacent elements from the vert index
  index_t get_adjacent_elements_from_vert(const index_t vert,
                                          const index_t* elems[]);

  //  Check for equality between two bounds/edges
  bool global_bound_equality(index_t na, const index_t a[], index_t nb,
                             const index_t b[]);
  bool global_edge_equality(const index_t a[], const index_t b[]);

  // Get elements associated with the given bound
  bool get_bound_elements(index_t bound, index_t* e1, index_t* e2);

  // Label the verts, edges and bounds that touch the list of vertices
  template <typename IdxType>
  void get_labels_from_verts(const index_t nv, const IdxType verts[],
                             index_t vert_labels[], index_t edge_labels[],
                             index_t bound_labels[]);

  // Get boundary bounds and labels
  index_t get_boundary_bounds(const index_t* boundary_bounds_[],
                              const index_t* boundary_labels_[]);

  // Return the number of boundary labels
  index_t get_num_boundary_labels() { return num_boundary_labels; }

  //  Count up the number of boundary bounds with the given label
  index_t get_num_boundary_bounds_with_label(const index_t label);

  // Add a boundary label from the vertices
  template <typename IdxType>
  index_t add_boundary_label_from_verts(index_t nv, const IdxType vert_list[]);

 protected:
  MeshConnectivityBase(index_t nverts, index_t nelems);
  ~MeshConnectivityBase();

  // Initialize all derived connectivity data and set up boundary, allocation
  // must be performed explicitly before calling this function
  void initialize();

  // Get a non-constant entities
  index_t get_element_bounds(index_t elem, index_t* bounds[]);
  index_t get_element_edges(index_t elem, index_t* edges[]);

  // Inputs: number of vertices, elements, bounds and edges
  index_t nelems, nbounds, nedges, nverts;

  // Number of bounds of different types
  index_t nline_bounds;               // for 2D mesh
  index_t ntri_bounds, nquad_bounds;  // for 3D mesh

  // Global data for pointing from verts to elements
  index_t* vert_element_ptr;
  index_t* vert_elements;

  // Bound -> element connectivity
  index_t* line_bound_elements;                       // for 2d elements
  index_t *tri_bound_elements, *quad_bound_elements;  // for 3d elements

  // Boundary labels for bounds that touch a boundary
  index_t num_boundary_labels;
  index_t num_boundary_bounds;
  index_t* boundary_labels;
  index_t* boundary_bounds;

 private:
  void init_vert_element_data();
  void init_bound_data();
  void init_edge_data();
};

// Mesh connectivity class for 2D meshes composed of triangle and quadrilateral
// elements
class MeshConnectivity2D final : public MeshConnectivityBase {
 public:
  template <typename I>
  MeshConnectivity2D(I nverts, I ntri, I* tri, I nquad, I* quad);
  ~MeshConnectivity2D();

  // Shift elem to get local index within its element type (tri or quad)
  // and return the meta data for this element type
  const inline ElemConnMetaData& get_local_elem_and_meta(index_t& elem);

 private:
  // Input counts of the triangle and quadrilateral elements
  index_t ntri, nquad;

  // Element -> vert connectivity
  index_t* tri_verts;
  index_t* quad_verts;

  // Element -> bound connectivity
  index_t* tri_bounds;
  index_t* quad_bounds;

  // Element meta data
  ElemConnMetaData meta_tri, meta_quad;
};

// Mesh connecivity class for 3D meshes composed of tetrahedral,
// hexahedral, wedge and pyramid elements
class MeshConnectivity3D final : public MeshConnectivityBase {
 public:
  template <typename I>
  MeshConnectivity3D(I nverts, I ntets, I* tets, I nhex, I* hex, I nwedge,
                     I* wedge, I npyrmd, I* pyrmd);
  ~MeshConnectivity3D();

  // Shift elem to get local index within its element type (tet, hex, etc.)
  // and return the meta data for this element type
  const inline ElemConnMetaData& get_local_elem_and_meta(index_t& elem);

 private:
  // Input counts of the tet, hex, wedge and pyramid elements
  index_t ntets, nhex, nwedge, npyrmd;

  // Element -> vert connectivity
  index_t* tet_verts;
  index_t* hex_verts;
  index_t* wedge_verts;
  index_t* pyrmd_verts;

  // Element -> bound connectivity
  index_t* tet_bounds;
  index_t* hex_bounds;
  index_t* wedge_bounds;
  index_t* pyrmd_bounds;

  // Element -> edge connectivity
  index_t* tet_edges;
  index_t* hex_edges;
  index_t* wedge_edges;
  index_t* pyrmd_edges;

  // Element meta data
  ElemConnMetaData meta_tet, meta_hex, meta_wedge, meta_pyrmd;
};

// This class provides a mapping between any boundary labels and
// associated boundary conditions
class DirichletBCInfo {
 public:
  static const index_t MAX_FIXED_FIELDS = INDEX_NBITS - 1;

  DirichletBCInfo() = default;

  //  Add a boundary condition
  void add_boundary_condition(index_t label, index_t basis = 0);
  void add_boundary_condition(index_t label, index_t basis, index_t nfixed,
                              const index_t fixed[]);

  // Search if a boundary condition is active for the specified label, basis
  // pair and dof
  bool active_for_label(index_t label);
  bool active_for_basis(index_t label, index_t basis);
  bool active_for_dof(index_t label, index_t basis, index_t dof);

 private:
  std::vector<std::tuple<index_t, index_t, index_t>> data;
};

/**
 * @brief ElementMesh base class
 *
 */
class ElementMeshBase {
 public:
  virtual ~ElementMeshBase() {}

  // Add all entries to the Jacobian matrix given a matrix with the specified
  // block size
  virtual void add_matrix_pairs(
      const index_t block_size,
      std::set<std::pair<index_t, index_t>>& pairs) const = 0;

  // Get the number of degrees of freedom
  virtual index_t get_num_elements() const = 0;
  virtual index_t get_num_dof() const = 0;

  static void create_block_csr(std::set<std::pair<index_t, index_t>>& pairs,
                               index_t& nrows, std::vector<index_t>& rowp,
                               std::vector<index_t>& cols) {
    nrows = 0;
    typename std::set<std::pair<index_t, index_t>>::iterator it;
    for (it = pairs.begin(); it != pairs.end(); it++) {
      if (it->first > nrows) {
        nrows = it->first;
      }
    }
    nrows++;

    // Find the number of nodes referenced by other nodes
    rowp.resize(nrows + 1);
    std::fill(rowp.begin(), rowp.end(), 0);

    for (it = pairs.begin(); it != pairs.end(); it++) {
      rowp[it->first + 1] += 1;
    }

    // Set the pointer into the rows
    rowp[0] = 0;
    for (index_t i = 0; i < nrows; i++) {
      rowp[i + 1] += rowp[i];
    }

    index_t nnz = rowp[nrows];
    cols.resize(nnz);

    for (it = pairs.begin(); it != pairs.end(); it++) {
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
};

// ElementMesh - Map from an element to the global to element local degrees
// of freedom
template <class Basis>
class ElementMesh : public ElementMeshBase {
 public:
  using ET = ElementTypes;
  static constexpr index_t dim = Basis::dim;

  // Number of degrees of freedom for each element
  static const index_t ndof_per_element = Basis::ndof;

  // Constructors
  ElementMesh(MeshConnectivityBase& conn);
  template <class InteriorBasis>
  ElementMesh(const index_t label, MeshConnectivityBase& conn,
              ElementMesh<InteriorBasis>& mesh);
  template <class HOrderBasis>
  ElementMesh(ElementMesh<HOrderBasis>& mesh);

  index_t get_num_elements() const { return nelems; }
  index_t get_num_dof() const { return num_dof; }
  index_t get_num_cumulative_dof(index_t basis) const {
    return num_dof_offset[basis];
  }

  template <index_t basis>
  int get_global_dof_sign(index_t elem, index_t index);
  template <index_t basis>
  index_t get_global_dof(index_t elem, index_t index);

  // Get the degrees of freedom associated with this element
  KOKKOS_FUNCTION void get_element_dof(const index_t elem,
                                       const index_t* dof[]) {
    *dof = &element_dof[ndof_per_element * elem];
  }

  // Get the signs associated with the degrees of freedom
  void get_element_signs(const index_t elem, const int* signs[]) {
    *signs = &element_sign[ndof_per_element * elem];
  }

  // Add (i, j) pairs from the element mesh
  void add_matrix_pairs(const index_t block_size,
                        std::set<std::pair<index_t, index_t>>& pairs) const;

 private:
  index_t nelems;                         // Total number of elements
  index_t num_dof;                        // Total number of degrees of freedom
  index_t num_dof_offset[Basis::nbasis];  // Cumulative number of degrees of
                                          // freedom each basis

  // Store the degrees of freedom for each element and the element sign
  index_t* element_dof;
  int* element_sign;
};

/*
  Base class for Dirichlet BCs
*/
template <typename T>
class DirichletBase {
 public:
  virtual ~DirichletBase() {}
  virtual index_t get_bcs(const index_t* bcs[], const T* values[]) const = 0;
};

/**
 * @brief Dirichlet boundary conditions for a given basis
 *
 * @tparam Basis The FEBasis type
 */
template <typename T, class Basis>
class DirichletBasis : public DirichletBase<T> {
 public:
  // Use the definitions from the element types
  using ET = ElementTypes;
  static constexpr index_t dim = Basis::dim;

  DirichletBasis(MeshConnectivityBase& conn, ElementMesh<Basis>& mesh,
                 DirichletBCInfo& bcinfo, T value = 0.0);
  ~DirichletBasis() {
    delete dof;
    delete[] vals;
  }

  index_t get_bcs(const index_t* bcs[], const T* values[]) const {
    if (bcs) {
      *bcs = dof;
    }
    if (values) {
      *values = vals;
    }
    return ndof;
  }

 private:
  index_t ndof;
  index_t* dof;
  T* vals;
};

/*
  A collection of Dirichlet BCs from different sources
*/
template <typename T>
class DirichletBCs {
 public:
  typedef std::shared_ptr<DirichletBase<T>> BCPtr;

  DirichletBCs() {}

  void add_bcs(BCPtr bc) {
    const index_t* dof;
    const T* vals;
    index_t count = bc->get_bcs(&dof, &vals);

    indices.reserve(indices.size() + count);
    indices.insert(indices.end(), dof, dof + count);

    values.reserve(values.size() + count);
    values.insert(values.end(), vals, vals + count);
  }

  template <class VecType>
  void zero_bcs(VecType& vec) {
    const index_t size = indices.size();
    for (index_t i = 0; i < size; i++) {
      vec[indices[i]] = T(0.0);
    }
  }

  template <class VecType>
  void set_bcs(VecType& vec) {
    const index_t size = indices.size();
    for (index_t i = 0; i < size; i++) {
      vec[indices[i]] = values[i];
    }
  }

  index_t get_bcs(const index_t* array[]) const {
    *array = indices.data();
    return indices.size();
  }
  index_t get_bcs(const index_t* array[], const T* vals[]) const {
    *array = indices.data();
    *vals = values.data();
    return indices.size();
  }

 private:
  std::vector<index_t> indices;
  std::vector<T> values;
};

}  // namespace A2D

#include "multiphysics/femesh-inl.h"

#endif  // A2D_FE_MESH_H
