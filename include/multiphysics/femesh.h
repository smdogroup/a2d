#ifndef A2D_FE_MESH_H
#define A2D_FE_MESH_H

#include <algorithm>
#include <set>

#include "a2dobjs.h"
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
                      (int)element);
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
                      (int)element);
        throw std::runtime_error(msg);
    }
  }

 private:
  static ElemConnMetaData create_tri(index_t* const* tri_bounds,
                                     const index_t* const* tri_verts) {
    return ElemConnMetaData(ET::TRI_NBOUNDS,             // NBOUNDS
                            0,                           // NEDGES
                            ET::TRI_NVERTS,              // NVERTS
                            ET::get_tri_bound_nverts(),  // BOUND_NVERTS
                            ET::get_tri_bound_verts(),   // BOUND_VERTS
                            nullptr,                     // BOUND_NEDGES
                            nullptr,                     // BOUND_EDGES
                            nullptr,                     // EDGE_VERTS
                            tri_bounds,                  // bounds
                            nullptr,                     // edges
                            tri_verts                    // verts
    );
  }

  static ElemConnMetaData create_quad(index_t* const* quad_bounds,
                                      const index_t* const* quad_verts) {
    return ElemConnMetaData(ET::QUAD_NBOUNDS,             // NBOUNDS
                            0,                            // NEDGES
                            ET::QUAD_NVERTS,              // NVERTS
                            ET::get_quad_bound_nverts(),  // BOUND_NVERTS
                            ET::get_quad_bound_verts(),   // BOUND_VERTS
                            nullptr,                      // BOUND_NEDGES
                            nullptr,                      // BOUND_EDGES
                            nullptr,                      // EDGE_VERTS
                            quad_bounds,                  // bounds
                            nullptr,                      // edges
                            quad_verts                    // verts
    );
  }

  static ElemConnMetaData create_tet(index_t* const* tet_bounds,
                                     index_t* const* tet_edges,
                                     const index_t* const* tet_verts) {
    return ElemConnMetaData(ET::TET_NBOUNDS,             // NBOUNDS
                            ET::TET_NEDGES,              // NEDGES
                            ET::TET_NVERTS,              // NVERTS
                            ET::get_tet_bound_nverts(),  // BOUND_NVERTS
                            ET::get_tet_bound_verts(),   // BOUND_VERTS
                            ET::get_tet_bound_nedges(),  // BOUND_NEDGES
                            ET::get_tet_bound_edges(),   // BOUND_EDGES
                            ET::get_tet_edge_verts(),    // EDGE_VERTS
                            tet_bounds,                  // bounds
                            tet_edges,                   // edges
                            tet_verts                    // verts
    );
  }
  static ElemConnMetaData create_hex(index_t* const* hex_bounds,
                                     index_t* const* hex_edges,
                                     const index_t* const* hex_verts) {
    return ElemConnMetaData(ET::HEX_NBOUNDS,             // NBOUNDS
                            ET::HEX_NEDGES,              // NEDGES
                            ET::HEX_NVERTS,              // NVERTS
                            ET::get_hex_bound_nverts(),  // BOUND_NVERTS
                            ET::get_hex_bound_verts(),   // BOUND_VERTS
                            ET::get_hex_bound_nedges(),  // BOUND_NEDGES
                            ET::get_hex_bound_edges(),   // BOUND_EDGES
                            ET::get_hex_edge_verts(),    // EDGE_VERTS
                            hex_bounds,                  // bounds
                            hex_edges,                   // edges
                            hex_verts                    // verts
    );
  }
  static ElemConnMetaData create_wedge(index_t* const* wedge_bounds,
                                       index_t* const* wedge_edges,
                                       const index_t* const* wedge_verts) {
    return ElemConnMetaData(ET::WEDGE_NBOUNDS,             // NBOUNDS
                            ET::WEDGE_NEDGES,              // NEDGES
                            ET::WEDGE_NVERTS,              // NVERTS
                            ET::get_wedge_bound_nverts(),  // BOUND_NVERTS
                            ET::get_wedge_bound_verts(),   // BOUND_VERTS
                            ET::get_wedge_bound_nedges(),  // BOUND_NEDGES
                            ET::get_wedge_bound_edges(),   // BOUND_EDGES
                            ET::get_wedge_edge_verts(),    // EDGE_VERTS
                            wedge_bounds,                  // bounds
                            wedge_edges,                   // edges
                            wedge_verts                    // verts
    );
  }
  static ElemConnMetaData create_pyrmd(index_t* const* pyrmd_bounds,
                                       index_t* const* pyrmd_edges,
                                       const index_t* const* pyrmd_verts) {
    return ElemConnMetaData(ET::PYRMD_NBOUNDS,             // NBOUNDS
                            ET::PYRMD_NEDGES,              // NEDGES
                            ET::PYRMD_NVERTS,              // NVERTS
                            ET::get_pyrmd_bound_nverts(),  // BOUND_NVERTS
                            ET::get_pyrmd_bound_verts(),   // BOUND_VERTS
                            ET::get_pyrmd_bound_nedges(),  // BOUND_NEDGES
                            ET::get_pyrmd_bound_edges(),   // BOUND_EDGES
                            ET::get_pyrmd_edge_verts(),    // EDGE_VERTS
                            pyrmd_bounds,                  // bounds
                            pyrmd_edges,                   // edges
                            pyrmd_verts                    // verts
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

// ElementMesh - Map from an element to the global to element local degrees
// of freedom
template <class Basis>
class ElementMesh {
 public:
  using ET = ElementTypes;
  static constexpr index_t dim = Basis::dim;

  // Number of degrees of freedom for each element
  static const index_t ndof_per_element = Basis::ndof;

  // Mult-dimensional array type
  using DofOffsetArray = MultiArrayNew<index_t[Basis::nbasis]>;
  using ElementDofArray = MultiArrayNew<index_t* [ndof_per_element]>;
  using ElementSignArray = MultiArrayNew<int* [ndof_per_element]>;

  // Constructors
  ElementMesh(MeshConnectivityBase& conn);
  template <class InteriorBasis>
  ElementMesh(const index_t label, MeshConnectivityBase& conn,
              ElementMesh<InteriorBasis>& mesh);
  template <class HOrderBasis>
  ElementMesh(const ElementMesh<HOrderBasis>& mesh);

  // Copy constructor, this is needed for parallel dispatch
  // Note, this is also a specialization of the third constructor above
  template <>
  KOKKOS_FUNCTION ElementMesh(const ElementMesh<Basis>& other)
      : nelems(other.nelems),
        num_dof(other.num_dof),
        num_dof_offset(other.num_dof_offset),
        element_dof(other.element_dof),
        element_sign(other.element_sign) {}

  index_t get_num_elements() const { return nelems; }
  index_t get_num_dof() const { return num_dof; }
  index_t get_num_cumulative_dof(index_t basis) const {
    return num_dof_offset[basis];
  }

  // Needed for parallel element execution
  template <index_t basis>
  KOKKOS_FUNCTION int get_global_dof_sign(const index_t elem,
                                          const index_t index) const;
  template <index_t basis>
  KOKKOS_FUNCTION index_t get_global_dof(const index_t elem,
                                         const index_t index) const;

  // Get the degrees of freedom associated with this element
  KOKKOS_FUNCTION auto get_element_dof(const index_t elem) {
    return Kokkos::subview(element_dof, elem, Kokkos::ALL);
  }

  // Get the signs associated with the degrees of freedom
  auto get_element_signs(const index_t elem) {
    return Kokkos::subview(element_sign, elem, Kokkos::ALL);
  }

  template <index_t M, index_t basis_offset = Basis::nbasis>
  void create_block_csr(index_t& nrows, std::vector<index_t>& rowp,
                        std::vector<index_t>& cols);

 private:
  index_t nelems;                 // Total number of elements
  index_t num_dof;                // Total number of degrees of freedom
  DofOffsetArray num_dof_offset;  // Cumulative number of degrees of
                                  // freedom each basis

  // Store the degrees of freedom for each element and the element sign
  ElementDofArray element_dof;
  ElementSignArray element_sign;
};

template <class Basis>
class DirichletBCs {
 public:
  // Use the definitions from the element types
  using ET = ElementTypes;
  static constexpr index_t dim = Basis::dim;

  DirichletBCs(MeshConnectivityBase& conn, ElementMesh<Basis>& mesh,
               DirichletBCInfo& bcinfo);

  index_t get_bcs(const index_t* bcs[]) const {
    if (bcs) {
      *bcs = dof;
    }
    return ndof;
  }

  template <class VecType, typename T>
  void set_bcs(VecType& vec, T value) {
    for (index_t i = 0; i < ndof; i++) {
      vec[dof[i]] = value;
    }
  }

 private:
  index_t ndof;
  index_t* dof;
};

}  // namespace A2D

#include "multiphysics/femesh-inl.h"

#endif  // A2D_FE_MESH_H
