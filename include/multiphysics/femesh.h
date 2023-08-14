#ifndef A2D_FE_MESH_H
#define A2D_FE_MESH_H

#include <algorithm>
#include <set>

#include "a2dmatops2d.h"
#include "a2dmatops3d.h"
#include "a2dobjs.h"
#include "multiphysics/feelementtypes.h"
#include "sparse/sparse_matrix.h"
#include "sparse/sparse_symbolic.h"
#include "utils/a2dprofiler.h"

namespace A2D {

// Meta data for a specific type of elements
struct ElemConnMetaData {
  using ET = ElementTypes;

  bool is_valid_element = false;

  // Number of quantities for a single element of this type
  index_t NVERTS = NO_INDEX;
  index_t NEDGES = NO_INDEX;
  index_t NFACES = NO_INDEX;

  // Face -> verts
  const index_t* FACE_NVERTS = nullptr;  // number of verts for each face
  const index_t (*FACE_VERTS)[ET::MAX_FACE_VERTS] =
      nullptr;  // verts of all faces of an element

  // Face -> edges
  const index_t* FACE_NEDGES = nullptr;
  const index_t (*FACE_EDGES)[ET::MAX_FACE_EDGES] = nullptr;

  // Edge -> verts
  const index_t (*EDGE_VERTS)[2] = nullptr;

  // Element -> vert connectivity
  index_t** verts = nullptr;

  // Element -> edge connectivity
  index_t** edges = nullptr;

  // Element -> face connectivity
  index_t** faces = nullptr;
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
  index_t get_num_verts() { return nverts; }
  index_t get_num_faces() { return nfaces; }
  index_t get_num_edges() { return nedges; }

  // Get quantities associated to an element
  index_t get_element_verts(index_t elem, const index_t* verts[]);
  index_t get_element_faces(index_t elem, const index_t* faces[]);
  index_t get_element_face_verts(index_t elem, index_t f, index_t verts[]);
  index_t get_element_face_vert_indices(index_t elem, index_t f,
                                        index_t verts[]);
  index_t get_element_face_edges(index_t elem, index_t f, index_t edges[]);
  index_t get_element_face_edge_indices(index_t elem, index_t f,
                                        index_t edges[]);
  index_t get_element_global_face_verts(index_t elem, index_t local_face,
                                        index_t verts[]);
  index_t get_element_edges(index_t elem, const index_t* edges[]);
  void get_element_edge_verts(index_t elem, index_t edge, index_t verts[]);

  // Get the adjacent elements from the vert index
  index_t get_adjacent_elements_from_vert(const index_t vert,
                                          const index_t* elems[]);

  //  Check for equality between two faces/edges
  bool global_face_equality(index_t na, const index_t a[], index_t nb,
                            const index_t b[]);
  bool global_edge_equality(const index_t a[], const index_t b[]);

  // Get elements associated with the given face
  bool get_face_elements(index_t face, index_t* e1, index_t* e2);

  // Label the verts, edges and faces that touch the list of vertices
  template <typename IdxType>
  void get_labels_from_verts(const index_t nv, const IdxType verts[],
                             index_t vert_labels[], index_t edge_labels[],
                             index_t face_labels[]);

  // Get boundary faces and labels
  index_t get_boundary_faces(const index_t* boundary_faces_[],
                             const index_t* boundary_labels_[]);

  // Return the number of boundary labels
  index_t get_num_boundary_labels() { return num_boundary_labels; }

  //  Count up the number of boundary faces with the given label
  index_t get_num_boundary_faces_with_label(const index_t label);

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
  index_t get_element_faces(index_t elem, index_t* faces[]);
  index_t get_element_edges(index_t elem, index_t* edges[]);

  // Inputs: number of vertices, elements, faces and edges
  index_t nverts, nelems, nfaces, nedges;

  // Derived count info
  index_t nline_faces, ntri_faces, nquad_faces;

  // Global data for pointing from verts to elements
  index_t* vert_element_ptr;
  index_t* vert_elements;

  // Face -> element connectivity
  index_t* line_face_elements;
  index_t* tri_face_elements;
  index_t* quad_face_elements;

  // Boundary labels for faces that touch a boundary
  index_t num_boundary_labels;
  index_t num_boundary_faces;
  index_t* boundary_labels;
  index_t* boundary_faces;

 private:
  void init_vert_element_data();
  void init_face_data();
  void init_edge_data();
};

// Mesh connecivity class for 2D meshes composed of triangle and quadrilateral
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

  // Element -> face connectivity
  index_t* tri_faces;
  index_t* quad_faces;

  // Element -> edge connectivity
  index_t* tri_edges;
  index_t* quad_edges;

  // Element meta data
  ElemConnMetaData meta_tri, meta_quad, meta_none;
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

  // Element -> face connectivity
  index_t* tet_faces;
  index_t* hex_faces;
  index_t* wedge_faces;
  index_t* pyrmd_faces;

  // Element -> edge connectivity
  index_t* tet_edges;
  index_t* hex_edges;
  index_t* wedge_edges;
  index_t* pyrmd_edges;

  // Element meta data
  ElemConnMetaData meta_tet, meta_hex, meta_wedge, meta_pyrmd, meta_none;
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

  // Constructors
  ElementMesh(MeshConnectivityBase& conn);
  template <class InteriorBasis>
  ElementMesh(const index_t label, MeshConnectivityBase& conn,
              ElementMesh<InteriorBasis>& mesh);
  template <class HOrderBasis>
  ElementMesh(ElementMesh<HOrderBasis>& mesh);

  index_t get_num_elements() { return nelems; }
  index_t get_num_dof() { return num_dof; }
  index_t get_num_cumulative_dof(index_t basis) {
    return num_dof_offset[basis];
  }

  template <index_t basis>
  int get_global_dof_sign(index_t elem, index_t index);
  template <index_t basis>
  index_t get_global_dof(index_t elem, index_t index);

  // Get the degrees of freedom associated with this element
  void get_element_dof(const index_t elem, const index_t* dof[]) {
    *dof = &element_dof[ndof_per_element * elem];
  }

  // Get the signs associated with the degrees of freedom
  void get_element_signs(const index_t elem, const int* signs[]) {
    *signs = &element_sign[ndof_per_element * elem];
  }

  template <index_t M, index_t basis_offset = Basis::nbasis>
  void create_block_csr(index_t& nrows, std::vector<index_t>& rowp,
                        std::vector<index_t>& cols);

 private:
  index_t nelems;                         // Total number of elements
  index_t num_dof;                        // Total number of degrees of freedom
  index_t num_dof_offset[Basis::nbasis];  // Cumulative number of degrees of
                                          // freedom each basis

  // Store the degrees of freedom for each element and the element sign
  index_t* element_dof;
  int* element_sign;
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
