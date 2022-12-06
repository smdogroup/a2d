#ifndef A2D_VTK_H
#define A2D_VTK_H

#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "a2dobjs.h"
#include "array.h"
#include "multiphysics/feelementtypes.h"

namespace A2D {

// Linear and nonlinear cell types in VTK
struct VTKID {
  static constexpr int VERTEX = 1;
  static constexpr int POLY_VERTEX = 2;
  static constexpr int LINE = 3;
  static constexpr int POLY_LINE = 4;
  static constexpr int TRIANGLE = 5;
  static constexpr int TRIANGLE_STRIP = 6;
  static constexpr int POLYGON = 7;
  static constexpr int PIXEL = 8;
  static constexpr int QUAD = 9;
  static constexpr int TETRA = 10;
  static constexpr int VOXEL = 11;
  static constexpr int HEXAHEDRON = 12;
  static constexpr int WEDGE = 13;
  static constexpr int PYRAMID = 14;
  static constexpr int QUADRATIC_EDGE = 21;
  static constexpr int QUADRATIC_TRIANGLE = 22;
  static constexpr int QUADRATIC_QUAD = 23;
  static constexpr int QUADRATIC_TETRA = 24;
  static constexpr int QUADRATIC_HEXAAHEDRON = 25;
};

struct VTK_NVERTS {
  static constexpr int VERTEX = 1;
  static constexpr int LINE = 2;
  static constexpr int TRIANGLE = 3;
  static constexpr int PIXEL = 4;
  static constexpr int QUAD = 4;
  static constexpr int TETRA = 4;
  static constexpr int VOXEL = 8;
  static constexpr int HEXAHEDRON = 8;
  static constexpr int WEDGE = 6;
  static constexpr int PYRAMID = 5;
  static constexpr int QUADRATIC_EDGE = 3;
  static constexpr int QUADRATIC_TRIANGLE = 6;
  static constexpr int QUADRATIC_QUAD = 8;
  static constexpr int QUADRATIC_TETRA = 10;
  static constexpr int QUADRATIC_HEXAAHEDRON = 20;
};

void write_real_val(std::FILE* fp, double val) {
  std::fprintf(fp, "%-20.15f", val);
}

void write_real_val(std::FILE* fp, A2D_complex_t<double> val) {
  std::fprintf(fp, "%-20.15f", val.real());
}

/**
 * @brief A VTK writer.
 *
 * @tparam ConnArray connectivity type, shape: (nelems, nnodes_per_elem)
 * @tparam NodeArray nodal location type, shape: (nnodes, spatial_dim)
 */
template <class ConnArray, class NodeArray>
class ToVTK {
 public:
  ToVTK(const ConnArray& conn, const NodeArray& X,
        const int _vtk_elem_type = -1,
        const std::string vtk_name = "result.vtk")
      : conn(conn), X(X), vtk_elem_type(_vtk_elem_type) {
    // Open file and destroy old contents
    fp = std::fopen(vtk_name.c_str(), "w+");

    // Get dimensions
    nnodes = X.extent(0);
    spatial_dim = X.extent(1);
    nelems = conn.extent(0);
    nnodes_per_elem = conn.extent(1);

    if (spatial_dim != 2 and spatial_dim != 3) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "Invalid spatial_dim, got %d, expect 2 or 3", spatial_dim);
      throw std::runtime_error(msg);
    }

    // If not provided, infer vtk type id from nnodes_per_elem
    if (vtk_elem_type == -1) {
      if (nnodes_per_elem == 4) {
        vtk_elem_type = VTKID::QUAD;
      } else if (nnodes_per_elem == 8) {
        vtk_elem_type = VTKID::HEXAHEDRON;
      }
    }

    vtk_has_sol_header = false;
  }

  ~ToVTK() {
    // Close file
    std::fclose(fp);
  }

  void write_mesh() {
    // Write header
    std::fprintf(fp, "# vtk DataFile Version 3.0\n");
    std::fprintf(fp, "my example\n");
    std::fprintf(fp, "ASCII\n");
    std::fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write nodes
    std::fprintf(fp, "POINTS %d double\n", nnodes);
    for (index_t i = 0; i < nnodes; i++) {
      if (spatial_dim == 2) {
        write_real_val(fp, X(i, 0));
        write_real_val(fp, X(i, 1));
        write_real_val(fp, 0.0);
        std::fprintf(fp, "\n");
      } else {
        write_real_val(fp, X(i, 0));
        write_real_val(fp, X(i, 1));
        write_real_val(fp, X(i, 2));
        std::fprintf(fp, "\n");
      }
    }

    // Write connectivity
    std::fprintf(fp, "CELLS %d %d\n", nelems, nelems * (1 + nnodes_per_elem));
    for (index_t i = 0; i < nelems; i++) {
      std::fprintf(fp, "%d ", nnodes_per_elem);
      for (index_t j = 0; j < nnodes_per_elem; j++) {
        std::fprintf(fp, "%d ", conn(i, j));
      }
      std::fprintf(fp, "\n");
    }

    // Write cell type
    std::fprintf(fp, "CELL_TYPES %d\n", nelems);
    for (index_t i = 0; i < nelems; i++) {
      std::fprintf(fp, "%d\n", vtk_elem_type);
    }
  }

  /**
   * @brief Write nodal solution to vtk.
   *
   * @tparam SolArray solution array type, shape: (nnodes, x)
   * @param sol_name solution name
   * @param sol_vec solution vector
   * @param second_dim which slice in the second dimension of sol_vec to write
   */
  template <class SolVector>
  void write_sol(const char sol_name[], const SolVector& sol_vec,
                 const index_t second_dim = 0) {
    // Check input
    if (sol_vec.extent(0) != nnodes) {
      char msg[256];
      std::snprintf(
          msg, sizeof(msg),
          "First dimension of sol_vec (%d) does not match nnodes (%d)",
          (int)sol_vec.extent(0), nnodes);
      throw std::runtime_error(msg);
    }

    // Write header
    if (!vtk_has_sol_header) {
      std::fprintf(fp, "POINT_DATA %d\n", nnodes);
      vtk_has_sol_header = true;
    }
    std::fprintf(fp, "SCALARS %s double 1\n", sol_name);
    std::fprintf(fp, "LOOKUP_TABLE default\n");

    // Write data
    for (index_t i = 0; i < nnodes; i++) {
      write_real_val(fp, sol_vec(i));  // , second_dim));
      std::fprintf(fp, "\n");
    }
  }

 private:
  const ConnArray& conn;
  const NodeArray& X;
  index_t nnodes;
  index_t spatial_dim;
  index_t nelems;
  index_t nnodes_per_elem;
  int vtk_elem_type;  // VTK type id
  std::FILE* fp;
  bool vtk_has_sol_header;
};

/**
 * @brief Write the new MeshConnectivity3D data to vtk
 */
class ToVTK3D {
 public:
  using ET = ElementTypes;
  ToVTK3D(index_t nverts, index_t ntets, index_t* tets, index_t nhex,
          index_t* hex, index_t nwedge, index_t* wedge, index_t npyrmd,
          index_t* pyrmd, double* xloc, const char vtk_name[] = "result.vtk") {
    // Open file and destroy old contents
    fp = std::fopen(vtk_name, "w");

    // Write header
    std::fprintf(fp, "# vtk DataFile Version 3.0\n");
    std::fprintf(fp, "my example\n");
    std::fprintf(fp, "ASCII\n");
    std::fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write vertices
    std::fprintf(fp, "POINTS %d double\n", nverts);
    for (index_t ivert = 0; ivert < nverts; ivert++) {
      for (index_t dim = 0; dim < 3; dim++) {
        write_real_val(fp, xloc[3 * ivert + dim]);
      }
      std::fprintf(fp, "\n");
    }

    // Write connectivities for tet, hex, wedge and pyrmd
    index_t ncells = ntets + nhex + nwedge + npyrmd;
    index_t data_size =
        (ET::TET_VERTS + 1) * ntets + (ET::HEX_VERTS + 1) * nhex +
        (ET::WEDGE_VERTS + 1) * nwedge + (ET::PYRMD_VERTS + 1) * npyrmd;
    std::fprintf(fp, "CELLS %d %d\n", ncells, data_size);
    write_cell_conn(fp, ntets, tets, ET::TET_VERTS);
    write_cell_conn(fp, nhex, hex, ET::HEX_VERTS);
    write_cell_conn(fp, nwedge, wedge, ET::WEDGE_VERTS);
    write_cell_conn(fp, npyrmd, pyrmd, ET::PYRMD_VERTS);

    // Write cell type
    std::fprintf(fp, "CELL_TYPES %d\n", ncells);
    write_cell_type(fp, ntets, VTKID::TETRA);
    write_cell_type(fp, nhex, VTKID::HEXAHEDRON);
    write_cell_type(fp, nwedge, VTKID::WEDGE);
    write_cell_type(fp, npyrmd, VTKID::PYRAMID);
  }

  ~ToVTK3D() {
    // Close file
    std::fclose(fp);
  }

 private:
  void write_cell_conn(std::FILE* fp, index_t ncells, index_t* cells,
                       index_t nverts_per_cell) {
    for (index_t icell = 0; icell < ncells; icell++) {
      std::fprintf(fp, "%d ", nverts_per_cell);
      for (index_t ivert = 0; ivert < nverts_per_cell; ivert++) {
        std::fprintf(fp, "%d ", cells[nverts_per_cell * icell + ivert]);
      }
      std::fprintf(fp, "\n");
    }
    return;
  }

  void write_cell_type(std::FILE* fp, index_t ncells, index_t label) {
    for (index_t i = 0; i < ncells; i++) {
      std::fprintf(fp, "%d\n", label);
    }
  }

  std::FILE* fp;
};

/**
 * @brief Given a VTK, extract information that MeshConnectivity3D needs.
 */
template <typename I, typename T>
class ReadVTK3D {
 public:
  ReadVTK3D(const std::string& vtk_name)
      : vtk_name(vtk_name),
        ntri(0),
        nquad(0),
        nverts(0),
        ntets(0),
        nhex(0),
        nwedge(0),
        npyramid(0),
        xyz_lower{1e20, 1e20, 1e20},
        xyz_upper{-1e20, -1e20, -1e20} {
    // If file doesn't exist, exit
    if (!std::filesystem::exists(vtk_name)) {
      char msg[256];
      std::snprintf(msg, sizeof(msg), "file %s does not exists!",
                    vtk_name.c_str());
      throw std::runtime_error(msg);
    }

    // Scan through the vtk file and save the information.
    // Note that vtk might contain element types that we don't support here,
    // we'll trim out those later.
    I nelems = 0;
    std::string dummy;  // store string to be discarded

    std::ifstream ifs(vtk_name);

    // Get number of points from vtk
    goto_line_with_str(ifs, "POINTS") >> dummy >> nverts;

    // Allocate and populate nodal location array, and record domain bounds
    T coord;
    Xloc.reserve(nverts * SPATIAL_DIM);
    std::string line;
    for (int i = 0; i < nverts; i++) {
      std::getline(ifs, line);
      std::istringstream iss(line);
      for (int j = 0; j < SPATIAL_DIM; j++) {
        iss >> coord;
        Xloc[SPATIAL_DIM * i + j] = coord;
        if (coord > xyz_upper[j]) {
          xyz_upper[j] = coord;
        }
        if (coord < xyz_lower[j]) {
          xyz_lower[j] = coord;
        }
      }
    }

    // Get number of all types of elements (including those that aren't
    // supported)
    goto_line_with_str(ifs, "CELLS") >> dummy >> nelems;

    // Populate connectivities for supported element types only.
    // Supported types are tet, hex, wedge, pyramid
    I nelems_celltype = 0, nelems_cellid = 0;

    std::ifstream ifs_celltype(vtk_name);
    std::ifstream ifs_cellid(vtk_name);
    goto_line_with_str(ifs_celltype, "CELL_TYPES") >> dummy >> nelems_celltype;
    goto_line_with_str(ifs_cellid, "CELL_DATA") >> dummy >> nelems_cellid;
    goto_line_with_str(ifs_cellid, "CellEntityIds");
    goto_line_with_str(ifs_cellid, "LOOKUP_TABLE");

    assert(nelems == nelems_celltype);
    assert(nelems == nelems_cellid);

    std::string line_celltype;
    std::string line_cellid;

    int cell_type = 0;
    int cell_id = 0;

    // Loop over all elements
    for (int i = 0; i < nelems; i++) {
      // Get element type
      std::getline(ifs_celltype, line_celltype);
      std::istringstream(line_celltype) >> cell_type;

      // Get element entity id
      std::getline(ifs_cellid, line_cellid);
      std::istringstream(line_cellid) >> cell_id;

      // Get element connectivity
      std::getline(ifs, line);

      // Populate cell id vectors and connectivity
      switch (cell_type) {
        case VTKID::TRIANGLE:
          id_tri.push_back(cell_id);
          insert_verts_to_conn(line, conn_tri);
          ntri++;
          break;
        case VTKID::QUAD:
          id_quad.push_back(cell_id);
          insert_verts_to_conn(line, conn_quad);
          nquad++;
          break;
        case VTKID::TETRA:
          id_tet.push_back(cell_id);
          insert_verts_to_conn(line, conn_tet);
          ntets++;
          break;
        case VTKID::HEXAHEDRON:
          id_hex.push_back(cell_id);
          insert_verts_to_conn(line, conn_hex);
          nhex++;
          break;
        case VTKID::WEDGE:
          id_wedge.push_back(cell_id);
          insert_verts_to_conn(line, conn_wedge);
          nwedge++;
          break;
        case VTKID::PYRAMID:
          id_pyramid.push_back(cell_id);
          insert_verts_to_conn(line, conn_pyramid);
          npyramid++;
          break;
      }
    }
  }

  /**
   * @brief Get vertex indices given cell entity id.
   *
   * In the given VTK, each cell is assigned an entity id, which is used to
   * group cells in the same physical group (e.g. boundary, etc.). Such id
   * labels are associated to cells, hence they need to be converted to
   * vertex indices for A2D to use.
   */
  std::vector<I> get_verts_given_cell_entity_id(const std::vector<int> id) {
    std::set<I> verts;
    insert_verts_to_set(VTK_NVERTS::TRIANGLE, id, verts, id_tri, conn_tri);
    insert_verts_to_set(VTK_NVERTS::QUAD, id, verts, id_quad, conn_quad);
    insert_verts_to_set(VTK_NVERTS::TETRA, id, verts, id_tet, conn_tet);
    insert_verts_to_set(VTK_NVERTS::HEXAHEDRON, id, verts, id_hex, conn_hex);
    insert_verts_to_set(VTK_NVERTS::WEDGE, id, verts, id_wedge, conn_wedge);
    insert_verts_to_set(VTK_NVERTS::PYRAMID, id, verts, id_pyramid,
                        conn_pyramid);
    std::shared_ptr<I[]> data(new I[verts.size()]);
    return std::vector<I>(verts.begin(), verts.end());
  }

  /**
   * @brief Get the vertex indices of those lie within the given box.
   *
   * This is handy to specify bc nodes directly in A2D.
   */
  std::vector<I> get_verts_within_box(double xmin, double xmax, double ymin,
                                      double ymax, double zmin, double zmax,
                                      double tol = 1e-6) {
    std::vector<I> verts;
    double x, y, z;
    for (I i = 0; i < nverts; i++) {
      x = Xloc[3 * i];
      y = Xloc[3 * i + 1];
      z = Xloc[3 * i + 2];
      if (x > xmin - tol && x < xmax + tol) {
        if (y > ymin - tol && y < ymax + tol) {
          if (z > zmin - tol && z < zmax + tol) {
            verts.push_back(i);
          }
        }
      }
    }
    return verts;
  }

  // Get number of vertices for each element type
  I get_nverts() { return nverts; }
  I get_ntri() { return ntri; }
  I get_nquad() { return nquad; }
  I get_ntets() { return ntets; }
  I get_nhex() { return nhex; }
  I get_nwedge() { return nwedge; }
  I get_npyrmd() { return npyramid; }

  // Get connectivity data for each element type
  I* get_tri() { return conn_tri.data(); }
  I* get_quad() { return conn_quad.data(); }
  I* get_tets() { return conn_tet.data(); }
  I* get_hex() { return conn_hex.data(); }
  I* get_wedge() { return conn_wedge.data(); }
  I* get_pyrmd() { return conn_pyramid.data(); }

  // Get nodal locations of vertices
  T* get_Xloc() { return Xloc.data(); }
  void get_bounds(T* lower, T* upper) {
    for (int i = 0; i < SPATIAL_DIM; i++) {
      lower[i] = xyz_lower[i];
      upper[i] = xyz_upper[i];
    }
  }

 private:
  /**
   * @brief Insert vertex indices of a given element type to a given set
   * based on the cell id.
   *
   * @param id cell ids
   * @param nverts_per_cell number of vertices per cell of this type
   * @param verts the set
   * @param id_vec id vector of this cell type
   * @param conn_vec connectivity vector of this cell type
   */
  void insert_verts_to_set(const int nverts_per_cell, const std::vector<int> id,
                           std::set<I>& verts, std::vector<int>& id_vec,
                           std::vector<I>& conn_vec) {
    for (int i = 0; i < id_vec.size(); i++) {
      if (std::count(id.begin(), id.end(), id_vec[i])) {
        for (int j = 0; j < nverts_per_cell; j++) {
          verts.insert(conn_vec[i * nverts_per_cell + j]);
        }
      }
    }
  }

  /**
   * @brief Insert the vertices of a given cell to corresponding
   * connectivity
   */
  void insert_verts_to_conn(std::string& conn_line, std::vector<I>& conn) {
    int elem_nverts;
    std::istringstream conn_iss = std::istringstream(conn_line);
    conn_iss >> elem_nverts;  // get number of vertices
    I vert_idx;
    for (int i = 0; i < elem_nverts; i++) {
      conn_iss >> vert_idx;
      conn.push_back(vert_idx);
    }
    return;
  }

  /**
   * @brief Find the first line that contains given string.
   *
   * This function search through the input text until hitting the first
   * line that contains the given string. Upon returning, the file handle
   * points to the next line.
   *
   * @param file the file input stream
   * @param str string to be matched
   * @return a std::istringstream type that can be used to extract data from
   *         the matching line. Usage: ret >> x >> y >> z ...
   */
  std::istringstream goto_line_with_str(std::ifstream& file,
                                        const std::string str) {
    std::string line;
    while (std::getline(file, line)) {
      if (line.find(str) != std::string::npos) {
        return std::istringstream(line);
      }
    }
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "\n[%s:%d]file %s does not contain given string %s!",
                  __FILE__, __LINE__, vtk_name.c_str(), str.c_str());
    throw std::runtime_error(msg);
  }

  index_t static constexpr SPATIAL_DIM = 3;
  std::string vtk_name;
  int ntri, nquad;
  int nverts, ntets, nhex, nwedge, npyramid;

  // Connectivity lists for 3d mesh cells
  std::vector<I> conn_tet, conn_hex, conn_wedge, conn_pyramid;
  std::vector<int> id_tet, id_hex, id_wedge, id_pyramid;

  // Connectivity lists for 2d mesh cells
  std::vector<I> conn_tri, conn_quad;
  std::vector<int> id_tri, id_quad;

  // Location coordinates for vertices
  std::vector<T> Xloc;

  // Domain bounds
  T xyz_lower[3], xyz_upper[3];
};

/**
 * @brief Generate a vtk containing a vector field in scattered data format.
 * No mesh contained.
 */
class VectorFieldToVTK {
 public:
  VectorFieldToVTK(index_t nverts, double* xloc, double* vec,
                   const std::string vtk_name = "vector_field.vtk") {
    fp = std::fopen(vtk_name.c_str(), "w");

    // Write header
    std::fprintf(fp, "# vtk DataFile Version 3.0\n");
    std::fprintf(fp, "my example\n");
    std::fprintf(fp, "ASCII\n");
    std::fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write vertices
    std::fprintf(fp, "POINTS %d double\n", nverts);
    for (index_t ivert = 0; ivert < nverts; ivert++) {
      for (index_t dim = 0; dim < 3; dim++) {
        write_real_val(fp, xloc[3 * ivert + dim]);
      }
      std::fprintf(fp, "\n");
    }

    // Write vector field data
    std::fprintf(fp, "POINT_DATA %d \n", nverts);
    std::fprintf(fp, "VECTORS vector double\n");
    for (index_t ivert = 0; ivert < nverts; ivert++) {
      for (index_t dim = 0; dim < 3; dim++) {
        write_real_val(fp, vec[3 * ivert + dim]);
      }
      std::fprintf(fp, "\n");
    }

    // Write coordinates as a second vector
    std::fprintf(fp, "VECTORS coord double\n");
    for (index_t ivert = 0; ivert < nverts; ivert++) {
      for (index_t dim = 0; dim < 3; dim++) {
        write_real_val(fp, xloc[3 * ivert + dim]);
      }
      std::fprintf(fp, "\n");
    }
  }

  ~VectorFieldToVTK() { std::fclose(fp); }

 private:
  std::FILE* fp;
};

/**
 * @brief Read conn and X from VTK.
 *
 * @tparam nnodes_per_elem second dimension of the connectivity array
 * @tparam spatial_dim second dimension of the nodal location array
 * @tparam T X entry data type
 * @tparam I connectivity entry data type
 */
template <int nnodes_per_elem, int spatial_dim, typename T, typename I>
class ReadVTK {
 public:
  static_assert(spatial_dim == 2 or spatial_dim == 3);
  using ConnArray_t = A2D::MultiArrayNew<I* [nnodes_per_elem]>;
  using NodeArray_t = A2D::MultiArrayNew<T* [spatial_dim]>;

  ReadVTK() {}

  /**
   * @brief Copy constructor
   */
  ReadVTK(const ReadVTK& src)
      : conn(src.conn),
        conn_(src.conn_),
        X(src.X),
        nnodes(src.nnodes),
        nelems(src.nelems),
        nelems_all(src.nelems_all),
        domain_lower(src.domain_lower),
        domain_upper(src.domain_upper),
        domain_sizes(src.domain_sizes) {}

  /**
   * @brief Assignment operator
   */
  ReadVTK& operator=(const ReadVTK& src) {
    conn = src.conn;
    conn_ = src.conn_;
    X = src.X;
    nnodes = src.nnodes;
    nelems = src.nelems;
    nelems_all = src.nelems_all;
    for (int i = 0; i < spatial_dim; i++) {
      domain_lower[i] = src.domain_lower[i];
      domain_upper[i] = src.domain_upper[i];
      domain_sizes[i] = src.domain_sizes[i];
    }
    return *this;
  }

  ReadVTK(const std::string& vtk_name) {
    // If file doesn't exist, exit
    if (!std::filesystem::exists(vtk_name)) {
      char msg[256];
      std::snprintf(msg, sizeof(msg), "file %s does not exists!",
                    vtk_name.c_str());
      throw std::runtime_error(msg);
    }

    nnodes = 0;
    nelems_all = 0;
    nelems = 0;

    for (int i = 0; i != spatial_dim; i++) {
      domain_lower[i] = 1e9;
      domain_upper[i] = -1e9;
      domain_sizes[i] = 0.0;
    }

    std::ifstream file(vtk_name);
    std::string line;

    // Get number of points
    while (std::getline(file, line)) {
      if (line.find("POINTS") != std::string::npos) {
        std::istringstream iss(line);
        std::string dummy;
        iss >> dummy >> nnodes;
        break;
      }
    }

    // Allocate nodal location array
    X = NodeArray_t("X", nnodes);

    // Populate nodal location array
    int count = 0;
    T pt[spatial_dim];
    while (std::getline(file, line)) {
      std::istringstream iss(line);

      for (int i = 0; i != spatial_dim; i++) {
        iss >> pt[i];
        if (pt[i] < domain_lower[i]) {
          domain_lower[i] = pt[i];
        }
        if (pt[i] > domain_upper[i]) {
          domain_upper[i] = pt[i];
        }
        X(count, i) = pt[i];
      }
      count++;
      if (count == nnodes) {
        break;
      }
    }

    // Compute domain sizes
    for (int i = 0; i != spatial_dim; i++) {
      domain_sizes[i] = domain_upper[i] - domain_lower[i];
    }

    // Get total number of elements (could have multiple types)
    while (std::getline(file, line)) {
      if (line.find("CELLS") != std::string::npos) {
        std::istringstream iss(line);
        std::string dummy;
        iss >> dummy >> nelems_all;
        break;
      }
    }

    // Allocate connectivity array
    conn_ = ConnArray_t("conn_", nelems_all);

    // Populate connectivity array
    int _nnodes = -1;
    count = 0;
    nelems = 0;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      iss >> _nnodes;
      if (_nnodes == nnodes_per_elem) {
        for (int i = 0; i != nnodes_per_elem; i++) {
          iss >> conn_(nelems, i);
        }
        nelems++;
      }
      count++;
      if (count == nelems_all) {
        break;
      }
    }

    // Trim connectivity
    conn = ConnArray_t("conn", nelems);
    BLAS::copy(conn, Kokkos::subview(conn_, Kokkos::make_pair(0, (int)nelems),
                                     Kokkos::ALL));

    // printf("X[%d]\n", nnodes);
    // printf("lower: %.2f %.2f %.2f\n", domain_lower[1], domain_lower[1],
    //        domain_lower[2]);
    // printf("upper: %.2f %.2f %.2f\n", domain_upper[1], domain_upper[1],
    //        domain_upper[2]);
    // printf("domain: %.2f %.2f %.2f\n", domain_sizes[1], domain_sizes[1],
    //        domain_sizes[2]);
    // for (int i = 0; i != nnodes; i++) {
    //   printf("%.2f %.2f %.2f\n", X(i, 0), X(i, 1), X(i, 2));
    // }

    // printf("conn[%d, 4]\n", nelems);
    // for (int i = 0; i != nelems; i++) {
    //   printf("%d %d %d %d\n", conn(i, 0), conn(i, 1), conn(i, 2), conn(i,
    //   3));
    // }
  }

  template <class XArray>
  void set_X(XArray& _X) {
    for (I i = 0; i != nnodes; i++) {
      for (I j = 0; j != spatial_dim; j++) {
        _X(i, j) = X(i, j);
      }
    }
  }

  template <class ConnArray>
  void set_conn(ConnArray& _conn) {
    for (I i = 0; i != nelems; i++) {
      for (I j = 0; j != nnodes_per_elem; j++) {
        _conn(i, j) = conn(i, j);
      }
    }
  }

  ConnArray_t conn_, conn;
  NodeArray_t X;

  I nnodes, nelems, nelems_all;

  T domain_lower[spatial_dim];  // smallest nodal coordinates
  T domain_upper[spatial_dim];  // largest nodal coordinates
  T domain_sizes[spatial_dim];  // domain length in each direction
};

}  // namespace A2D

#endif  // A2D_VTK_H