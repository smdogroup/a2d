#ifndef A2D_VTK_H
#define A2D_VTK_H

#include <complex>
#include <cstdio>
#include <map>
#include <stdexcept>

namespace A2D {

void write_real_val(std::FILE* fp, double val) {
  std::fprintf(fp, "%-20.15f", val);
}

void write_real_val(std::FILE* fp, std::complex<double> val) {
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
        const char vtk_name[] = "result.vtk")
      : conn(conn), X(X) {
    // Open file and destroy old contents
    fp = std::fopen(vtk_name, "w+");

    // Get dimensions
    nnodes = X.extent(0);
    spatial_dim = X.extent(1);
    nelems = conn.extent(0);
    nnodes_per_elem = conn.extent(1);

    if (spatial_dim != 2 and spatial_dim != 3) {
      char msg[256];
      std::sprintf(msg, "Invalid spatial_dim, got %d, expect 2 or 3",
                   spatial_dim);
      throw std::runtime_error(msg);
    }

    elem_types[4] = 9;   // 4-node quadrilateral element
    elem_types[8] = 12;  // 8-node brick element

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
      std::fprintf(fp, "%d\n", elem_types[nnodes_per_elem]);
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
      std::sprintf(msg,
                   "First dimension of sol_vec (%d) does not match nnodes (%d)",
                   sol_vec.extent(0), nnodes);
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
      write_real_val(fp, sol_vec(i, second_dim));
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
  std::map<int, int> elem_types;  // VTK has a type id for each element type
  std::FILE* fp;
  bool vtk_has_sol_header;
};

}  // namespace A2D

#endif  // A2D_VTK_H