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

namespace A2D {

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
        const int _vtk_elem_type = -1, const char vtk_name[] = "result.vtk")
      : conn(conn), X(X), vtk_elem_type(_vtk_elem_type) {
    // Open file and destroy old contents
    fp = std::fopen(vtk_name, "w+");

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
        vtk_elem_type = 9;
      } else if (nnodes_per_elem == 8) {
        vtk_elem_type = 12;
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