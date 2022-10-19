/**
 * @file a2dmesh.h
 * @author Yicong Fu (fuyicong1996@gmail.com)
 * @brief A collection of basic meshes for demonstrative purpose.
 * @version 0.1
 * @date 2022-07-13
 */

#ifndef A2D_MESH_H
#define A2D_MESH_H

#include <string>
#include <vector>

#include "a2dprofiler.h"
#include "a2dvtk.h"
#include "array.h"

namespace A2D {

/**
 * @brief A structured mesh generator for 2-dimensional rectangular domain.
 */
class MesherRect2D {
 public:
  /**
   * @brief Constructor.
   *
   * @param nx number of elements along x direction
   * @param ny number of elements along y direction
   * @param lx length along x direction
   * @param ly length along y direction
   */
  MesherRect2D(int nx, int ny, double lx, double ly)
      : nx(nx), ny(ny), lx(lx), ly(ly) {}

  /**
   * @brief populate X and conn
   *
   * @param X nodal location multiarray
   * @param conn connectivity multiarray
   * @param lx length along x direction
   * @param ly length along y direction
   */
  template <class ConnArray, class XArray>
  void set_X_conn(XArray& X, ConnArray& conn) {
    // Set X
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * j;

        X(node, 0) = lx * i / nx;
        X(node, 1) = ly * j / ny;
      }
    }

    // Set connectivity
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int elem = i + nx * j;

        int conn_coord[4];
        for (int jj = 0, index = 0; jj < 2; jj++) {
          for (int ii = 0; ii < 2; ii++, index++) {
            conn_coord[index] = (i + ii) + (nx + 1) * (j + jj);
          }
        }

        // Convert to the correct connectivity
        conn(elem, 0) = conn_coord[0];
        conn(elem, 1) = conn_coord[1];
        conn(elem, 2) = conn_coord[3];
        conn(elem, 3) = conn_coord[2];
      }
    }
  }

  /**
   * @brief populate boundary conditions: fix all dofs along edge x = 0
   *
   * @param bcs boundary condition multiarray
   */
  template <class BcsArray>
  void set_bcs(BcsArray& bcs) {
    index_t index = 0;
    for (int j = 0; j < ny + 1; j++) {
      int i = 0;
      int node = i + (nx + 1) * j;

      // Set the boundary conditions
      bcs(index, 0) = node;
      for (int ii = 0; ii < 2; ii++) {
        bcs(index, 1) |= 1U << ii;
      }
      index++;
    }
  }

  /**
   * @brief uniformly filled design variables to 1.0
   *
   * @param x design variable multiarray
   */
  template <class DvArray>
  void set_dv(DvArray& x) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * j;
        x(node, 0) = 1.0;
      }
    }
  }

  /**
   * @brief For testing helmholtz solver - set lower 1/4 rectangle to be 1.
   *
   * @param x design variable multiarray
   */
  template <class DvArray>
  void set_helm_dv(DvArray& x) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * j;
        if (i < nx / 2 and j < ny / 2) {
          x(node, 0) = 0.1;
        } else {
          x(node, 0) = 1.0;
        }
      }
    }
  }

  /**
   * @brief Set the system residual given a PDE model with boundary conditions
   *
   * @param model the PDE model
   * @param residual the residual multiarray
   */
  template <class Model, class RhsArray>
  void set_force(Model& model, RhsArray& residual) {
    A2D::BLAS::zero(*residual);
    (*residual)(nx, 1) = -1e2;
    model->zero_bcs(residual);
  }

 private:
  int nx, ny;
  double lx, ly;
};

/**
 * @brief A structured mesh generator for 3-dimensional brick domain.
 */
class MesherBrick3D {
 public:
  /**
   * @brief Constructor.
   *
   * @param nx number of elements along x direction
   * @param ny number of elements along y direction
   * @param nz number of elements along z direction
   * @param lx length along x direction
   * @param ly length along y direction
   * @param lz length along z direction
   */
  MesherBrick3D(int _nx, int _ny, int _nz, double _lx, double _ly, double _lz)
      : nx(_nx), ny(_ny), nz(_nz), lx(_lx), ly(_ly), lz(_lz) {}

  /**
   * @brief populate X and conn
   *
   * @param X nodal location multiarray
   * @param conn connectivity multiarray
   */
  template <class ConnArray, class XArray>
  void set_X_conn(XArray& X, ConnArray& conn) {
    Timer t("MesherBrick3D::set_X_conn()");
    // Set X
    for (int k = 0; k < nz + 1; k++) {
      for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx + 1; i++) {
          int node = i + (nx + 1) * (j + (ny + 1) * k);

          X(node, 0) = lx * i / nx;
          X(node, 1) = ly * j / ny;
          X(node, 2) = lz * k / nz;
        }
      }
    }

    // Set connectivity
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          int elem = i + nx * (j + ny * k);

          int conn_coord[8];
          for (int kk = 0, index = 0; kk < 2; kk++) {
            for (int jj = 0; jj < 2; jj++) {
              for (int ii = 0; ii < 2; ii++, index++) {
                conn_coord[index] =
                    (i + ii) + (nx + 1) * ((j + jj) + (ny + 1) * (k + kk));
              }
            }
          }

          // Convert to the correct connectivity
          conn(elem, 0) = conn_coord[0];
          conn(elem, 1) = conn_coord[1];
          conn(elem, 2) = conn_coord[3];
          conn(elem, 3) = conn_coord[2];

          conn(elem, 4) = conn_coord[4];
          conn(elem, 5) = conn_coord[5];
          conn(elem, 6) = conn_coord[7];
          conn(elem, 7) = conn_coord[6];
        }
      }
    }
  }

  /**
   * @brief populate boundary conditions: fix all dofs along face x = 0
   *
   * @param bcs boundary condition multiarray
   */
  template <class BcsArray>
  void set_bcs(BcsArray& bcs) {
    Timer t("MesherBrick3D::set_bcs()");
    index_t index = 0;
    for (int k = 0; k < nz + 1; k++) {
      for (int j = 0; j < ny + 1; j++) {
        int i = 0;
        int node = i + (nx + 1) * (j + (ny + 1) * k);

        // Set the boundary conditions
        bcs(index, 0) = node;
        for (int ii = 0; ii < 3; ii++) {
          bcs(index, 1) |= 1U << ii;
        }
        index++;
      }
    }
  }

  /**
   * @brief uniformly filled design variables to 1.0
   *
   * @param x design variable multiarray
   */
  template <class DvArray>
  void set_dv(DvArray& x) {
    Timer t("MesherBrick3D::set_dv()");
    for (int k = 0; k < nz + 1; k++) {
      for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx + 1; i++) {
          int node = i + (nx + 1) * (j + (ny + 1) * k);
          x(node, 0) = 1.0;
        }
      }
    }
  }

  /**
   * @brief Set the system residual given a PDE model with boundary conditions
   *
   * @param model the PDE model
   * @param residual the residual multiarray
   */
  template <class Model, class RhsArray>
  void set_force(Model model, RhsArray& residual) {
    Timer t("MesherBrick3D::set_force()");
    A2D::BLAS::zero(*residual);
    for (int k = nz / 4; k < 3 * nz / 4; k++) {
      int node = nx + (nx + 1) * (0 + (ny + 1) * k);
      (*residual)(node, 1) = -1e2;
    }
    model->zero_bcs(residual);
  }

 private:
  int nx, ny, nz;
  double lx, ly, lz;
};

/**
 * @brief Load the mesh (connectivity and X only for now) from vtk.
 *
 * @tparam nnodes_per_elem second dimension of the connectivity array
 * @tparam T X entry data type
 * @tparam I connectivity entry data type
 */
template <int nnodes_per_elem, typename T, typename I>
class MesherFromVTK3D {
 public:
  MesherFromVTK3D(const std::string& vtk_name) {
    vtk_reader = VTK_t(vtk_name);
    nnodes = vtk_reader.nnodes;
    nelems = vtk_reader.nelems;

    // Loop over X to find bc and loaded nodes
    // We set those nodes with minimum x-coordinates to be bc nodes
    nbcs = 0;
    nforces_x = 0;
    nforces_y = 0;
    nforces_z = 0;
    T tol = 1e-6;
    for (I i = 0; i != nnodes; i++) {
      // Fix xmin
      if (vtk_reader.X(i, 0) - vtk_reader.domain_lower[0] < tol) {
        bc_nodes.push_back(i);
        nbcs++;
      }

      // // Fix zmax
      // if (vtk_reader.domain_upper[2] - vtk_reader.X(i, 2) < tol) {
      //   bc_nodes.push_back(i);
      //   nbcs++;
      // }

      // Apply force to xmax
      if (vtk_reader.domain_upper[0] - vtk_reader.X(i, 0) < tol) {
        force_nodes_y.push_back(i);
        nforces_y++;
      }

      // // Apply force to ymax
      // if (vtk_reader.domain_upper[1] - vtk_reader.X(i, 1) < tol) {
      //   force_nodes_z.push_back(i);
      //   nforces_z++;
      // }
    }
  }

  I get_nnodes() { return nnodes; }
  I get_nelems() { return nelems; }
  I get_nbcs() { return nbcs; }

  template <class ConnArray, class XArray>
  void set_X_conn(XArray& X, ConnArray& conn) {
    vtk_reader.set_X(X);
    vtk_reader.set_conn(conn);
  }

  template <class BcsArray>
  void set_bcs(BcsArray& bcs) {
    I index = 0;
    for (auto it = bc_nodes.begin(); it != bc_nodes.end(); it++) {
      bcs(index, 0) = *it;
      for (int ii = 0; ii < 3; ii++) {
        bcs(index, 1) |= 1U << ii;
      }
      index++;
    }
  }

  template <class DvArray>
  void set_dv(DvArray& x) {
    for (I i = 0; i != nnodes; i++) {
      x(i, 0) = 1.0;
    }
  }

  template <class Type, class Model, class RhsArray>
  void set_force(Model& model, RhsArray& residual, const Type force) {
    A2D::BLAS::zero(residual);

    for (auto it = force_nodes_y.begin(); it != force_nodes_y.end(); it++) {
      (*residual)(*it, 1) = -force / T(nforces_y);
    }

    for (auto it = force_nodes_z.begin(); it != force_nodes_z.end(); it++) {
      (*residual)(*it, 2) = -force / T(nforces_z);
    }

    model->zero_bcs(residual);
  }

 private:
  using VTK_t = ReadVTK<nnodes_per_elem, 3, T, I>;

  VTK_t vtk_reader;
  I nnodes, nelems, nbcs;
  I nforces_x, nforces_y, nforces_z;

  std::vector<I> bc_nodes;
  std::vector<I> force_nodes_x, force_nodes_y, force_nodes_z;
};

}  // namespace A2D

#endif  // A2D_MESH_H