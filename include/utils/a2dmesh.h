/**
 * @file a2dmesh.h
 * @author Yicong Fu (fuyicong1996@gmail.com)
 * @brief A collection of basic meshes for demonstrative purpose.
 * @version 0.1
 * @date 2022-07-13
 */

#ifndef A2D_MESH_H
#define A2D_MESH_H

namespace A2D {

/**
 * @brief A structured mesh generator for 2-dimensional rectangular domain.
 *
 * @tparam nx number of elements along x direction
 * @tparam ny number of elements along y direction
 */
template <int nx, int ny>
class MesherRect2D {
 public:
  /**
   * @brief populate X and conn
   *
   * @param X nodal location multiarray
   * @param conn connectivity multiarray
   * @param lx length along x direction
   * @param ly length along y direction
   */
  template <class ConnArray, class XArray>
  static void set_X_conn(XArray& X, ConnArray& conn, const double lx = 1.0,
                         const double ly = 1.0) {
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
  static void set_bcs(BcsArray& bcs) {
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
   * @param nx number of elements along x direction
   * @param ny number of elements along y direction
   * @param x design variable multiarray
   */
  template <class DvArray>
  static void set_dv(DvArray& x) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * j;
        x(node, 0) = 1.0;
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
  static void set_force(Model& model, RhsArray& residual) {
    residual->zero();
    (*residual)(nx, 1) = -1e2;
    model->zero_bcs(residual);
  }
};

/**
 * @brief A structured mesh generator for 3-dimensional brick domain.
 *
 * @tparam nx number of elements along x direction
 * @tparam ny number of elements along y direction
 * @tparam nz number of elements along z direction
 */
template <int nx, int ny, int nz>
class MesherBrick3D {
 public:
  /**
   * @brief populate X and conn
   *
   * @param X nodal location multiarray
   * @param conn connectivity multiarray
   * @param lx length along x direction
   * @param ly length along y direction
   * @param lz length along z direction
   */
  template <class ConnArray, class XArray>
  static void set_X_conn(XArray& X, ConnArray& conn, const double lx = 1.0,
                         const double ly = 1.0, const double lz = 1.0) {
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
  static void set_bcs(BcsArray& bcs) {
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
  static void set_dv(DvArray& x) {
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
  static void set_force(Model model, RhsArray& residual) {
    residual->zero();
    for (int k = nz / 4; k < 3 * nz / 4; k++) {
      int node = nx + (nx + 1) * (0 + (ny + 1) * k);
      (*residual)(node, 1) = -1e2;
    }
    model->zero_bcs(residual);
  }
};

}  // namespace A2D

#endif  // A2D_MESH_H