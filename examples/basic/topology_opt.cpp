#include "topology.h"
#include "utils/a2dvtk.h"

int main(int argc, char *argv[]) {
  using I = A2D::index_t;
  using T = double;

  A2D::Timer timer("main()");
  Kokkos::initialize();

  std::cout << "Topology linear elasticity\n";
  A2D::TopoLinearElasticity<T, 3> elasticity(70e3, 0.3, 5.0);
  A2D::TestPDEImplementation<T>(elasticity);

  A2D::ReadVTK3D<I, T> readvtk("3d_hex_twisted.vtk");
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();

  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  std::vector<int> ids{100, 101};
  auto verts = readvtk.get_verts_given_cell_entity_id(ids);

  A2D::index_t end_label =
      conn.add_boundary_label_from_verts(verts.size(), verts.data());

  A2D::ToVTK3D tovtk3d(nverts, ntets, tets, nhex, hex, nwedge, wedge, npyrmd,
                       pyrmd, Xloc, "output.vtk");

#if 0
  // Number of elements in each dimension
  const int nx = 8, ny = 4, nz = 4;
  auto node_num = [](int i, int j, int k) {
    return i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  };

  // Number of edges
  const int nverts = (nx + 1) * (ny + 1) * (nz + 1);
  int ntets = 0, nwedge = 0, npyrmd = 0;
  const int nhex = nx * ny * nz;

  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  int hex[8 * nhex];

  using ET = A2D::ElementTypes;

  for (int k = 0, e = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, e++) {
        for (int ii = 0; ii < ET::HEX_VERTS; ii++) {
          hex[8 * e + ii] = node_num(i + ET::HEX_VERTS_CART[ii][0],
                                     j + ET::HEX_VERTS_CART[ii][1],
                                     k + ET::HEX_VERTS_CART[ii][2]);
        }
      }
    }
  }

  double Xloc[3 * nverts];
  for (int k = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        Xloc[3 * node_num(i, j, k)] = (2.0 * i) / nx;
        Xloc[3 * node_num(i, j, k) + 1] = (1.0 * j) / ny;
        Xloc[3 * node_num(i, j, k) + 2] = (1.0 * k) / nz;
      }
    }
  }

  // Constrain the nodes at either end of the block
  const int num_boundary_verts = 2 * (ny + 1) * (nz + 1);
  int boundary_verts[num_boundary_verts];
  for (int k = 0, index = 0; k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(0, j, k);
    }
  }

  for (int k = 0, index = (ny + 1) * (nz + 1); k < nz + 1; k++) {
    for (int j = 0; j < ny + 1; j++, index++) {
      boundary_verts[index] = node_num(nx, j, k);
    }
  }

  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  A2D::index_t end_label =
      conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);
#endif

  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(end_label, basis);

  // Create the finite-element model
  const A2D::index_t degree = 6;
  T E = 70.0e3, nu = 0.3, q = 5.0;
  T design_stress = 200.0, ks_penalty = 50.0;
  TopoOpt<T, degree> topo(conn, bcinfo, E, nu, q);

  // Set the geometry from the node locations
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<TopoOpt<T, degree>::GeoBasis>(nhex, hex, Xloc,
                                                            elem_geo);
  topo.reset_geometry();
  topo.solve();

  // Write the problem to a vtk file
  topo.tovtk("filename.vtk");

  return 0;
}
