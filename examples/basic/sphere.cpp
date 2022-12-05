
#include "topology.h"

template <class GeoBasis, typename T, class GeoElemVec>
void set_geo_spherical(const T alpha, const T R, GeoElemVec &elem_geo) {
  using ET = A2D::ElementTypes;

  double a = alpha * R;

  T Xloc[3 * ET::HEX_VERTS];
  for (A2D::index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
    double u = 1.0 * ET::HEX_VERTS_CART[ii][0];
    double v = 1.0 * ET::HEX_VERTS_CART[ii][1];
    double w = 1.0 * ET::HEX_VERTS_CART[ii][2];

    Xloc[3 * ii] = (2.0 * u - 1.0) * a;
    Xloc[3 * ii + 1] = (2.0 * v - 1.0) * a;
    Xloc[3 * ii + 2] = (2.0 * w - 1.0) * a;
  }

  // Get the geometry values
  typename GeoElemVec::FEDof geo_dof(0, elem_geo);

  for (A2D::index_t ii = 0; ii < GeoBasis::ndof; ii++) {
    double pt[3];
    GeoBasis::get_dof_point(ii, pt);

    T x = pt[0] * a;
    T y = pt[1] * a;
    T z = pt[2] * a;

    // Interpolate to find the basis
    if (ii % 3 == 0) {
      geo_dof[ii] = x;
    } else if (ii % 3 == 1) {
      geo_dof[ii] = y;
    } else if (ii % 3 == 2) {
      geo_dof[ii] = z;
    }
  }

  elem_geo.set_element_values(0, geo_dof);

  // We know the 2-direction will be radial...
  for (A2D::index_t face = 0; face < ET::HEX_FACES; face++) {
    // Get the geometry values
    typename GeoElemVec::FEDof geo_dof(face + 1, elem_geo);

    for (A2D::index_t ii = 0; ii < GeoBasis::ndof; ii++) {
      double pt[3];
      GeoBasis::get_dof_point(ii, pt);

      double N[4];
      N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
      N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
      N[2] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);
      N[3] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);

      // Find the x-y-z coordinate on the surface
      T x0 = 0.0;
      T y0 = 0.0;
      T z0 = 0.0;

      for (A2D::index_t i = 0; i < 4; i++) {
        x0 += N[i] * Xloc[3 * ET::HEX_FACE_VERTS[face][i]];
        y0 += N[i] * Xloc[3 * ET::HEX_FACE_VERTS[face][i] + 1];
        z0 += N[i] * Xloc[3 * ET::HEX_FACE_VERTS[face][i] + 2];
      }
      T norm = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);

      T x1 = x0 * R / norm;
      T y1 = y0 * R / norm;
      T z1 = z0 * R / norm;

      // Project
      T x = 0.5 * (1.0 - pt[2]) * x0 + 0.5 * (1.0 + pt[2]) * x1;
      T y = 0.5 * (1.0 - pt[2]) * y0 + 0.5 * (1.0 + pt[2]) * y1;
      T z = 0.5 * (1.0 - pt[2]) * z0 + 0.5 * (1.0 + pt[2]) * z1;

      // Interpolate to find the basis
      if (ii % 3 == 0) {
        geo_dof[ii] = x;
      } else if (ii % 3 == 1) {
        geo_dof[ii] = y;
      } else if (ii % 3 == 2) {
        geo_dof[ii] = z;
      }
    }

    elem_geo.set_element_values(face + 1, geo_dof);
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();

  const A2D::index_t dim = 3;
  const A2D::index_t degree = 20;
  const A2D::index_t filter_degree = 2;

  using ET = A2D::ElementTypes;
  using T = double;

  // Set up the connectivity between
  double alpha = 0.25;  // Value < 1
  double R = 1.0;       // Radius of the sphere
  double Ra = alpha * R;
  double R45 = R / std::sqrt(2.0);

  const int nverts = 16;
  const int nhex = 7;

  // Set up the connectivity between elements
  int hex[8 * nhex];
  int hex_ext[8];

  // Set up the interior connectivity for the first face
  for (A2D::index_t i = 0; i < 8; i++) {
    hex[i] = i;
    hex_ext[i] = 8 + i;
  }

  // For each face add the faces
  for (A2D::index_t face = 0; face < ET::HEX_FACES; face++) {
    for (A2D::index_t i = 0; i < 4; i++) {
      hex[(face + 1) * ET::HEX_VERTS + i] = hex[ET::HEX_FACE_VERTS[face][i]];
      hex[(face + 1) * ET::HEX_VERTS + i + 4] =
          hex_ext[ET::HEX_FACE_VERTS[face][i]];
    }
  }

  int ntets = 0, nwedge = 0, npyrmd = 0;
  int *tets = NULL, *wedge = NULL, *pyrmd = NULL;
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  A2D::index_t bc_label = 0;
  A2D::index_t traction_label = 1;

  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label, basis);

  // Set the traction components
  T t[3] = {0.0, 0.0, 0.0};

  // Create the finite-element model
  T E = 70.0e3, nu = 0.3, q = 5.0;
  T design_stress = 200.0, ks_penalty = 50.0;
  TopoOpt<T, degree, filter_degree> topo(conn, bcinfo, E, nu, q, traction_label,
                                         t);

  // Set the geometry from the node locations
  auto elem_geo = topo.get_geometry();
  set_geo_spherical<TopoOpt<T, degree, filter_degree>::GeoBasis>(alpha, R,
                                                                 elem_geo);
  topo.reset_geometry();

  // Write the problem to a vtk file
  topo.tovtk("filename.vtk");

  return 0;
}