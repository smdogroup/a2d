#ifndef A2D_HEX_TOOLS_H
#define A2D_HEX_TOOLS_H

#include <string>

#include "feelementvector.h"
#include "utils/a2dprofiler.h"
#include "utils/a2dvtk.h"

namespace A2D {

template <class GeoBasis, typename I, typename T, class GeoElemVec,
          ElemVecType evtype>
void set_geo_from_hex_nodes(const index_t nhex, const I hex[], const T Xloc[],
                            ElementVectorBase<evtype, GeoElemVec> &elem_geo) {
  for (int e = 0; e < nhex; e++) {
    // Get the geometry values
    typename GeoElemVec::FEDof geo_dof(e, elem_geo);

    for (int ii = 0; ii < GeoBasis::ndof; ii++) {
      double pt[3];
      GeoBasis::get_dof_point(ii, pt);

      double N[8];
      N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
      N[2] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[3] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
      N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
      N[6] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
      N[7] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

      // Interpolate to find the basis
      if (ii % 3 == 0) {
        T x = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          x += N[kk] * Xloc[3 * hex[8 * e + kk]];
        }
        geo_dof[ii] = x;
      } else if (ii % 3 == 1) {
        T y = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          y += N[kk] * Xloc[3 * hex[8 * e + kk] + 1];
        }
        geo_dof[ii] = y;
      } else if (ii % 3 == 2) {
        T z = 0.0;
        for (index_t kk = 0; kk < 8; kk++) {
          z += N[kk] * Xloc[3 * hex[8 * e + kk] + 2];
        }
        geo_dof[ii] = z;
      }
    }

    if constexpr (evtype == ElemVecType::Serial) {
      elem_geo.set_element_values(e, geo_dof);
    }
  }
  if constexpr (evtype == ElemVecType::Parallel) {
    elem_geo.set_values();
  }
}

template <index_t outputs, index_t degree, typename T, class DataBasis,
          class GeoBasis, class Basis, class PDE, class DataElemVec,
          class GeoElemVec, class ElemVec, class FunctorType>
void write_hex_to_vtk(PDE &pde, DataElemVec &elem_data, GeoElemVec &elem_geo,
                      ElemVec &elem_sol, const std::string filename,
                      const FunctorType &func) {
  Timer timer("write_hex_to_vtk()");
  using ET = ElementTypes;
  const index_t nex = degree;
  using QuadPts = HexGaussLobattoQuadrature<nex + 1>;

  const index_t nhex = elem_sol.get_num_elements();
  index_t nvtk_elems = nhex * nex * nex * nex;
  index_t nvtk_nodes = nhex * (nex + 1) * (nex + 1) * (nex + 1);

  auto vtk_node_num = [](index_t i, index_t j, index_t k) {
    return i + (nex + 1) * (j + (nex + 1) * k);
  };

  MultiArrayNew<int *[8]> vtk_conn("vtk_elems", nvtk_elems);
  MultiArrayNew<T *[3]> vtk_nodes("vtk_nodes", nvtk_nodes);
  MultiArrayNew<T *[outputs]> vtk_outputs("vtk_outputs", nvtk_nodes);

  for (index_t n = 0, counter = 0; n < nhex; n++) {
    // Get the data values
    typename DataElemVec::FEDof data_dof(n, elem_data);
    elem_data.get_element_values(n, data_dof);
    QptSpace<QuadPts, typename PDE::DataSpace> data;
    DataBasis::template interp(data_dof, data);

    // Get the geometry values
    typename GeoElemVec::FEDof geo_dof(n, elem_geo);
    elem_geo.get_element_values(n, geo_dof);
    QptSpace<QuadPts, typename PDE::FiniteElementGeometry> geo;
    GeoBasis::template interp(geo_dof, geo);

    // Get the degrees of freedom for the element
    typename ElemVec::FEDof sol_dof(n, elem_sol);
    elem_sol.get_element_values(n, sol_dof);
    QptSpace<QuadPts, typename PDE::FiniteElementSpace> sol;
    Basis::template interp(sol_dof, sol);

    const index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

    for (index_t k = 0; k < nex + 1; k++) {
      for (index_t j = 0; j < nex + 1; j++) {
        for (index_t i = 0; i < nex + 1; i++) {
          const index_t index = vtk_node_num(i, j, k);
          const index_t node = off + index;

          typename PDE::FiniteElementSpace &sref = sol.get(index);
          typename PDE::FiniteElementGeometry &gref = geo.get(index);

          // Initialize the transform object
          T detJ;
          typename PDE::SolutionMapping transform(gref, detJ);

          // Transform the solution the physical element
          typename PDE::FiniteElementSpace x, s;
          transform.transform(sref, s);

          auto X = gref.template get<0>().get_value();
          vtk_nodes(node, 0) = X(0);
          vtk_nodes(node, 1) = X(1);
          vtk_nodes(node, 2) = X(2);

          for (index_t kk = 0; kk < outputs; kk++) {
            vtk_outputs(node, kk) = func(kk, data.get(index), gref, s);
          }
        }
      }
    }

    for (index_t k = 0; k < nex; k++) {
      for (index_t j = 0; j < nex; j++) {
        for (index_t i = 0; i < nex; i++, counter++) {
          index_t off = n * (nex + 1) * (nex + 1) * (nex + 1);

          for (index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
            vtk_conn(counter, ii) =
                off + vtk_node_num(i + ET::HEX_VERTS_CART[ii][0],
                                   j + ET::HEX_VERTS_CART[ii][1],
                                   k + ET::HEX_VERTS_CART[ii][2]);
          }
        }
      }
    }
  }

  ToVTK vtk(vtk_conn, vtk_nodes, -1, filename);
  vtk.write_mesh();

  MultiArrayNew<T *> vtk_vec("vtk_vec", nvtk_nodes);
  for (index_t i = 0; i < outputs; i++) {
    for (index_t j = 0; j < nvtk_nodes; j++) {
      vtk_vec(j) = vtk_outputs(j, i);
    }
    char name[256];
    std::snprintf(name, sizeof(name), "solution%d", i + 1);
    vtk.write_sol(name, vtk_vec);
  }
}

}  // namespace A2D

#endif  // A2D_HEX_TOOLS_H