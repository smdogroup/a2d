#include "a2d.h"

using namespace A2D;

typedef index_t I;
typedef double T;

template <class XArray>
void get_bc_nodes(XArray& X, std::vector<I>& nodes_left_cylinder,
                  std::vector<I>& nodes_right_cylinder) {
  const double tol = 1e-10;
  double d2, r2 = 0.3 * 0.3;
  double x0 = 0.6, y0 = 0.5;
  double x1 = 5.4, y1 = 0.5;

  // Loop over nodes
  I nnodes = X.layout.get_extent(0);
  for (I i = 0; i != nnodes; i++) {
    d2 = (X(i, 0) - x0) * (X(i, 0) - x0) + (X(i, 1) - y0) * (X(i, 1) - y0);
    if (absfunc(d2 - r2) < tol) {
      nodes_left_cylinder.push_back(i);
    }
    d2 = (X(i, 0) - x1) * (X(i, 0) - x1) + (X(i, 1) - y1) * (X(i, 1) - y1);
    if (absfunc(d2 - r2) < tol) {
      nodes_right_cylinder.push_back(i);
    }
  }
}

int main(int argc, char* argv[]) {
  // Define the problem
  static const int SPATIAL_DIM = 3;
  static const int nnodes_per_elem = 4;
  static const int vtk_type_id = 10;  // 3D tetrahedral element

  // Define basis and PDE info
  using Basis = BasisOps<SPATIAL_DIM, TetraLinearBasisFunc, Tetra4ptQuadrature>;
  using ElasticityPDE = ElasticityPDEInfo<SPATIAL_DIM, I, T>;

  // Create the vtk loader
  MesherFromVTK3D<nnodes_per_elem, T, I> mesher("rigid_body.vtk");

  // Set PDE model and element
  I nnodes = mesher.get_nnodes();
  I nelems = mesher.get_nelems();
  I nbcs = mesher.get_nbcs();
  auto model = std::make_shared<FEModel<I, T, ElasticityPDE>>(nnodes, nbcs);
  auto element = std::make_shared<LinElasticityElement<I, T, Basis>>(nelems);
  model->add_element(element);

  // Set the node locations and connectivity
  auto X = model->get_nodes();
  auto conn = element->get_conn();
  mesher.set_X_conn(X, conn);

  // Set the node locations - Note: This must be done after setting the
  // connectivity!
  model->init();

  // Find indices of the boundary nodes
  std::vector<I> nodes_left, nodes_right;
  get_bc_nodes(X, nodes_left, nodes_right);

  // Check
  auto left = MultiArray<T, CLayout<1>>(nnodes);
  auto right = MultiArray<T, CLayout<1>>(nnodes);
  for (auto it = nodes_left.begin(); it != nodes_left.end(); it++) {
    left(*it, 0) = 1.0;
  }
  for (auto it = nodes_right.begin(); it != nodes_right.end(); it++) {
    right(*it, 0) = 1.0;
  }

  // Write to vtk
  ToVTK<decltype(element->get_conn()), decltype(model->get_nodes())> vtk(
      conn, X, vtk_type_id);
  vtk.write_mesh();
  vtk.write_sol("left", left, 0);
  vtk.write_sol("right", right, 0);

  return 0;
}
