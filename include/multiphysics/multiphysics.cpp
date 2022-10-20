#include <iostream>
#include <set>
#include <tuple>
#include <vector>

#include "a2dobjs.h"
#include "a2dtmp2d.h"
#include "a2dtmp3d.h"

using namespace A2D;

template <typename T, index_t D>
class L1ScalarSpace {
 public:
  L1ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const index_t components = 1;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const index_t ndof = 1;

  // Set the input seed
  void set_seed(const index_t seed) { u = 1.0; }
  T get_value(const index_t seed) { return u; }

  T& get_value() { return u; }
  const T& get_value() const { return u; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, L1ScalarSpace<T, D>& s) {
    s.u = u;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, L1ScalarSpace<T, D>& s) {
    s.u = u;
  }

 private:
  T u;
};

template <typename T, index_t D>
class H1ScalarSpace {
 public:
  H1ScalarSpace() { u = 0.0; }

  // Number of solution components
  static const index_t components = 1;

  // Spatial dimension
  static const index_t dim = D;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const index_t ndof = 1 + D;

  // Set the seed value
  void set_seed(const index_t seed) {
    if (seed == 0) {
      u = 1.0;
    } else {
      grad(seed - 1) = 1.0;
    }
  }
  T get_value(const index_t seed) {
    if (seed == 0) {
      return u;
    } else {
      return grad(seed - 1);
    }
  }

  T& get_value() { return u; }
  A2D::Vec<T, D>& get_grad() { return grad; }

  const T& get_value() const { return u; }
  const A2D::Vec<T, D>& get_grad() const { return grad; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }

 private:
  T u;
  A2D::Vec<T, D> grad;
};

template <typename T, index_t C, index_t D>
class H1Space {
 public:
  H1Space() {}

  // Number of solution components
  static const index_t components = C;

  // Spatial dimension
  static const index_t dim = D;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const index_t ndof = (D + 1) * C;

  // Set the seed value
  void set_seed(const index_t seed) {
    if (seed < C) {
      u(seed) = 1.0;
    } else {
      grad((seed - C) / D, (seed - C) % D) = 1.0;
    }
  }
  T get_value(const index_t seed) {
    if (seed < C) {
      return u(seed);
    } else {
      return grad((seed - C) / D, (seed - C) % D);
    }
  }

  A2D::Vec<T, C>& get_value() { return u; }
  A2D::Mat<T, C, D>& get_grad() { return grad; }

  const A2D::Vec<T, C>& get_value() const { return u; }
  const A2D::Mat<T, C, D>& get_grad() const { return grad; }

  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, H1ScalarSpace<T, D>& s) {
    s.u = u;
    s.grad = grad;
  }

 private:
  A2D::Vec<T, C> u;
  A2D::Mat<T, C, D> grad;
};

template <typename T>
class Hdiv2DSpace {
 public:
  Hdiv2DSpace() { div = 0.0; }
  // Number of solution components
  static const index_t components = 2;

  // Spatial dimension
  static const index_t dim = 2;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const index_t ndof = 3;

  void set_seed(const index_t seed) {
    if (seed < 2) {
      u(seed) = 1.0;
    } else {
      div = 1.0;
    }
  }
  T get_value(const index_t seed) {
    if (seed < 2) {
      return u(seed);
    } else {
      return div;
    }
  }

  A2D::Vec<T, 2>& get_value() { return u; }
  T& get_div() { return div; }

  const A2D::Vec<T, 2>& get_value() const { return u; }
  const T& get_div() const { return div; }

  void transform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                 const A2D::Mat<T, 2, 2>& Jinv, Hdiv2DSpace<T>& s) {
    s.u = u;
    s.div = div;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                  const A2D::Mat<T, 2, 2>& Jinv, Hdiv2DSpace<T>& s) {
    s.u = u;
    s.div = div;
  }

 private:
  A2D::Vec<T, 2> u;
  T div;
};

template <typename T>
class Hcurl2DSpace {
 public:
  Hcurl2DSpace() { curl = 0.0; }

  // Number of solution components
  static const index_t components = 2;

  // Spatial dimension
  static const index_t dim = 2;

  // Number of degrees of freedom. This count treats each element of the
  // solution/gradient or other resultant as independent
  static const index_t ndof = 3;

  void set_seed(const index_t seed) {
    if (seed < 2) {
      u(seed) = 1.0;
    } else {
      curl = 1.0;
    }
  }
  T get_value(const index_t seed) {
    if (seed < 2) {
      return u(seed);
    } else {
      return curl;
    }
  }

  A2D::Vec<T, 2>& get_value() { return u; }
  T& get_curl() { return curl; }

  const A2D::Vec<T, 2>& get_value() const { return u; }
  const T& get_curl() const { return curl; }

  void transform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                 const A2D::Mat<T, 2, 2>& Jinv, Hcurl2DSpace<T>& s) {
    s.u = u;
    s.curl = curl;
  }
  void rtransform(const T& detJ, const A2D::Mat<T, 2, 2>& J,
                  const A2D::Mat<T, 2, 2>& Jinv, Hcurl2DSpace<T>& s) {
    s.u = u;
    s.curl = curl;
  }

 private:
  A2D::Vec<T, 2> u;
  T curl;
};

template <class... Spaces>
struct __count_solution_ndof;

template <>
struct __count_solution_ndof<> {
  static const index_t ndof = 0;
};

template <class First, class... Remain>
struct __count_solution_ndof<First, Remain...> {
  static const index_t ndof =
      First::ndof + __count_solution_ndof<Remain...>::ndof;
};

/*
  A collection of spaces used in the finite-element problem
*/
template <typename T, index_t D, class... Spaces>
class FESpace {
 public:
  typedef std::tuple<Spaces...> SolutionSpace;
  static constexpr index_t nspaces = std::tuple_size<std::tuple<Spaces...>>();

  /*
    Count up the total number of degrees of freedom
  */
  static constexpr index_t ndof = __count_solution_ndof<Spaces...>::ndof;

  /*
    Set a seed on a given solution field
  */
  void set_seed(const index_t seed) { set_seed_<0, Spaces...>(seed); }

  /*
    Get a solution value based on the index
  */
  T get_value(const index_t seed) { return get_value_<0, Spaces...>(seed); }

  /*
    Extract the specified solution space from the FESpace object
  */
  template <index_t index>
  typename std::tuple_element<index, SolutionSpace>::type& get() noexcept {
    return std::get<index>(u);
  }

  template <index_t index>
  const typename std::tuple_element<index, SolutionSpace>::type& get()
      const noexcept {
    return std::get<index>(u);
  }

  /*
    Transform the solution space from the reference element to the physical
    element
  */
  void transform(const T& detJ, const A2D::Mat<T, D, D>& J,
                 const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    transform_<0, Spaces...>(detJ, J, Jinv, s);
  }

  /*
    Perform the reverse of the transform - transfer the derivative from the
    physical element to the reference element
  */
  void rtransform(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    rtransform_<0, Spaces...>(detJ, J, Jinv, s);
  }

 private:
  // Solution space tuple object
  SolutionSpace u;

  template <index_t index>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {}

  template <index_t index, class First, class... Remain>
  void transform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                  const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    std::get<index>(u).transform(detJ, J, Jinv, std::get<index>(s.u));
    transform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <index_t index>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
  }

  template <index_t index, class First, class... Remain>
  void rtransform_(const T& detJ, const A2D::Mat<T, D, D>& J,
                   const A2D::Mat<T, D, D>& Jinv, FESpace<T, D, Spaces...>& s) {
    std::get<index>(u).rtransform(detJ, J, Jinv, std::get<index>(s.u));
    rtransform_<index + 1, Remain...>(detJ, J, Jinv, s);
  }

  template <index_t index>
  void set_seed_(const index_t seed) {}

  template <index_t index, class First, class... Remain>
  void set_seed_(const index_t seed) {
    if (seed < First::ndof) {
      std::get<index>(u).set_seed(seed);
    } else {
      set_seed_<index + 1, Remain...>(seed - First::ndof);
    }
  }

  template <index_t index>
  T get_value_(const index_t seed) {
    return 0.0;
  }

  template <index_t index, class First, class... Remain>
  T get_value_(const index_t seed) {
    if (seed < First::ndof) {
      return std::get<index>(u).get_value(seed);
    } else {
      return get_value_<index + 1, Remain...>(seed - First::ndof);
    }
  }
};

/*
  Quadrature class for triangles
*/
const double TriangleWts3[] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
const double TrianglePts3[] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0,
                               1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

class TriQuadrature3 {
 public:
  static index_t get_num_points() { return 3; }
  static void get_point(const index_t n, double pt[]) {
    pt[0] = TrianglePts3[2 * n];
    pt[1] = TrianglePts3[2 * n + 1];
  }
  static double get_weight(const index_t n) { return TriangleWts3[n]; }
};

/*
  Lagrange basis for a triangle
*/
template <typename T>
class LagrangeTri0 {
 public:
  static const int ndof = 1;

  template <class Quadrature, class SolnType>
  static void interp(index_t n, const SolnType& sol, L1ScalarSpace<T, 2>& out) {
    T& u = out.get_value();
    u = sol[0];
  }

  template <class Quadrature, class SolnType>
  static void add(index_t n, const L1ScalarSpace<T, 2>& in, SolnType& res) {
    const T& u = in.get_value();
    res[0] += u;
  }
};

/*
  Lagrange basis for a triangle
*/
template <typename T, index_t C>
class LagrangeTri1 {
 public:
  static const index_t ndof = 3 * C;

  template <class Quadrature, class SolnType>
  static void interp(index_t n, const SolnType& sol, H1ScalarSpace<T, 2>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    T& u = out.get_value();
    A2D::Vec<T, 2>& grad = out.get_grad();

    u = N[0] * sol[0] + N[1] * sol[1] + N[2] * sol[2];
    grad(0) = sol[1] - sol[0];
    grad(1) = sol[2] - sol[0];
  }

  template <class Quadrature, class SolnType>
  static void add(index_t n, const H1ScalarSpace<T, 2>& in, SolnType& res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    const T& u = in.get_value();
    const A2D::Vec<T, 2>& grad = in.get_grad();

    res[0] += N[0] * u - grad(0) - grad(1);
    res[1] += N[1] * u + grad(0);
    res[2] += N[2] * u + grad(1);
  }

  template <class Quadrature, class SolnType>
  static void interp(index_t n, const SolnType& sol, H1Space<T, C, 2>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    A2D::Vec<T, C>& u = out.get_value();
    A2D::Mat<T, C, 2>& grad = out.get_grad();

    for (index_t i = 0; i < C; i++) {
      u(i) = N[0] * sol[i] + N[1] * sol[C + i] + N[2] * sol[2 * C + i];
      grad(i, 0) = sol[C + i] - sol[i];
      grad(i, 1) = sol[2 * C + i] - sol[i];
    }
  }

  template <class Quadrature, class SolnType>
  static void add(index_t n, const H1Space<T, C, 2>& in, SolnType& res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    // Compute the shape functions
    double N[3];
    N[0] = 1.0 - pt[0] - pt[1];
    N[1] = pt[0];
    N[2] = pt[1];

    const A2D::Vec<T, C>& u = in.get_value();
    const A2D::Mat<T, C, 2>& grad = in.get_grad();

    for (index_t i = 0; i < C; i++) {
      res[i] += N[0] * u(i) - grad(i, 0) - grad(i, 1);
      res[C + i] += N[1] * u(i) + grad(i, 0);
      res[2 * C + i] += N[2] * u(i) + grad(i, 1);
    }
  }
};

/*
  Raviart-Thomas element for H(div) in 2D
*/
template <typename T>
class RT2DTri1 {
 public:
  static const int ndof = 3;

  template <class Quadrature, class SolnType>
  static void interp(index_t n, const SolnType& sol, Hdiv2DSpace<T>& out) {
    double pt[2];
    Quadrature::get_point(n, pt);

    A2D::Vec<T, 2>& u = out.get_value();
    T& div = out.get_div();

    u(0) = pt[0] * sol[0] + (pt[0] - 1.0) * sol[1] + pt[0] * sol[2];
    u(1) = pt[1] * sol[0] + pt[1] * sol[1] + (pt[1] - 1.0) * sol[2];
    div = 2.0 * (sol[0] + sol[1] + sol[2]);
  }

  template <class Quadrature, class SolnType>
  static void add(index_t n, const Hdiv2DSpace<T>& in, SolnType& res) {
    double pt[2];
    Quadrature::get_point(n, pt);

    const A2D::Vec<T, 2>& u = in.get_value();
    const T& div = in.get_div();

    res[0] += pt[0] * u(0) + pt[1] * u(1) + 2.0 * div;
    res[1] += (pt[0] - 1.0) * u(0) + pt[1] * u(1) + 2.0 * div;
    res[2] += pt[0] * u(0) + (pt[1] - 1.0) * u(1) + 2.0 * div;
  }
};

const int tri_edge_nodes[3][2] = {{2, 1}, {2, 0}, {0, 1}};

/*
  Compute a node to triangle or node to quad data structure
*/
void compute_nodes_to_elem(index_t nnodes, index_t nelems,
                           index_t num_elem_nodes, const index_t conn[],
                           index_t** _ptr, index_t** _node_to_elems) {
  // Set the pointer
  index_t* ptr = new index_t[nnodes + 1];
  std::fill(ptr, ptr + nnodes + 1, 0);

  // Count up the references
  const index_t conn_size = nelems * num_elem_nodes;
  for (index_t i = 0; i < conn_size; i++) {
    ptr[conn[i] + 1]++;
  }

  // Set the pointer into the array
  for (index_t i = 0; i < nnodes; i++) {
    ptr[i + 1] += ptr[i];
  }

  // Compute the node to quads
  index_t* node_to_elems = new index_t[ptr[nnodes]];
  const index_t* conn_ptr = conn;
  for (index_t i = 0; i < nelems; i++) {
    for (index_t j = 0; j < num_elem_nodes; j++) {
      index_t node = conn_ptr[0];
      if (node >= 0) {
        node_to_elems[ptr[node]] = i;
        ptr[node]++;
        conn_ptr++;
      }
    }
  }

  // Reset the pointer array
  for (index_t i = nnodes; i > 0; i--) {
    ptr[i] = ptr[i - 1];
  }
  ptr[0] = 0;

  // Set the output points
  *_ptr = ptr;
  *_node_to_elems = node_to_elems;
}

/*
  Compute all of the edges within the triangular mesh
*/
index_t compute_planar_edges(index_t nnodes, index_t ntris,
                             const index_t tris[], index_t* tri_edge_nums,
                             index_t* tri_orient) {
  // Compute the edges in the triangular mesh
  index_t* ptr;
  index_t* node_to_tris;
  compute_nodes_to_elem(nnodes, ntris, 3, tris, &ptr, &node_to_tris);

  // Set the no-label
  const index_t no_label = std::numeric_limits<index_t>::max();

  // Now compute the neighbors for each triangle
  for (index_t i = 0; i < 3 * ntris; i++) {
    tri_edge_nums[i] = no_label;
  }

  index_t ne = 0;
  for (index_t i = 0; i < ntris; i++) {
    // Search through each edge of the each triangle
    for (index_t j = 0; j < 3; j++) {
      if (tri_edge_nums[3 * i + j] == no_label) {
        tri_edge_nums[3 * i + j] = ne;
        tri_orient[3 * i + j] = 1;

        index_t e0[2];
        e0[0] = tris[3 * i + tri_edge_nodes[j][0]];
        e0[1] = tris[3 * i + tri_edge_nodes[j][1]];

        // Search for the neighboring that shares this edge
        index_t kp = ptr[e0[0]];
        index_t kpend = ptr[e0[0] + 1];
        for (; kp < kpend; kp++) {
          // Find the potential neighbor
          index_t n = node_to_tris[kp];

          // Don't count the same edge twice
          if (n == i) {
            continue;
          }

          // Flag to indicate that we have found the other edge (there
          // will only be at most one other match since this is
          // planar in parameter space)
          index_t quit = 0;

          // Search over all the edges on this quad, and see
          // if they match
          for (index_t e = 0; e < 3; e++) {
            index_t e1[2];
            e1[0] = tris[3 * n + tri_edge_nodes[e][0]];
            e1[1] = tris[3 * n + tri_edge_nodes[e][1]];

            // Check if the adjacent edge matches in either direction
            if ((e0[0] == e1[0] && e0[1] == e1[1]) ||
                (e0[0] == e1[1] && e0[1] == e1[0])) {
              // Label the other edge that shares this same node
              tri_edge_nums[3 * n + e] = ne;

              // Opposite orientation
              tri_orient[3 * i + j] = 0;

              quit = 1;
            }
          }
          if (quit) {
            break;
          }
        }

        // Increment the edge number
        ne++;
      }
    }
  }

  // Free the data that is no longer required
  delete[] ptr;
  delete[] node_to_tris;

  return ne;
}

/*
  The element connectivity class

  This class will need a substantial overhaul
  1. The connectivity is fixed at this point
  2.
*/
class ElementConnectivity {
 public:
  ElementConnectivity(index_t nnodes, index_t nelems, index_t* conn_)
      : nnodes(nnodes), nelems(nelems) {
    conn = new index_t[3 * nelems];
    for (index_t i = 0; i < 3 * nelems; i++) {
      conn[i] = conn_[i];
    }

    edges = new index_t[3 * nelems];
    orient = new index_t[3 * nelems];
    nedges = compute_planar_edges(nnodes, nelems, conn, edges, orient);
  }

  index_t get_face_dof(index_t elem, index_t index) { return elem; }

  index_t get_edge_dof(index_t elem, index_t index, int& ort) {
    if (orient[3 * elem + index]) {
      ort = 1;
    } else {
      ort = -1;
    }

    return edges[3 * elem + index];
  }

  index_t get_node_dof(index_t elem, index_t index) {
    return conn[3 * elem + index];
  }

  index_t get_num_elements() { return nelems; }
  index_t get_num_edges() { return nedges; }
  index_t get_num_nodes() { return nnodes; }

 private:
  index_t nnodes;
  index_t nedges;
  index_t nelems;
  index_t* conn;
  index_t* edges;
  index_t* orient;  // edge orientation 1 for +ve, 0 for -ve
};

/*
  ElementMesh - Map from an element to the global to element local degrees
  of freedom
*/
enum SpaceType { L2, H1, EDGE };

template <class... Basis>
class ElementMesh {
 public:
  ElementMesh(ElementConnectivity& conn, SpaceType spaces_[],
              index_t dim_[] = NULL)
      : conn(conn) {
    for (index_t i = 0; i < nspaces; i++) {
      spaces[i] = spaces_[i];
    }
    if (dim_) {
      for (index_t i = 0; i < nspaces; i++) {
        dim[i] = dim_[i];
      }
    } else {
      for (index_t i = 0; i < nspaces; i++) {
        dim[i] = 1;
      }
    }
    init_<0, Basis...>();
  }

  static const index_t nspaces = sizeof...(Basis);

  index_t get_num_elements() { return conn.get_num_elements(); }
  index_t get_num_dof() { return offset[nspaces]; }

  int get_global_dof_sign(index_t elem, index_t space, index_t index) {
    if (spaces[space] == EDGE) {
      int orient;
      conn.get_edge_dof(elem, index, orient);
      return orient;
    }
    return 1;
  }
  index_t get_global_dof(index_t elem, index_t space, index_t index) {
    if (spaces[space] == L2) {
      return offset[space] +
             dim[space] * conn.get_face_dof(elem, index / dim[space]) +
             index % dim[space];
    } else if (spaces[space] == H1) {
      return offset[space] +
             dim[space] * conn.get_node_dof(elem, index / dim[space]) +
             index % dim[space];
    } else if (spaces[space] == EDGE) {
      int orient;
      return offset[space] +
             dim[space] * conn.get_edge_dof(elem, index / dim[space], orient) +
             index % dim[space];
    }
    return 0;
  }

 private:
  template <index_t index>
  void init_() {}

  template <index_t index, class First, class... Remain>
  void init_() {
    if (index == 0) {
      offset[index] = 0;
    }

    if (spaces[index] == L2) {
      counts[index] = dim[index] * conn.get_num_elements();
    } else if (spaces[index] == H1) {
      counts[index] = dim[index] * conn.get_num_nodes();
    } else if (spaces[index] == EDGE) {
      counts[index] = dim[index] * conn.get_num_edges();
    }

    offset[index + 1] = offset[index] + counts[index];

    init_<index + 1, Remain...>();
  }

  ElementConnectivity conn;
  index_t counts[nspaces];
  index_t offset[nspaces + 1];
  SpaceType spaces[nspaces];
  index_t dim[nspaces];
};

/*
  The solution vector
*/
template <typename T>
class SolutionVector {
 public:
  SolutionVector(index_t ndof) : ndof(ndof), x(ndof) {}
  T& operator[](index_t index) { return x[index]; }
  const T& operator[](index_t index) const { return x[index]; }

 private:
  const index_t ndof;
  std::vector<T> x;
};

/*
  In-place element vector implementation
*/
template <typename T, class FiniteElementSpace, class... Basis>
class ElementVectorInPlace {
 public:
  ElementVectorInPlace(ElementMesh<Basis...>& mesh, SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {}

  SolutionVector<T>& get_vector() { return vec; }
  void set_solution() {}
  void get_residual() {}

  template <class Quadrature>
  void interp(index_t elem, index_t pt, FiniteElementSpace& s) {
    interp_<Quadrature, 0, Basis...>(elem, pt, s);
  }

  template <class Quadrature>
  void add(index_t elem, index_t pt, const FiniteElementSpace& s) {
    add_<Quadrature, 0, Basis...>(elem, pt, s);
  }

 private:
  template <class Quadrature, index_t index, class First, class... Remain>
  void interp_(index_t elem, index_t pt, FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    for (index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const index_t dof = mesh.get_global_dof(elem, index, i);
      values[i] = sign * vec[dof];
    }

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, index_t index>
  void interp_(index_t elem, index_t pt, FiniteElementSpace& s) {}

  template <class Quadrature, index_t index, class First, class... Remain>
  void add_(index_t elem, index_t pt, const FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    std::fill(values, values + First::ndof, T(0.0));

    // Add the interpolation
    First::template add<Quadrature>(pt, s.template get<index>(), values);

    for (index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const index_t dof = mesh.get_global_dof(elem, index, i);
      vec[dof] += sign * values[i];
    }

    // Do the next solution space, if any...
    add_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, index_t index>
  void add_(index_t elem, index_t pt, const FiniteElementSpace& s) {}

  ElementMesh<Basis...>& mesh;
  SolutionVector<T>& vec;
};

/*
template <typename T, class FiniteElementSpace, class... Basis>
class ElementVector2 {
 public:
  static const index_t ndof_per_element = dof_in_basis<Basis...>;

  ElementVector2(ElementMesh<Basis...>& mesh, SolutionVector<T>& vec)
      : mesh(mesh), vec(vec) {
    index_t nelems = mesh.get_num_elements();
  }

  SolutionVector<T>& get_vector() { return vec; }

  void set_solution() {
    // Go set values into element_dof and transfer them to the device
    for (index_t elem = 0; elem < num_elements; elem++) {
      for (index_t i = 0; i < First::ndof; i++) {
        const int sign = mesh.get_global_dof_sign(elem, index, i);
        const index_t dof = mesh.get_global_dof(elem, index, i);
        element_dof[elem, i] = sign * vec[dof];
      }
    }
  }

  void get_residual() {
    for (index_t elem = 0; elem < num_elements; elem++) {
      for (index_t i = 0; i < First::ndof; i++) {
        const int sign = mesh.get_global_dof_sign(elem, index, i);
        const index_t dof = mesh.get_global_dof(elem, index, i);
        vec[dof] += sign * element_dof[elem, i];
      }
    }
  }

  template <class Quadrature>
  void interp(index_t elem, index_t pt, FiniteElementSpace& s) {
    interp_<Quadrature, 0, Basis...>(elem, pt, s);
  }

  template <class Quadrature, index_t index>
  void interp_(index_t elem, index_t pt, FiniteElementSpace& s) {}

  template <class Quadrature, index_t index, class First, class... Remain>
  void interp_(index_t elem, index_t pt, FiniteElementSpace& s) {
    const values = subview(element_dof, elem, Kokkos::ALL);

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    interp_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, index_t index, class First, class... Remain>
  void add_(index_t elem, index_t pt, const FiniteElementSpace& s) {
    values = subview(element_dof, elem, Kokkos::ALL);

    // Add the interpolation
    First::template add<Quadrature>(pt, s.template get<index>(), values);

    // Do the next solution space, if any...
    add_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, index_t index>
  void add_(index_t elem, index_t pt, const FiniteElementSpace& s) {}

 private:
  Kokkos::View<T* [ndof_per_element]> element_dof;
};
*/

template <typename T, class FiniteElementSpace, class... Basis>
class ElementMatrix {
 public:
  ElementMatrix(ElementMesh<Basis...>& mesh, SolutionMatrix<T>& vec)
      : mesh(mesh), vec(vec) {}

  template <class Quadrature>
  void outer(index_t elem, index_t pt, FiniteElementSpace& s) {
    outer_<Quadrature, 0, Basis...>(elem, pt, s);
  }

 private:
  template <class Quadrature, index_t index, class First, class... Remain>
  void outer_(index_t elem, index_t pt, FiniteElementSpace& s) {
    // Un-pack to a local array
    T values[First::ndof];
    for (index_t i = 0; i < First::ndof; i++) {
      const int sign = mesh.get_global_dof_sign(elem, index, i);
      const index_t dof = mesh.get_global_dof(elem, index, i);
      values[i] = sign * vec[dof];
    }

    // Interpolate
    First::template interp<Quadrature>(pt, values, s.template get<index>());

    // Do the next solution space, if any...
    outer_<Quadrature, index + 1, Remain...>(elem, pt, s);
  }

  template <class Quadrature, index_t index>
  void outer_(index_t elem, index_t pt, FiniteElementSpace& s) {}

  ElementMesh<Basis...>& mesh;
  SolutionMatrix<T>& mat;
};

/*
  The PDE object
*/
template <typename T>
class MixedPoisson2D {
 public:
  // Spatial dimension
  static const index_t dim = 2;

  // Finite element space
  typedef FESpace<T, dim, L1ScalarSpace<T, dim>, Hdiv2DSpace<T>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1 in 2D
  typedef FESpace<T, dim, H1Space<T, dim, dim>> FiniteElementGeometry;

  /**
   * Evaluate the variational form at a quadrature point
   *
   * <tau, sigma> - <div(tau), u> + <v, div(sigma)>
   */
  static T eval_weak_form(T wdetJ, FiniteElementSpace& s,
                          FiniteElementSpace& t) {
    // Field objects for solution functions
    L1ScalarSpace<T, 2>& u = s.template get<0>();
    Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    L1ScalarSpace<T, 2>& v = t.template get<0>();
    Hdiv2DSpace<T>& tau = t.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // return A2D::VecDot(tau_val, sigma_val) - tau_div * u + v * sigma_div;
    return wdetJ * (tau_val(0) * sigma_val(0) + tau_val(0) * sigma_val(0) -
                    tau_div * u + v * sigma_div);
  }

  static void eval_weak_coef(T wdetJ, FiniteElementSpace& s,
                             FiniteElementSpace& coef) {
    // Field objects for solution functions
    L1ScalarSpace<T, 2>& u = s.template get<0>();
    Hdiv2DSpace<T>& sigma = s.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    L1ScalarSpace<T, 2>& v = coef.template get<0>();
    Hdiv2DSpace<T>& tau = coef.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    tau_val(0) = wdetJ * sigma_val(0);
    tau_val(1) = wdetJ * sigma_val(1);

    tau_div = -wdetJ * u_val;

    v_val = wdetJ * sigma_div;
  }

  static void eval_weak_jacobian_vec_product(T wdetJ, FiniteElementSpace& s,
                                             FiniteElementSpace& p,
                                             FiniteElementSpace& coef) {
    // Field objects for solution functions
    L1ScalarSpace<T, 2>& u = p.template get<0>();
    Hdiv2DSpace<T>& sigma = p.template get<1>();

    // Solution function values
    A2D::Vec<T, 2>& sigma_val = sigma.get_value();
    T& sigma_div = sigma.get_div();
    T& u_val = u.get_value();

    // Test function values
    L1ScalarSpace<T, 2>& v = coef.template get<0>();
    Hdiv2DSpace<T>& tau = coef.template get<1>();

    // Test function values
    A2D::Vec<T, 2>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    tau_val(0) = wdetJ * sigma_val(0);
    tau_val(1) = wdetJ * sigma_val(1);

    tau_div = -wdetJ * u_val;

    v_val = wdetJ * sigma_div;
  }
};

template <typename T, class Quadrature, class PDE, class GeoBasis,
          class... Basis>
class FiniteElement {
 public:
  FiniteElement(ElementMesh<GeoBasis>& geomesh, SolutionVector<T>& nodes,
                ElementMesh<Basis...>& mesh, SolutionVector<T>& solvec,
                SolutionVector<T>& resvec)
      : geomesh(geomesh),
        mesh(mesh),
        geo(geomesh, nodes),
        sol(mesh, solvec),
        res(mesh, resvec) {}

  /*
    Copy values to the solution vector
  */
  void set_solution(SolutionVector<T>& sol) {}

  void add_residual() {
    const index_t num_elements = mesh.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    for (index_t i = 0; i < num_elements; i++) {
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace s;
        sref.transform(detJ, J, Jinv, s);

        // Compute the coefficients for the weak form of the PDE
        typename PDE::FiniteElementSpace coef;
        double weight = Quadrature::get_weight(j);
        PDE::eval_weak_coef(weight * detJ, s, coef);

        // Compute the coefficients in the reference element
        typename PDE::FiniteElementSpace cref;
        coef.rtransform(detJ, J, Jinv, cref);

        // Add the contributions back to the residual
        res.template add<Quadrature>(i, j, cref);
      }
    }
  }

  /*
    Compute y = y + J * x
  */
  void add_jacobian_vector_product(SolutionVector<T>& xvec,
                                   SolutionVector<T>& yvec) {
    ElementVector<T, typename PDE::FiniteElementSpace, Basis...> x(mesh, xvec);
    ElementVector<T, typename PDE::FiniteElementSpace, Basis...> y(mesh, yvec);

    const index_t num_elements = mesh.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    for (index_t i = 0; i < num_elements; i++) {
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace xref;
        x.template interp<Quadrature>(i, j, xref);

        // Transform to the local coordinate system
        typename PDE::FiniteElementSpace s, p;
        sref.transform(detJ, J, Jinv, s);
        xref.transform(detJ, J, Jinv, p);

        // Compute the coefficients for the weak form of the PDE
        typename PDE::FiniteElementSpace coef;
        double weight = Quadrature::get_weight(j);
        PDE::eval_weak_jacobian_vec_product(weight * detJ, s, p, coef);

        // Compute the coefficients in the reference element
        typename PDE::FiniteElementSpace cref;
        coef.rtransform(detJ, J, Jinv, cref);

        // Add the contributions back to the residual
        y.template add<Quadrature>(i, j, cref);
      }
    }
  }

  void add_jacobian() {
    const index_t ndof = PDE::FiniteElementSpace::ndof;
    const index_t num_elements = mesh.get_num_elements();
    const index_t num_quadrature_points = Quadrature::get_num_points();

    for (index_t i = 0; i < num_elements; i++) {
      for (index_t j = 0; j < num_quadrature_points; j++) {
        // Extract the Jacobian of the element transformation
        typename PDE::FiniteElementGeometry gk;
        geo.template interp<Quadrature>(i, j, gk);
        A2D::Mat<T, PDE::dim, PDE::dim>& J = (gk.template get<0>()).get_grad();

        // Compute the inverse of the transformation
        A2D::Mat<T, PDE::dim, PDE::dim> Jinv;
        A2D::MatInverse(J, Jinv);

        // Compute the determinant of the Jacobian matrix
        T detJ;
        A2D::MatDet(J, detJ);

        // Interpolate the solution vector using the basis
        typename PDE::FiniteElementSpace sref;
        sol.template interp<Quadrature>(i, j, sref);

        typename PDE::FiniteElementSpace cref[ndof];
        for (index_t k = 0; k < ndof; k++) {
          typename PDE::FiniteElementSpace pref;
          pref.set_seed(k);

          // Transform to the local coordinate system
          typename PDE::FiniteElementSpace s, p;
          sref.transform(detJ, J, Jinv, s);
          pref.transform(detJ, J, Jinv, p);

          // Compute the coefficients for the weak form of the PDE
          typename PDE::FiniteElementSpace coef;
          double weight = Quadrature::get_weight(j);
          PDE::eval_weak_jacobian_vec_product(weight * detJ, s, p, coef);

          // Compute the coefficients in the reference element
          coef.rtransform(detJ, J, Jinv, cref[k]);
        }
      }
    }
  }

 private:
  // The element mesh
  ElementMesh<Basis...>& mesh;
  ElementMesh<GeoBasis>& geomesh;

  // Element-wise views of the solution and residual vector
  ElementVector<T, typename PDE::FiniteElementGeometry, GeoBasis> geo;
  ElementVector<T, typename PDE::FiniteElementSpace, Basis...> sol;
  ElementVector<T, typename PDE::FiniteElementSpace, Basis...> res;
};

int main(int argc, char* argv[]) {
  typedef double T;
  typedef MixedPoisson2D<T> PDE;
  typedef TriQuadrature3 Quadrature;
  typedef FiniteElement<T, Quadrature, PDE, LagrangeTri1<T, 2>, LagrangeTri0<T>,
                        RT2DTri1<T>>
      PoissonFE;

  // Set the node locations
  index_t nx = 10, ny = 10;
  index_t nnodes = (nx + 1) * (ny + 1);
  index_t nelems = 2 * nx * ny;

  // Create the node vector
  SolutionVector<T> nodes(2 * nnodes);

  for (index_t j = 0; j < ny + 1; j++) {
    for (index_t i = 0; i < nx + 1; i++) {
      index_t node = i + j * (nx + 1);
      nodes[2 * node] = 1.0 * i;
      nodes[2 * node + 1] = 1.0 * j;
    }
  }

  // Set the connectivity
  index_t* conn = new index_t[3 * nelems];

  for (index_t j = 0; j < ny + 1; j++) {
    for (index_t i = 0; i < nx + 1; i++) {
      index_t elem = 2 * (i + j * nx);
      conn[3 * elem] = i + j * (nx + 1);
      conn[3 * elem + 1] = i + 1 + j * (nx + 1);
      conn[3 * elem + 2] = i + 1 + (j + 1) * (nx + 1);

      elem += 1;
      conn[3 * elem] = i + j * (nx + 1);
      conn[3 * elem + 1] = i + 1 + (j + 1) * (nx + 1);
      conn[3 * elem + 2] = i + (j + 1) * (nx + 1);
    }
  }

  ElementConnectivity connect(nnodes, nelems, conn);
  delete[] conn;

  // Create the mesh for the geometry
  SpaceType geo_space[] = {H1};
  index_t dims[] = {2};
  ElementMesh<LagrangeTri1<T, 2>> geomesh(connect, geo_space, dims);

  SpaceType sol_space[] = {L2, EDGE};
  ElementMesh<LagrangeTri0<T>, RT2DTri1<T>> mesh(connect, sol_space);

  // Get the total number of degrees of freedom
  index_t ndof = mesh.get_num_dof();

  SolutionVector<T> sol(ndof), res(ndof), pert(ndof);
  PoissonFE poisson(geomesh, nodes, mesh, sol, res);

  for (index_t i = 0; i < ndof; i++) {
    pert[i] = 1.0;
  }

  poisson.add_residual();
  poisson.add_jacobian_vector_product(pert, res);
  poisson.add_jacobian();

  return (0);
}