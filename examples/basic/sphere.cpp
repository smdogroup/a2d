#include <iostream>
#include <memory>
#include <string>

#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/lagrange_quad_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dprofiler.h"

template <A2D::index_t nx, class GeoBasis, typename T, class GeoElemVec>
void set_geo_spherical(const T alpha, const T R, GeoElemVec& elem_geo) {
  using ET = A2D::ElementTypes;

  T a = alpha * R;
  T len = 2.0 * a / nx;
  T ulen = 2.0 / nx;

  T Xloc[3 * ET::HEX_VERTS];
  for (A2D::index_t ii = 0; ii < ET::HEX_VERTS; ii++) {
    T u = 1.0 * ET::HEX_VERTS_CART[ii][0];
    T v = 1.0 * ET::HEX_VERTS_CART[ii][1];
    T w = 1.0 * ET::HEX_VERTS_CART[ii][2];

    Xloc[3 * ii] = (2.0 * u - 1.0) * a;
    Xloc[3 * ii + 1] = (2.0 * v - 1.0) * a;
    Xloc[3 * ii + 2] = (2.0 * w - 1.0) * a;
  }

  for (A2D::index_t k = 0; k < nx; k++) {
    for (A2D::index_t j = 0; j < nx; j++) {
      for (A2D::index_t i = 0; i < nx; i++) {
        A2D::index_t elem = i + nx * (j + nx * k);
        typename GeoElemVec::FEDof geo_dof(elem, elem_geo);

        for (A2D::index_t ii = 0; ii < GeoBasis::ndof; ii++) {
          double pt[3];
          GeoBasis::get_dof_point(ii, pt);

          T x = (0.5 * pt[0] + 0.5) * len + i * len - a;
          T y = (0.5 * pt[1] + 0.5) * len + j * len - a;
          T z = (0.5 * pt[2] + 0.5) * len + k * len - a;

          // Interpolate to find the basis
          if (ii % 3 == 0) {
            geo_dof[ii] = x;
          } else if (ii % 3 == 1) {
            geo_dof[ii] = y;
          } else if (ii % 3 == 2) {
            geo_dof[ii] = z;
          }
        }

        elem_geo.set_element_values(elem, geo_dof);
      }
    }
  }

  // We know the 2-direction will be radial...
  for (A2D::index_t face = 0; face < ET::HEX_FACES; face++) {
    for (A2D::index_t k = 0; k < nx; k++) {
      for (A2D::index_t j = 0; j < nx; j++) {
        for (A2D::index_t i = 0; i < nx; i++) {
          A2D::index_t elem = i + nx * (j + nx * k) + (face + 1) * nx * nx * nx;
          typename GeoElemVec::FEDof geo_dof(elem, elem_geo);

          for (A2D::index_t ii = 0; ii < GeoBasis::ndof; ii++) {
            double pt[3];
            GeoBasis::get_dof_point(ii, pt);

            T u = ulen * (0.5 * pt[0] + 0.5) + ulen * i - 1.0;
            T v = ulen * (0.5 * pt[1] + 0.5) + ulen * j - 1.0;
            T w = ulen * (0.5 * pt[2] + 0.5) + ulen * k - 1.0;

            T N[4];
            N[0] = 0.25 * (1.0 - u) * (1.0 - v);
            N[1] = 0.25 * (1.0 + u) * (1.0 - v);
            N[2] = 0.25 * (1.0 + u) * (1.0 + v);
            N[3] = 0.25 * (1.0 - u) * (1.0 + v);

            // Find the x-y-z coordinate on the surface
            T x0 = 0.0;
            T y0 = 0.0;
            T z0 = 0.0;

            for (A2D::index_t jj = 0; jj < 4; jj++) {
              x0 += N[jj] * Xloc[3 * ET::HEX_FACE_VERTS[face][jj]];
              y0 += N[jj] * Xloc[3 * ET::HEX_FACE_VERTS[face][jj] + 1];
              z0 += N[jj] * Xloc[3 * ET::HEX_FACE_VERTS[face][jj] + 2];
            }
            T norm = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);

            T x1 = x0 * R / norm;
            T y1 = y0 * R / norm;
            T z1 = z0 * R / norm;

            // Project
            T x = 0.5 * (1.0 - w) * x0 + 0.5 * (1.0 + w) * x1;
            T y = 0.5 * (1.0 - w) * y0 + 0.5 * (1.0 + w) * y1;
            T z = 0.5 * (1.0 - w) * z0 + 0.5 * (1.0 + w) * z1;

            // Interpolate to find the basis
            if (ii % 3 == 0) {
              geo_dof[ii] = x;
            } else if (ii % 3 == 1) {
              geo_dof[ii] = y;
            } else if (ii % 3 == 2) {
              geo_dof[ii] = z;
            }
          }

          elem_geo.set_element_values(elem, geo_dof);
        }
      }
    }
  }
}

template <typename T, A2D::index_t D>
class PoissonForSphere {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;
  static const A2D::index_t data_dim = 0;

  using DataSpace = A2D::FESpace<T, data_dim>;
  using FiniteElementSpace = A2D::FESpace<T, dim, A2D::H1Space<T, 1, dim>>;
  using FiniteElementGeometry = A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>>;
  using QMatType = A2D::SymmMat<T, FiniteElementSpace::ncomp>;
  using SolutionMapping = A2D::InteriorMapping<T, dim>;

  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& dobj,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    const A2D::H1Space<T, 1, dim>& u = s.template get<0>();
    const A2D::Vec<T, dim>& u_grad = u.get_grad();

    A2D::H1Space<T, 1, dim>& v = coef.template get<0>();
    A2D::Vec<T, dim>& v_grad = v.get_grad();

    // Compute the right-hand-side
    T& v_value = v.get_value();
    const A2D::Vec<T, dim>& x = (geo.template get<0>()).get_value();
    T r = std::sqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2));
    v_value = -wdetJ * r * r * r * r;

    // Set the terms from the variational statement
    for (A2D::index_t k = 0; k < dim; k++) {
      v_grad(k) = wdetJ * u_grad(k);
    }
  }

  /**
   * @brief Construct the JacVecProduct functor
   *
   * This functor computes a Jacobian-vector product of the weak form
   *
   * @param pde The PDE object for this class
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param geo The geometry at the quadrature point
   * @param s The solution at the quadrature point
   */
  class JacVecProduct {
   public:
    A2D_INLINE_FUNCTION JacVecProduct(const PoissonForSphere<T, D>& pde,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        : wdetJ(wdetJ) {}

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      const A2D::H1Space<T, 1, dim>& u = p.template get<0>();
      const A2D::Vec<T, dim>& u_grad = u.get_grad();

      A2D::H1Space<T, 1, dim>& v = Jp.template get<0>();
      A2D::Vec<T, dim>& v_grad = v.get_grad();

      // Set the terms from the variational statement
      for (A2D::index_t k = 0; k < dim; k++) {
        v_grad(k) = wdetJ * u_grad(k);
      }
    }

   private:
    T wdetJ;
  };
};

template <typename T, A2D::index_t D>
class PoissonForSphereError {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;
  static const A2D::index_t data_dim = 0;

  using DataSpace = typename PoissonForSphere<T, D>::DataSpace;
  using FiniteElementSpace =
      typename PoissonForSphere<T, D>::FiniteElementSpace;
  using FiniteElementGeometry =
      typename PoissonForSphere<T, D>::FiniteElementGeometry;
  using SolutionMapping = typename PoissonForSphere<T, D>::SolutionMapping;

  PoissonForSphereError(T R = 1.0) : R(R) {}

  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) {
    T u = (s.template get<0>()).get_value();

    // Compute the right-hand-side
    const A2D::Vec<T, dim>& x = (geo.template get<0>()).get_value();
    T r = std::sqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2));
    T u0 = (r * r * r * r * r * r - R * R * R * R * R * R) / 42.0;

    T err = (u - u0);

    return wdetJ * err * err;
  }

  double R;
};

template <typename T, A2D::index_t degree>
class PoissonSphere {
 public:
  // Alias templates
  template <class... Args>
  using ElementVector = A2D::ElementVector_Serial<Args...>;

  // Basic types
  using I = A2D::index_t;

  // Magic integers
  static constexpr I spatial_dim = 3;  // spatial dimension
  static constexpr I var_dim = 1;      // dimension of the PDE solution variable
  static constexpr I data_dim = 1;     // dimension of material data
  static constexpr I low_degree = 1;   // low order preconditinoer mesh degree
  static constexpr I block_size = var_dim;  // block size for BSR matrix

  // Problem PDE
  using PDE = PoissonForSphere<T, spatial_dim>;

  // The type of solution vector to use
  using BasisVecType = A2D::SolutionVector<T>;

  // Quadrature, basis and element views for original mesh
  using Quadrature = A2D::HexGaussQuadrature<degree + 1>;
  using DataBasis = A2D::FEBasis<T>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using DataElemVec = A2D::EmptyElementVector;
  using GeoElemVec = ElementVector<T, GeoBasis, BasisVecType>;
  using ElemVec = ElementVector<T, Basis, BasisVecType>;

  // Quadrature, basis and element views for low order preconditioner mesh
  using LOrderQuadrature = A2D::HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = A2D::FEBasis<T>;
  using LOrderGeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, low_degree>>;
  using LOrderBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, low_degree>>;
  using LOrderDataElemVec = A2D::EmptyElementVector;
  using LOrderGeoElemVec = ElementVector<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = ElementVector<T, LOrderBasis, BasisVecType>;

  // FE type
  using FE_PDE =
      A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Finite element functional for low order preconditioner mesh
  using LOrderFE = A2D::FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                      LOrderGeoBasis, LOrderBasis>;

  // Block compressed row sparse matrix
  using BSRMatType = A2D::BSRMat<I, T, block_size, block_size>;

  // Matrix-free operator
  using MatFree =
      A2D::MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Algebraic multigrid solver
  static constexpr I null_size = 1;
  using BSRMatAmgType = A2D::BSRMatAmg<I, T, block_size, null_size>;

  PoissonSphere(A2D::MeshConnectivity3D& conn, A2D::DirichletBCInfo& bcinfo,
                int amg_nlevels, int cg_it, double cg_rtol, double cg_atol,
                bool verbose = false)
      :  // Meshes for the solution, geometry and data
        mesh(conn),
        geomesh(conn),

        // Boundary conditions
        bcs(conn, mesh, bcinfo),

        // Project the meshes onto the low-order meshes
        lorder_mesh(mesh),
        lorder_geomesh(geomesh),

        // Solution, geometry and data vectors
        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),

        // Element-level views of the solution geometry and data
        elem_sol(mesh, sol),
        elem_geo(geomesh, geo),

        // Low-order views of the solution geometry and data
        lorder_elem_sol(lorder_mesh, sol),
        lorder_elem_geo(lorder_geomesh, geo),

        B("B", sol.get_num_dof() / block_size),

        amg_nlevels(amg_nlevels),
        cg_it(cg_it),
        cg_rtol(cg_rtol),
        cg_atol(cg_atol),
        verbose(verbose) {
    // Create the matrix for the low-order mesh
    I nrows;
    std::vector<I> rowp, cols;
    lorder_mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    // Create the shared pointer
    mat = std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);

    // Initialize the near null-space to an appropriate vector
    for (I i = 0; i < B.extent(0); i++) {
      B(i, 0, 0) = 1.0;
    }

    // Zero out the boundary conditions
    const I* bc_dofs;
    I nbcs = bcs.get_bcs(&bc_dofs);
    for (I i = 0; i < nbcs; i++) {
      I dof = bc_dofs[i];
      for (I j = 0; j < null_size; j++) {
        B(dof / block_size, dof % block_size, j) = 0.0;
      }
    }
  }

  A2D::index_t get_num_dof() { return sol.get_num_dof(); }

  GeoElemVec& get_geometry() { return elem_geo; }
  void reset_geometry() {}

  /**
   * @brief Solve the governing equations and set the new solution vector
   *
   */
  void solve() {
    // Create a view of the low-order element matrix
    A2D::ElementMat_Serial<T, LOrderBasis, BSRMatType> elem_mat(lorder_mesh,
                                                                *mat);

    // Initialie the Jacobian matrix
    lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                           lorder_elem_sol, elem_mat);

    // Apply the boundary conditions
    const I* bc_dofs;
    I nbcs = bcs.get_bcs(&bc_dofs);
    mat->zero_rows(nbcs, bc_dofs);

    // Initialize the matrix-free data
    matfree.initialize(pde, elem_data, elem_geo, elem_sol);

    // Allocate space for temporary variables with the matrix-vector code
    A2D::SolutionVector<T> xvec(mesh.get_num_dof());
    A2D::SolutionVector<T> yvec(mesh.get_num_dof());
    ElemVec elem_xvec(mesh, xvec);
    ElemVec elem_yvec(mesh, yvec);

    auto mat_vec = [&](A2D::MultiArrayNew<T* [block_size]>& in,
                       A2D::MultiArrayNew<T* [block_size]>& out) -> void {
      xvec.zero();
      yvec.zero();
      for (I i = 0; i < xvec.get_num_dof(); i++) {
        xvec[i] = in(i / block_size, i % block_size);
      }
      matfree.add_jacobian_vector_product(elem_xvec, elem_yvec);

      for (I i = 0; i < yvec.get_num_dof(); i++) {
        out(i / block_size, i % block_size) = yvec[i];
      }

      // Set the boundary conditions as equal to the inputs
      const I* bc_dofs;
      I nbcs = bcs.get_bcs(&bc_dofs);
      for (I i = 0; i < nbcs; i++) {
        I dof = bc_dofs[i];

        out(dof / block_size, dof % block_size) =
            in(dof / block_size, dof % block_size);
      }
    };

    // Allocate the solver - we should add some of these as solver options
    double omega = 4.0 / 3.0;
    double epsilon = 0.0;
    bool print_info = false;
    if (verbose) {
      print_info = true;
    }
    BSRMatAmgType amg(amg_nlevels, omega, epsilon, mat, B, print_info);

    // Create the solution and right-hand-side vectors
    I size = sol.get_num_dof() / block_size;
    A2D::MultiArrayNew<T* [block_size]> sol_vec("sol_vec", size);
    A2D::MultiArrayNew<T* [block_size]> rhs_vec("rhs_vec", size);

    // Zero the solution
    sol.zero();

    // Assemble the body force contribution
    A2D::SolutionVector<T> res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

    // Set the right-hand-side
    for (I i = 0; i < sol.get_num_dof(); i++) {
      rhs_vec(i / block_size, i % block_size) = -res[i];
    }

    // Zero out the boundary conditions
    for (A2D::index_t i = 0; i < nbcs; i++) {
      A2D::index_t dof = bc_dofs[i];
      rhs_vec(dof / block_size, dof % block_size) = 0.0;
    }

    // Solve the problem
    I monitor = 0;
    if (verbose) {
      monitor = 5;
    }
    bool succ =
        amg.cg(mat_vec, rhs_vec, sol_vec, monitor, cg_it, cg_rtol, cg_atol);
    if (!succ) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "%s:%d: CG failed to converge after %d iterations given "
                    "rtol=%.1e, atol=%.1e",
                    __FILE__, __LINE__, cg_it, cg_rtol, cg_atol);
      throw std::runtime_error(msg);
    }

    // Record the solution
    for (I i = 0; i < sol.get_num_dof(); i++) {
      sol[i] = -sol_vec(i / block_size, i % block_size);
    }
  }

  T compute_solution_error() {
    // Compute and print out the solution error
    double R = 1.0;
    PoissonForSphereError<T, spatial_dim> error(R);
    A2D::FiniteElement<T, PoissonForSphereError<T, spatial_dim>,
                       A2D::HexGaussQuadrature<degree + 10>, DataBasis,
                       GeoBasis, Basis>
        functional;

    T err = functional.integrate(error, elem_data, elem_geo, elem_sol);
    return std::sqrt(err);
  }

  void tovtk(const std::string filename) {
    A2D::write_hex_to_vtk<1, degree, T, DataBasis, GeoBasis, Basis>(
        pde, elem_data, elem_geo, elem_sol, filename,
        [](I k, typename PDE::DataSpace& d,
           typename PDE::FiniteElementGeometry& g,
           typename PDE::FiniteElementSpace& s) { return s[0]; });
  }

 private:
  double alpha, R;

  T E, nu, q;

  A2D::ElementMesh<Basis> mesh;
  A2D::ElementMesh<GeoBasis> geomesh;
  A2D::DirichletBCs<Basis> bcs;

  A2D::ElementMesh<LOrderBasis> lorder_mesh;
  A2D::ElementMesh<LOrderGeoBasis> lorder_geomesh;

  A2D::SolutionVector<T> sol;
  A2D::SolutionVector<T> geo;

  ElemVec elem_sol;
  GeoElemVec elem_geo;
  DataElemVec elem_data;

  LOrderElemVec lorder_elem_sol;
  LOrderGeoElemVec lorder_elem_geo;
  LOrderDataElemVec lorder_elem_data;

  PDE pde;
  FE_PDE fe;

  LOrderFE lorder_fe;
  MatFree matfree;

  // The near null-space to an appropriate vector
  A2D::MultiArrayNew<T* [block_size][null_size]> B;

  // System matrix
  std::shared_ptr<BSRMatType> mat;

  // AMG settings
  int amg_nlevels;

  // CG settings
  int cg_it;
  double cg_rtol, cg_atol;
  bool verbose;
};

template <A2D::index_t nx, A2D::index_t degree>
void find_spherical_error(bool write_sphere = false) {
  using ET = A2D::ElementTypes;
  using T = double;

  // Set up the connectivity between
  double alpha = 0.3;  // Value < 1
  double R = 1.0;      // Radius of the sphere

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
  A2D::MeshConnectivity3D conn_coarse(nverts, ntets, tets, nhex, hex, nwedge,
                                      wedge, npyrmd, pyrmd);

  // Set up a fake basis that will never be evaluated to create the mesh from
  // the spherical elements
  using FakeBasis =
      A2D::FEBasis<double, A2D::LagrangeH1HexBasis<double, 1, nx>>;
  using FakeLOrderBasis = A2D::FEBasis<
      double, typename A2D::LagrangeH1HexBasis<double, 1, nx>::LOrderBasis>;
  A2D::ElementMesh<FakeBasis> fake_mesh(conn_coarse);

  // Extract the low-order element connectivity
  A2D::ElementMesh<FakeLOrderBasis> lorder_mesh(fake_mesh);

  int nverts_refine = lorder_mesh.get_num_dof();
  const int coord_to_hex[] = {0, 1, 3, 2, 4, 5, 7, 6};

  int nhex_refine = lorder_mesh.get_num_elements();
  int* hex_refine = new int[8 * nhex_refine];
  for (int i = 0; i < nhex_refine; i++) {
    for (int j = 0; j < 8; j++) {
      hex_refine[8 * i + coord_to_hex[j]] =
          lorder_mesh.template get_global_dof<0>(i, j);
    }
  }

  // Create the full mesh connectivity
  A2D::MeshConnectivity3D conn(nverts_refine, ntets, tets, nhex_refine,
                               hex_refine, nwedge, wedge, npyrmd, pyrmd);

  // Set the label for the entire surface
  A2D::index_t bc_label = 0;
  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(bc_label, basis);

  // Add the boundary
  int amg_nlevels = 3;
  if (degree > 20) {
    amg_nlevels = 4;
  }
  int cg_it = 300;
  double cg_rtol = 1e-14;
  double cg_atol = 1e-30;

  PoissonSphere<T, degree> sphere(conn, bcinfo, amg_nlevels, cg_it, cg_rtol,
                                  cg_atol);

  // Set the geometry from the node locations
  auto elem_geo = sphere.get_geometry();
  set_geo_spherical<nx, typename PoissonSphere<T, degree>::GeoBasis>(alpha, R,
                                                                     elem_geo);
  sphere.reset_geometry();

  // Find best p and nx to give similar dof at the end
  // std::printf("degree: %2d, nx: %2d, ndof: %10d\n", degree, nx,
  //             sphere.get_num_dof());
  // return;

  // Solve the spherical problem
  Kokkos::Timer timer;
  sphere.solve();
  double elapsed_time = timer.seconds();

  T error = sphere.compute_solution_error();

  if (nx == 1 && degree == 1) {
    std::printf("%10s%10s%10s%20s%25s\n", "p", "nx", "ndof", "elapsed_time(s)",
                "error");
  }
  A2D::index_t ndof = sphere.get_num_dof();
  std::printf("%10d%10d%10d%20.5e%25.15e\n", degree, nx, ndof, elapsed_time,
              error);

  if (write_sphere) {
    sphere.tovtk("filename.vtk");
  }
}

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  find_spherical_error<1, 1>();
  find_spherical_error<5, 1>();
  find_spherical_error<10, 1>();
  find_spherical_error<15, 1>();
  find_spherical_error<20, 1>();
  find_spherical_error<25, 1>();
  find_spherical_error<30, 1>();
  find_spherical_error<35, 1>();
  find_spherical_error<40, 1>();
  find_spherical_error<45, 1>();
  find_spherical_error<50, 1>();
  find_spherical_error<54, 1>();

  find_spherical_error<1, 2>();
  find_spherical_error<4, 2>();
  find_spherical_error<8, 2>();
  find_spherical_error<12, 2>();
  find_spherical_error<16, 2>();
  find_spherical_error<20, 2>();
  find_spherical_error<24, 2>();
  find_spherical_error<27, 2>();

  find_spherical_error<1, 4>();
  find_spherical_error<3, 4>();
  find_spherical_error<5, 4>();
  find_spherical_error<7, 4>();
  find_spherical_error<9, 4>();
  find_spherical_error<11, 4>();
  find_spherical_error<13, 4>();

  find_spherical_error<1, 6>();
  find_spherical_error<3, 6>();
  find_spherical_error<5, 6>();
  find_spherical_error<7, 6>();
  find_spherical_error<9, 6>();

  find_spherical_error<1, 8>();
  find_spherical_error<2, 8>();
  find_spherical_error<3, 8>();
  find_spherical_error<4, 8>();
  find_spherical_error<5, 8>();
  find_spherical_error<6, 8>();
  find_spherical_error<7, 8>();

  find_spherical_error<1, 10>();
  find_spherical_error<2, 10>();
  find_spherical_error<3, 10>();
  find_spherical_error<4, 10>();
  find_spherical_error<5, 10>(true);

  return 0;
}