#include <iostream>
#include <memory>

#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/integrand_heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/integrand_poisson.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "multiphysics/static_condensation.h"
#include "sparse/sparse_amg.h"

using namespace A2D;

/**
 * @brief Mixed Helmholtz problem
 *
 * @tparam T Scalar type for the calculation
 * @tparam D Dimension of the problem
 */
template <typename T, A2D::index_t D>
class MixedHelmholtz {
 public:
  MixedHelmholtz(T r) : r(r) {}

  T r;  // The Helmholtz radius for the filter

  // Spatial dimension
  static const A2D::index_t dim = D;

  // No data associated with this element
  static const A2D::index_t data_dim = 0;

  // Space for the finite-element data
  typedef A2D::FESpace<T, data_dim> DataSpace;

  // Finite element space
  typedef A2D::FESpace<T, dim, A2D::HdivSpace<T, dim>, A2D::L2Space<T, 1, dim>>
      FiniteElementSpace;

  // Space for the element geometry - parametrized by H1 in 2D
  typedef A2D::FESpace<T, dim, A2D::H1Space<T, dim, dim>> FiniteElementGeometry;

  // The type of matrix used to store data at each quadrature point
  typedef A2D::SymMat<T, FiniteElementSpace::ncomp> QMatType;

  // Mapping of the solution from the reference element to the physical element
  using SolutionMapping = A2D::InteriorMapping<T, dim>;

  /**
   * @brief Evaluate the weak form of the coefficients
   *
   * v * u + <tau, sigma> - r * v * div(sigma) - r * div(tau) * u - v * f
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param dobj The data at the quadrature point
   * @param geo The geometry evaluated at the current point
   * @param s The trial solution
   * @param coef Derivative of the weak form w.r.t. coefficients
   */
  A2D_INLINE_FUNCTION void weak(T wdetJ, const DataSpace& dobj,
                                const FiniteElementGeometry& geo,
                                const FiniteElementSpace& s,
                                FiniteElementSpace& coef) {
    // Field objects for solution functions
    const A2D::HdivSpace<T, dim>& sigma = s.template get<0>();
    const A2D::L2Space<T, 1, dim>& u = s.template get<1>();

    // Solution function values
    const A2D::Vec<T, dim>& sigma_val = sigma.get_value();
    const T& sigma_div = sigma.get_div();
    const T& u_val = u.get_value();

    // Test function values
    A2D::HdivSpace<T, dim>& tau = coef.template get<0>();
    A2D::L2Space<T, 1, dim>& v = coef.template get<1>();

    // Test function values
    A2D::Vec<T, dim>& tau_val = tau.get_value();
    T& tau_div = tau.get_div();
    T& v_val = v.get_value();

    // Set the terms from the variational statement
    for (A2D::index_t k = 0; k < dim; k++) {
      tau_val(k) = wdetJ * sigma_val(k);
    }

    const A2D::Vec<T, dim>& x = (geo.template get<0>()).get_value();
    T r0 = std::sqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2));
    T f = 1.0 - r0;

    // Use f = 1.0 for now...
    v_val = -wdetJ * (u_val + r * sigma_div - f);
    tau_div = -r * wdetJ * u_val;
  }

  /**
   * @brief Construct a JacVecProduct functor
   *
   * This functor computes a Jacobian-vector product of the
   *
   * @param wdetJ The quadrature weight times determinant of the Jacobian
   * @param data The data at the quadrature point
   * @param s The solution at the quadrature point
   */
  class JacVecProduct {
   public:
    A2D_INLINE_FUNCTION JacVecProduct(const MixedHelmholtz<T, D>& integrand,
                                      T wdetJ, const DataSpace& data,
                                      const FiniteElementGeometry& geo,
                                      const FiniteElementSpace& s)
        : wdetJ(wdetJ), r(integrand.r) {}

    A2D_INLINE_FUNCTION void operator()(const FiniteElementSpace& p,
                                        FiniteElementSpace& Jp) {
      // Field objects for solution functions
      const A2D::HdivSpace<T, dim>& sigma = p.template get<0>();
      const A2D::L2Space<T, 1, dim>& u = p.template get<1>();

      // Solution function values
      const A2D::Vec<T, dim>& sigma_val = sigma.get_value();
      const T& sigma_div = sigma.get_div();
      const T& u_val = u.get_value();

      // Test function values
      A2D::HdivSpace<T, dim>& tau = Jp.template get<0>();
      A2D::L2Space<T, 1, dim>& v = Jp.template get<1>();

      // Test function values
      A2D::Vec<T, dim>& tau_val = tau.get_value();
      T& tau_div = tau.get_div();
      T& v_val = v.get_value();

      // Set the terms from the variational statement
      for (A2D::index_t k = 0; k < dim; k++) {
        tau_val(k) = wdetJ * sigma_val(k);
      }

      v_val = -wdetJ * (u_val + r * sigma_div);
      tau_div = -r * wdetJ * u_val;
    }

   private:
    T wdetJ, r;
  };
};

// - (2 *x^2 - 2 * x - e^(1-x) - 4 * e^x + e^(1 + x) + 4)/(2 * x)
template <typename T, A2D::index_t D>
class HelmholtzForSphereError {
 public:
  // Spatial dimension
  static const A2D::index_t dim = D;
  static const A2D::index_t data_dim = 0;

  using DataSpace = typename MixedHelmholtz<T, D>::DataSpace;
  using FiniteElementSpace = typename MixedHelmholtz<T, D>::FiniteElementSpace;
  using FiniteElementGeometry =
      typename MixedHelmholtz<T, D>::FiniteElementGeometry;
  using SolutionMapping = typename MixedHelmholtz<T, D>::SolutionMapping;

  HelmholtzForSphereError() {}

  T integrand(T wdetJ, const DataSpace& data, const FiniteElementGeometry& geo,
              const FiniteElementSpace& s) {
    T u = (s.template get<1>()).get_value();

    // Compute the right-hand-side
    const A2D::Vec<T, dim>& xpt = (geo.template get<0>()).get_value();
    T x = std::sqrt(xpt(0) * xpt(0) + xpt(1) * xpt(1) + xpt(2) * xpt(2));
    if (std::real(x) == 0) {
      x = 1e-15;
    }
    // T inv = 0.5 / r;
    // T u0 = -inv * (2.0 * r * r - 2.0 * r - std::exp(1.0 - r) -
    //                4.0 * std::exp(r) + std::exp(r + 1.0) + 4.0);

    // y(x) = (e^(-32 x) (-33 e^(32 x) (512 x^2 - 512 x + 1) - 31 e^(32 (x + 2))
    // (512 x^2 - 512 x + 1) + 33 e^(64 x) + 511 e^(64 x + 32) + 31 e^64 - 511
    // e^32))/(512 (33 + 31 e^64) x)

    T u0 = (std::exp(-32.0 * x) *
            (-33.0 * std::exp(32 * x) * (512.0 * x * x - 512.0 * x + 1.0) -
             31.0 * std::exp(32 * (x + 2)) * (512.0 * x * x - 512.0 * x + 1.0) +
             33.0 * std::exp(64 * x) + 511 * std::exp(64 * x + 32) +
             31.0 * std::exp(64.0) - 511 * std::exp(32.0))) /
           (512.0 * (33.0 + 31.0 * std::exp(64.0) * x));

    T err = (u - u0);

    return wdetJ * err * err;
  }

  double R;
};

template <A2D::index_t nx, class GeoBasis, typename T, class GeoElemVec>
void set_geo_spherical(const T alpha, const T R, GeoElemVec& elem_geo) {
  using ET = A2D::ElementTypes;

  T a = alpha * R;
  T len = 2.0 * a / nx;
  T ulen = 2.0 / nx;

  T Xloc[3 * ET::HEX_NVERTS];
  for (A2D::index_t ii = 0; ii < ET::HEX_NVERTS; ii++) {
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
  for (A2D::index_t face = 0; face < ET::HEX_NBOUNDS; face++) {
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
              x0 += N[jj] * Xloc[3 * ET::HEX_BOUND_VERTS[face][jj]];
              y0 += N[jj] * Xloc[3 * ET::HEX_BOUND_VERTS[face][jj] + 1];
              z0 += N[jj] * Xloc[3 * ET::HEX_BOUND_VERTS[face][jj] + 2];
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

template <typename T, A2D::index_t degree>
class HelmholtzSphere {
 public:
  // Alias templates
  template <class... Args>
  using ElementVector = A2D::ElementVector_Serial<Args...>;
  using ElementVectorEmpty = A2D::ElementVector_Empty<A2D::ElemVecType::Serial>;

  // Basic types
  using I = A2D::index_t;

  // Magic integers
  static constexpr I spatial_dim = 3;  // spatial dimension
  static constexpr I var_dim =
      1;  // dimension of the Integrand solution variable
  static constexpr I data_dim = 1;    // dimension of material data
  static constexpr I low_degree = 1;  // low order preconditinoer mesh degree
  static constexpr I block_size = var_dim;  // block size for BSR matrix

  // Offset for static condenseation
  static constexpr I static_condensation_basis_offset = 1;

  // Null-space size
  static constexpr I null_size = 1;

  // Problem Integrand
  using Integrand = MixedHelmholtz<T, spatial_dim>;

  // The type of solution vector to use
  using BasisVecType = A2D::SolutionVector<T>;

  // Quadrature, basis and element views for original mesh
  using Quadrature = A2D::HexGaussQuadrature<degree + 1>;
  using DataBasis = A2D::FEBasis<T>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = FEBasis<T, QHdivHexBasis<T, degree>,
                        LagrangeL2HexBasis<T, 1, degree - 1>>;
  using DataElemVec = ElementVectorEmpty;
  using GeoElemVec = ElementVector<T, GeoBasis, BasisVecType>;
  using ElemVec = ElementVector<T, Basis, BasisVecType>;

  // Quadrature, basis and element views for low order preconditioner mesh
  using LOrderQuadrature = A2D::HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis = A2D::FEBasis<T>;
  using LOrderGeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, low_degree>>;
  using LOrderBasis = FEBasis<T, QHdivHexBasis<T, low_degree>,
                              LagrangeL2HexBasis<T, 1, low_degree - 1>>;
  using LOrderDataElemVec = ElementVectorEmpty;
  using LOrderGeoElemVec = ElementVector<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = ElementVector<T, LOrderBasis, BasisVecType>;

  // FE type
  using FE_PDE =
      A2D::FiniteElement<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis>;

  // Finite element functional for low order preconditioner mesh
  using LOrderFE =
      A2D::FiniteElement<T, Integrand, LOrderQuadrature, LOrderDataBasis,
                         LOrderGeoBasis, LOrderBasis>;

  // Matrix-free operator
  using MatFree =
      A2D::MatrixFree<T, Integrand, Quadrature, DataBasis, GeoBasis, Basis>;

  // Static condensation matrix information
  using SDMatType =
      StaticCondensationMat<T, block_size, static_condensation_basis_offset,
                            LOrderBasis>;

  HelmholtzSphere(T r, A2D::MeshConnectivity3D& conn,
                  A2D::DirichletBCInfo& bcinfo, int amg_nlevels, int cg_it,
                  double cg_rtol, double cg_atol,
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

        // Store the static condensation matrix
        elem_mat(lorder_mesh),

        integrand(r),

        amg_nlevels(amg_nlevels),
        cg_it(cg_it),
        cg_rtol(cg_rtol),
        cg_atol(cg_atol),
        verbose(verbose) {}

  A2D::index_t get_num_dof() { return sol.get_num_dof(); }

  GeoElemVec& get_geometry() { return elem_geo; }
  void reset_geometry() {}

  /**
   * @brief Solve the governing equations and set the new solution vector
   *
   */
  void solve() {
    // Initialie the Jacobian matrix
    lorder_fe.add_jacobian(integrand, lorder_elem_data, lorder_elem_geo,
                           lorder_elem_sol, elem_mat);

    // Apply Boundary conditions to the matrix
    elem_mat.zero_bcs(bcs);

    // Factor the matrix
    elem_mat.factor();

    // Initialize the matrix-free data
    matfree.initialize(integrand, elem_data, elem_geo, elem_sol);

    // Allocate space for temporary variables with the matrix-vector code
    A2D::SolutionVector<T> xvec(mesh.get_num_dof());
    A2D::SolutionVector<T> yvec(mesh.get_num_dof());
    ElemVec elem_xvec(mesh, xvec);
    ElemVec elem_yvec(mesh, yvec);

    // Apply the boundary conditions
    const I* bc_dofs;
    I nbcs = bcs.get_bcs(&bc_dofs);

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

    // Create the solution and right-hand-side vectors
    I size = sol.get_num_dof() / block_size;
    A2D::MultiArrayNew<T* [block_size]> sol_vec("sol_vec", size);
    A2D::MultiArrayNew<T* [block_size]> rhs_vec("rhs_vec", size);

    // Zero the solution
    sol.zero();

    // Assemble the body force contribution
    A2D::SolutionVector<T> res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    fe.add_residual(integrand, elem_data, elem_geo, elem_sol, elem_res);

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
    const I gmres_size = 50;
    I nrestart = 10;
    bool succ = A2D::fgmres<T, block_size, gmres_size>(
        mat_vec,
        [&](A2D::MultiArrayNew<T* [block_size]>& in,
            A2D::MultiArrayNew<T* [block_size]>& out) -> void {
          elem_mat.apply_factor(in, out);
        },
        rhs_vec, sol_vec, monitor, nrestart, cg_rtol, cg_atol);
    if (!succ) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "%s:%d: GMRES failed to converge after %d iterations given "
                    "rtol=%.1e, atol=%.1e",
                    __FILE__, __LINE__, cg_it, cg_rtol, cg_atol);
      std::cout << msg << std::endl;
      // throw std::runtime_error(msg);
    }

    // Record the solution
    for (I i = 0; i < sol.get_num_dof(); i++) {
      sol[i] = sol_vec(i / block_size, i % block_size);
    }
  }

  T compute_solution_error() {
    // Compute and print out the solution error
    HelmholtzForSphereError<T, spatial_dim> error;
    A2D::FiniteElement<T, HelmholtzForSphereError<T, spatial_dim>,
                       A2D::HexGaussQuadrature<degree + 10>, DataBasis,
                       GeoBasis, Basis>
        functional;

    T err = functional.integrate(error, elem_data, elem_geo, elem_sol);
    return std::sqrt(err);
  }

  void tovtk(const std::string filename) {
    A2D::write_hex_to_vtk<5, degree, T, DataBasis, GeoBasis, Basis>(
        integrand, elem_data, elem_geo, elem_sol, filename,
        [](I k, typename Integrand::DataSpace& d,
           typename Integrand::FiniteElementGeometry& g,
           typename Integrand::FiniteElementSpace& s) { return s[k]; });
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

  // Matrix type
  SDMatType elem_mat;

  Integrand integrand;
  FE_PDE fe;

  LOrderFE lorder_fe;
  MatFree matfree;

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
  for (A2D::index_t face = 0; face < ET::HEX_NBOUNDS; face++) {
    for (A2D::index_t i = 0; i < 4; i++) {
      hex[(face + 1) * ET::HEX_NVERTS + i] = hex[ET::HEX_BOUND_VERTS[face][i]];
      hex[(face + 1) * ET::HEX_NVERTS + i + 4] =
          hex_ext[ET::HEX_BOUND_VERTS[face][i]];
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

  double r = 1.0 / 32.0;  // Set r = 1.0 for the error computation
  HelmholtzSphere<T, degree> sphere(r, conn, bcinfo, amg_nlevels, cg_it,
                                    cg_rtol, cg_atol);

  // Set the geometry from the node locations
  auto elem_geo = sphere.get_geometry();
  set_geo_spherical<nx, typename HelmholtzSphere<T, degree>::GeoBasis>(
      alpha, R, elem_geo);
  sphere.reset_geometry();

  // Solve the spherical problem
  sphere.solve();

  T error = sphere.compute_solution_error();

  if (nx == 1 && degree == 1) {
    std::cout << std::setw(10) << "p" << std::setw(10) << "nx" << std::setw(10)
              << "ndof" << std::setw(25) << "error" << std::endl;
  }
  A2D::index_t ndof = sphere.get_num_dof();
  std::cout << std::setw(10) << degree << std::setw(10) << nx << std::setw(10)
            << ndof << std::setw(25) << std::setprecision(16) << error
            << std::endl;

  if (write_sphere) {
    sphere.tovtk("filename.vtk");
  }
}

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  find_spherical_error<2, 1>();
  find_spherical_error<3, 1>();
  find_spherical_error<4, 1>();
  find_spherical_error<5, 1>();
  find_spherical_error<6, 1>();
  find_spherical_error<7, 1>();
  find_spherical_error<8, 1>();
  find_spherical_error<9, 1>();
  find_spherical_error<10, 1>();
  find_spherical_error<11, 1>();
  find_spherical_error<12, 1>();
  find_spherical_error<13, 1>();
  find_spherical_error<14, 1>();
  find_spherical_error<15, 1>();
  find_spherical_error<16, 1>();
  find_spherical_error<17, 1>();
  find_spherical_error<18, 1>();
  find_spherical_error<19, 1>();
  find_spherical_error<20, 1>();

  find_spherical_error<1, 2>();
  find_spherical_error<2, 2>();
  find_spherical_error<3, 2>();
  find_spherical_error<4, 2>();
  find_spherical_error<5, 2>();
  find_spherical_error<6, 2>();
  find_spherical_error<7, 2>();
  find_spherical_error<8, 2>();
  find_spherical_error<9, 2>();
  find_spherical_error<10, 2>();

  find_spherical_error<1, 3>();
  find_spherical_error<2, 3>();
  find_spherical_error<3, 3>();
  find_spherical_error<4, 3>();
  find_spherical_error<5, 3>();
  find_spherical_error<6, 3>();
  find_spherical_error<7, 3>();

  find_spherical_error<1, 4>();
  find_spherical_error<2, 4>();
  find_spherical_error<3, 4>();
  find_spherical_error<4, 4>();
  find_spherical_error<5, 4>(true);

  // find_spherical_error<1, 5>();
  // find_spherical_error<2, 5>();
  // find_spherical_error<3, 5>();
  // find_spherical_error<4, 5>();

  // find_spherical_error<5, 4>();
  // find_spherical_error<6, 4>(true);

  return (0);
}