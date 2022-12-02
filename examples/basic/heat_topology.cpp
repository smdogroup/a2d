#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "multiphysics/elasticity.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/heat_conduction.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/lagrange_hex_basis.h"
#include "multiphysics/poisson.h"
#include "multiphysics/qhdiv_hex_basis.h"
#include "sparse/sparse_amg.h"
#include "utils/a2dprofiler.h"

template <typename T, A2D::index_t degree>
class HeatTopoOpt {
 public:
  static const A2D::index_t spatial_dim = 3;
  static const A2D::index_t var_dim = 1;
  static const A2D::index_t data_dim = 1;
  using PDE = A2D::HeatConduction<T, spatial_dim>;

  template <class... Args>
  using ElementVector = A2D::ElementVector_Serial<Args...>;

  // The type of solution vector to use
  using BasisVecType = A2D::SolutionVector<T>;

  using Quadrature = A2D::HexGaussQuadrature<degree + 1>;
  using DataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, degree - 1>>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using DataElemVec = ElementVector<T, DataBasis, BasisVecType>;
  using GeoElemVec = ElementVector<T, GeoBasis, BasisVecType>;
  using ElemVec = ElementVector<T, Basis, BasisVecType>;
  using FE = A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  // Matrix-free operator for the problem
  using MatFree =
      A2D::MatrixFree<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  static const A2D::index_t low_degree = 1;
  using LOrderQuadrature = A2D::HexGaussQuadrature<low_degree + 1>;
  using LOrderDataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, low_degree - 1>>;
  using LOrderGeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, low_degree>>;
  using LOrderBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, low_degree>>;
  using LOrderDataElemVec = ElementVector<T, LOrderDataBasis, BasisVecType>;
  using LOrderGeoElemVec = ElementVector<T, LOrderGeoBasis, BasisVecType>;
  using LOrderElemVec = ElementVector<T, LOrderBasis, BasisVecType>;
  using LOrderFE = A2D::FiniteElement<T, PDE, LOrderQuadrature, LOrderDataBasis,
                                      LOrderGeoBasis, LOrderBasis>;

  // Block-size for the finite-element problem
  static const A2D::index_t block_size = var_dim;
  static const A2D::index_t null_size = 1;
  using BSRMatType = A2D::BSRMat<A2D::index_t, T, block_size, block_size>;
  using BSRMatAmgType = A2D::BSRMatAmg<A2D::index_t, T, block_size, null_size>;

  // Functional definitions
  using VolumeFunctional =
      A2D::FiniteElement<T, A2D::TopoVolume<T, var_dim, spatial_dim>,
                         Quadrature, DataBasis, GeoBasis, Basis>;

  HeatTopoOpt(A2D::MeshConnectivity3D &conn, A2D::DirichletBCInfo &bcinfo,
              T kappa, T q)
      :  // Material parameters and penalization
        kappa(kappa),
        q(q),

        // Meshes for the solution, geometry and data
        mesh(conn),
        geomesh(conn),
        datamesh(conn),

        bcs(conn, mesh, bcinfo),

        // Project the meshes onto the low-order meshes
        lorder_mesh(mesh),
        lorder_geomesh(geomesh),
        lorder_datamesh(datamesh),

        // Solution, geometry and data vectors
        sol(mesh.get_num_dof()),
        geo(geomesh.get_num_dof()),
        data(datamesh.get_num_dof()),

        // Element-level views of the solution geometry and data
        elem_data(datamesh, data),
        elem_geo(geomesh, geo),
        elem_sol(mesh, sol),

        // Low-order views of the solution geometry and data
        lorder_elem_data(lorder_datamesh, data),
        lorder_elem_geo(lorder_geomesh, geo),
        lorder_elem_sol(lorder_mesh, sol),

        B("B", sol.get_num_dof() / block_size) {
    // Initialize the design variable
    for (A2D::index_t i = 0; i < data.get_num_dof(); i++) {
      data[i] = 1.0;
    }

    // Create the matrix for the low-order mesh
    A2D::index_t nrows;
    std::vector<A2D::index_t> rowp, cols;
    lorder_mesh.template create_block_csr<block_size>(nrows, rowp, cols);

    // Create the shared pointer
    mat = std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);
  }

  GeoElemVec &get_geometry() { return elem_geo; }

  void reset_geometry() {
    A2D::SolutionVector<T> x(mesh.get_num_dof());
    A2D::SolutionVector<T> y(mesh.get_num_dof());
    A2D::SolutionVector<T> z(mesh.get_num_dof());

    ElemVec elem_x(mesh, x);
    ElemVec elem_y(mesh, y);
    ElemVec elem_z(mesh, z);

    A2D::DOFCoordinates<T, PDE, GeoBasis, Basis> coords;
    coords.get_dof_coordinates(elem_geo, elem_x, elem_y, elem_z);

    // Initialize the near null-space to an appropriate vector
    A2D::BLAS::fill(B, 1.0);

    // Zero out the boundary conditions
    const A2D::index_t *bc_dofs;
    A2D::index_t nbcs = bcs.get_bcs(&bc_dofs);
    for (A2D::index_t i = 0; i < nbcs; i++) {
      A2D::index_t dof = bc_dofs[i];
      for (A2D::index_t j = 0; j < null_size; j++) {
        B(dof / block_size, dof % block_size, j) = 0.0;
      }
    }
  }

  /**
   * @brief Solve the governing equations and set the new solution vector
   *
   */
  void solve() {
    A2D::Timer timer("TopoOpt::solve()");
    // Create a view of the low-order element matrix
    A2D::ElementMat_Serial<T, LOrderBasis, BSRMatType> elem_mat(lorder_mesh,
                                                                *mat);

    // Initialie the Jacobian matrix
    lorder_fe.add_jacobian(pde, lorder_elem_data, lorder_elem_geo,
                           lorder_elem_sol, elem_mat);

    // Apply the boundary conditions
    const A2D::index_t *bc_dofs;
    A2D::index_t nbcs = bcs.get_bcs(&bc_dofs);
    mat->zero_rows(nbcs, bc_dofs);

    // Initialize the matrix-free data
    matfree.initialize(pde, elem_data, elem_geo, elem_sol);

    // Allocate space for temporary variables with the matrix-vector code
    A2D::SolutionVector<T> xvec(mesh.get_num_dof());
    A2D::SolutionVector<T> yvec(mesh.get_num_dof());
    ElemVec elem_xvec(mesh, xvec);
    ElemVec elem_yvec(mesh, yvec);

    auto mat_vec = [&](A2D::MultiArrayNew<T *[block_size]> &in,
                       A2D::MultiArrayNew<T *[block_size]> &out) -> void {
      xvec.zero();
      yvec.zero();
      for (A2D::index_t i = 0; i < xvec.get_num_dof(); i++) {
        xvec[i] = in(i / block_size, i % block_size);
      }
      matfree.add_jacobian_vector_product(elem_xvec, elem_yvec);

      for (A2D::index_t i = 0; i < yvec.get_num_dof(); i++) {
        out(i / block_size, i % block_size) = yvec[i];
      }

      // Set the boundary conditions as equal to the inputs
      const A2D::index_t *bc_dofs;
      A2D::index_t nbcs = bcs.get_bcs(&bc_dofs);
      for (A2D::index_t i = 0; i < nbcs; i++) {
        A2D::index_t dof = bc_dofs[i];

        out(dof / block_size, dof % block_size) =
            in(dof / block_size, dof % block_size);
      }
    };

    // Allocate the solver - we should add some of these as solver options
    A2D::index_t num_levels = 3;
    double omega = 4.0 / 3.0;
    double epsilon = 0.0;
    bool print_info = true;
    BSRMatAmgType amg(num_levels, omega, epsilon, mat, B, print_info);

    // Create the solution and right-hand-side vectors
    A2D::index_t size = sol.get_num_dof() / block_size;
    A2D::MultiArrayNew<T *[block_size]> sol_vec("sol_vec", size);
    A2D::MultiArrayNew<T *[block_size]> rhs_vec("rhs_vec", size);

    A2D::SolutionVector<T> res(mesh.get_num_dof());
    ElemVec elem_res(mesh, res);
    sol.zero();
    fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

    for (A2D::index_t i = 0; i < sol.get_num_dof(); i++) {
      rhs_vec(i / block_size, i % block_size) = -res[i];
    }

    // Zero out the boundary conditions
    for (A2D::index_t i = 0; i < nbcs; i++) {
      A2D::index_t dof = bc_dofs[i];
      rhs_vec(dof / block_size, dof % block_size) = 0.0;
    }

    // Solve the problem
    // amg.applyFactor(rhs_vec, sol_vec);
    amg.cg(mat_vec, rhs_vec, sol_vec, 5, 100);

    // Record the solution
    for (A2D::index_t i = 0; i < sol.get_num_dof(); i++) {
      sol[i] = sol_vec(i / block_size, i % block_size);
    }
  }

  /**
   * @brief Get the number of design variables
   *
   * @return The number of design variables
   */
  A2D::index_t get_num_design_vars() { return datamesh.get_num_dof(); }

  /**
   * @brief Set the design variable values into the topology
   *
   * @tparam VecType The type of design vector
   * @param xvec The vector of design variable values
   */
  template <class VecType>
  void set_design_vars(const VecType &xvec) {
    for (A2D::index_t i = 0; i < datamesh.get_num_dof(); i++) {
      data[i] = xvec[i];
    }
  }

  /**
   * @brief Evaluate the volume
   *
   * @return The volume of the parametrized topology
   */
  T eval_volume() {
    A2D::Timer timer("TopoOpt::eval_volume()");
    VolumeFunctional functional;
    A2D::TopoVolume<T, var_dim, spatial_dim> volume;
    T vol = functional.integrate(volume, elem_data, elem_geo, elem_sol);

    return vol;
  }

  /**
   * @brief Add the volume gradient to the specified functional
   *
   * @tparam VecType The derivative vector type
   * @param dfdx The derivative value
   */
  template <class VecType>
  void add_volume_gradient(VecType &dfdx) {
    A2D::Timer timer("TopoOpt::add_volume_gradient()");
    VolumeFunctional functional;
    A2D::TopoVolume<T, var_dim, spatial_dim> volume;

    // Create the element-view for the derivative
    ElementVector<T, DataBasis, VecType> elem_dfdx(datamesh, dfdx);

    functional.add_data_derivative(volume, elem_data, elem_geo, elem_sol,
                                   elem_dfdx);
  }

  void tovtk(const std::string filename) {
    A2D::write_hex_to_vtk<2, degree, T, DataBasis, GeoBasis, Basis>(
        pde, elem_data, elem_geo, elem_sol, filename,
        [](A2D::index_t k, typename PDE::DataSpace &d,
           typename PDE::FiniteElementGeometry &g,
           typename PDE::FiniteElementSpace &s) {
          if (k == 0) {
            return (s.template get<0>()).get_value();  // state
          } else {
            return (d.template get<0>()).get_value();  // design
          }
        });
  }

 private:
  T kappa, q;

  A2D::ElementMesh<Basis> mesh;
  A2D::ElementMesh<GeoBasis> geomesh;
  A2D::ElementMesh<DataBasis> datamesh;

  A2D::DirichletBCs<Basis> bcs;

  A2D::ElementMesh<LOrderBasis> lorder_mesh;
  A2D::ElementMesh<LOrderGeoBasis> lorder_geomesh;
  A2D::ElementMesh<LOrderDataBasis> lorder_datamesh;

  A2D::SolutionVector<T> sol;
  A2D::SolutionVector<T> geo;
  A2D::SolutionVector<T> data;

  DataElemVec elem_data;
  GeoElemVec elem_geo;
  ElemVec elem_sol;

  LOrderDataElemVec lorder_elem_data;
  LOrderGeoElemVec lorder_elem_geo;
  LOrderElemVec lorder_elem_sol;

  PDE pde;
  FE fe;
  LOrderFE lorder_fe;
  MatFree matfree;

  // The near null-space to an appropriate vector
  A2D::MultiArrayNew<T *[block_size][null_size]> B;

  // System matrix
  std::shared_ptr<BSRMatType> mat;
};

void test_analysis(int argc, char *argv[]) {
  // Magic numbers
  const int degree = 1;            // polynomial degree, 1 = linear
  const int order = 2;             // order = degree + 1
  const int spatial_dim = 3;       // spatial dimension
  const int var_dim = 1;           // solution variable dimension
  const int data_dim = 1;          // design variable dimension
  const int block_size = var_dim;  // block size for BSR matrix

  /* Types */

  // Basic type
  using T = double;
  using I = A2D::index_t;

  // Quadrature and basis
  using Quadrature = A2D::HexGaussQuadrature<order>;
  using Basis = A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, var_dim, degree>>;
  using GeoBasis =
      A2D::FEBasis<T, A2D::LagrangeH1HexBasis<T, spatial_dim, degree>>;
  using DataBasis =
      A2D::FEBasis<T, A2D::LagrangeL2HexBasis<T, data_dim, degree - 1>>;

  // Block sparse compressed row matrix
  using BSRMatType = A2D::BSRMat<I, T, block_size, block_size>;

  // Element-centric views
  using GlobalVecType = A2D::SolutionVector<T>;
  using ElemVec = A2D::ElementVector_Serial<T, Basis, GlobalVecType>;
  using GeoElemVec = A2D::ElementVector_Serial<T, GeoBasis, GlobalVecType>;
  using DataElemVec = A2D::ElementVector_Serial<T, DataBasis, GlobalVecType>;
  using ElemMat = A2D::ElementMat_Serial<T, Basis, BSRMatType>;

  // Physics and functional
  using PDE = A2D::HeatConduction<T, spatial_dim>;
  using FE = A2D::FiniteElement<T, PDE, Quadrature, DataBasis, GeoBasis, Basis>;

  /* Load mesh and boundary vertices from vtk */

  // Load vtk
  std::string vtk_name = "3d_hex_small.vtk";
  if (argc > 1) {
    vtk_name = argv[1];
  }
  A2D::ReadVTK3D<I, T> readvtk(vtk_name);

  // Get connectivity for each element type
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;

  // Construct connectivity
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Extract boundary vertices
  std::vector<int> ids{100, 101};
  std::vector<I> bc_verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Construct A2D mesh objects
  A2D::ElementMesh<Basis> mesh(conn);
  A2D::ElementMesh<GeoBasis> geomesh(conn);
  A2D::ElementMesh<DataBasis> datamesh(conn);

  /* Construct solution, data and X vectors and their element views */

  // Allocate vectors
  GlobalVecType sol(mesh.get_num_dof());
  GlobalVecType res(mesh.get_num_dof());
  GlobalVecType data(datamesh.get_num_dof());
  GlobalVecType geo(geomesh.get_num_dof());

  // Allocate elem views
  ElemVec elem_sol(mesh, sol);
  ElemVec elem_res(mesh, res);
  DataElemVec elem_data(datamesh, data);
  GeoElemVec elem_geo(geomesh, geo);

  // Populate data and geo
  data.fill(1.0);
  A2D::set_geo_from_hex_nodes<GeoBasis>(nhex, hex, Xloc, elem_geo);

  /* Construct the Jacobian matrix */

  // Construct the nonzero pattern
  I nrows;
  std::vector<I> rowp, cols;
  mesh.create_block_csr<block_size>(nrows, rowp, cols);
  std::shared_ptr<BSRMatType> mat =
      std::make_shared<BSRMatType>(nrows, nrows, cols.size(), rowp, cols);

  // Populate the stiffness matrix via element view
  ElemMat elem_mat(mesh, *mat);
  FE fe;
  PDE pde;
  fe.add_jacobian(pde, elem_data, elem_geo, elem_sol, elem_mat);

  mat->write_mtx("heat_jacobian_nobc.mtx");

  // Apply boundary condition
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(bc_verts.size(), bc_verts.data()));
  A2D::DirichletBCs<Basis> bcs(conn, mesh, bcinfo);
  const I *bc_dofs;
  I nbcs = bcs.get_bcs(&bc_dofs);
  mat->zero_rows(nbcs, bc_dofs);

  // Write the Jacobian matrix to mtx for visualization
  mat->write_mtx("heat_jacobian_withbc.mtx");

  /* Construct the residual and rhs */

  fe.add_residual(pde, elem_data, elem_geo, elem_sol, elem_res);

  I size = sol.get_num_dof() / block_size;
  A2D::MultiArrayNew<T *[block_size]> sol_vec("sol_vec", size);
  A2D::MultiArrayNew<T *[block_size]> rhs_vec("rhs_vec", size);

  printf("rhs, before applying bcs:\n");
  for (I i = 0; i < sol.get_num_dof(); i++) {
    rhs_vec(i / block_size, i % block_size) = -res[i];
    printf("rhs[%d] = %20.10f\n", i, rhs_vec(i / block_size, i % block_size));
  }

  // Set the temperature boundary conditions
  for (I i = 0; i < nbcs; i++) {
    I dof = bc_dofs[i];
    rhs_vec(dof / block_size, dof % block_size) = 0.5;
  }

  printf("rhs, after applying bcs:\n");
  for (I i = 0; i < sol.get_num_dof(); i++) {
    printf("rhs[%d] = %20.10f\n", i, rhs_vec(i / block_size, i % block_size));
  }

  /* Construct the AMG solver and solve the problem */

  // Create amg
  A2D::index_t num_levels = 3;
  double omega = 4.0 / 3.0;
  double epsilon = 0.0;
  bool print_info = false;
  const int null_size = 1;
  A2D::MultiArrayNew<T *[block_size][null_size]> B(
      "B", sol.get_num_dof() / block_size);
  A2D::BLAS::fill(B, 1.0);
  A2D::BSRMatAmg<I, T, block_size, null_size> amg(num_levels, omega, epsilon,
                                                  mat, B, print_info);

  // Solve
  auto mat_vec = [&](A2D::MultiArrayNew<T *[block_size]> &in,
                     A2D::MultiArrayNew<T *[block_size]> &out) -> void {
    A2D::BSRMatVecMult<I, T, block_size, block_size>(*mat, in, out);
  };
  amg.cg(mat_vec, rhs_vec, sol_vec, 5, 100);

  // Record the solution
  for (I i = 0; i < sol.get_num_dof(); i++) {
    sol[i] = sol_vec(i / block_size, i % block_size);
  }

  // Write result to vtk
  A2D::write_hex_to_vtk<2, degree, T, DataBasis, GeoBasis, Basis>(
      pde, elem_data, elem_geo, elem_sol, "heat_analysis.vtk",
      [](A2D::index_t k, typename PDE::DataSpace &d,
         typename PDE::FiniteElementGeometry &g,
         typename PDE::FiniteElementSpace &s) {
        if (k == 0) {
          return (s.template get<0>()).get_value();  // state
        } else {
          return (d.template get<0>()).get_value();  // design
        }
      });
}

void main_body(int argc, char *argv[]) {
  using I = A2D::index_t;
  using T = double;
  constexpr int DEGREE = 8;

  // Load connectivity and vertex coordinates from vtk
  std::string vtk_name = "3d_hex_small.vtk";
  if (argc > 1) {
    vtk_name = argv[1];
  }
  A2D::ReadVTK3D<I, T> readvtk(vtk_name);
  T *Xloc = readvtk.get_Xloc();
  I nverts = readvtk.get_nverts();
  I nhex = readvtk.get_nhex();
  I *hex = readvtk.get_hex();
  I ntets = 0, nwedge = 0, npyrmd = 0;
  I *tets = nullptr, *wedge = nullptr, *pyrmd = nullptr;

  // Construct connectivity
  A2D::MeshConnectivity3D conn(nverts, ntets, tets, nhex, hex, nwedge, wedge,
                               npyrmd, pyrmd);

  // Extract boundary vertices
  std::vector<int> ids{100, 101};
  std::vector<I> verts = readvtk.get_verts_given_cell_entity_id(ids);

  // Set up boundary condition information
  A2D::index_t basis = 0;
  A2D::DirichletBCInfo bcinfo;
  bcinfo.add_boundary_condition(
      conn.add_boundary_label_from_verts(verts.size(), verts.data()), basis);

  // Initialize the analysis instance
  T kappa = 1.0, q = 5.0;
  HeatTopoOpt<T, DEGREE> topo(conn, bcinfo, kappa, q);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<HeatTopoOpt<T, DEGREE>::GeoBasis>(nhex, hex, Xloc,
                                                                elem_geo);
  topo.reset_geometry();

  // Initialize the optimization object
  int nvars = topo.get_num_design_vars();
  int ncon = 1;
  int nineq = 1;
  topo.set_design_vars(std::vector<T>(nvars, T(1.0)));
  topo.solve();

  topo.tovtk("heat_topo.vtk");
}

int main(int argc, char *argv[]) {
  Kokkos::initialize();
  {
    A2D::HeatConduction<double, 3> heat;
    A2D::TestPDEImplementation<double>(heat);
    // test_analysis(argc, argv);
    main_body(argc, argv);
  }
  Kokkos::finalize();
  return 0;
}