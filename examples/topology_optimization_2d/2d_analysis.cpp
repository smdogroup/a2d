#include <vector>

#include "ParOptOptimizer.h"
#include "a2ddefs.h"
#include "multiphysics/feanalysis.h"
#include "multiphysics/febasis.h"
#include "multiphysics/feelement.h"
#include "multiphysics/feelementmat.h"
#include "multiphysics/femesh.h"
#include "multiphysics/fequadrature.h"
#include "multiphysics/hex_tools.h"
#include "multiphysics/integrand_elasticity.h"
#include "multiphysics/integrand_helmholtz.h"
#include "multiphysics/lagrange_hypercube_basis.h"
#include "sparse/sparse_cholesky.h"
#include "sparse/sparse_utils.h"
#include "utils/a2dmesh.h"

using namespace A2D;

/**
 * @brief Make the analysis look like it's using the
 *
 */
template <class FltrImpl, class AnlyImpl>
class TopoFilterAnalysis : public Analysis<AnlyImpl> {
 public:
  using T = typename AnlyImpl::type;
  using Vec_t = typename AnlyImpl::Vec_t;

  static_assert(std::is_same<Vec_t, typename FltrImpl::Vec_t>::value,
                "Vector types must be the same");

  TopoFilterAnalysis(std::shared_ptr<Analysis<FltrImpl>> filter,
                     std::shared_ptr<Analysis<AnlyImpl>> analysis)
      : filter(filter), analysis(analysis) {}

  std::shared_ptr<Vec_t> get_data() { return filter->get_data(); }
  std::shared_ptr<Vec_t> get_geo() { return analysis->get_geo(); }
  std::shared_ptr<Vec_t> get_sol() { return analysis->get_sol(); }
  std::shared_ptr<Vec_t> get_res() { return analysis->get_res(); }

  // x -> rho -> topology optimization
  void linear_solve() {
    filter->linear_solve();
    analysis->linear_solve();
  }

  void nonlinear_solve() {
    filter->linear_solve();
    analysis->nonlinear_solve();
  }

  // Evaluate a function
  T evaluate(FunctionalBase<AnlyImpl> &func) {
    return func.evaluate(*this->get_data(), *this->get_geo(), *this->get_sol());
  }

  void add_derivative(FunctionalBase<AnlyImpl> &func, FEVarType wrt, T alpha,
                      Vec_t &dfdx) {
    func.add_derivative(wrt, alpha, *this->get_data(), *this->get_geo(),
                        *this->get_sol(), dfdx);
  }

  void eval_adjoint_derivative(FunctionalBase<AnlyImpl> &func, FEVarType wrt,
                               Vec_t &dfdx) {
    if (wrt == FEVarType::GEOMETRY) {
      analysis->eval_adjoint_derivative(func, wrt, dfdx);
    } else if (wrt == FEVarType::DATA) {
      analysis->eval_adjoint_derivative(func, wrt, *filter->get_res());
      filter->eval_adjoint_derivative(wrt, dfdx);
    }
  }

  void eval_adjoint_derivative(FEVarType wrt, Vec_t &dfdx) {}

  void to_vtk(const std::string filename) {
    filter->to_vtk(filename + std::string("-filter.vtk"));
    analysis->to_vtk(filename + std::string("-elasticity.vtk"));
  }

 private:
  std::shared_ptr<Analysis<FltrImpl>> filter;
  std::shared_ptr<Analysis<AnlyImpl>> analysis;
};

template <class FltrImpl, class AnlyImpl>
class TopOptProb : public ParOptProblem {
 public:
  TopOptProb(std::string prefix, MPI_Comm comm,
             std::shared_ptr<TopoFilterAnalysis<FltrImpl, AnlyImpl>> analysis,
             std::shared_ptr<FunctionalBase<AnlyImpl>> obj,
             std::shared_ptr<FunctionalBase<AnlyImpl>> volume,
             double target_volume,
             std::shared_ptr<typename AnlyImpl::Vec_t> dfdx)
      : ParOptProblem(comm),
        prefix(prefix),
        analysis(analysis),
        obj(obj),
        volume(volume),
        dfdx(dfdx),
        target_volume(target_volume) {
    nvars = analysis->get_data()->size();
    ncon = 1;
    nineq = 1;

    if (!std::filesystem::is_directory(prefix)) {
      std::filesystem::create_directory(prefix);
    }
    setProblemSizes(nvars, ncon, 0);
    setNumInequalities(nineq, 0);

    fobj_ref = 1.0;
  }

  //! Create the quasi-def matrix associated with this problem
  ParOptQuasiDefMat *createQuasiDefMat() {
    int nwblock = 0;
    return new ParOptQuasiDefBlockMat(this, nwblock);
  }

  void getVarsAndBounds(ParOptVec *xvec, ParOptVec *lbvec, ParOptVec *ubvec) {
    ParOptScalar *x, *lb, *ub;
    xvec->getArray(&x);
    lbvec->getArray(&lb);
    ubvec->getArray(&ub);

    // Set the design variable bounds
    for (int i = 0; i < nvars; i++) {
      x[i] = 0.95;
      lb[i] = 1e-3;
      ub[i] = 1.0;
    }
  }

  // Note: constraints > 0
  int evalObjCon(ParOptVec *xvec, ParOptScalar *fobj, ParOptScalar *cons) {
    Timer timer("TopOptProb::evalObjCon()");
    // Set design variables by paropt
    ParOptScalar *x;
    xvec->getArray(&x);

    auto data = analysis->get_data();

    // Set the design variables
    for (int i = 0; i < nvars; i++) {
      (*data)[i] = x[i];
    }

    // Solve the problem
    analysis->linear_solve();

    // Evaluate objective
    *fobj = analysis->evaluate(*obj) / fobj_ref;

    // Evaluate constraint
    cons[0] = 1.0 - analysis->evaluate(*volume) / target_volume;

    // Write design to vtk
    if (opt_iter % vtk_freq == 0) {
      char vtk_name[256];
      std::snprintf(vtk_name, sizeof(vtk_name), "result_%d.vtk", opt_iter);
      std::filesystem::path path =
          std::filesystem::path(prefix) / std::filesystem::path(vtk_name);
      analysis->to_vtk(path);
    }

    opt_iter++;

    std::printf("[%2d] obj: %20.10e, con: %20.10e\n", opt_iter, *fobj, cons[0]);

    return 0;
  }

  int evalObjConGradient(ParOptVec *xvec, ParOptVec *gvec, ParOptVec **Ac) {
    Timer timer("TopOptProb::evalObjConGradient()");
    ParOptScalar *g, *c;

    // Evaluate the gradient
    dfdx->zero();
    analysis->eval_adjoint_derivative(*obj, FEVarType::DATA, *dfdx);

    // Set the gradient objective gradient
    gvec->getArray(&g);
    for (int i = 0; i < nvars; i++) {
      g[i] = (*dfdx)[i] / fobj_ref;
    }

    // Evaluate the gradient of the volume constraint
    dfdx->zero();
    analysis->add_derivative(*volume, FEVarType::DATA, -1.0 / target_volume,
                             *dfdx);

    // Set the gradient objective gradient
    Ac[0]->getArray(&c);
    for (int i = 0; i < nvars; i++) {
      c[i] = (*dfdx)[i];
    }

    return 0;
  }

 private:
  std::string prefix;
  MPI_Comm comm;
  std::shared_ptr<TopoFilterAnalysis<FltrImpl, AnlyImpl>> analysis;
  std::shared_ptr<FunctionalBase<AnlyImpl>> obj;
  std::shared_ptr<FunctionalBase<AnlyImpl>> volume;
  std::shared_ptr<typename AnlyImpl::Vec_t> dfdx;

  int nvars;
  int ncon;
  int nineq;
  ParOptScalar fobj_ref;
  ParOptScalar target_volume;

  int opt_iter;
  bool verbose;
  int vtk_freq;
};

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    MPI_Comm comm = MPI_COMM_SELF;
    using T = double;

    // Number of elements in each dimension
    const index_t degree = 1;
    const index_t nx = 8, ny = 8;
    const double lx = 2.0, ly = 2.0;
    const double target_volume = 0.4 * lx * ly;

    // Set up mesh
    const index_t nverts = (nx + 1) * (ny + 1);
    const index_t ntri = 0, nquad = nx * ny;
    index_t *tri = nullptr;
    std::vector<index_t> quad(4 * nquad);
    std::vector<double> Xloc(2 * nverts);
    MesherRect2D mesher(nx, ny, lx, ly);
    bool randomize = false;
    unsigned int seed = 0;
    double fraction = 0.2;
    mesher.set_X_conn<index_t, T>(Xloc.data(), quad.data(), randomize, seed,
                                  fraction);

    MeshConnectivity2D conn(nverts, ntri, tri, nquad, quad.data());
    // Set up bcs
    auto node_num = [](index_t i, index_t j) { return i + j * (nx + 1); };

    const index_t num_boundary_verts = (ny + 1);
    index_t boundary_verts[num_boundary_verts];

    for (index_t j = 0; j < ny + 1; j++) {
      boundary_verts[j] = node_num(0, j);
    }

    index_t bc_label =
        conn.add_boundary_label_from_verts(num_boundary_verts, boundary_verts);

    DirichletBCInfo bcinfo;
    bcinfo.add_boundary_condition(bc_label);

    // Set up the type of implementation we're using
    const index_t dim = 2;
    const index_t data_size = 1;
    const index_t block_size = dim;
    using FltrImpl_t = typename DirectCholeskyAnalysis<T, data_size>::Impl_t;
    using AnlyImpl_t = typename DirectCholeskyAnalysis<T, block_size>::Impl_t;
    using Vec_t = typename AnlyImpl_t::Vec_t;

    // Set the types of elements
    using Filter_t = QuadHelmholtzFilterElement<FltrImpl_t, degree>;

    constexpr GreenStrainType etype = GreenStrainType::LINEAR;
    using Elem_t = QuadTopoElement<AnlyImpl_t, etype, degree>;
    using BodyForce_t = QuadBodyForceTopoElement<AnlyImpl_t, degree>;

    // Create the meshes for the Hex elements
    auto data_mesh = std::make_shared<ElementMesh<Elem_t::DataBasis>>(conn);
    auto geo_mesh = std::make_shared<ElementMesh<Elem_t::GeoBasis>>(conn);
    auto sol_mesh = std::make_shared<ElementMesh<Elem_t::Basis>>(conn);

    // Get the number of different dof for the analysis
    index_t ndata = data_mesh->get_num_dof();
    index_t ngeo = geo_mesh->get_num_dof();
    index_t ndof = sol_mesh->get_num_dof();

    // Create the data vector for the filter
    auto filter_data = std::make_shared<Vec_t>(ndata);

    // Create the solution vector for the filter - same as the analysis data
    // vector
    auto filter_sol = std::make_shared<Vec_t>(ndata);
    auto filter_res = std::make_shared<Vec_t>(ndata);

    // Derivative of the function of interest
    auto dfdx = std::make_shared<Vec_t>(ndata);

    // Create the geometry vector - same for both filter/analysis
    auto geo = std::make_shared<Vec_t>(ngeo);

    // Create the solution/residual vector for the analysis
    auto sol = std::make_shared<Vec_t>(ndof);
    auto res = std::make_shared<Vec_t>(ndof);

    // Set the data
    for (index_t i = 0; i < ndata; i++) {
      (*filter_data)[i] = 1.0;
    }

    // Set the geometry
    typename AnlyImpl_t::ElementVector<Elem_t::GeoBasis> elem_geo(*geo_mesh,
                                                                  *geo);
    set_geo_from_quad_nodes<Elem_t::GeoBasis>(nquad, quad.data(), Xloc.data(),
                                              elem_geo);

    // Create the filter
    T length = 1.0;
    T r0 = 0.05 * length;
    HelmholtzFilter<T, dim> filter_integrand(r0);

    auto filer_assembler = std::make_shared<ElementAssembler<FltrImpl_t>>();
    filer_assembler->add_element(std::make_shared<Filter_t>(
        filter_integrand, data_mesh, geo_mesh, data_mesh));

    auto filter = std::make_shared<DirectCholeskyAnalysis<T, 1>>(
        filter_data, geo, filter_sol, filter_res, filer_assembler);

    // Create the element integrand
    T E = 70.0, nu = 0.3, q = 5.0;
    TopoElasticityIntegrand<T, dim, etype> elem_integrand(E, nu, q);

    // Create the body force integrand
    T tx[] = {0.0, 10.0};
    TopoBodyForceIntegrand<T, dim> body_integrand(q, tx);

    auto assembler = std::make_shared<ElementAssembler<AnlyImpl_t>>();
    assembler->add_element(std::make_shared<Elem_t>(elem_integrand, data_mesh,
                                                    geo_mesh, sol_mesh));
    assembler->add_element(std::make_shared<BodyForce_t>(
        body_integrand, data_mesh, geo_mesh, sol_mesh));

    // Set up the boundary conditions
    auto bcs = std::make_shared<DirichletBCs<T>>();
    bcs->add_bcs(std::make_shared<DirichletBasis<T, Elem_t::Basis>>(
        conn, *sol_mesh, bcinfo, 0.0));

    // Create the assembler object
    auto analysis = std::make_shared<DirectCholeskyAnalysis<T, block_size>>(
        filter_sol, geo, sol, res, assembler, bcs);

    T design_stress = 100.0, ks_param = 0.01;
    using Func_t = QuadTopoVonMises<AnlyImpl_t, etype, degree>;
    TopoVonMisesKS<T, dim, etype> func_integrand(E, nu, q, design_stress,
                                                 ks_param);
    auto functional =
        std::make_shared<Func_t>(func_integrand, data_mesh, geo_mesh, sol_mesh);

    using Volume_t = QuadTopoVolume<AnlyImpl_t, degree>;
    TopoVolume<T, dim> vol_integrand;
    auto volume = std::make_shared<Volume_t>(vol_integrand, data_mesh, geo_mesh,
                                             sol_mesh);

    auto topo = std::make_shared<TopoFilterAnalysis<FltrImpl_t, AnlyImpl_t>>(
        filter, analysis);

    // Set up the topology optimization problem
    std::string prefix("./");
    TopOptProb<FltrImpl_t, AnlyImpl_t> prob(prefix, comm, topo, functional,
                                            volume, target_volume, dfdx);

    prob.checkGradients(1e-6);

    // // Evaluate the function and its derivative
    // topo->linear_solve();
    // T f0 = topo->evaluate(*functional);
    // topo->eval_adjoint_derivative(*functional, FEVarType::DATA, *dfdx);

    // // Perturb x and test the gradient
    // double dh = 1e-6;
    // T dfdp = 0.0;
    // for (int i = 0; i < ndata; i++) {
    //   dfdp += (*dfdx)[i];
    //   (*filter_data)[i] += dh;
    // }

    // topo->linear_solve();
    // T f1 = topo->evaluate(*functional);
    // T fd = (f1 - f0) / dh;

    // std::cout << "dfdp:  " << std::setw(20) << std::setprecision(12) << dfdp
    //           << std::endl;
    // std::cout << "fd:    " << std::setw(20) << std::setprecision(12) << fd
    //           << std::endl;
    // std::cout << "err:   " << std::setw(20) << std::setprecision(12)
    //           << (fd - dfdp) / fd << std::endl;
  }

  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}