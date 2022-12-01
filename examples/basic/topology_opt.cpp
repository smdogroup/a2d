#include "ParOptOptimizer.h"
#include "topology.h"
#include "utils/a2dvtk.h"

using I = A2D::index_t;
using T = ParOptScalar;
constexpr int DEGREE = 6;

class TopOptProb : public ParOptProblem {
 public:
  TopOptProb(MPI_Comm comm, int nvars, int ncon, int nineq, T domain_vol,
             T volume_frac, TopoOpt<T, DEGREE> &topo)
      : ParOptProblem(comm, nvars, ncon, nineq, 0, 0),
        comm(comm),
        nvars(nvars),
        ncon(ncon),
        nineq(nineq),
        domain_vol(domain_vol),
        volume_frac(volume_frac),
        topo(topo),
        opt_iter(0) {}

  void getVarsAndBounds(ParOptVec *xvec, ParOptVec *lbvec, ParOptVec *ubvec) {
    T *x, *lb, *ub;
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
  int evalObjCon(ParOptVec *xvec, T *fobj, T *cons) {
    // Set design variables by paropt
    T *x;
    xvec->getArray(&x);
    topo.set_design_vars(x);

    topo.solve();

    // Evaluate objective
    *fobj = topo.eval_compliance();

    // Evaluate constraints
    cons[0] = 1.0 - volume_frac * topo.eval_volume() / domain_vol;

    // Write design to vtk
    char name[256];
    std::snprintf(name, sizeof(name), "result_%d.vtk", opt_iter);
    topo.tovtk(name);

    opt_iter++;

    return 0;
  }

  int evalObjConGradient(ParOptVec *xvec, ParOptVec *gvec, ParOptVec **Ac) {
    T *x, *g, *c;

    // Set design variables
    xvec->getArray(&x);

    // topo.set_design_vars(x);
    // topo.solve();

    // Evaluate objective gradient
    gvec->zeroEntries();
    gvec->getArray(&g);

    topo.add_compliance_gradient(g);

    // Evaluate constraint gradient
    Ac[0]->zeroEntries();
    Ac[0]->getArray(&c);
    topo.add_volume_gradient(c);

    // Scale the volume constraint
    for (int i = 0; i < nvars; i++) {
      c[i] = -volume_frac * c[i] / domain_vol;
    }

    return 0;
  }

 private:
  MPI_Comm comm;
  int nvars;
  int ncon;
  int nineq;
  T domain_vol;
  T volume_frac;
  TopoOpt<T, DEGREE> &topo;
  int opt_iter;
};

void main_body() {
  // Load connectivity and vertex coordinates from vtk
  A2D::ReadVTK3D<I, T> readvtk("3d_hex_twisted.vtk");
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
  T E = 70.0e3, nu = 0.3, q = 5.0;
  TopoOpt<T, DEGREE> topo(conn, bcinfo, E, nu, q);
  auto elem_geo = topo.get_geometry();
  A2D::set_geo_from_hex_nodes<TopoOpt<T, DEGREE>::GeoBasis>(nhex, hex, Xloc,
                                                            elem_geo);
  topo.reset_geometry();

  // Initialize the optimization object
  int nvars = topo.get_num_design_vars();
  int ncon = 1;
  int nineq = 1;
  topo.set_design_vars(std::vector<T>(nvars, T(1.0)));
  T domain_vol = topo.eval_volume();
  T volume_frac = 0.4;
  TopOptProb *prob = new TopOptProb(MPI_COMM_WORLD, nvars, ncon, nineq,
                                    domain_vol, volume_frac, topo);
  prob->incref();

  // Sanity check
  // prob->checkGradients(1e-6);

  // Define paropt optimizer
  ParOptOptions *options = new ParOptOptions;
  options->incref();
  ParOptOptimizer::addDefaultOptions(options);

  options->setOption("algorithm", "mma");
  options->setOption("mma_asymptote_contract", 0.7);
  options->setOption("mma_asymptote_relax", 1.2);
  options->setOption("mma_bound_relax", 0);
  options->setOption("mma_delta_regularization", 1e-05);
  options->setOption("mma_eps_regularization", 0.001);
  options->setOption("mma_infeas_tol", 1e-05);
  options->setOption("mma_init_asymptote_offset", 0.25);
  options->setOption("mma_l1_tol", 1e-06);
  options->setOption("mma_linfty_tol", 1e-06);
  options->setOption("mma_max_asymptote_offset", 10);
  options->setOption("mma_max_iterations", 50);
  options->setOption("mma_min_asymptote_offset", 0.01);
  options->setOption("mma_use_constraint_linearization", true);

  ParOptOptimizer *opt = new ParOptOptimizer(prob, options);
  opt->incref();
  opt->optimize();

  // Delete objects
  prob->decref();
  options->incref();
  opt->decref();
  return;
}

int main(int argc, char *argv[]) {
  A2D::Timer timer("main()");
  MPI_Init(&argc, &argv);
  Kokkos::initialize();
  { main_body(); }
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}