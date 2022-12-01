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
        topo(topo) {}

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
    cons[0] = volume_frac * topo.eval_volume() / domain_vol - 1.0;

    std::printf("evalObjCon()\n");
    for (int i = 0; i < 3; i++) {
      std::printf("x[%2d] = %20.10e\n", i, x[i]);
    }

    std::printf("obj: %20.10f\n", *fobj);
    std::printf("con: %20.10f\n", cons[0]);

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
      c[i] = volume_frac * c[i] / domain_vol;
    }

    std::printf("evalObjConGradient()\n");
    for (int i = 0; i < 3; i++) {
      std::printf("x[%2d] = %20.10e\n", i, x[i]);
    }

    for (int i = 0; i < 3; i++) {
      std::printf("g[%2d] = %20.10e\n", i, g[i]);
    }

    for (int i = 0; i < 3; i++) {
      std::printf("c[%2d] = %20.10e\n", i, c[i]);
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
  TopoOpt<T, DEGREE> topo;
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
  TopOptProb prob(MPI_COMM_WORLD, nvars, ncon, nineq, domain_vol, volume_frac,
                  topo);

  prob.checkGradients(1e-3);

  // Everything's fine here, but not within checkGradients()
  topo.set_design_vars(std::vector<T>(nvars, 0.951));
  topo.solve();
  topo.solve();
  topo.solve();
  T obj = topo.eval_compliance();
  std::vector<T> g(nvars, 0.0);
  topo.add_compliance_gradient(g);

  std::printf("Outside TopOptProb:\n");
  std::printf("obj: %20.10f\n", obj);
  for (int i = 0; i < 3; i++) {
    std::printf("g[%2d] = %20.10e\n", i, g[i]);
  }

  // Write the problem to a vtk file
  topo.tovtk("filename.vtk");
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