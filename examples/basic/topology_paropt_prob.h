#ifndef TOPOLOGY_PAROPT_PROB_H
#define TOPOLOGY_PAROPT_PROB_H

#include "ParOptOptimizer.h"

template <A2D::index_t degree>
class TopOptProb : public ParOptProblem {
 public:
  TopOptProb(MPI_Comm comm, int nvars, int ncon, int nineq,
             ParOptScalar ref_comp, ParOptScalar domain_vol,
             ParOptScalar volume_frac, TopoAnalysis<ParOptScalar, degree> &topo)
      : ParOptProblem(comm, nvars, ncon, nineq, 0, 0),
        comm(comm),
        nvars(nvars),
        ncon(ncon),
        nineq(nineq),
        ref_comp(ref_comp),
        domain_vol(domain_vol),
        volume_frac(volume_frac),
        topo(topo),
        opt_iter(0) {}

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
    // Set design variables by paropt
    ParOptScalar *x;
    xvec->getArray(&x);
    topo.set_design_vars(x);

    topo.solve();

    // Evaluate objective
    *fobj = topo.eval_compliance() / ref_comp;

    // Evaluate constraidegree
    cons[0] = volume_frac - topo.eval_volume() / domain_vol;

    // Write design to vtk
    char name[256];
    std::snprintf(name, sizeof(name), "result_%d.vtk", opt_iter);
    topo.tovtk(name);

    opt_iter++;

    std::printf("[%2d]obj: %20.10e, con: %20.10e\n", opt_iter, *fobj, cons[0]);

    return 0;
  }

  int evalObjConGradient(ParOptVec *xvec, ParOptVec *gvec, ParOptVec **Ac) {
    ParOptScalar *x, *g, *c;

    // Set design variables
    xvec->getArray(&x);

    // Evaluate objective gradient
    gvec->zeroEntries();
    gvec->getArray(&g);

    topo.add_compliance_gradient(g);

    // Evaluate constraint gradient
    Ac[0]->zeroEntries();
    Ac[0]->getArray(&c);
    topo.add_volume_gradient(c);

    // Scale gradients
    for (int i = 0; i < nvars; i++) {
      g[i] = g[i] / ref_comp;
      c[i] = -c[i] / domain_vol;
    }

    return 0;
  }

 private:
  MPI_Comm comm;
  int nvars;
  int ncon;
  int nineq;
  ParOptScalar ref_comp;
  ParOptScalar domain_vol;
  ParOptScalar volume_frac;
  TopoAnalysis<ParOptScalar, degree> &topo;
  int opt_iter;
};

#endif