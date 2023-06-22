#ifndef TOPOLOGY_PAROPT_PROB_H
#define TOPOLOGY_PAROPT_PROB_H

#include <string>

#include "ParOptOptimizer.h"
#include "utils/a2dprofiler.h"

template <class Analysis>
class TopOptProb : public ParOptProblem {
 public:
  TopOptProb(std::string prefix, MPI_Comm comm, int nvars, int ncon, int nineq,
             ParOptScalar ref_comp, ParOptScalar domain_vol,
             ParOptScalar volume_frac, Analysis &topo, bool verbose,
             int vtk_freq)
      : ParOptProblem(comm),
        prefix(prefix),
        comm(comm),
        nvars(nvars),
        ncon(ncon),
        nineq(nineq),
        ref_comp(ref_comp),
        domain_vol(domain_vol),
        volume_frac(volume_frac),
        topo(topo),
        opt_iter(0),
        verbose(verbose),
        vtk_freq(vtk_freq) {
    setProblemSizes(nvars, ncon, 0);
    setNumInequalities(nineq, 0);
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
    A2D::Timer timer("TopOptProb::evalObjCon()");
    // Set design variables by paropt
    ParOptScalar *x;
    xvec->getArray(&x);

    topo.set_design_vars(x);
    topo.solve();

    // Evaluate objective
    *fobj = topo.eval_compliance() / ref_comp;

    // Evaluate constraint
    ParOptScalar vol = topo.eval_volume();
    cons[0] = volume_frac - vol / domain_vol;

    // Write design to vtk
    if (opt_iter % vtk_freq == 0) {
      char vtk_name[256];
      std::snprintf(vtk_name, sizeof(vtk_name), "result_%d.vtk", opt_iter);
      std::filesystem::path path =
          std::filesystem::path(prefix) / std::filesystem::path(vtk_name);
      topo.tovtk(path);
    }

    opt_iter++;

    std::printf("[%2d] obj: %20.10e, con: %20.10e\n", opt_iter, *fobj, cons[0]);

    return 0;
  }

  int evalObjConGradient(ParOptVec *xvec, ParOptVec *gvec, ParOptVec **Ac) {
    A2D::Timer timer("TopOptProb::evalObjConGradient()");
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
  std::string prefix;
  MPI_Comm comm;
  int nvars;
  int ncon;
  int nineq;
  ParOptScalar ref_comp;
  ParOptScalar domain_vol;
  ParOptScalar volume_frac;
  Analysis &topo;
  int opt_iter;
  bool verbose;
  int vtk_freq;
};

#endif