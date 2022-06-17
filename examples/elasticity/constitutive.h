#ifndef A2D_CONSTITUTIVE_H
#define A2D_CONSTITUTIVE_H

namespace A2D {

/*
  Base class for a constitutive object.

  These objects provide a way to set the design variables into the elements.
  The mapping from design variables to element data is not defined and depends
  on the nature of the design variable.
*/
template <typename I, typename T, class PDE>
class Constitutive {
 public:
  virtual void set_design_vars(typename PDE::DesignArray& x) {}
  virtual void add_adjoint_dfdx(typename PDE::SolutionArray& psi,
                                typename PDE::DesignArray& dfdx) {}
};

}  // namespace A2D

#endif  // A2D_CONSTITUTIVE_H