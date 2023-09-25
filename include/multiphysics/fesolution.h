#ifndef A2D_FE_SOLUTION_H
#define A2D_FE_SOLUTION_H

#include <vector>

#include "a2ddefs.h"

namespace A2D {

/*
  The solution vector, TODO: we actually don't need this, drop this and use
  Kokkos::view dierctly
*/
template <typename T>
class SolutionVector {
 public:
  SolutionVector(index_t ndof) : ndof(ndof), array("array", ndof) { zero(); }
  KOKKOS_FUNCTION T& operator[](index_t index) const { return array(index); }
  // KOKKOS_FUNCTION const T& operator[](index_t index) const {
  //   return array(index);
  // }

  KOKKOS_FUNCTION index_t get_num_dof() const { return ndof; }

  KOKKOS_FUNCTION void zero() { BLAS::zero(array); }
  KOKKOS_FUNCTION void fill(T val) { BLAS::fill(array, val); }
  KOKKOS_FUNCTION T* data() { return array.data(); }

 private:
  const index_t ndof;
  MultiArrayNew<T*> array;
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_H