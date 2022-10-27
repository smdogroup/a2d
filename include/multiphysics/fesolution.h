#ifndef A2D_FE_SOLUTION_H
#define A2D_FE_SOLUTION_H

#include <vector>

namespace A2D {

/*
  The solution vector
*/
template <typename T>
class SolutionVector {
 public:
  SolutionVector(A2D::index_t ndof) : ndof(ndof), x(ndof) {}
  T& operator[](A2D::index_t index) { return x[index]; }
  const T& operator[](A2D::index_t index) const { return x[index]; }

 private:
  const A2D::index_t ndof;
  std::vector<T> x;
};

}  // namespace A2D

#endif  // A2D_FE_SOLUTION_H