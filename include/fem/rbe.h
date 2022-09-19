#ifndef A2D_RBE_H
#define A2D_RBE_H

namespace A2D {

template <typename I, typename T, class PDEInfo>
class RBE3Element : public ElementBase<I, T, PDEInfo> {
 public:
  template <typename IdxType>
  RBE3Element(const index_t nnodes, const IdxType nodes_[], ) : nnodes(nnodes) {
    for (index_t i = 0; i < nnodes; i++) {
      nodes(i) = nodes_[i];
    }
  }

  ~RBE3Element() {}

  // Set the node locations
  void set_nodes(typename PDEInfo::NodeArray& X) {
    for (index_t i = 0; i < nnodes; i++) {
      for (index_t j = 0; j < 3; j++) {
        Xn(i, j) = X(i, j);
      }
    }
  }

  // Set the non-zero pattern
  void add_node_set(std::set<std::pair<I, I>>& node_set) {}

  // Set the solution
  void set_solution(typename PDEInfo::SolutionArray& U) {}

  // Add the residual
  void add_residual(typename PDEInfo::SolutionArray& res) {}

  // Add the Jacobian matrix
  void add_jacobian(typename PDEInfo::SparseMat& jac) {}

 private:
  const index_t nnodes;
};

}  // namespace A2D

#endif  // A2D_RBE_H