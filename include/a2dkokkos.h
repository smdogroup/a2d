#include "Kokkos_Core.hpp"
#include "a2dobjs.h"

namespace A2D {

// template <typename T, index_t N0, index_t N1 = 0, index_t N2 = 0,
//           index_t N3 = 0>
// struct KokkosLayoutWrapper;

// template <typename T, index_t N0>
// struct KokkosLayoutWrapper<T, N0, 0, 0, 0> {
//   static T* dummy_array[N0];
//   KokkosLayoutWrapper(const index_t _leading_dim) : leading_dim(_leading_dim)
//   {} const index_t leading_dim;
// };

// template <typename T, index_t N0, index_t N1>
// struct KokkosLayoutWrapper<T, N0, N1, 0, 0> {
//   static T* dummy_array[N0][N1];
//   KokkosLayoutWrapper(const index_t _leading_dim) : leading_dim(_leading_dim)
//   {} const index_t leading_dim;
// };

}  // namespace A2D
