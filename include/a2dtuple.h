#ifndef A2D_TUPLE_H
#define A2D_TUPLE_H

#include <type_traits>

#define USE_CUSTOM_TUPLE

#ifdef USE_CUSTOM_TUPLE
template <class _Tp>
inline constexpr _Tp&& a2d_forward(
    typename std::remove_reference<_Tp>::type& __t) noexcept {
  return static_cast<_Tp&&>(__t);
}
#else
#include <utility>
#define a2d_forward std::forward
#endif

#endif  //  A2D_TUPLE_H
