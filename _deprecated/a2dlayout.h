#ifndef A2D_LAYOUT_H
#define A2D_LAYOUT_H

#include "multiarray.h"

namespace A2D {

/**
 * @brief Control the multiarray layout
 */
template <index_t... dims>
using A2D_Layout = CLayout<dims...>;

/**
 * @brief Check if a Layout is CLayout or FLayout
 *
 * Usage:
 *   if Layout is CLayout<...>, then
 *     is_layout<Layout, CLayout>::value == true
 *   if Layout is FLayout<...>, then
 *     is_layout<Layout, FLayout>::value == true
 */
template <class, template <index_t...> class>
struct is_layout : public std::false_type {};

template <index_t... dims, template <index_t...> class U>
struct is_layout<U<dims...>, U> : public std::true_type {};

}  // namespace A2D

#endif  // A2D_LAYOUT_H