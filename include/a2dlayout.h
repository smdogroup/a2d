#ifndef A2D_LAYOUT_H
#define A2D_LAYOUT_H

#include "multiarray.h"

namespace A2D {
template <index_t... dims>
using A2D_Layout = CLayout<dims...>;
}

#endif  // A2D_LAYOUT_H