#ifndef A2D_TYPES_H
#define A2D_TYPES_H

// #include "a2dobjs.h"

namespace A2D {

template <class A>
class ADExpression {
 public:
  KOKKOS_FUNCTION A& cast() { return static_cast<A&>(*this); }
  KOKKOS_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  KOKKOS_FUNCTION void forward() { cast().forward(); }
  KOKKOS_FUNCTION void reverse() { cast().reverse(); }
};

template <class A>
class A2DExpression {
 public:
  KOKKOS_FUNCTION A& cast() { return static_cast<A&>(*this); }
  KOKKOS_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  KOKKOS_FUNCTION void reverse() { cast().reverse(); }
  KOKKOS_FUNCTION void hforward() { cast().hforward(); }
  KOKKOS_FUNCTION void hreverse() { cast().hreverse(); }
};

}  // namespace A2D

#endif  // A2D_TYPES_H
