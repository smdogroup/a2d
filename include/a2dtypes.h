#ifndef A2D_TYPES_H
#define A2D_TYPES_H

// #include "a2dobjs.h"

namespace A2D {

template <class A>
class ADExpression {
 public:
  A2D_INLINE_FUNCTION A& cast() { return static_cast<A&>(*this); }
  A2D_INLINE_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  A2D_INLINE_FUNCTION void forward() { cast().forward(); }
  A2D_INLINE_FUNCTION void reverse() { cast().reverse(); }
};

template <class A>
class A2DExpression {
 public:
  A2D_INLINE_FUNCTION A& cast() { return static_cast<A&>(*this); }
  A2D_INLINE_FUNCTION const A& cast() const {
    return static_cast<const A&>(*this);
  }

  A2D_INLINE_FUNCTION void reverse() { cast().reverse(); }
  A2D_INLINE_FUNCTION void hforward() { cast().hforward(); }
  A2D_INLINE_FUNCTION void hreverse() { cast().hreverse(); }
};

}  // namespace A2D

#endif  // A2D_TYPES_H
