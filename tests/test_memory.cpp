/*
  This is a set of unit tests for a2dtmp3d.h using Google Test framework.
 */

#include "a2dmemory.h"
#include "cstddef"
#include "test_commons.h"

class MemoryTest : public ::testing::Test {};

// Test suite: memory allocation
class Malloc : public MemoryTest {
 public:
  template <typename T>
  void populate_array_and_expect_eq(T* ptr, std::size_t size, T val) {
    for (std::size_t i = 0; i != size; i++) {
      ptr[i] = val;
    }
    for (std::size_t i = 0; i != size; i++) {
      EXPECT_FLOAT_EQ(ptr[i], val);
    }
    return;
  }
};

TEST_F(Malloc, a2d_malloc_free) {
  double* ptr = nullptr;
  std::size_t array_size = 1L << 20;  // 1M values
  A2D::a2d_malloc(&ptr, array_size);
  populate_array_and_expect_eq(ptr, array_size, 4.2);
  A2D::a2d_free(ptr);
}
