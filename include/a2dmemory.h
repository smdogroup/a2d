#ifndef A2D_MEMORY_H
#define A2D_MEMORY_H

#include <cstddef>
#include <iostream>
#include <stdexcept>

#ifdef A2D_USE_CUDA
#include "a2dcuda.h"
#endif

namespace A2D {

/**
 * @brief malloc wrapper for a2d
 *
 * @tparam T data type
 * @param ptr_addr address of the pointer to memory to be allocated
 * @param n size of data in the array
 */
template <typename T>
void a2d_malloc(T** ptr_addr, std::size_t n) {
#ifdef A2D_USE_CUDA
  std::printf("a2dmemory.h:calling cuda_malloc()\n");

  // Allocate CUDA unified memory
  int fail = cuda_malloc((void**)(void*)ptr_addr, sizeof(T) * n);
  if (fail) {
    *ptr_addr = nullptr;
    char msg[256];
    std::sprintf(msg, "cudaMallocManaged() failed with exit code %d.", fail);
    throw std::runtime_error(msg);
  }
#else
  std::printf("a2dmemory.h:using operator new[]\n");

  // Allocate CPU memory
  *ptr_addr = new T[n];
#endif
  return;
}

/**
 * @brief free wrapper for a2d
 *
 * @tparam T data type
 * @param ptr pointer to the allocated memory
 */
template <typename T>
void a2d_free(T* ptr) {
#ifdef A2D_USE_CUDA
  std::printf("a2dmemory.h:calling cuda_free()\n");
  int fail = cuda_free(ptr);
  if (fail) {
    char msg[256];
    std::sprintf(msg, "cudaFree() failed with exit code %d.", fail);
    throw std::runtime_error(msg);
  }
#else
  std::printf("a2dmemory.h:using operator delete[]\n");
  delete[] ptr;
#endif
  return;
}

}  // namespace A2D

#endif  // A2D_MEMORY_H