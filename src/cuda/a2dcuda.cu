/**
 * CUDA functionalities: unified memory allocation, kernel wrapper, etc.
 */

#include "a2dcuda.h"

int cuda_malloc(void **ptr_addr, std::size_t size) {
  cudaError_t ret = cudaMallocManaged(ptr_addr, size);
  return static_cast<int>(ret);
}

int cuda_free(void *ptr) {
  cudaError_t ret = cudaFree(ptr);
  return static_cast<int>(ret);
}
