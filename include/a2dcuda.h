/**
 * CUDA functionalities: unified memory allocation, kernel wrapper, etc.
 */

#ifndef A2D_CUDA_H
#define A2D_CUDA_H

#include <cstddef>
#include <iostream>
#include <stdexcept>

#include "a2dobjs.h"

/**
 * @brief Check device execution
 *
 *  Regularly call this after execution of device code to make sure it returns
 * successfully
 */
#define cuda_check_error() \
  { __cuda_check_error(__FILE__, __LINE__); }

void __cuda_check_error(const char *file, int line);

/**
 * @brief Allocate memory in CUDA's unified memory space.
 *
 * @param ptr_addr address of the pointer to unified memory to be allocated
 * @param size memory size in byte
 * @return int error code, 0 = success
 */
int cuda_malloc(void **ptr_addr, std::size_t size);

/**
 * @brief Free memory in CUDA's unified memory space
 *
 * @param ptr pointer to allocated unified memory
 * @return int error code, 0 = success
 */
int cuda_free(void *ptr);

void cuda_device_synchronize();

/**
 * @brief CUDA kernel to perform parallel for
 *
 * @tparam FunctorType the functor type
 * @param N number of items for the parallel for loop
 * @param func the unary function -> void
 */
template <class FunctorType>
__global__ void __kernel_parallel_for(A2D::index_t N, const FunctorType func) {
  A2D::index_t index = blockIdx.x * blockDim.x + threadIdx.x;
  A2D::index_t stride = blockDim.x * gridDim.x;
  for (A2D::index_t i = index; i != N; i += stride) {
    func(i);
  }
}

/**
 * @brief Execute parallel for on CUDA device
 *
 * @tparam FunctorType the functor type
 * @param N number of items for the parallel for loop
 * @param func the unary function -> void
 * @param blockSize number of threads in each threadblock, default = 256
 */
template <class FunctorType>
void cuda_parallel_for(A2D::index_t N, const FunctorType func,
                       int blockSize = 256) {
  int numBlocks = (N + blockSize - 1) / blockSize;
  __kernel_parallel_for<<<numBlocks, blockSize>>>(N, func);
  cuda_device_synchronize();
  cuda_check_error();
}

#endif  // A2D_CUDA_H
