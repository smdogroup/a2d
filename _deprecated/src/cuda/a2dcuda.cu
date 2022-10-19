/**
 * CUDA functionalities: unified memory allocation, kernel wrapper, etc.
 */

#include "a2dcuda.h"

void __cuda_check_error(const char *file, int line) {
  // Get error code
  cudaError_t err = cudaGetLastError();

  if (err != cudaSuccess) {
    char msg[256];
    std::sprintf(msg, "%s:%d:CUDA kernel launch failed with exit code %d (%s)",
                 file, line, err, cudaGetErrorString(err));
    throw std::runtime_error(msg);
  }
}

int cuda_malloc(void **ptr_addr, std::size_t size) {
  cudaError_t ret = cudaMallocManaged(ptr_addr, size);
  return static_cast<int>(ret);
}

int cuda_free(void *ptr) {
  cudaError_t ret = cudaFree(ptr);
  return static_cast<int>(ret);
}

void cuda_device_synchronize() { cudaDeviceSynchronize(); }

// template <class FunctorType>
// __global__ void __kernel_parallel_for(A2D::index_t N, const FunctorType func)
// {
//   A2D::index_t index = blockIdx.x * blockDim.x + threadIdx.x;
//   A2D::index_t stride = blockDim.x * gridDim.x;
//   for (A2D::index_t i = index; i != N; i += stride) {
//     func(i);
//   }
// }

// template <class FunctorType>
// void cuda_parallel_for(A2D::index_t N, const FunctorType func,
//                        int blockSize = 256) {
//   int numBlocks = (N + blockSize - 1) / blockSize;
//   __kernel_parallel_for<<<numBlocks, blockSize>>>(N, func);
//   cudaDeviceSynchronize();
//   CUDACheckError();
// }