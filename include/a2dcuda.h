/**
 * CUDA functionalities: unified memory allocation, kernel wrapper, etc.
 */

#ifndef A2D_CUDA_H
#define A2D_CUDA_H

#include <cstddef>

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

#endif  // A2D_CUDA_H
