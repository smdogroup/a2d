#ifndef TOOLKIT_H
#define TOOLKIT_H

#include <bits/stdc++.h>
#include <stdlib.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include "cuda_runtime.h"

#if defined(__GNUC__) || defined(__clang__)
#include <cxxabi.h>
#endif
#include "helper_cuda.h"

// using for showing the type of a variable
#define SHOW_TYPE(T) std::cout << type_name<T>() << std::endl;

// debug
#define __CHECK__ printf("Line: %i => pass\n", __LINE__);

// using for timing
#define TIMER(key)                                                             \
  {                                                                            \
    static void timing(std::string const &key) {                               \
      static std::map<std::string, std::chrono::steady_clock::time_point>      \
          saves;                                                               \
      auto it = saves.find(key);                                               \
      if (it == saves.end()) {                                                 \
        saves.emplace(key, std::chrono::steady_clock::now());                  \
      } else {                                                                 \
        double dt = std::chrono::duration_cast<std::chrono::duration<double>>( \
                        std::chrono::steady_clock::now() - it->second)         \
                        .count();                                              \
        std::cout << key << ": " << dt << " s" << std::endl;                   \
      }                                                                        \
    }                                                                          \
  }

// using for generating random numbers
#define RAND_NUM                                         \
  {                                                      \
    static float frand() {                               \
      static std::mt19937 gen;                           \
      static std::uniform_real_distribution<float> unif; \
      return unif(gen);                                  \
    }                                                    \
  }

// // using for checking the GPU error
#define CHECK(ans) \
  { gpuErrorCheck((ans), __FILE__, __LINE__); }

/* CUDA error macro */
#define CUDA_SAFE_CALL(call) \
  { cuda_call_sell(call); }

template <class T>
std::string type_name() {
  const char *name = typeid(T).name();
#if defined(__GNUC__) || defined(__clang__)
  int status;
  char *p = abi::__cxa_demangle(name, 0, 0, &status);
  std::string s = p;
  std::free(p);
#else
  std::string s = name;
#endif
  if (std::is_const_v<std::remove_reference_t<T>>) s += " const";
  if (std::is_volatile_v<std::remove_reference_t<T>>) s += " volatile";
  if (std::is_lvalue_reference_v<T>) s += " &";
  if (std::is_rvalue_reference_v<T>) s += " &&";
  return s;
}

void gpuErrorCheck(cudaError_t code, const char *file, int line,
                   bool abort = true) {
  if (code != cudaSuccess) {
    std::cerr << "GPU Error: " << cudaGetErrorString(code) << " at " << file
              << ":" << line << std::endl;
    // checkCudaErrors(code);
    // printf("GPU: %s\n", cudaGetErrorString(cudaGetLastError()));
    if (abort) exit(code);
  }
}

void cuda_call_sell(cudaError_t err) {
  if (cudaSuccess != err) {
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", __FILE__,
            __LINE__, cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
  while (0)
    ;
}

#endif  // TOOLKIT_H