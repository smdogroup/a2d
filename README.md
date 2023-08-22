# A2D - a PDE solver using Almost Automatic Differentiation

A toolkit for almost automatic differentiation of vector and matrix expressions.

This code relies heavily on the approach for deriving auto-diff expressions by
M. B. Giles, "Collected matrix derivative results for forward and reverse mode
AD".

A2D is a header only c++ templated library with python binding created with
[pybind11](https://pybind11.readthedocs.io/en/stable/).

## Theory
[So, how does A2D solve PDEs?](docs/theory.md)

## Dependencies
A2D requires following dependencies
- OpenMP
- LAPACK
- [Kokkos](https://github.com/kokkos/kokkos)
- [Kokkos-kernels](https://github.com/kokkos/kokkos-kernels)

### Install Kokkos
To build Kokkos with OpenMP and CUDA backend, use:
```
cd extern/kokkos
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../../installs/kokkos -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -G Ninja
ninja install
```
For a complete instruction on installing Kokkos, see [Kokkos
documentation](https://kokkos.github.io/kokkos-core-wiki/ProgrammingGuide/Compiling.html).

### Install Kokkos-kernels
To build Kokkos-kernels with Kokkos installation, use:
```
cd extern/kokkos-kernels
mkdir build
cd build
cmake .. -DKokkos_ROOT=../../../installs/kokkos -DCMAKE_INSTALL_PREFIX=../../../installs/kokkos-kernels -G Ninja
ninja install
```
For a complete instruction on installing Kokkos-kernels, see [Kokkos-kernels
documentation](https://github.com/kokkos/kokkos-kernels/wiki/Building).

## Build examples, python bindings and tests

Build system [CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html) is used.
A2D requires CMAKE_BUILD_TYPE to be set explicitly to either Debug or Release.
It is recommended to use out-of-source builds and separate different build types.
For example, to build examples with **debug**/**optimization** flags, do:

```
mkdir build-debug
cd build-debug
cmake .. -DCMAKE_BUILD_TYPE=[Debug,Release]
make -j <nproc>
````

To install python extensions (copy binaries from ```./build_<BUILD_TYPE>``` to ```./a2d```),
execute
```
cd build_<BUILD_TYPE>
cmake --install .
```

To access python extension as well as python utilities, add root directory to
```PYTHONPATH```. For example, add this to your shell startup.
```
export PYTHONPATH=${PYTHONPATH}:~/git/a2d
```

To see a full list of CMake options and their values for the current build, execute
```
ccmake .
```
in ```build_<BUILD_TYPE>``` folder.


## Testing
Unit tests are implemented using [Google
Test](https://google.github.io/googletest/primer.html) framework, which is
automatically downloaded when building tests.

To run unit tests, need to `cmake` with `-DA2D_BUILD_UNIT_TESTS=ON`. After build, execute
```cd build_<BUILD_TYPE>/tests && ctest```.

## Coding style
```clangFormat``` is used as the auto-formatter, with style ```Google```. If you would
like to contribute to the project, please make sure you set up the auto-formatter accordingly.

