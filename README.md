# A2D - **A**lmost **a**utomatic **d**ifferentiation toolkit

A toolkit for almost automatic differentiation of vector and matrix expressions.

This code relies heavily on the approach for deriving auto-diff expressions by
M. B. Giles, "Collected matrix derivative results for forward and reverse mode
AD".

A2D is a header only c++ templated library with python binding created with
[pybind11](https://pybind11.readthedocs.io/en/stable/).

## Dependency
To use A2D in your c++ application, the following libraries need to be linked:
- OpenMP
- LAPACK

## Build examples, python bindings and tests

Build system [CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html) is used.
To build all examples and python binding, execute the following commands in
this directory:

```
mkdir build
cd build
cmake .. -DA2D_BUILD_EXAMPLES_BASIC=ON -DA2D_BUILD_EXAMPLES_AMGX=ON -DA2D_BUILD_EXAMPLES_KOKKOS=ON -DA2D_BUILD_EXTENSION=ON
make -j <nproc>
````

To install python extensions (copy binaries from ```./build``` to ```./a2d```),
execute
```
cd build
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
in ```build``` folder.


## Testing
Unit tests are implemented using [Google
Test](https://google.github.io/googletest/primer.html) framework, which is
automatically downloaded when building tests.

To run unit tests, need to `cmake` with `-DA2D_BUILD_UNIT_TESTS=ON`. After build, execute
```cd build/tests && ctest```.

## Coding style
```clangFormat``` is used as the auto-formatter, with style ```Google```. If you would
like to contribute to the project, please make sure you set up the auto-formatter accordingly.
