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

## Testing
Unit tests are implemented using [Google
Test](https://google.github.io/googletest/primer.html) framework, which is
automatically downloaded when building tests with
[CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html).

## Build examples, python bindings and tests

It is recommended to build the binaries with
[CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html).
To build all examples and python binding, execute the following commands in
this directory:

```
mkdir build
cd build
cmake .. -DBUILD_EXTENSION=ON -DBUILD_EXAMPLES=ON
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
