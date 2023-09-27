
<h2 align="center">
    <img src="docs/a2d_logo.svg" width="400" />
</h2>

# A2D - a PDE discretization library using almost automatic differentiation

A toolkit for almost automatic differentiation of vector and matrix expressions.

This code relies heavily on the approach for deriving auto-diff expressions by
M. B. Giles, "Collected matrix derivative results for forward and reverse mode
AD".

A2D is a header only c++ templated library.
<!-- with python binding created with
[pybind11](https://pybind11.readthedocs.io/en/stable/). -->

## TODO list
- GPU support
    - [ ] replace plain C-array with Kokkos array for all data that are needed from GPU
    - [ ] decorate all functions that are needed from device
    - [ ] make all element operations in ```feelement.h``` support GPU parallelization
- AD
    - [x] complete core operations
    - [x] complete AD/A2D expressions
    - [x] unit test for core operations and AD/A2D expressions
- General
    - [ ] ```feelementmat.h``` add a parallel implementation of ElementMat (similar to ElementVector)
    - [ ] ```array.h``` replace reference implementations under ```BLAS``` namespace with Kokkos kernel functions
    - [ ] ```a2dobjs.h``` drop ```A2D::fill()``` and use Kokkos kernel
    - [ ] ```integrand_[poisson/elasticity/heat_conduction].h``` rename ```weak()``` -> ```residual()```
    - [ ] ```integrand_[poisson/elasticity/heat_conduction].h``` drop Jacobian vector product functor and Adjoint jacobian product functor, compute second order derivative matrices directly using ```jtransform()``` (see ```IntegrandTopoLinearElasticity::jacobian()``` for example)
    - [ ] ```feelement.h``` change ```add_jacobian()```, ```add_jacobian_vector_product()```, etc. accordingly
- Matrix
    - [ ] External library for parallel matrix solution on GPU
    - [ ] External library for serial non-symmetric matrix factorization



## Install

A2D is currently header-only so no need to build itself.
The following command can be used to build examples and unit tests:
```
# You must start from the root directory of a2d
mkdir build &&
cd build &&
cmake .. -DCMAKE_BUILD_TYPE=[Debug,Release] -DA2D_KOKKOS_DIR=<your Kokkos install dir> -DA2D_METIS_DIR=<your metis install dir> -DA2D_BUILD_EXAMPLES=ON -DA2D_BUILD_UNIT_TESTS=ON &&
make -j # parallel make using maximum number of processors
```

Note: metis and Kokkos are assumed to be installed in ```a2d/installs/metis``` and
```a2d/installs/kokkos``` if corresponding CMake variables are not specified.
See [CMake variables](#cmake-variables) for a complete list of A2D CMake variables
and defaults.
See [Install Kokkos](#install-kokkos) and [Install METIS](#install-metis) for
instructions on installing Kokkos and METIS.

## CMake variables

Below is the complete table of CMake variables that A2D accepts to control the
compilation.

_Recall that to give the variable VARIABLE value VAL, use the following syntax
int the command line:_
```
cmake ... -DVARIABLE=VAL ...
```

| Variable | Description | Default | Choices |
|----------|-------------|---------|---------|
| CMAKE_BUILD_TYPE | whether this is a release (optimized) build or debug (containing debug info) build | No default | Debug/Release |
| A2D_KOKKOS_DIR | directory of kokkos installation | a2d/installs/kokkos | a valid path |
| A2D_METIS_DIR | directory of metis installation | a2d/installs/metis | a valid path |
| A2D_BUILD_EXAMPLES | build examples if set to ON | ON | ON/OFF |
| A2D_BUILD_UNIT_TESTS | build unit tests if set to ON | OFF | ON/OFF |



## Theory
[So, how does A2D solve PDEs?](docs/theory.md)

## Dependencies
A2D requires following dependencies
- OpenMP
- LAPACK
- [Kokkos](https://github.com/kokkos/kokkos)
- [METIS (5.1.0)](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
<!-- - [Kokkos-kernels](https://github.com/kokkos/kokkos-kernels) -->

### Install Kokkos
To build Kokkos with OpenMP and CUDA backend, use the following commands:
```
# You must start from the root directory of a2d
cd extern &&
git clone https://github.com/kokkos/kokkos.git &&
cd kokkos &&
mkdir build &&
cd build &&
cmake .. -DCMAKE_INSTALL_PREFIX=../../../installs/kokkos -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -G Ninja &&
ninja install
```
For a complete instruction on installing Kokkos, see [Kokkos
documentation](https://kokkos.github.io/kokkos-core-wiki/ProgrammingGuide/Compiling.html).

### Install METIS

Obtain the tarball from [here](https://src.fedoraproject.org/lookaside/pkgs/metis/metis-5.1.0.tar.gz/5465e67079419a69e0116de24fce58fe/).
The following commands can be used to install METIS to a2d/extern/metis:
```
# You must start from the root directory of a2d
METIS_SOURCE_DIR=$(pwd)/extern &&
METIS_INSTALL_DIR=$(pwd)/installs/metis &&
cd $METIS_SOURCE_DIR &&
wget https://src.fedoraproject.org/lookaside/pkgs/metis/metis-5.1.0.tar.gz/5465e67079419a69e0116de24fce58fe/metis-5.1.0.tar.gz &&
tar -zxvf metis-5.1.0.tar.gz &&
cd metis-5.1.0 &&
make config prefix=$METIS_INSTALL_DIR &&
make &&
make install
```


<!-- ### Install Kokkos-kernels
To build Kokkos-kernels with Kokkos installation, use:
```
cd extern/kokkos-kernels
mkdir build
cd build
cmake .. -DKokkos_ROOT=../../../installs/kokkos -DCMAKE_INSTALL_PREFIX=../../../installs/kokkos-kernels -G Ninja
ninja install
```
For a complete instruction on installing Kokkos-kernels, see [Kokkos-kernels
documentation](https://github.com/kokkos/kokkos-kernels/wiki/Building). -->

<!-- ## Build examples, python bindings and tests

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
in ```build_<BUILD_TYPE>``` folder. -->


## Testing
Unit tests are implemented using [Google
Test](https://google.github.io/googletest/primer.html) framework, which is
automatically downloaded when building tests.

CTest (bundled with CMake) is used to execute tests, simply go to```<build
dir>/tests``` and execute
```
ctest
```

## Coding style
```clangFormat``` is used as the auto-formatter, with style ```--style=Google```.
 If you would like to contribute to the project, please make sure you set up the
 auto-formatter accordingly.

