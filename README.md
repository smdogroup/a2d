
<h2 align="center">
    <img src="docs/a2d_logo.svg" width="400" />
</h2>

# A2D - a PDE discretization library using almost automatic differentiation

A toolkit for almost automatic differentiation of vector and matrix expressions.

This code relies heavily on the approach for deriving auto-diff expressions by
M. B. Giles, "Collected matrix derivative results for forward and reverse mode
AD".

## Dependencies

A2D is a header-only c++ templated library. The only requirement for using A2D
is a C++17 supported compiler.

## Install and use

### CMake

CMake is preferred to install (copy over) A2D. For basic installation, use the
following command:
```
mkdir build && cd build && cmake .. && make install
```
This installs A2D (headers and CMake files) into ```${HOME}/installs/a2d```.
Then in the application, add
```
find_package(A2D REQUIRED PATHS <path-to-sparse-utils-installation>)
```
to ```CMakeLists.txt```, and use
```
target_link_libraries(<app-target> A2D::A2D)
```
in the ```CMakeLists.txt``` for the application executables. See
[examples/CMakeLists.txt](examples/CMakeLists.txt) for example.


### Manual
Alternatively, you can directly include ```include/a2dcore.h``` and manually
manage the include path for the compiler.

## Examples
Use the following to build examples:
```
cd examples &&
mkdir build &&
cd build &&
cmake .. &&
make -j
```

## Test
Unit tests are implemented using [Google
Test](https://google.github.io/googletest/primer.html) framework, which is
automatically downloaded when building tests. Use the following snippet to build
and run unit tests.
```
mkdir build &&
cd build &&
cmake .. -DA2D_BUILD_TESTS=ON &&
make -j &&
ctest
```

## Code style
```clangFormat``` is used as the auto-formatter, with style ```--style=Google```.
 If you would like to contribute to the project, please make sure you set up the
 auto-formatter accordingly.

