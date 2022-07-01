## Compile the code and tests with CMake

To build a2d with [CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html), 
run the following commands in this directory:

```
mkdir build
cd build
cmake ..
make -j <nproc>
````

This builds the executables, pybind11 extensions as well as the unit tests implemented using 
[Google Test](https://google.github.io/googletest/primer.html) framework.
Unit tests are in ```build/tests``` and are executable.
