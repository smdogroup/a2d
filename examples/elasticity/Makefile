
include Makefile.in

default: test elasticity extension

elasticity: elasticity.cpp
	${CXX} ${OPENMP_ARGS} -O3 -I../../src/ -std=c++11 -o elasticity elasticity.cpp ${OPENMP_LIB}

test: test.cpp
	${CXX} ${OPENMP_ARGS} -O3 -I../../src/ -std=c++11 -o test test.cpp ${OPENMP_LIB}

extension:
	${CXX} -O3 ${OPENMP_ARGS} ${LIB_ARGS} ${shell python -m pybind11 --includes} example.cpp -o example${shell python3-config --extension-suffix} ${OPENMP_LIB}

clean:
	rm -rf test elasticity ./build *.so
