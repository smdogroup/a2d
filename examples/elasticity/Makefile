
include Makefile.in

default: test_numeric test_grad_hess elasticity extension

elasticity: elasticity.cpp
	${CXX} ${OPENMP_ARGS} ${CFLAGS} -std=c++11 -o elasticity elasticity.cpp ${OPENMP_LIB}

test_grad_hess: test_grad_hess.cpp
	${CXX} ${OPENMP_ARGS} ${CFLAGS} -std=c++11 -o test_grad_hess test_grad_hess.cpp ${OPENMP_LIB}

test_numeric: test_numeric.cpp
	${CXX} ${OPENMP_ARGS} ${CFLAGS} -std=c++11 -o test_numeric test_numeric.cpp ${OPENMP_LIB}

extension:
	${CXX} ${CFLAGS} -std=c++11 ${OPENMP_ARGS} ${LIB_ARGS} ${shell python3 -m pybind11 --includes} example.cpp -o example${shell python3-config --extension-suffix} ${OPENMP_LIB}

clean:
	rm -rf test_numeric test_grad_hess elasticity  *.so
