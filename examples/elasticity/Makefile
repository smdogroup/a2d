
default: test elasticity elasticity_impl

elasticity_impl:
	${CXX} -O3 -g -I../../src/ -std=c++11 -c elasticity_impl.cpp

elasticity:
	${CXX} -O3 -g -I../../src/ -std=c++11 -o elasticity elasticity.cpp

test:
	${CXX} -03 -g -I../../src/ -std=c++11 -o test test.cpp

interface:
	python3 setup.py build_ext --inplace
