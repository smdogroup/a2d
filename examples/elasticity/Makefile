
default: test elasticity interface

CXX = mpicxx

elasticity: elasticity.cpp
	${CXX} -O3 -g -I../../src/ -std=c++11 -o elasticity elasticity.cpp

test: test.cpp
	${CXX} -O2 -g -I../../src/ -std=c++11 -o test test.cpp

interface: elasticity_impl.cpp
	python3 setup.py build_ext --inplace

clean:
	rm -rf test elasticity ./build model.*.so
