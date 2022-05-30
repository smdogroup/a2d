
default: test elasticity interface

CXX = mpicxx

elasticity: elasticity.cpp
	${CXX} -O3 -g -I../../src/ -std=c++11 -o elasticity elasticity.cpp

test: test.cpp
	${CXX} -O3 -g -I../../src/ -std=c++11 -o test test.cpp

interface: elasticity_impl.cpp
	python3 setup.py build_ext --inplace

clean:
	rm -rf test elasticity ./build model.*.so


extension:
	${CXX} -g -Wall -dynamiclib -undefined dynamic_lookup -std=c++11 -fPIC ${shell python -m pybind11 --includes} example.cpp -o example${shell python-config --extension-suffix}
