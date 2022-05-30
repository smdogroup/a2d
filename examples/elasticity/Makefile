
default: test elasticity extension

CXX = mpicxx

elasticity: elasticity.cpp
	${CXX} -Xpreprocessor -fopenmp -O3 -I../../src/ -std=c++11 -o elasticity elasticity.cpp -lomp

test: test.cpp
	${CXX} -O3 -I../../src/ -std=c++11 -o test test.cpp

extension:
	${CXX} -O3 -Wall -dynamiclib -undefined dynamic_lookup -std=c++11 -fopenmp -fPIC ${shell python -m pybind11 --includes} example.cpp -o example${shell python-config --extension-suffix}

clean:
	rm -rf test elasticity ./build *.so
