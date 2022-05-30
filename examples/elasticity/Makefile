
default: test elasticity extension

CXX = mpicxx

elasticity: elasticity.cpp
	${CXX} -Xpreprocessor -fopenmp -O3 -I../../src/ -std=c++11 -o elasticity elasticity.cpp -lomp

test: test.cpp
	${CXX} -Xpreprocessor -fopenmp -O3 -I../../src/ -std=c++11 -o test test.cpp -lomp

extension:
	${CXX} -O3 -Wall -dynamiclib -undefined dynamic_lookup -std=c++11 -Xpreprocessor -fopenmp -fPIC ${shell python -m pybind11 --includes} example.cpp -lomp -o example${shell python-config --extension-suffix}

clean:
	rm -rf test elasticity ./build *.so
