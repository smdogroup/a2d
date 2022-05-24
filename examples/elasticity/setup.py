from setuptools import Extension, setup
import numpy as np
from Cython.Build import cythonize

setup(ext_modules=cythonize(Extension(
           "model",                                # the extension name
           sources=["linear_elasticity.pyx", "elasticity_impl.cpp"], # the Cython source and
                                                  # additional C++ source files
           language="c++",                        # generate and compile C++ code
           include_dirs=[np.get_include()],
           extra_compile_args=['-std=c++11']
      )))