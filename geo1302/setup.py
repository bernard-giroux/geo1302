#
#  Compiler avec la commande     python setup.py build_ext --inplace
#
#
import os
from distutils.extension import Extension
from distutils.core import setup

import numpy as np

from Cython.Build import cythonize

os.environ['CC'] = 'clang'
os.environ['CXX'] = 'clang++'

extra_compile_args = ['-std=c++11', '-stdlib=libc++', '-O3']
include_dirs = [np.get_include(), ]


extensions = [
    Extension('em_c',
              sources=['em_c.pyx'],
              include_dirs=include_dirs,
              language='c++',
              extra_compile_args=extra_compile_args,
              ),
]

setup(
    ext_modules=cythonize(extensions, include_path=['.', ], language_level=3),
)
