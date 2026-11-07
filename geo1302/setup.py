#
#  Compiler avec la commande     python setup.py build_ext --inplace
#
import os
import platform
import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

if sys.platform == 'win32':
    extra_compile_args = ['/O2']
else:
    extra_compile_args = ['-std=c++11', '-O3']
    if sys.platform == 'darwin':
        extra_compile_args.append('-stdlib=libc++')
        # Python (python.org) est compilé pour une vieille version de macOS,
        # qui n'est plus supportée par libc++ : on cible la version courante
        os.environ.setdefault('MACOSX_DEPLOYMENT_TARGET',
                              platform.mac_ver()[0].split('.')[0] + '.0')

extensions = [
    Extension('em_c',
              sources=['em_c.pyx'],
              include_dirs=[np.get_include()],
              language='c++',
              extra_compile_args=extra_compile_args,
              ),
]

setup(
    ext_modules=cythonize(extensions, include_path=['.'], language_level=3),
)
