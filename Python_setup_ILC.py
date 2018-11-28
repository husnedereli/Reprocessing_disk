

## Compile the code with :
# python Python_setup_ILC.py build_ext -i -lcfitsio --inplace

# setup.py file
import sys
import os
import shutil

from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

## AD stands for accretion disk model

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("AD_module",
		sources= ["AD_module.pyx", "illuminating_LC.c"],
		language = "c",
        include_dirs=[numpy.get_include()])]
)
