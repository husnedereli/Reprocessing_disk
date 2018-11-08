# python setup.py build_ext -i -lcfitsio --inplace

# setup.py file
import sys
import os
import shutil

from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("CC_module",
		sources= ["CC_module.pyx", "CC.c"],
		language = "c",
                libraries=["matrix"],
		extra_link_args=["-L/home/cayley/Desktop/Thesis/Programmation/AAA_Library"],
                include_dirs=[numpy.get_include(), "/home/cayley/Desktop/Thesis/Programmation/AAA_Library" ],
		runtime_library_dirs=["/home/cayley/Desktop/Thesis/Programmation/AAA_Library"],
	)]
)
#/home/cayley/Desktop/Thesis/Programmation/Kinetic/Asaf/C_library/matrix.h

#library_dirs=['/home/cayley/Desktop/Thesis/Programmation/AAA_Library/'],
# extra_link_args=["-L./home/cayley/Desktop/Thesis/Programmation/AAA_Library/"],
#	        extra_compile_args=["-I./home/cayley/Desktop/Thesis/Programmation/AAA_Library/"],
#                runtime_library_dirs=['/home/cayley/Desktop/Thesis/Programmation/AAA_Library/'],

#libraries=["cfitsio", "quadmath", "matrix"],

