#!/usr/bin/env python
import numpy, os, setuptools

setuptools.setup(name='pyfalcon',
    description='Python interface for the fast-multipole gravity solver falcON from NEMO',
    ext_modules=[setuptools.Extension(
    name='pyfalcon',
    sources=['pyfalcon.cpp',
    'src/src/exception.cc',
    'src/src/numerics.cc',
    'src/src/basic.cc',
    'src/src/body.cc',
    'src/src/gravity.cc',
    'src/src/kernel.cc',
    'src/src/tree.cc'],
    include_dirs=[numpy.get_include(), 'src/inc', 'src/inc/utils'],
    define_macros=[('falcON_SINGLE', None)],
    ) ],
    requires=['numpy'] )
