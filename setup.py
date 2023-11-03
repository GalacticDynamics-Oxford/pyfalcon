#!/usr/bin/env python
import sys, numpy, os, setuptools

if sys.version_info[0]==3 and sys.version_info[1]>=10 and 'install' in sys.argv:
    print('If you are scared by a deprecation warning about running "setup.py install", try "pip install ." instead\n')

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
