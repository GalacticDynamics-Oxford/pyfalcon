import numpy, os, setuptools

nemopath = os.getenv('NEMO')
if nemopath is None:
    raise RuntimeError('NEMO environment variable is not found')
setuptools.setup(name='pyfalcon',
    description='Python interface for the fast-multipole gravity solver falcON from NEMO',
    ext_modules=[setuptools.Extension(
    name='pyfalcon',
    sources=['pyfalcon.cpp'],
    include_dirs=[numpy.get_include(), nemopath+'/usr/dehnen/falcON/inc', nemopath+'/usr/dehnen/utils/inc'],
    define_macros=[('falcON_SINGLE', None), ('falcON_NEMO', None)],
    extra_objects=[
        nemopath+'/usr/dehnen/falcON/lib/libfalcON.a',
        nemopath+'/usr/dehnen/utils/lib/libWDutils.a',
        nemopath+'/lib/libnemo.a'],
    ) ],
    requires=['numpy'] )
