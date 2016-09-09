from distutils.core import setup, Extension
import os
import glob

assert 0 == os.system('make -j pyobjs')

susi_native = Extension(
    'susi_native',
    sources = [],
    extra_objects = (
        glob.glob('build/src_python/*.o') +
        glob.glob('build/dss/*.o') +
        glob.glob('build/sdsl/*.o')
    ),
    libraries = ['stdc++'],
)

setup(
    name = 'susi',
    version = '1.33.7',
    description = 'A frontend for susi',
    ext_modules = [susi_native],
    packages = ['susi'],
    )
