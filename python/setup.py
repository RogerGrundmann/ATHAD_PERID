from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ( "pyperid",
              ['pyperid.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args=["-std=c++17"],
              libraries = ['cond'],
              include_dirs = ['../atmosphere', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyperid',
    ext_modules = cythonize(extensions, language_level = "2"),
    )
