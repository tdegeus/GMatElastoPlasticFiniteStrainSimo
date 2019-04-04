desc = '''
Elasto-plastic material model.
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import pyxtensor

header = open('include/GMatElastoPlasticFiniteStrainSimo/config.h','r').read()
world  = re.split(r'(.*)(\#define GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split(r'(.*)(\#define GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split(r'(.*)(\#define GMATELASTOPLASTICFINITESTRAINSIMO_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GMatElastoPlasticFiniteStrainSimo',
    ['python/main.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      pyxtensor.get_include(False),
      pyxtensor.get_include(True ),
      pyxtensor.find_xtensor(),
      pyxtensor.find_xtl(),
      pyxtensor.find_xsimd(),
    ],
    language='c++'
  ),
]

setup(
  name             = 'GMatElastoPlasticFiniteStrainSimo',
  description      = 'Elasto-plastic material model',
  long_description = desc,
  version          = __version__,
  license          = 'MIT',
  author           = 'Tom de Geus',
  author_email     = 'tom@geus.me',
  url              = 'https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo',
  ext_modules      = ext_modules,
  install_requires = ['pybind11>=2.2.0'],
  cmdclass         = {'build_ext': pyxtensor.BuildExt},
  zip_safe         = False,
)
