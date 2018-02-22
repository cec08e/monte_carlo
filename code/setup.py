from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Heisenberg 3D',
  ext_modules = cythonize("heisenberg3d_opt_cy.pyx"),
)
