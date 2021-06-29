'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#!/usr/bin/env python3
# =============================================================================
# @file    setup.py
# @brief   SpatialPy setup file
# @license Please see the file named LICENSE in the project directory
# @website https://github.com/StochSS/SpatialPy
#
# Note: how to do a PyPI release
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run the following commands:
#
#   python3 setup.py sdist bdist_wheel
#   twine upload dist/*
#
# =============================================================================


from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install
import os
from os import path


# Read the contents of auxiliary files.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SETUP_DIR = path.dirname(os.path.abspath(__file__))

with open(path.join(SETUP_DIR, 'requirements.txt')) as f:
    reqs = f.read().rstrip().splitlines()

with open(path.join(SETUP_DIR, 'README.md'), 'r', errors = 'ignore') as f:
    readme = f.read()

# The following reads the variables without doing an "import handprint",
# because the latter will cause the python execution environment to fail if
# any dependencies are not already installed -- negating most of the reason
# we're using setup() in the first place.  This code avoids eval, for security.

version = {}
with open(path.join(SETUP_DIR, 'spatialpy/__version__.py')) as f:
    text = f.read().rstrip().splitlines()
    vars = [line for line in text if line.startswith('__') and '=' in line]
    for v in vars:
        setting = v.split('=')
        version[setting[0].strip()] = setting[1].strip().replace("'", '')


# Finally, define our namesake.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


setup(name                 = version['__title__'].lower(),
      version              = version['__version__'],
      author               = version['__author__'],
      author_email         = version['__email__'],
      maintainer           = version['__author__'],
      maintainer_email     = version['__email__'],
      license              = version['__license__'],
      url                  = version['__url__'],
      download_url         = version['__download_url__'],
      description          = version['__description__'],
      long_description     = readme,
      long_description_content_type = "text/markdown",
      keywords             = "biochemical simulation, Gillespie algorithm, stochastic simulation, biology, spatial simulation, RDME",
      project_urls         = {
          "Tracker": "https://github.com/StochSS/SpatialPy/issues",
          "Source" : "https://github.com/StochSS/SpatialPy",
      },
      packages             = find_packages('.'),
      include_package_data = True,
      install_requires     = reqs,

      classifiers      = [
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research'
      ],
)
