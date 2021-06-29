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

#__all__=['model','spatialpy']
#from .spatialpy import *

import sys
if (sys.version_info < (3,0)):
    raise Exception("SpatialPy only works in Python 3.0 and higher")

from spatialpy.Model import *
from spatialpy.Solver import *
from spatialpy.Geometry import *
from spatialpy.Domain import *
from spatialpy.DataFunction import DataFunction
from spatialpy.InitialCondition import *
from spatialpy.BoundaryCondition import BoundaryCondition
from spatialpy.VTKReader import *
