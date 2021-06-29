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

from models.mincde import mincde
from models.mincde_5r import MinCDE5R
import scipy.fftpack
import numpy
import unittest
import spatialpy

class Membrane(spatialpy.Geometry):
    def inside(self,x,on_boundary):
        return on_boundary
class Cytosol(spatialpy.Geometry):
    def inside(self,x,on_boundary):
        return not on_boundary
class MeshSize(spatialpy.DataFunction):
    def __init__(self,mesh):
        spatialpy.DataFunction.__init__(self, name="MeshSize")
        self.mesh = mesh
        self.h = mesh.get_mesh_size()

    def map(self, x):
        ret = self.h[self.mesh.closest_vertex(x)]
        return ret

class TestMinCDE(unittest.TestCase):

    def test_mincde_oscillation_period(self):
        """ Check that the MinCDE model is producing oscillation of the right period. """
        model = mincde()
        result = model.run()
        mindm = result.get_species("MinD_m")
        y_vals = model.mesh.coordinates()[:, 1]
        idx = (y_vals < 1e-6)
        mindmsum = numpy.sum(mindm[:,idx],axis=1)
        mindfft = scipy.fftpack.fft(mindmsum[200:]-numpy.mean(mindmsum[200:]))
        N = len(mindfft)
        T = model.tspan[-1] - model.tspan[200]
        mindpsd = numpy.abs(mindfft[:numpy.floor((N-1)/2)])
        mindfreq = numpy.arange(len(mindpsd), dtype=float)/T
        mind_max_period = 1/mindfreq[1+numpy.argmax(mindpsd[1:])]
        print(mind_max_period)
        self.assertTrue(mind_max_period > 50)
        self.assertTrue(mind_max_period < 70)

    def test_mincde_oscillation_period_fange(self):
        """ Check that the MinCDE model according to Fange et. al. is producing oscillation of the right period (approximately 28s).
            This is the same test as used in the URDME paper based on the legacy URDME Matlab interface. """
        model = MinCDE5R()
        result = model.run()
        mindm = result.get_species("MinD_m")
        y_vals = model.mesh.coordinates()[:, 1]
        idx = (y_vals < 1e-6)
        mindmsum = numpy.sum(mindm[:,idx],axis=1)
        mindfft = scipy.fftpack.fft(mindmsum[200:]-numpy.mean(mindmsum[200:]))
        N = len(mindfft)
        T = model.tspan[-1] - model.tspan[200]
        mindpsd = numpy.abs(mindfft[:numpy.floor((N-1)/2)])
        mindfreq = numpy.arange(len(mindpsd), dtype=float)/T
        mind_max_period = 1/mindfreq[1+numpy.argmax(mindpsd[1:])]
        print(mind_max_period)
        self.assertTrue(mind_max_period > 26)
        self.assertTrue(mind_max_period < 30)


if __name__ == '__main__':
    unittest.main()
