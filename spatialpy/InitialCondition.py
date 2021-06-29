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

import numpy
from spatialpy.Model import ModelError


class InitialCondition():
    """ Class used to defined initial conditions in SpatialPy.
        SubClasses must implement the 'apply(model)' method, which
        direction modifies the model.u0[species,voxel] matrix.
    """

    def apply(self, model):
        raise ModelError("spatialpy.InitialCondition subclasses must implement apply()")

#TODO: implement InitialConditionFromResult()

class PlaceInitialCondition(InitialCondition):
    def __init__(self, species, count, location):
        self.species = species
        self.count = count
        self.location = location

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, at: {str(self.location)}"
        return print_string

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break
        vtx = model.domain.closest_vertex(self.location)
        model.u0[spec_ndx, vtx] += self.count

class UniformInitialCondition(InitialCondition):
    def __init__(self, species, count, types=None):
        self.species = species
        self.count = count
        self.types = types

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, Uniformly distrbuted in: {str(self.types)}"
        return print_string

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break

        if self.types is None:
            nvox = model.domain.get_num_voxels()
            for vtx in range(nvox):
                model.u0[spec_ndx, vtx] += self.count
        else:
            for i in range(model.domain.get_num_voxels()):
                if model.domain.type[i] in self.types:
                    model.u0[spec_ndx, i] += self.count


class ScatterInitialCondition(InitialCondition):

    def __init__(self, species, count, types=None):
        """ Scatter 'count' of 'species' randomly over the list of types
            (all types if None)."""
        self.species = species
        self.count = count
        self.types = types

    def __str__(self):
        print_string = f"{self.species.name}: {str(self.count)}, Scatter over: {str(self.types)}"
        return print_string

    def apply(self, model):
        spec_name = self.species.name
        for spec_ndx, spec_name in enumerate(model.listOfSpecies.keys()):
            if model.listOfSpecies[spec_name] == self.species: break

        if self.types is None:
            nvox = model.domain.get_num_voxels()
            for mol in range(self.count):
                vtx = numpy.random.randint(0, nvox)
                model.u0[spec_ndx, vtx] += 1
        else:
            allowed_voxels = []
            for i in range(model.domain.get_num_voxels()):
                if model.domain.type[i] in self.types:
                    allowed_voxels.append(i)
            nvox = len(allowed_voxels)
            if nvox==0: raise ModelError("ScatterInitialCondition has zero voxels to scatter in. Species={0} count={1} types={2}".format(self.species.name, self.count, self.types))
            for mol in range(self.count):
                v_ndx = numpy.random.randint(0, nvox)
                vtx = allowed_voxels[v_ndx]
                model.u0[spec_ndx, vtx] += 1
