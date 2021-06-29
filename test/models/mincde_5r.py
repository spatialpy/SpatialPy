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
""" spatialpy model file for the MinCDE example. """

import os.path
from spatialpy.InitialCondition import InitialCondition, ScatterInitialCondition
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

class MinCDE5R(spatialpy.Model):
    """ Model of MinCDE oscillations in E. Coli based on the model by Fange and Elf. """

    def __init__(self, model_name="mincde"):
        spatialpy.URDMEModel.__init__(self,model_name)

        # Species
        MinD_m     = spatialpy.Species(name="MinD_m", diffusion_constant=1e-14, D=2)
        MinD_c_adp = spatialpy.Species(name="MinD_c_adp", diffusion_constant=2.5e-12, D=3)
        MinD_c_atp = spatialpy.Species(name="MinD_c_atp", diffusion_constant=2.5e-12, D=3)
        MinD_e     = spatialpy.Species(name="MinD_e", diffusion_constant=2.5e-12, D=3)
        MinDE      = spatialpy.Species(name="MinDE", diffusion_constant=1e-14, D=2)

        self.add_species([MinD_m,MinD_c_atp,MinD_c_adp,MinD_e,MinDE])

        # Make sure that we have the correct path to the mesh file even if we are not executing from the basedir.
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.mesh = spatialpy.Mesh.read_xml_mesh(basedir+"/data/coli.xml")

        interior = dolfin.CellFunction("size_t",self.mesh)
        interior.set_all(1)
        boundary = dolfin.FacetFunction("size_t",self.mesh)
        boundary.set_all(0)

        # Mark the boundary points
        membrane = Membrane()
        membrane.mark(boundary,2)

        self.add_subdomain(interior)
        self.add_subdomain(boundary)

        # Average mesh size to feed into the propensity functions
        h = self.mesh.get_mesh_size()
        self.add_data_function(MeshSize(self.mesh))

        # Parameters
        NA = spatialpy.Parameter(name="NA", expression=6.022e23)
        sigma_d  = spatialpy.Parameter(name="sigma_d", expression=1.25e-8)
        sigma_dD = spatialpy.Parameter(name="sigma_dD", expression="9.0e6/(1000.0*NA)")
        sigma_e  = spatialpy.Parameter(name="sigma_e", expression="5.58e7/(1000.0*NA)")
        sigma_de = spatialpy.Parameter(name="sigma_de", expression=0.7)
        sigma_dt = spatialpy.Parameter(name="sigma_dt", expression=0.5)

        self.add_parameter([NA,sigma_d,sigma_dD,sigma_e,sigma_de,sigma_dt])

        # List of Physical domain markers that match those in the  Gmsh .geo file.
        interior = [1]
        boundary = [2]

        # Reactions
        R1 = spatialpy.Reaction(name="R1", reactants={MinD_c_atp:1}, products={MinD_m:1}, propensity_function="MinD_c_atp*sigma_d/MeshSize", restrict_to=boundary)
        R2 = spatialpy.Reaction(name="R2", reactants={MinD_c_atp:1,MinD_m:1}, products={MinD_m:2}, massaction=True, rate=sigma_dD)
        R3 = spatialpy.Reaction(name="R3", reactants={MinD_m:1,MinD_e:1}, products={MinDE:1}, massaction=True, rate=sigma_e)
        R4 = spatialpy.Reaction(name="R4", reactants={MinDE:1}, products={MinD_c_adp:1,MinD_e:1}, massaction=True, rate=sigma_de)
        R5 = spatialpy.Reaction(name="R5", reactants={MinD_c_adp:1}, products={MinD_c_atp:1}, massaction=True, rate=sigma_dt)

        self.add_reaction([R1,R2,R3,R4,R5])

        # Restrict to boundary
        self.restrict(MinD_m,boundary)
        self.restrict(MinDE,boundary)

        # Distribute molecules over the mesh according to their initial values
        self.add_initial_condition(spatialpy.ScatterInitialCondition(MinD_c_adp, 4000))
        self.add_initial_condition(spatialpy.ScatterInitialCondition(MinD_e, 1000))

        self.timespan(range(900))
