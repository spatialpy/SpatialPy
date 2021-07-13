import os
import numpy
import sys
from collections import OrderedDict
sys.path.append("../..")
import spatialpy
# Global Constants
MAX_X_DIM = 5.0
MIN_X_DIM = -5.0
TOL = 1e-9

class Edge1(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MAX_X_DIM) < 0.05
class Edge2(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MIN_X_DIM) < 0.05
class Middle(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return abs(x[0] - MIN_X_DIM) >= 0.05
class cylinderDemo3D(spatialpy.Model):
    def __init__(self, model_name="cylinder_demo3d"):
        spatialpy.Model.__init__(self, model_name)

        self.timestep_size = 1

        # System constants
        D_const = 0.1

        # Define Species
        A = spatialpy.Species(name="A", diffusion_constant=D_const)
        B = spatialpy.Species(name="B", diffusion_constant=D_const)
        self.add_species([A, B])

        # Define Geometry
        self.domain = spatialpy.Domain.read_xml_mesh('cylinder.xml')

        # Define Types
        self.set_type(Middle(), 1)
        self.set_type(Edge1(), 2)
        self.set_type(Edge2(), 3)

        # Restrict the movement of Chemical Species
        self.restrict(A,[1,2])
        self.restrict(B,[1,3])

        vol = self.domain.get_vol()
        type_id = self.domain.type
        left = numpy.sum(vol[type_id == 2])
        right = numpy.sum(vol[type_id == 3])
        print("left "+str(left)+" right "+str(right))

        k_react = spatialpy.Parameter(name="k_react", expression=1.0)
        k_creat1 = spatialpy.Parameter(name="k_creat1",
                                     expression=100/left)
        k_creat2 = spatialpy.Parameter(name="k_creat2",
                                     expression=100/right)
        self.add_parameter([k_react, k_creat1,k_creat2])


        # Define Reactions
        R1 = spatialpy.Reaction(reactants=None, products={A:1},
                                rate=k_creat1, restrict_to=2)
        R2 = spatialpy.Reaction(reactants=None, products={B:1},
                              rate=k_creat2, restrict_to=3)
        R3 = spatialpy.Reaction(reactants={A:1, B:1}, products=None,
                              rate=k_react)
        self.add_reaction([R1, R2, R3])

        # Define simulation timespan
        #self.set_timesteps(1, 200)
        self.timespan(range(500))

if __name__=="__main__":
    model = cylinderDemo3D()
    result = model.run()
    A_sum = numpy.sum(result.get_species("A"), axis=1)
    B_sum = numpy.sum(result.get_species("B"), axis=1)
    print(A_sum[-1])
    print(B_sum[-1])
