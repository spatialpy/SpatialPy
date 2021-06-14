#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt
import numpy
import sys
sys.path.append('../..')
import spatialpy
print(spatialpy.__file__)


class All(spatialpy.Geometry):
    def inside(self, x, on_boundary):
        return True


class Walls(spatialpy.Geometry):
    ''' Outside of the unit square'''

    def inside(self, x, on_boundary):
        if x[0] < 0.0 or x[0] > 1.0 or x[1] < 0.0 or x[1] > 1.0:
            return True
        return False


class ChemicalGradient(spatialpy.Model):
    def __init__(self):
        spatialpy.Model.__init__(self, "ChemicalGradient")

        # System constants
        D = 0.1              # diffusion constant of the chemical species
        cLow = 50            # low value for boundary condition concentration
        cHigh = 100          # high value for boundary condition concentration

        nxF, nyF = 50, 50      # number of fluid particles in x and y-direction
        nu = 1.             # fluid viscosity
        # characteristic lenght of the cavity (= width = height)
        L = 1.
        nW = 2              # number of wall points
        rho = 1.             # fluid density

        # Discretization
        # total number of particles in x-direction (including walls)
        nxTot = nxF + 2*nW
        # total number of particles in y-direction (including walls)
        nyTot = nyF + 2*nW

        # Compute domain bounds (including the boundary)
        dx, dy = L/(nxF-1), L/(nyF-1)
        xLim = ((0-(nW-1)*dx), 1+(nW-1)*dx)
        yLim = ((0-(nW-1)*dy), 1+(nW-1)*dy)

        # Compute volume and mass per particle
        # in 2D simulations, consider z-lenght = 1
        vol = (xLim[1]-xLim[0])*(yLim[1]-yLim[0])*1.0
        # density * total volume / total number of particles
        mPP = rho*vol/(nxTot*nyTot)

        # Create mesh object
        print("xlim={0} ylim={1} nx={2} ny={3}".format(
            xLim, yLim, nxTot, nyTot))
        self.domain = spatialpy.Domain.create_2D_domain(
            xLim, yLim, nxTot, nyTot,
            type_id=1, mass=mPP, nu=nu, fixed=False)

        # Define Subdomains
        self.set_type(All(), 1)
        self.set_type(Walls(), 2, fixed=True)

        # Define a chemical Species
        self.add_species(spatialpy.Species('A', diffusion_constant=D))

        # Boundary conditions
        self.add_boundary_condition(spatialpy.BoundaryCondition(
            xmax=0.0,
            species='A',
            deterministic=True,
            value=cLow
        ))
        self.add_boundary_condition(spatialpy.BoundaryCondition(
            xmin=1.0,
            species='A',
            deterministic=True,
            value=cHigh
        ))

        # Time span
        # self.timespan(numpy.linspace(0,5,1000))
        self.timespan(numpy.linspace(0, 5, 6))


model = ChemicalGradient()

model.timestep_size = 1e-3   # timestep size
model.num_timesteps = math.ceil(
    model.tspan[-1]/model.timestep_size)    # number of timesteps
# frequency of outputting results
model.output_freq = math.ceil(
    (model.tspan[1]-model.tspan[0]) / model.timestep_size)


print("model.tspan", model.tspan)
print("model.timestep_size", model.timestep_size)
print("model.num_timesteps", model.num_timesteps)
print("model.output_freq", model.output_freq)


sol = spatialpy.Solver(model, debug_level=0)
sol.compile()

result = sol.run()


def plot_step(key='type', index=0, time=0):
    pts, data = result.read_step(time)
    plt.figure(figsize=(15, 10))
    if (key == 'v'):
        d = data[key]
        d = [d[i][index] for i in range(0, len(d))]
    else:
        d = data[key]
    plt.scatter(pts[:, 0], pts[:, 1], c=d)
    plt.axis('equal')
    plt.colorbar()
    plt.title('t={0}'.format(time))
    #plt.xticks(numpy.arange(-0.6, 0.7, 0.1))
    #plt.yticks(numpy.arange(-0.6, 0.7, 0.1))
    plt.grid(linestyle='--', linewidth=1)


def plot_all(key='type'):
    for i, t in enumerate(result.get_timespan()):
        plot_step(key, time=i)


result.get_timespan()

plot_step()

plot_all('C[A]')

plot_step('C[A]', time=1)

result.plot_species("A", -1, deterministic=True)

model.tspan

len(result.get_timespan())
