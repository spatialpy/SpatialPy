/**
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
**/

/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef model_h
#define model_h

#include "particle.hpp"
#include "particle_system.hpp"

namespace Spatialpy{
    void filterDensity(Particle* me, ParticleSystem* system);

    void pairwiseForce(Particle* me, ParticleSystem* system);

    void computeBoundaryVolumeFraction(Particle* me, ParticleSystem* system);

    void applyBoundaryVolumeFraction(Particle* me, ParticleSystem* system);
}

#endif //model_h
