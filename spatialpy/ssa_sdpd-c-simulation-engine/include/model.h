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
