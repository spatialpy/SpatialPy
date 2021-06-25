/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef output_h
#define output_h

#include "particle_system.hpp"

namespace Spatialpy{
    void output_csv(ParticleSystem*system, int current_step);
    void output_vtk__sync_step(ParticleSystem*system, int current_step);
    void output_vtk__async_step(ParticleSystem *system);
}
#endif // output_h
