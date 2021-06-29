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

/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef simulate_rdme_h
#define simulate_rdme_h

#include <random>

#include "particle_system.hpp"
#include "propensities.hpp"

namespace Spatialpy{

void initialize_rdme(ParticleSystem *system, size_t *irN, size_t *jcN,int *prN,
                     size_t *irG,size_t *jcG, unsigned int*u0);
void simulate_rdme(ParticleSystem*system, unsigned int step);
void destroy_rdme(ParticleSystem*system);


/******************************************************************/

void nsm_core__create(ParticleSystem*system, size_t *irN, size_t *jcN,int *prN,
                      size_t *irG, size_t *jcG);
void nsm_core__destroy(ParticleSystem*system);

void nsm_core__initialize_chem_populations(ParticleSystem*system, unsigned int*u0);

void nsm_core__initialize_rxn_propensities(ParticleSystem*system);
void nsm_core__initialize_diff_propensities(ParticleSystem*system);
void nsm_core__initialize_heap(ParticleSystem*system);


//void nsm_core__build_diffusion_matrix(ParticleSystem*system);
//void nsm_core__destroy_diffusion_matrix(ParticleSystem*system);

void nsm_core__take_step(ParticleSystem*system, double current_time,
                         double step_size);
}

#endif /* simulate_rdme_h */
