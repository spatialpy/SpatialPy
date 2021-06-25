/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef simulate_rdme_h
#define simulate_rdme_h

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
