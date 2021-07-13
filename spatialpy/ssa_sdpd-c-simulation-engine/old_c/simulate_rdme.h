/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef simulate_rdme_h
#define simulate_rdme_h
#include "linked_list.h"
#include "particle.h"
#include "propensities.h"
#include  "dSFMT/dSFMT.h"

// Global flags
extern dsfmt_t dsfmt;

// typedef struct __rdme_data_t rdme_t;

// Moved to ParticleSystem
// struct __rdme_data_t {
//     //size_t *irD; // replaced by neighbor_node_t->D_i_j
//     //size_t *jcD;
//     //double *prD;
//     const size_t *irN;
//     const size_t *jcN;
//     const int *prN;
//     const size_t *irG;
//     const size_t *jcG;
//     //int num_subdomains;  // system->num_types
//     //const double*subdomain_diffusion_matrix;  // in system
//     //size_t Ncells; // system->particle_list->count
//     //size_t Mspecies;  //   system->num_stoch_species
//     //size_t Mreactions; //  system->num_stoch_rxn
//     //size_t Ndofs;  // deprecated
//     // unsigned int *xx; // moved to particle_t->xx
//     int initialized;
//     //PropensityFun *rfun; // system->stoch_rxn_propensity_functions
//     //double *srrate; // moved to rdme_voxel_t
//     //double *rrate;
//     //double *sdrate;
//     //double *Ddiag;
//     //double *rtimes; // using  ordered_list_t*heap
//     //int *node;
//     //int *heap;
//     ordered_list_t*heap; // ordered list of reaction/diffusion times


//     long int total_reactions;
//     long int total_diffusion;
//     //char** species_names; //  system->species_names
// };

// // Moved to Particle
// typedef struct __rdme_voxel_t rdme_voxel_t;
// struct __rdme_voxel_t {
//     double srrate;
//     double* rrate;
//     double sdrate;
//     double* Ddiag;
//     ordered_node_t*heap_index;
// };



void initialize_rdme(ParticleSystem*system, size_t *irN, size_t *jcN,int *prN,size_t *irG,size_t *jcG,
                        unsigned int*u0);
void simulate_rdme(ParticleSystem*system, unsigned int step);
void destroy_rdme(ParticleSystem*system);


/******************************************************************/

void nsm_core__create(ParticleSystem*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG);
void nsm_core__destroy(ParticleSystem*system);

void nsm_core__initialize_chem_populations(ParticleSystem*system, unsigned int*u0);

void nsm_core__initialize_rxn_propensities(ParticleSystem*system);
void nsm_core__initialize_diff_propensities(ParticleSystem*system);
void nsm_core__initialize_heap(ParticleSystem*system);


void nsm_core__build_diffusion_matrix(ParticleSystem*system);
void nsm_core__destroy_diffusion_matrix(ParticleSystem*system);

void nsm_core__take_step(ParticleSystem*system, double current_time, double step_size);



#endif /* simulate_rdme_h */

