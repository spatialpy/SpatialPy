/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef simulate_rdme_h
#define simulate_rdme_h
#include "particle.h"
#include "propensities.h"


typedef struct __rdme_data_t rdme_t;

struct __rdme_data_t {
    //size_t *irD; // replaced by neighbor_node_t->D_i_j
    //size_t *jcD;
    //double *prD;
    const size_t *irN;
    const size_t *jcN;
    const int *prN;
    const size_t *irG;
    const size_t *jcG;
    //int num_subdomains;  // system->num_types
    //const double*subdomain_diffusion_matrix;  // in system
    //size_t Ncells; // system->particle_list->count
    //size_t Mspecies;  //   system->num_stoch_species
    //size_t Mreactions; //  system->num_stoch_rxn
    //size_t Ndofs;  // deprecated
    // unsigned int *xx; // moved to particle_t->xx
    int initialized;
    //PropensityFun *rfun; // system->stoch_rxn_propensity_functions
    //double *srrate; // moved to rdme_voxel_t
    //double *rrate;
    //double *sdrate;
    //double *Ddiag;
    //double *rtimes; // using  ordered_list_t*heap 
    //int *node;
    //int *heap;
    ordered_list_t*heap; //TODO: check initialization of this


    long int total_reactions;
    long int total_diffusion;
    //char** species_names; //  system->species_names
};

typedef struct __rdme_voxel_t rdme_voxel_t;
struct __rdme_voxel_t {
    double srrate;
    double* rrate;
    double sdrate;
    double* Ddiag;
}



void initialize_rdme(system_t*system, const int Ncells, const int Mspecies, const int Mreactions, 
                        size_t *irN, size_t *jcN,int *prN,size_t *irG,size_t *jcG,
                        const char* const species_names[], const unsigned int*u0,
                        const int num_subdomains, const double*subdomain_diffusion_matrix);
void simulate_rdme(system_t*system, unsigned int step);
void destroy_rdme(system_t*system);


/******************************************************************/

void nsm_core__create(system_t*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG);
void nsm_core__destroy(rdme_t*rdme);

void nsm_core__initialize_chem_populations(rdme_t* rdme, const unsigned int*u0);

void nsm_core__initialize_rxn_propensities(rdme_t* rdme);
void nsm_core__initialize_diff_propensities(rdme_t* rdme);
void nsm_core__initialize_heap(rdme_t* rdme);


void nsm_core__build_diffusion_matrix(rdme_t*rdme,system_t*system);
void nsm_core__destroy_diffusion_matrix(rdme_t*rdme);

void nsm_core__take_step(system_t*system, double current_time, double step_size);



#endif /* simulate_rdme_h */

