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
#ifndef simulate_rdme_h
#define simulate_rdme_h
#include "linked_list.h"
#include "particle.h"
#include "propensities.h"
#include  "dSFMT/dSFMT.h"

// Global flags
extern dsfmt_t dsfmt;

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
    ordered_list_t*heap; // ordered list of reaction/diffusion times


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
    ordered_node_t*heap_index;
};



void initialize_rdme(system_t*system, size_t *irN, size_t *jcN,int *prN,size_t *irG,size_t *jcG,
                        unsigned int*u0);
void simulate_rdme(system_t*system, unsigned int step);
void destroy_rdme(system_t*system);


/******************************************************************/

void nsm_core__create(system_t*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG);
void nsm_core__destroy(rdme_t*rdme);

void nsm_core__initialize_chem_populations(system_t*system, unsigned int*u0);

void nsm_core__initialize_rxn_propensities(system_t*system);
void nsm_core__initialize_diff_propensities(system_t*system);
void nsm_core__initialize_heap(system_t*system);


void nsm_core__build_diffusion_matrix(rdme_t*rdme,system_t*system);
void nsm_core__destroy_diffusion_matrix(rdme_t*rdme);

void nsm_core__take_step(system_t*system, double current_time, double step_size);



#endif /* simulate_rdme_h */

