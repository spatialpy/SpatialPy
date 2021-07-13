/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef particle_h
#define particle_h

typedef struct __particle_t particle_t;
typedef struct __system_t system_t;
typedef struct __bond_t bond_t;

#include "linked_list.h"
#include "simulate_rdme.h"

struct __particle_t {
    unsigned int id;
    int type;
    double x[3];
    double v[3];
    double vt[3];
    double mass;
    double rho;
    double nu;
    int solidTag;
    double bvf_phi;
    double normal[3];
    double F[3];
    double Frho;
    double Fbp[3];
    // Data Function
    double* data_fn;
    // chem_rxn_system
    unsigned int*xx; // populaion of discrete/stochastic species
    double *C;  // concentration of chem species
    double *Q;  // flux of chem species
    // below here for simulation
    node_t* x_index;  // 
    neighbor_list_t*neighbors;
    rdme_voxel_t*rdme;  // RDME solver data
};

struct __system_t {
    int dimension;
    double dt;
    unsigned int nt; 
    unsigned int current_step;
    unsigned int output_freq;
    double xlo, xhi, ylo, yhi, zlo, zhi;
    double h;
    double c0;
    double rho0;
    double P0;
    linked_list_t* particle_list;
    //linked_list_t* bond_list;
    linked_list_t* x_index;

    char boundary_conditions[3];
    int static_domain;
    size_t num_types;
    size_t num_data_fn;

    ChemRxnFun* chem_rxn_rhs_functions;
    size_t num_chem_species;
    size_t num_chem_rxns;

    rdme_t*rdme;
    PropensityFun* stoch_rxn_propensity_functions;
    size_t num_stoch_species;
    size_t num_stoch_rxns;

    const char * const* species_names; 

    const double *subdomain_diffusion_matrix;
    //int *stochic_matrix;
    int *stoichiometric_matrix;
    double* gravity;

};

//struct __bond_t {
//    unsigned int id;
//    particle_t* p1;
//    particle_t* p2;
//    double param_k;
//    double rest_distance;
//};


void find_neighbors(particle_t* me, system_t* system);
system_t* create_system(size_t num_types, size_t num_chem_species, size_t num_chem_rxns, 
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn);
particle_t* create_particle(int id);
void add_particle(particle_t* me, system_t* system);
double particle_dist(particle_t* p1, particle_t*p2);
double particle_dist_sqrd(particle_t* p1, particle_t*p2);

bond_t* create_bond(particle_t*p1, particle_t*p2, double k, double rest_distance);


// Global flags
extern int debug_flag;


#endif //particle_h

