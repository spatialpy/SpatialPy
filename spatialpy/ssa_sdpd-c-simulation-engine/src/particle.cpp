/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>

// Include ANN KD Tree
#include "ANN/ANN.h"
#include "particle.hpp"
#include "particle_system.hpp"

namespace Spatialpy{
    
    std::mutex mtx ;

    ParticleSystem::ParticleSystem(size_t num_types, size_t num_chem_species, size_t num_chem_rxns,
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn):
                            num_types(num_types), num_chem_species(num_chem_species),
                            num_chem_rxns(num_chem_rxns), num_stoch_species(num_stoch_species),
                            num_stoch_rxns(num_stoch_rxns), num_data_fn(num_data_fn){
        dimension = 3;
        boundary_conditions[0] = 'n';
        boundary_conditions[1] = 'n';
        boundary_conditions[2] = 'n';
        static_domain = 0;
        gravity = (double*) calloc(3,sizeof(double));
        kdTree_initialized = false;
    }

    void ParticleSystem::add_particle(Particle *me){
    	// x_index.push(me) ;
    	particles.push_back(*me) ;
    }
    ParticleSystem::~ParticleSystem(){}

    Particle::Particle(ParticleSystem *sys, unsigned int id):sys(sys), id(id){
    	nu = 0.01;
    	mass = 1;
    	rho = 1;
    	solidTag = 0;
    	x[0] = x[1] = x[2] = 0.0;
    	v[0] = v[1] = v[2] = 0.0;
    	Q = (double*) calloc(sys->num_chem_species, sizeof(double));
    	C = (double*) calloc(sys->num_chem_species, sizeof(double));
    	data_fn = (double*) calloc(sys->num_data_fn, sizeof(double));
        //Q = std::make_shared<double[]>(sys.num_chem_species) ;
        //C = std::make_shared<double[]>(sys.num_chem_species) ;
        //data_fn = std::make_shared<double[]>(sys.num_data_fn) ;
    }

    NeighborNode::NeighborNode(Particle *data, double dist, double dWdr, double D_i_j):data(data), dist(dist), dWdr(dWdr), D_i_j(D_i_j){}
    //EventNode::EventNode(Particle *data, double tt):data(data), tt(tt){}

    void Particle::check_particle_nan(){
        if(
            isnan(x[0]) || !isfinite(x[0]) ||
            isnan(x[1]) || !isfinite(x[1]) ||
            isnan(x[2]) || !isfinite(x[2]) ||
            isnan(v[0]) || !isfinite(v[0]) ||
            isnan(v[1]) || !isfinite(v[1]) ||
            isnan(v[2]) || !isfinite(v[2]) ||
            isnan(rho)  || !isfinite(rho) ){
            printf("ERROR: nan/inf detected!!!\n");
            printf("number of neighbors: %li\n", neighbors.size()) ;
            printf("id=%i\n",id);
            printf("x[0]=%e\n",x[0]);
            printf("x[1]=%e\n",x[1]);
            printf("x[2]=%e\n",x[2]);
            printf("v[0]=%e\n",v[0]);
            printf("v[1]=%e\n",v[1]);
            printf("v[2]=%e\n",v[2]);
            printf("vt[0]=%e\n",vt[0]);
            printf("vt[1]=%e\n",vt[1]);
            printf("vt[2]=%e\n",vt[2]);
            printf("F[0]=%e\n",F[0]);
            printf("F[1]=%e\n",F[1]);
            printf("F[2]=%e\n",F[2]);
            printf("Fbp[0]=%e\n",Fbp[0]);
            printf("Fbp[1]=%e\n",Fbp[1]);
            printf("Fbp[2]=%e\n",Fbp[2]);
            printf("old_x[0]=%e\n",old_x[0]);
            printf("old_x[1]=%e\n",old_x[1]);
            printf("old_x[2]=%e\n",old_x[2]);
            printf("old_v[0]=%e\n",old_v[0]);
            printf("old_v[1]=%e\n",old_v[1]);
            printf("old_v[2]=%e\n",old_v[2]);
            printf("sys->dimension: %i\n", sys->dimension) ;
            printf("sys->current_step=%i\n",sys->current_step);
            printf("sys->dt=%e\n",sys->dt);
            exit(1);

        }
        old_x[0] =  x[0];
        old_x[1] =  x[1];
        old_x[2] =  x[2];
        old_v[0] =  v[0];
        old_v[1] =  v[1];
        old_v[2] =  v[2];
        old_rho =  rho;
    }

    double Particle::particle_dist(Particle *p2){
        double a = x[0] - p2->x[0];
        double b = x[1] - p2->x[1];
        double c = x[2] - p2->x[2];
        return sqrt( a*a + b*b + c*c);
    }

    double Particle::particle_dist_sqrd(Particle *p2){
        double a = x[0] - p2->x[0];
        double b = x[1] - p2->x[1];
        double c = x[2] - p2->x[2];
        return ( a*a + b*b + c*c);
    }

    int Particle::add_to_neighbor_list(Particle *neighbor, ParticleSystem *system, double r2){
     //    double a = x[0] - neighbor.x[0];
    	// double b = x[1] - neighbor.x[1];
    	// double c = x[2] - neighbor.x[2];
    	// double r2 =  ( a*a + b*b + c*c);
    	// Make sure the distance was actually set by the annkFRSearch
        //double r2_old = r2;
        if(r2 == ANN_DIST_INF) {
            r2 = particle_dist_sqrd(neighbor) ;
        }
        double r = sqrt(r2);

    	if(r > system->h){ return 0; } // do not add, out side support radius

    	// calculate dWdr
    	double h = system->h;
    	double R = r / h;
    	double alpha = 105 / (16 * M_PI * h * h * h); // 3D
    	double dWdr = alpha * (-12 * r / (h * h)) * ((1 - R)* (1 - R));
    	// calculate D_i_j

    	// Eq (13-14), Drawert et al 2019
    	double ih = 1.0 / h;
    	double ihsq = ih * ih;
    	double dhr = h - r;
    	double wfd = -25.066903536973515383e0 * dhr * dhr * ihsq * ihsq * ihsq * ih; //3D
    	// Eq 28 of Drawert et al 2019, Tartakovsky et. al., 2007, JCP
    	double D_i_j = -2.0*(mass*neighbor->mass)/(mass+neighbor->mass)*(rho+neighbor->rho)/(rho*neighbor->rho) * r2 * wfd / (r2+0.01*h*h);

    	if(isnan(D_i_j)){
    	    printf("Got NaN calculating D_i_j for me=%i, neighbor=%i\n",id, neighbor->id);
    	    printf("r2=%e ",r2);
    	    printf("r2_old=%e ",r2);
    	    printf("r=%e ",r);
    	    printf("h=%e ",h);
    	    printf("alpha=%e ",alpha);
    	    printf("dWdr=%e ",dWdr);
    	    printf("mass=%e ",mass);
    	    printf("rho=%e ",rho);
    	    printf("n->mass=%e ", neighbor->mass);
    	    printf("n->rho=%e ", neighbor->rho);
            printf("x=[%e,%e,%e] ",this->x[0],this->x[1],this->x[2]);
            printf("n->x=[%e,%e,%e] ",neighbor->x[0],neighbor->x[1],neighbor->x[2]);

    	    exit(1);
    	}
    	neighbors.emplace_back(neighbor, r, dWdr, D_i_j) ;

    	return 1;
    }

    void Particle::get_k_cleanup(ANNidxArray nn_idx, ANNdistArray dists) {
        delete [] nn_idx;
        delete [] dists;
    }

    void Particle::search_cleanup(ANNpoint queryPt, ANNidxArray nn_idx, ANNdistArray dists) {
        //printf("***TOP OF SEARCH CLEANUP***\n") ;
        get_k_cleanup(nn_idx, dists);
        annDeallocPt(queryPt);
        //printf("***BOTTOM OF SEARCH CLEANUP***\n") ;
    }

    int Particle::get_k__approx(ParticleSystem *system) {
        return system->kdTree->nPoints();
    }

    int Particle::get_k__exact(ANNpoint queryPt, ANNdist dist, ANNkd_tree *tree) {
        //printf("***TOP OF GET_K__EXACT***\n") ;
        //printf("USING DIST %f...", dist) ;
        for(int i = 0; i < 3; i++){
        }
        mtx.lock();
        int k = tree->annkFRSearch(queryPt, dist, 0);
        mtx.unlock(); 
        //printf("Found %i neighbors!\n", k) ;
        return k;
    }

    void Particle::find_neighbors(ParticleSystem *system, bool use_exact_k){
    	//clean out previous neighbors
    	//printf("find_neighbors.empty_linked_list\n");
    	//TODO: empty_neighbor_list(neighbors);
    	// search for points forward: (assumes the list is sorted ascending)
    	//printf("find_neighbors.search forward\n");

        // ANN KD Tree k fixed radius nearest neighbor search
        ANNpoint queryPt = annAllocPt(system->dimension);
        for(int i = 0; i < system->dimension; i++) {
            queryPt[i] = x[i];
        }
        // Squared radius
        ANNdist dist = system->h * system->h;
        // Number of neighbors to search for
        int k;
        if(use_exact_k) {
            k = get_k__exact(queryPt, dist, system->kdTree);
        }else{
            k = get_k__approx(system);
        }
        // Holds indicies that identify the neighbor in system.kdTree_pts
        ANNidxArray nn_idx = new ANNidx[k];
        // Holds squared distances to the neighbor (r2)
        ANNdistArray dists = new ANNdist[k];
        // Search for k nearest neighbors in the fixed squared radius dist from queryPt
        /*
        printf("Searching for Neighbors...\n") ;
        printf("Point at (%f, %f, %f), dist %f, finding %i neighbors out of %i indices and %i dists\n",
                    queryPt[0], queryPt[1], queryPt[2], dist, k, sizeof(nn_idx), sizeof(dists)) ;
        printf("CURRENT SIZE OF PARTICLES: %i\n", system->particles.size()) ;
        printf("NPOINTS IN KDTREE: %i\n", system->kdTree->nPoints()) ;
        printf("SEARCH RETURNED %i\n", system->kdTree->annkFRSearch(queryPt, dist, k, nn_idx, dists));
        printf("Neighbor search complete!\n") ;
        */
        mtx.lock();
        system->kdTree->annkFRSearch(queryPt, dist, k, nn_idx, dists);
        mtx.unlock();
        neighbors.clear() ;
        for(int i = 0; i < k && nn_idx[i] != ANN_NULL_IDX; i++) {
            //printf("Adding neighbor %i to list.  Index of neighbor: %i with dist: %f\n", i, nn_idx[i], dists[i]) ;
            Particle *neighbor = &system->particles[nn_idx[i]];
            add_to_neighbor_list(neighbor, system, dists[i]);
            //printf("ADDED %i with dist %f!\n", nn_idx[i], dists[i]) ;
            if(debug_flag > 2) {
                printf("find_neighbors(%i) forward found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
                    id, neighbor->id, sqrt(dists[i]),
                    x[0] - neighbor->x[0],
                    x[1] - neighbor->x[1],
                    x[2] - neighbor->x[2]);
                }
            }
        // Cleanup after the search
        search_cleanup(queryPt, nn_idx, dists);
        }

        // Brute force neighbor look-up
        //  for(Particle n : system.x_index){
    	//     if(n.data.x[0] > (x[0] + system.h)) break; //stop searching
    	//     if( (n.data.x[1] > (x[1] + system.h)) || (n.data.x[1] < (x[1] - system.h) ) ) continue;
    	// 	//TODO: Should these be breaks instead of continue?
    	//     if( (n.data.x[2] > (x[2] + system.h)) || (n.data.x[2] < (x[2] - system.h) ) ) continue;
    	//     add_to_neighbor_list(n, system);
    	//     if(debug_flag>2){
    	// 	printf("find_neighbors(%i) forward found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
    	// 	    id,n.data.id, particle_dist(n.data),
    	// 	    x[0] - n.data.x[0],
    	// 	    x[1] - n.data.x[1],
    	// 	    x[2] - n.data.x[2]);
    	//     }
    	// }

    	// //search for points backward
    	// for(std::vector<Particle>::iterator n = x_index.end(); n!=x_index.begin(); --n){
    	//     if(n->data->x[0] < (x[0] - system->h)){
    	// 	break; //stop searching backward
    	//     }
    	//     if( (n->data.x[1] > (x[1] + system.h)) || (n->data.x[1] < (x[1] - system.h) ) ){
    	// 	continue;
    	//     }
    	//     if( (n->data.x[2] > (x[2] + system.h)) || (n->data.x[2] < (x[2] - system.h) ) ){
    	// 	continue;
    	//     }
    	//     add_to_neighbor_list(n->data, system);
    	//     if(debug_flag>2){
    	// 	printf("find_neighbors(%i) backwards found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
    	// 	    id,n->data.id, particle_dist(n->data),
    	// 	    x[0] - n->data.x[0],
    	// 	    x[1] - n->data.x[1],
    	// 	    x[2] - n->data.x[2]);
    	//     }
    	// }

}
