
/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef particle_hpp
#define particle_hpp

#include <vector>

#include "ANN/ANN.h" // ANN KD Tree

extern int debug_flag;

namespace Spatialpy{

    struct Particle;
    struct ParticleSystem;
    struct NeighborNode;
    struct EventNode;

    struct Particle{
        Particle(ParticleSystem *sys, unsigned int id=0);
        ParticleSystem *sys;
        unsigned int id;
        int type;
        double old_x[3];
        double old_v[3];
        double old_rho;
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
        double * data_fn;
        //std::shared_ptr<double[]> data_fn;
        // chem_rxn_system
        unsigned int*xx; // populaion of discrete/stochastic species
        double *C;
        double *Q;
        //std::shared_ptr<double[]> C;  // concentration of chem species
        //std::shared_ptr<double[]> Q;  // flux of chem species
        // below here for simulation
        std::vector<NeighborNode> neighbors;

        // Moved from rdme_voxel_t
        double srrate;
        double* rrate;
        double sdrate;
        double* Ddiag;
        //EventNode*heap_index;
        std::size_t particle_index;

        void check_particle_nan();

        double particle_dist(Particle *p2);
        double particle_dist_sqrd(Particle *p2);
        int add_to_neighbor_list(Particle *neighbor, ParticleSystem *system, double r2);

        // KD TREE FUNCTIONS
        void get_k_cleanup(ANNidxArray nn_idx, ANNdistArray dists);
        void search_cleanup(ANNpoint queryPt, ANNidxArray nn_idx, ANNdistArray dists);
        int get_k__approx(ParticleSystem *system);
        int get_k__exact(ANNpoint queryPt, ANNdist dist, ANNkd_tree *tree);
        void find_neighbors(ParticleSystem *system, bool use_exact_k=true);

        bool operator<(const Particle& p2){
            return x[0] > p2.x[0];
        }
    };
}

#endif
