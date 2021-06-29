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

#include <time.h>
#include <random>
#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//#include "binheap.h"
#include "output.h"
#include "particle_system.hpp"
#include "propensities.hpp"
#include "simulate_rdme.hpp"


/**************************************************************************/
namespace Spatialpy{

    void initialize_rdme(ParticleSystem *system, size_t *irN, size_t *jcN,int *prN,
                        size_t *irG,size_t *jcG, unsigned int*u0){

        if(debug_flag){printf("initialize_rdme()\n");fflush(stdout);}
        //printf("nsm_core__create() BEGIN\n");fflush(stdout);
        nsm_core__create(system,irN,jcN,prN,irG,jcG);
        //printf("nsm_core__create() END\n");fflush(stdout);

        //printf("nsm_core__initialize_chem_populations() BEGIN\n");fflush(stdout);
        nsm_core__initialize_chem_populations(system, u0);
        //printf("nsm_core__initialize_chem_populations() END\n");fflush(stdout);

    }

    /**************************************************************************/
    // This function get called by the main simulation loop.  It advances the
    // state the of RDME system by dt
    void simulate_rdme(ParticleSystem*system,unsigned int step){
        // rdme_t*rdme = system->rdme;
        // if(rdme == NULL){
        //     return;
        // }
        if(!system->static_domain || !system->initialized){
    // All of below is replaced by code in find_neighbor()
    //            // if the  domain is not static, rebuild the diffusion matrix after movement
    //        if(!rdme->initialized){
    //          if(debug_flag) {
    //            printf("\tnsm_core__build_diffusion_matrix\n");
    //            nsm_core__build_diffusion_matrix(rdme,system);
    //          }
              system->initialized=1;
    //        }else{
    //            if(debug_flag) printf("Rebuilding diffusion matrix\n");
    //            if(debug_flag) printf("\tnsm_core__destroy_diffusion_matrix\n");
    //            nsm_core__destroy_diffusion_matrix(rdme);
    //            if(debug_flag) printf("\tnsm_core__build_diffusion_matrix\n");
    //            nsm_core__build_diffusion_matrix(rdme,system);
    //        }
            if(debug_flag) printf("\tnsm_core__initialize_rxn_propensities\n");
            nsm_core__initialize_rxn_propensities(system);
            if(debug_flag) printf("\tnsm_core__initialize_diff_propensities\n");
            nsm_core__initialize_diff_propensities(system);
            if(debug_flag) printf("\tnsm_core__initialize_heap\n");
            nsm_core__initialize_heap(system);
        }
        if(debug_flag) printf("Simulating RDME for %e seconds\n",system->dt);
        nsm_core__take_step(system, system->dt*step, system->dt);
    }
    /**************************************************************************/
    void destroy_rdme(ParticleSystem*system){
        if(debug_flag) printf("NSM: total # reacton events %lu\n",system->total_reactions);
        if(debug_flag) printf("NSM: total # diffusion events %lu\n",system->total_diffusion);
    }



    //===================================================
    // Adapted from PyURDME's nsmcore.c
    //===================================================

    void print_current_state(Particle*subvol, ParticleSystem*system){
        unsigned long int i;
        printf("Current state in voxel %i:\n",subvol->id);
        for(i=0;i<system->num_stoch_species;i++){
            printf("xx[%li] = %i\n",i,subvol->xx[i]);
        }
        printf("Neighbors:\n");
        Particle *p2;
        for(auto nn: subvol->neighbors){
            p2 = nn.data;
            printf("%i: nn->D_i_j=%e \n",p2->id,nn.D_i_j);
        }

    }

    /*void nsm_core(const size_t *irD,const size_t *jcD,const double *prD,
                  const size_t *irN,const size_t *jcN,const int *prN,
                  const size_t *irG,const size_t *jcG,
              const double *vol,
              const int *sd,
              const double *data,
                  const size_t Ncells,
                  const size_t Mspecies,
              const size_t Mreactions,
              const size_t dsize,
              int report_level
                  int *xx,
              double end_time)
    */
    /* Specification of the inputs:

     Ncells
     Number of subvolumes.

     Mspecies
     Number of species.

     Hence Ndofs = Ncells*Mspecies.

     Mreactions
     Total number of reactions.

     dsize
     Size of data vector sent to propensities.

     end_time
     length of the simulation.  This function simulates the RDME on the
     time interval [0, end_time].

     report_level
     The desired degree of feedback during simulations. 0, 1, and 2 are
     currently supported options.

     Diffusion matrix D. Double sparse (Ndofs X Ndofs).
     Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to
     dof #i. This matrix uses the CSR-format and not CSC because fast access to
     rows is needed.

     State vector 'xx'. Integer (Mspecies X Ncells).
     Gives the initial copy number of the species in each subvolume.
     !!! This is also the output.

     Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
     N(:,j) describes how reaction j changes the number of species.

     Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
     G(i,Mspecies+j) is non-zero if executing reaction j means that reaction i
     needs to be re-evaluated. The first Mspecies columns of G similarily cover
     diffusion events.
     vol. Double vector (length Ncells).
     vol[i] gives the volume of cell #i.

     data. Double matrix (dsize X Ncells).
     Generalized data matrix, data(:,j) gives a data vector for cell #j.

     sd. Integer vector (length Ncells).
     Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
     be used to separate boundaries, line segments and points.

     Format of sparse matrices:
     G, N and S are sparse matrices in compressed column format (CCS). D is sparse
     but in compressed row format (CRS), or equivalently, a transposed matrix in
     CCS format.
     jcD, irD, prD (double *)
     jcN, irN, prN (int *)
     jcG, irG (int *)

     Propensities:
     a vector of function pointers (length Mreactions) is input by
     linking with the prototypes in propensities.h and function
     definitions in a user-specified .c-file. The type of this vector is
     PropensityFun which defines the input to a property function. See
     propensities.h for more details.

     Ordering of the dofs:
     Dof #i is located in cell #(i/Mspecies), and the dofs located in
     cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
     Mspecies X Ncells.
     */

    /**************************************************************************/
    void nsm_core__create(ParticleSystem*system, size_t *irN, size_t *jcN,int *prN, size_t *irG, size_t *jcG){
        /* Create the RDME object */
        // rdme_t* rdme = (rdme_t*) malloc(sizeof(rdme_t));

        system->irN = irN;
        system->jcN = jcN;
        system->prN = prN;
        system->irG = irG;
        system->jcG = jcG;
        system->total_reactions = 0;
        system->total_diffusion = 0;
        system->initialized = 0;

        //NO MORE ORDERED LIST
        //system->heap = create_ordered_list();



        Particle*p;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i];
            // p->rdme = (rdme_voxel_t*) malloc(sizeof(rdme_voxel_t));
            p->srrate = 0;
            p->rrate = (double*) malloc(system->num_stoch_rxns * sizeof(double));
            p->sdrate = 0;
            p->Ddiag = (double*) malloc(system->num_stoch_species * sizeof(double));

            //p->heap_index = ordered_list_add(system->heap, p );

        }



        /* Create reaction rate matrix (Mreactions X Ncells) and total rate
         vector. In rrate we store all propensities for chemical rections,
         and in srrate the sum of propensities in every subvolume. */
        //rdme->rrate = (double *)malloc(system->num_stoch_rxn*rdme->Ncells*sizeof(double));
        //rdme->srrate = (double *)malloc(rdme->Ncells*sizeof(double));

        //nsm_core__initialize_rxn_propensities(rdme);

        /* Total diffusion rate vector (length Mcells). It will hold
         the total diffusion rates in each subvolume. */
        //rdme->sdrate = (double *)malloc(rdme->Ncells*sizeof(double));
        /* The diagonal value of the D-matrix is used frequently. For
         efficiency, we store the negative of D's diagonal in Ddiag. */
        //rdme->Ddiag = (double *)malloc(rdme->Ndofs*sizeof(double));

        //nsm_core__initialize_diff_propensities(rdme);

        /* Create binary (min)heap. */
        //rdme->rtimes = (double *)malloc(rdme->Ncells*sizeof(double));
        //rdme->node = (int *)malloc(rdme->Ncells*sizeof(int));
        //rdme->heap = (int *)malloc(rdme->Ncells*sizeof(int));


        /* return rdme structure */
       //return rdme;
       // system->rdme = rdme;
    }

    /**************************************************************************/
    void nsm_core__initialize_rxn_propensities(ParticleSystem*system){
        long unsigned int j;
        /* Calculate the propensity for every reaction and every
         subvolume. Store the sum of the reaction intensities in each
         subvolume in srrate. */
        Particle*p;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p = &system->particles[i];
            p->srrate = 0.0;
            for (j = 0; j < system->num_stoch_rxns; j++) {
                //rrate[i*Mreactions+j] =
                //(*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i],i,xx,irK,jcK,prK);
                //srrate[i] += rrate[i*Mreactions+j];
                //rdme->rrate[i*system->num_stoch_rxn+j] = (*rdme->rfun[j])(&rdme->xx[i*system->num_stoch_species],0.0,rdme->vol[i],&rdme->data[i*rdme->dsize],rdme->sd[i]);
                //rdme->srrate[i] += rdme->rrate[i*system->num_stoch_rxn+j];
                double vol = (p->mass / p->rho);
                p->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(p->xx,0.0,vol,p->data_fn,p->type);
                p->srrate += p->rrate[j];
            }
        }
    }

    /**************************************************************************/
    void nsm_core__initialize_diff_propensities(ParticleSystem* system){
        long unsigned int s_ndx;

    //    for (i = 0; i < rdme->Ndofs; i++) {
    //        rdme->Ddiag[i] = 0.0;
    //        for (j = rdme->jcD[i]; j < rdme->jcD[i+1]; j++)
    //        if (rdme->irD[j] == i) rdme->Ddiag[i] = -1*rdme->prD[j];
    //    }
    //
    //    /* Calculate the total diffusion rate for each subvolume. */
    //    for(i = 0; i < rdme->Ncells; i++) {
    //        rdme->sdrate[i] = 0.0;
    //        for(j = 0; j < system->num_stoch_species; j++)
    //        rdme->sdrate[i] += rdme->Ddiag[i*system->num_stoch_species+j]*rdme->xx[i*system->num_stoch_species+j];
    //    }

        Particle *p,*p2;
        double diff_const;
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p= &system->particles[i];
            if(p->neighbors.size() == 0){
                p->find_neighbors(system);
            }
            p->sdrate = 0.0;
            for(s_ndx=0; s_ndx<system->num_stoch_species; s_ndx++){
                p->Ddiag[s_ndx] = 0.0;  //Ddiag is sum of (diff_const*n2->D_i_j)
                for(auto n2=p->neighbors.begin(); n2!=p->neighbors.end(); ++n2){
                    p2 = n2->data;
                    diff_const = system->subdomain_diffusion_matrix[s_ndx*system->num_types + (p2->type-1)];
                    p->Ddiag[s_ndx] += diff_const*n2->D_i_j;
                }
                p->sdrate += p->Ddiag[s_ndx] * p->xx[s_ndx];
            }
        }
    }

    /**************************************************************************/
    void nsm_core__initialize_heap(ParticleSystem*system){
        /** TODO: UPDATE THIS FOR NEW EVENT DATA STRUCTURE **/
        /** TODO: INSERT DR. SANFT'S CODE HERE **/

        // rdme_t* rdme = system->rdme;
        /* Calculate times to next event (reaction or diffusion)
         in each subvolume and initialize heap. */

        std::vector<double> propensities(system->particles.size());
        double propensitySum=0.0;
        std::size_t activeChannels=0; // initial number of nonzero propensities (typically improves initial performance)

        long unsigned int p_i = 0 ;
        for(auto p = system->particles.begin(); p!=system->particles.end(); p++){
            propensities[p_i] = p->srrate + p->sdrate;
            propensitySum += propensities[p_i];
            if(propensities[p_i] > 0){
                activeChannels++ ;
            }
            if(p->particle_index != p_i){
                // particles can move around in the sys->particles vector,
                // this ensures that we know where in the vector each particle is
                p->particle_index = p_i;
            }
            p_i++ ;
        }
        /**
        for(long unsigned int i=0; i<system->particles.size(); i++){
        //for (i = 0; i < rdme->Ncells; i++) {
            //rdme->rtimes[i] = -log(1.0-dsfmt_genrand_close_open(&dsfmt))/(rdme->srrate[i]+rdme->sdrate[i]);
            //rdme->heap[i] = rdme->node[i] = i;
            Particle *p ;
            p = &system->particles[i] ;
            //p = e->data;

            //long unsigned int srng = rng() ;
            //long unsigned int rng_max = rng.max() ;
            //double tt = -log(1.0-(rng() * 1.0 / rng.max())) / (p->srrate+p->sdrate);
            propensities[i] = p->srrate + p->sdrate;
            propensitySum += propensities[i];
            if(propensities[i] > 0){
                activeChannels++
            }
            if(p->particle_index != i){
                // particles can move around in the sys->particles vector,
                // this ensures that we know where in the vector each particle is
                p->particle_index = i;
            }

            //system->event_v.emplace_back(p, tt) ;
        }
    **/
        //initialize_heap(rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
        // ordered_list_sort(system->heap);

        double timeOffset = system->current_step * system->dt;

        // TODO: does this deallocate memory on the 2nd call?  No -> check with Kevin, this might be ok
        // TODO: make a deallocation function
        if(system->num_stoch_rxns > 0){
            system->rdme_event_q.build(propensities, propensitySum, activeChannels,
                                   rng, timeOffset );
        }

    }

    /**************************************************************************/
    /**
    void nsm_core__destroy(rdme_t*rdme){
        // free(rdme);
    }
    **/
    /**************************************************************************/
    void nsm_core__initialize_chem_populations(ParticleSystem*system, unsigned int*u0){
        /* Set xx to the initial state. xx will always hold the current solution. */
        //printf("malloc Ndofs = %li\n",rdme->Ndofs);
        //rdme->xx = (unsigned int *)malloc(rdme->Ndofs*sizeof(unsigned int));
        //memcpy(rdme->xx,u0,rdme->Ndofs*sizeof(unsigned int));

        //printf("       Ndofs = %li\n",rdme->Ndofs);
        //printf("xx = [ ");
        //int i;
        //for(i=0;i<rdme->Ndofs;i++){
        //    printf("%u ",rdme->xx[i]);
        //}
        //printf("]\n");

        int num_s = system->num_stoch_species;
        int i = 0 ;
        for(auto p = system->particles.begin(); p!=system->particles.end(); p++){
        //for(unsigned long int i = 0; i < system->particles.size(); i++){
            //p1 = &system->particles[i];
            //memcpy(p1->xx,&u0[i*num_s],system->num_stoch_species*sizeof(unsigned int));
            p->xx = &u0[i*num_s];
            i++;
        }
    }


    /**************************************************************************/
    /**void nsm_core__build_diffusion_matrix(ParticleSystem*system){
        printf("*************** build_diffusion_matrix ***************\n");fflush(stdout);
        double off_diag_sum,diff_const,dist2;
        //NeighborNode *n2;
        Particle *p1,*p2;
        unsigned long int s_ndx;
        //double D_i_j;
        double ih,ihsq,wfd;
        double h = system->h;
        printf("System->h = %e\n",system->h);
        size_t Ncells = system->particles.size();
        size_t jcD_length = Ncells + 1;
        size_t irD_length = 0;
        size_t prD_length = 0;
        // find total length of jc & pr arrays: O(n)
        for(long unsigned int i = 0; i < system->particles.size(); i++){
            p1 = &system->particles[i];
            if(p1->neighbors.size() == 0){
                if(debug_flag){printf("find_neighbors(%i)\n",p1->id);}
                p1->find_neighbors(system);
            }
            //printf("node %i # neighbors %lu\n",p1->id,p1->neighbors->count);
            irD_length += (p1->neighbors.size() + 1);
        }
        prD_length = irD_length;
        printf("irD_length= %li\n",irD_length);fflush(stdout);
        printf("jcD_length= %li\n",jcD_length);fflush(stdout);
        // allocate space for each array
        printf("MALLOC irD [%li]\n",irD_length*system->num_stoch_species);
        size_t*irD = (size_t*) malloc(sizeof(size_t)*irD_length*system->num_stoch_species);
        size_t irD_ndx = 0;
        printf("MALLOC jcD [%li]\n",jcD_length*system->num_stoch_species);
        size_t*jcD = (size_t*) malloc(sizeof(size_t)*jcD_length*system->num_stoch_species);
        size_t jcD_ndx = 0;
        jcD[jcD_ndx++] = 0;
        printf("MALLOC prD [%li]\n",prD_length*system->num_stoch_species);
        double *prD = (double*) malloc(sizeof(double)*prD_length*system->num_stoch_species);
        size_t prD_ndx = 0;
        // for each particle p, look at each neighbor p2
        for(unsigned long int i = 0; i < system->particles.size(); i++){
            p1 = &system->particles[i];
            //if(p1->neighbors->count == 0){
                p1->find_neighbors(system);
            //}
            //Ordering is very inefficient here. Should do all species for a single p2. Requires reordering.
            for(s_ndx=0; s_ndx<system->num_stoch_species; s_ndx++){
                off_diag_sum = 0.0; // keep track of the total off diagonal sum
                for(auto n2=p1->neighbors.begin(); n2!=p1->neighbors.end(); ++n2){
                    p2 = n2->data;
                    //diff_const = system->subdomain_diffusion_matrix[s_ndx*system->num_subdomains + (p2->type-1)];
                    diff_const = system->subdomain_diffusion_matrix[s_ndx*system->num_types + (p2->type-1)];
                    dist2 = p1->particle_dist_sqrd(p2);
                    // Eq (13-14), Drawert et al 2019
                    ih = 1.0 / h;
                    ihsq = ih * ih;
                    wfd = h - sqrt(dist2);
                    if(wfd <= 0.0){
                        continue; // outside support of basis function
                    }
                    wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih; //3D
                    printf("wfd=%e\n",wfd);
                    // Eq 28 of Drawert et al 2019, Tartakovsky et. al., 2007, JCP
                    double tmp__D_i_j = -2.0*(p1->mass*p2->mass)/(p1->mass+p2->mass)*(p1->rho+p2->rho)/(p1->rho*p2->rho) * dist2 * wfd / (dist2+0.01*h*h);


                    printf("p1=%i p2=%i n2->D_i_j=%e D_i_j=%e\n",p1->id, p2->id, n2->D_i_j, tmp__D_i_j);
                    printf("n2->dist=%e\n",n2->dist);
                    printf("n2->dWdr=%e\n",n2->dWdr);
                    printf("dist2=%e sqr(dist2)=%e\n",dist2,sqrt(dist2));
                    printf("s_ndx=%li\n",s_ndx);
                    printf("p1->mass=%e p2->mass=%e\n",p1->mass,p2->mass);
                    printf("p1->rho=%e p2->rho=%e\n",p1->rho,p2->rho);

                    if(tmp__D_i_j != n2->D_i_j){
                        printf("error! D_i_j not right\n");
                        fflush(stdout);
                        exit(1);
                    }

                    if(diff_const > 0.0){
                        irD[irD_ndx++] = p2->id*system->num_stoch_species + s_ndx;
                        prD[prD_ndx++] = diff_const * n2->D_i_j;
                        off_diag_sum += diff_const * n2->D_i_j;
                    }
                }
                irD[irD_ndx++] = p1->id*system->num_stoch_species + s_ndx;
                prD[prD_ndx++] = -1*off_diag_sum;
                jcD[jcD_ndx++] = prD_ndx;
            }
        }
        printf("irD_ndx (%li) length irD (%li)\n",irD_ndx,irD_length*system->num_stoch_species);
        printf("jcD_ndx (%li) length jcD (%li)\n",jcD_ndx,jcD_length*system->num_stoch_species);
        printf("prD_ndx (%li) length prD (%li)\n",prD_ndx,prD_length*system->num_stoch_species);
        if( prD_ndx != irD_ndx){
                printf("Assembly: prD_ndx (%zu) != irD_ndx (%zu)\n",prD_ndx,irD_ndx);
        }
        if( irD_ndx != irD_length*system->num_stoch_species ){
            printf("Assembly: irD_ndx (%zu) != irD_length*Mspecies (%li)\n", irD_ndx, irD_length*system->num_stoch_species);
        }
        char filename[256];
        time_t seconds;
        size_t i;
        seconds = time(NULL);
        sprintf(filename,"diffusion_matrix_%ld", seconds);
        printf("Writing out diffusion matrix to '%s'\n",filename);
        FILE*fp = fopen(filename,"w+");
        fprintf(fp, "irD = [");
        for(i=0;i<irD_ndx;i++){
            if(i>0){ fprintf(fp,",");}
            fprintf(fp, "%zu",irD[i]);
        }
        fprintf(fp, "]\n");
        fprintf(fp, "jcD = [");
        for(i=0;i<jcD_ndx;i++){
            if(i>0){ fprintf(fp,",");}
            fprintf(fp, "%zu",jcD[i]);
        }
        fprintf(fp, "]\n");
        fprintf(fp, "prD = [");
        for(i=0;i<prD_ndx;i++){
            if(i>0){ fprintf(fp,",");}
            fprintf(fp, "%e",prD[i]);
        }
        fprintf(fp, "]\n");
        fprintf(fp, "D = scipy.sparse.csc_matrix(prD,irD,jcD)\n");
        fclose(fp);

    }**/

    /**************************************************************************/
    /**
    void nsm_core__destroy_diffusion_matrix(rdme_t*rdme){
        //free(rdme->irD);
        //free(rdme->jcD);
        //free(rdme->prD);
    }
    **/
    /**************************************************************************/
    // Update to use priority queue
    void nsm_core__take_step(ParticleSystem*system, double current_time, double step_size){

        if( system->num_stoch_rxns == 0 ){ return;}
        // rdme_t*rdme = system->rdme;
        double tt = current_time;
        double end_time = current_time + step_size;
        double totrate,cum,rdelta,rrdelta;
        int event,errcode = 0;
        long unsigned int re, spec = 0 ;
        Particle* subvol;
        size_t i,j = 0;
        //double old_rrate = 0.0,old_drate = 0.0;
        double rand1,rand2,cum2,old;
        double vol,diff_const;
        Particle *dest_subvol = NULL;
        NeighborNode *nn = NULL;


        // Check the integrety of the heap
    //    ordered_node_t*on,*bon,*aon;
    //    printf("========================================================\n");
    //    on= system->rdme->heap->head;
    //    bon=NULL;
    //    aon=system->rdme->heap->head->next;
    //    while(on!=NULL){
    //        if( on->prev != bon ||
    //            on->next != aon ){
    //            printf("next/prev mismatch\n");exit(1);
    //        }
    //        on = on->next;
    //        if(bon==NULL){bon=system->rdme->heap->head;}else{bon=bon->next;}
    //        if(aon!=NULL){aon=aon->next;}
    //    }
    //    if(system->rdme->heap->tail != bon){
    //        printf("tail mismatch\n");exit(1);
    //    }
    //    printf("========================================================\n");


        std::pair<double,int> timeRxnPair;
        int reactionIndex, subvol_index;
        /* Main loop. */
        while(tt <= end_time){

            /* Get the subvolume in which the next event occurred.
             This subvolume is on top of the heap. */
            //told = tt;
            //tt   = rdme->rtimes[0];
            //subvol = rdme->node[0];
            //subvol = system->event_v.front().data;
            //tt = system->event_v.front().tt;

            timeRxnPair = system->rdme_event_q.selectReaction();
            tt = timeRxnPair.first;
            subvol_index = timeRxnPair.second;
            subvol = &system->particles[subvol_index];
            vol = (subvol->mass / subvol->rho);
            if(debug_flag){printf("nsm: tt=%e subvol=%i\n",tt,subvol->id);}
            /* First check if it is a reaction or a diffusion event. */
            totrate = subvol->srrate + subvol->sdrate;

    //        if(totrate <= 0){ // Sanity check, is there a non-zero reaction and diffusion propensity
    //            totrate = rdme->srrate[subvol]+rdme->sdrate[subvol];
    //            if (totrate > 0.0)
    //                rdme->rtimes[0] = -log(1.0-dsfmt_genrand_close_open(&dsfmt))/totrate+tt;
    //            else
    //                rdme->rtimes[0] = INFINITY;
    //            /* Update the heap. */
    //            update(0,rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
    //            // go to the next element
    //            continue;
    //        }

            rand1 = rng() * 1.0 / rng.max();

            if (rand1 <= subvol->srrate/totrate) { // use normalized floating point comparision
                /* Reaction event. */
                event = 0;

                /* a) Determine the reaction re that did occur (direct SSA). */
                double rand_rval = rand1 * subvol->srrate;
                for (re = 0, cum = subvol->rrate[0]; re < system->num_stoch_rxns && rand_rval > cum; re++, cum += subvol->rrate[re]);
                if(re >= system->num_stoch_rxns){
                    if(cum != subvol->srrate){
                        printf("Reaction propensity mismatch in voxel %i. re=%li, srrate[subvol]=%e cum=%e rand_rval=%e\n",subvol->id,re,subvol->srrate,cum,rand_rval);
                        rdelta = 0.0;
                        for (j=0;j<system->num_stoch_rxns; j++) {
                            rdelta += (subvol->rrate[j] = (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type));
                        }
                        subvol->srrate = rdelta;
                    }
                    if(subvol->srrate == 0.0){ continue; }


                    double rand_rval2 = rand1 * subvol->srrate; // sum of propensitiess is not propensity sum, re-roll

                    for (re = 0, cum = subvol->rrate[0]; re < system->num_stoch_rxns && rand_rval2 > cum; re++, cum += subvol->rrate[re]);
                    if(re >= system->num_stoch_rxns){ // failed twice, problems!
                        printf("Propensity sum overflow, rand=%e rand_rval=%e rand_rval2=%e srrate[%i]=%e cum=%e\n",rand1,rand_rval,rand_rval2,subvol->id,subvol->srrate,cum);
                        exit(1);
                    }
                }
                if(debug_flag){printf("nsm: tt=%e subvol=%i type=%i ",tt,subvol->id,subvol->type);}
                if(debug_flag){printf("Rxn %li \n",re);}
                /* b) Update the state of the subvolume subvol and sdrate[subvol]. */
                for (i = system->jcN[re]; i < system->jcN[re+1]; i++) {
                    int prev_val = subvol->xx[system->irN[i]];
                    subvol->xx[system->irN[i]] += system->prN[i];
                    if (subvol->xx[system->irN[i]] < 0){
                        errcode = 1;
                        printf("Negative state detected after reaction %li, subvol %i, species %li at time %e (was %i now %i)\n",re,subvol->id,system->irN[i],tt,prev_val,subvol->xx[system->irN[i]]);

                        print_current_state(subvol,system);
                        exit(1);
                    }
                    subvol->sdrate += subvol->Ddiag[system->irN[i]]*system->prN[i];
                }

                /* c) Recalculate srrate[subvol] using dependency graph. */
                for (i = system->jcG[system->num_stoch_species+re], rdelta = 0.0; i < system->jcG[system->num_stoch_species+re+1]; i++) {
                    old = subvol->rrate[system->irG[i]];
                    j = system->irG[i];
                    rdelta +=
                    (subvol->rrate[j] =
                     //(*rdme->rfun[j])(&rdme->xx[subvol*system->num_stoch_species],tt,rdme->vol[subvol],&rdme->data[subvol*rdme->dsize],rdme->sd[subvol])

                     (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                     )-old;

                }
                subvol->srrate += rdelta;

                system->total_reactions++; /* counter */
            }
            else {
                /* Diffusion event. */
                event = 1;

                /* a) Determine which species... */
                double diff_rand = rand1 * subvol->sdrate;

                //for (spec = 0, dof = subvol*system->num_stoch_species, cum = rdme->Ddiag[dof]*rdme->xx[dof];
                //     spec < system->num_stoch_species && diff_rand > cum;
                //     spec++, cum += rdme->Ddiag[dof+spec]*rdme->xx[dof+spec]);
                for (spec = 0, cum = subvol->Ddiag[spec]*subvol->xx[spec];
                     spec < system->num_stoch_species && diff_rand > cum;
                     spec++, cum += subvol->Ddiag[spec]*subvol->xx[spec]);
                if(spec >= system->num_stoch_species){
                    //printf("Diffusion species overflow\n");
                    // try again, 'cum' is a better estimate of the propensity sum
                    if(cum != subvol->sdrate){
                        printf("Diffusion propensity mismatch in voxel %i. spec=%li, sdrate[subvol]=%e cum=%e diff_rand=%e\n",subvol->id,spec,subvol->sdrate,cum,diff_rand);
                        rdelta = 0.0;
                        for(j = 0; j < system->num_stoch_species; j++){
                            rdelta += subvol->Ddiag[j]*subvol->xx[j];
                        }
                        subvol->sdrate = rdelta;
                    }
                    if(subvol->sdrate == 0.0){ continue; }

                    diff_rand = cum *rand1;
                    for (spec = 0, cum = subvol->Ddiag[spec]*subvol->xx[spec];
                         spec < system->num_stoch_species && diff_rand > cum;
                         spec++, cum += subvol->Ddiag[spec]*subvol->xx[spec]);
                    if(spec >= system->num_stoch_species){
                        spec--;
                        while(subvol->xx[spec] <= 0){
                            spec--;
                            if(spec <=0){
                                printf("Error: diffusion event in voxel %i was selected, but no molecues to move\n",subvol->id);
                                print_current_state(subvol,system);
                                exit(1);
                            }
                        }
                    }
                }


                /* b) and then the direction of diffusion. */
                double r2 = rng() * 1.0 / rng.max();
                rand2 = r2 * subvol->Ddiag[spec];

                /* Search for diffusion direction. */
    //            for (i = rdme->jcD[col], cum2 = 0.0; i < rdme->jcD[col+1]; i++)
    //                if (rdme->irD[i] != col && (cum2 += rdme->prD[i]) > rand2)
    //                    break;
    //
    //            /* paranoia fix: */
    //            // This paranoia fix creates errors if the final rate has a zero propensity.  It can cause negative populations.
    //            if (i >= rdme->jcD[col+1]){
    //                //printf("Diffusion direction overflow\n");
    //                // try again, 'cum2' is a better estimate of propensity sum
    //                rand2 = r2*cum2;
    //                for (i = rdme->jcD[col], cum2 = 0.0; i < rdme->jcD[col+1]; i++)
    //                    if (rdme->irD[i] != col && (cum2 += rdme->prD[i]) > rand2)
    //                        break;
    //                if (i >= rdme->jcD[col+1]){
    //                    i--;
    //                }
    //            }
    //
    //            to_node = rdme->irD[i];
    //            to_vol = to_node/system->num_stoch_species;

                Particle*p2;
                cum2 = 0.0;
                for(auto nn=subvol->neighbors.begin(); nn!=subvol->neighbors.end(); ++nn){
                    p2 = nn->data;
                    diff_const = system->subdomain_diffusion_matrix[spec*system->num_types + (p2->type-1)];
                    cum2 += nn->D_i_j * diff_const;
                    if(cum2 > rand2){
                        dest_subvol = p2;
                        break;
                    }
                }
                if(dest_subvol==NULL){ //Diffusion direction overflow
                    // try again, 'cum2' is a better estimate of propensity sum
                    rand2 = r2*cum2;
                    cum2 = 0.0;
                    for(auto nn=subvol->neighbors.begin(); nn!=subvol->neighbors.end(); ++nn){
                        p2 = nn->data;
                        diff_const = system->subdomain_diffusion_matrix[spec*system->num_types + (p2->type-1)];
                        cum2 += nn->D_i_j * diff_const;
                        if(cum2 > rand2){
                            dest_subvol = p2;
                            break;
                        }
                    }
                    if(dest_subvol==NULL){
                        printf("Error: overflow in trying to determine which destination voxel subvol=%i\n",subvol->id);
                        printf("subvol->id=%u ",subvol->id);
                        printf("rand2=%e ",rand2);
                        printf("cum2=%e ",cum2);
                        fflush(stdout);
                        print_current_state(subvol,system);
                        exit(1);
                    }
                }


                /* c) Execute the diffusion event (check for negative elements). */
                subvol->xx[spec]--;
                if (subvol->xx[spec] < 0){
                        errcode = 1;
                        printf("Negative state detected after diffusion, voxel %i -> %i, species %li at time %e\n",subvol->id,dest_subvol->id,spec,tt);
                        print_current_state(subvol,system);
                        exit(1);
                }

                dest_subvol->xx[spec]++;


                if(debug_flag){printf("nsm: tt=%e subvol=%i type=%i ",tt,subvol->id,subvol->type);}
                if(debug_flag){printf("Diff %i->%i\n",subvol->id,dest_subvol->id);}

                /* Save reaction and diffusion rates. */
                //old_rrate = dest_subvol->srrate;
                //old_drate = dest_subvol->sdrate;
                /* Recalculate the reaction rates using dependency graph G. */
                if (system->num_stoch_rxns > 0){
                    for (i = system->jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < system->jcG[spec+1]; i++) {

    // herehere
                        j = system->irG[i];
                        // update this voxel's reactions
                        old = subvol->rrate[j];
                        rdelta +=
                          (subvol->rrate[j] =
                            (*system->stoch_rxn_propensity_functions[j])(subvol->xx,tt,vol,subvol->data_fn,subvol->type)
                          )-old;


                        // update dest voxel's reactions
                        old = dest_subvol->rrate[j];
                        rrdelta += (dest_subvol->rrate[j] =
                            (*system->stoch_rxn_propensity_functions[j])(dest_subvol->xx,tt,vol,dest_subvol->data_fn,dest_subvol->type)
                          )-old;
                    }

                    subvol->srrate += rdelta;
                    dest_subvol->srrate += rrdelta;
                }

                /* Adjust diffusion rates. */
                subvol->sdrate -= subvol->Ddiag[spec];
                dest_subvol->sdrate += dest_subvol->Ddiag[spec];

                system->total_diffusion++; /* counter */

            }

            /* Compute time to new event for this subvolume. */
            totrate = subvol->srrate+subvol->sdrate;
            //if(totrate > 0.0){
            //    system->event_v.front().tt = -log(1.0-(rng() * 1.0 / rng.max()))/totrate+tt;
            //}else{
            //    system->event_v.front().tt = INFINITY;
            //}

            /* Update the heap. */
            //update(0,rdme->rtimes,rdme->node,rdme->heap,rdme->Ncells);
            //ordered_list_bubble_up_down(system->event_v, subvol->heap_index);

            system->rdme_event_q.update(subvol_index, totrate, tt, rng);



            /* If it was a diffusion event, also update the other affected
             node. */
            if(event) {

                totrate = dest_subvol->srrate+dest_subvol->sdrate;
//
//                if(totrate > 0.0) {
//                    printf("dest_subvol->heap_index->tt: %f\n", dest_subvol->heap_index->tt) ;
//                    if(!isinf(dest_subvol->heap_index->tt)){
//                        dest_subvol->heap_index->tt =
//                        (old_rrate+old_drate)/totrate*(dest_subvol->heap_index->tt - tt)+tt;
//                    }else{
//                        // generate a new waiting time
//                        dest_subvol->heap_index->tt = -log(1.0-(rng() * 1.0 / rng.max()))/totrate+tt;
//                    }
//
//                }else{
//                    dest_subvol->heap_index->tt = INFINITY;
//                }
                //ordered_list_bubble_up_down(system->heap, dest_subvol->heap_index);

                system->rdme_event_q.update(dest_subvol->particle_index, totrate, tt, rng);
            }

            // re-sort the heap
            //ordered_list_sort(system->rdme->heap);  // this resorts the whole list

            /* Check for error codes. */
            if (errcode) {
                /* Cannot continue. Clear this solution and exit. */
                printf("Exiting due to errcode %i\n",errcode);
                print_current_state(subvol,system);
                exit(1);
            }
        }


    }
}
/**************************************************************************/
