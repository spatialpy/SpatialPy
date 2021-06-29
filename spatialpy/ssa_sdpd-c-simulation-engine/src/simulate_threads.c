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
#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Include ANN KD Tree
#include "ANN/ANN.h"
#include "model.h"
#include "output.h"
#include "particle_system.hpp"
#include "simulate.hpp"
#include "simulate_rdme.hpp"
#include "pthread_barrier.h"

namespace Spatialpy{
    struct arg {
        ParticleSystem* system;
        unsigned int thread_id;
        unsigned int num_threads;
        unsigned int num_my_particles;
        //unsigned int num_my_bonds;
        int my_first_particle;
        //bond*my_first_bond;
    };

    pthread_barrier_t begin_step_barrier;
    pthread_barrier_t end_step_barrier;
    pthread_barrier_t begin_sort_barrier;
    pthread_barrier_t end_sort_barrier;
    pthread_barrier_t begin_output_barrier;
    pthread_barrier_t end_output_barrier;

    void* output_system_thread(void* targ){
        ParticleSystem* system = (ParticleSystem*) targ;
        while(1){
            pthread_barrier_wait(&begin_output_barrier);
            output_csv(system, system->current_step);
            if(debug_flag) printf("[OUT] start output_vtk__sync_step()\n");
            output_vtk__sync_step(system, system->current_step);
            if(debug_flag) printf("[OUT] done output_vtk__sync_step()\n");
            pthread_barrier_wait(&end_output_barrier);
            if(debug_flag) printf("[OUT] start output_vtk__async_step()\n");
            output_vtk__async_step(system);
            if(debug_flag) printf("[OUT] done output_vtk__async_step()\n");
        }
    }

    struct sarg {
        ParticleSystem *system;
        int sort_ndx;
    };

    void buildKDTree(ParticleSystem *system) {
        // cleanup KD Tree
        if(debug_flag) printf("\tstarting buildKDTree()\n");

        if(system->kdTree_initialized) {
            if(debug_flag){ printf("\tsystem->kdTree_initialized=True\n");}
            if(system->static_domain) {
                return;} // do not rebuild tree for static domains
            if(debug_flag){ printf("\tannDeallocPts()\n");}
            annDeallocPts(system->kdTree_pts);
            if(debug_flag){ printf("\tdeleting system->kdTree (no close)\n"); fflush(stdout);}
            delete system->kdTree; // NO [] on your delete!!!!
            if(debug_flag){ printf("\tdone deleting system->kdTree (no close)\n"); fflush(stdout);}
            //if(debug_flag) printf("\tannClose()\n");
            //annClose();
        }
        if(debug_flag){ printf("\tsystem->particles.size()\n");fflush(stdout);}
        int nPts = system->particles.size();
        if(debug_flag){ printf("\tannAllocPts()\n");fflush(stdout);}
        system->kdTree_pts = annAllocPts(nPts, system->dimension);
        for(int i = 0; i < nPts; i++) {
            for(int j = 0; j < system->dimension; j++) {
                system->kdTree_pts[i][j] = system->particles[i].x[j];
            }
        }
        if(debug_flag){ printf("\tsystem->kdTree = new ANNkd_tree()\n");}
        system->kdTree = new ANNkd_tree(system->kdTree_pts, nPts, system->dimension);
        system->kdTree_initialized = true;
    }

    void* sort_index_thread(void* targ_in){
        struct sarg *targ = (struct sarg*) targ_in;
        while(1){
            pthread_barrier_wait(&begin_sort_barrier);
            if(debug_flag) printf("[SORT] begin sort\n");
            // ANN KD Tree
            buildKDTree(targ->system);
            if(debug_flag) printf("[SORT] sort complete\n");
            pthread_barrier_wait(&end_sort_barrier);
        }
    }

    void* run_simulation_thread(void *targ_in){
        struct arg* targ = (struct arg*)targ_in;
        ParticleSystem *system = targ->system;
        unsigned int step;
        Particle* p;
        unsigned int i;
        int count = 0;
        // each thread will take a step with each of it's particles
        for(step=0; step < system->nt; step++){
            // block on the begin barrier
            if(debug_flag) printf("[WORKER %i] waiting to begin step %i\n",targ->thread_id,step);
            pthread_barrier_wait(&begin_step_barrier);
            //---------------------------------------
            unsigned int nsubsteps = get_number_of_substeps();
            for(unsigned int substep=0;substep < nsubsteps; substep++){
                // take_step
                count = 0;
                for(i=0; i<targ->num_my_particles; i++){
                    p = &system->particles[i+targ->my_first_particle] ;
                    take_step(p,system,step,substep);
                    count++;
                }
                // block on the end barrier
                if(debug_flag)printf("[WORKER %i] completed step %i, substep %i/%i, processed %i particles\n",targ->thread_id,step,substep,nsubsteps,count);
                pthread_barrier_wait(&end_step_barrier);
            }
            //---------------------------------------
            // compute_bond_forces
            /*
            count = 0;
            n=targ->my_first_bond;
            for(i=0; i<targ->num_my_bonds; i++){
                if(n==NULL) break;
                compute_bond_forces(n->data,system,step);
                count++;
                n=n->next;
            }
            // block on the end barrier
            if(debug_flag)printf("[WORKER %i] completed compute_bond_forces %i, processed %i particles\n",targ->thread_id,step,count);
            pthread_barrier_wait(&end_step_barrier);
            */
            //---------------------------------------
        } //end for(step)
        return NULL;
    }

    void run_simulation(int num_threads, ParticleSystem *system){
        //
        int i,j;
        int num_particles_per_thread = system->particles.size() / num_threads;
        //int num_bonds_per_thread = system->bond_list->count / num_threads;
        // start worked threads
        struct arg*targs = (struct arg*) malloc(sizeof(struct arg)*num_threads);
        pthread_t* thread_handles = (pthread_t*) malloc(sizeof(pthread_t)*num_threads);
        // create barrier that unblocks when "num_threads+1" call wait() on it
        pthread_barrier_init(&begin_step_barrier, NULL, num_threads+1);
        pthread_barrier_init(&end_step_barrier, NULL, num_threads+1);
        // create all the worker threads
        int num_particles_left = system->particles.size();
        long unsigned int particle_list_ittr = 0;
        if(debug_flag){ printf("Creating %i simulation threads\n",num_threads);}
        for (i=0; i < num_threads; i++) {
            targs[i].system = system;
            targs[i].num_threads = num_threads;
            targs[i].thread_id = i;
            targs[i].my_first_particle = particle_list_ittr;
            if(i==num_threads-1){
                targs[i].num_my_particles = num_particles_left;
            }else{
                targs[i].num_my_particles = num_particles_per_thread;
                num_particles_left -= num_particles_per_thread;
                for(j=0;j<num_particles_per_thread;j++){
                    particle_list_ittr++;
                    if(particle_list_ittr >= system->particles.size()){
                        printf("ERROR: particle_list_ittr was unexpectly NULL\n");
                        exit(1);
                    }
                }
            }
            if(debug_flag) printf("Creating thread to process %i particles\n", targs[i].num_my_particles);
            pthread_create(&thread_handles[i], NULL, run_simulation_thread, &targs[i]);
        }
        // Create threads to sort indexes
        if(debug_flag) printf("Creating thread to update sort index\n");
        pthread_t* sort_thread_handles = (pthread_t*) malloc(sizeof(pthread_t)*3);
        pthread_barrier_init(&begin_sort_barrier, NULL, 2);
        pthread_barrier_init(&end_sort_barrier, NULL, 2);
        struct sarg sort_args[1];
        sort_args[0].system = system ;
        sort_args[0].sort_ndx = 0;
        pthread_create(&sort_thread_handles[0], NULL, sort_index_thread, &sort_args[0]);
        /*sort_args[1].ll = system->y_index;
        sort_args[1].sort_index = 1;
        pthread_create(&sort_thread_handles[1], NULL, sort_index_thread, &sort_args[1]);
        sort_args[2].ll = system->z_index;
        sort_args[2].sort_index = 2;
        pthread_create(&sort_thread_handles[2], NULL, sort_index_thread, &sort_args[2]);*/

        // Create thread to handle output
        pthread_t* output_thread_handles = (pthread_t*) malloc(sizeof(pthread_t)*1);
        pthread_barrier_init(&begin_output_barrier, NULL, 2);
        pthread_barrier_init(&end_output_barrier, NULL, 2);
        pthread_create(&output_thread_handles[0], NULL, output_system_thread, system);
        if(debug_flag) printf("Creating thread to create output files\n");

        // Start simulation, coordinate simulation
        unsigned int next_output_step = 0;
        for(system->current_step=0; system->current_step < system->nt; system->current_step++){

            // Release the Sort Index threads
            if(debug_flag) printf("[%i] Starting the Sort Index threads\n",system->current_step);
            pthread_barrier_wait(&begin_sort_barrier);
            // Output state
            if(system->current_step >= next_output_step){
                // Release output thread
                if(debug_flag) printf("[%i] Starting the Output threads\n",system->current_step);
                pthread_barrier_wait(&begin_output_barrier);
                next_output_step += system->output_freq;
                // Wait until output threads are done
                pthread_barrier_wait(&end_output_barrier);
                if(debug_flag) printf("[%i] Output threads finished\n",system->current_step);

            }
            // Wait until Sort Index threads are done
            pthread_barrier_wait(&end_sort_barrier);
            if(debug_flag) printf("[%i] Sort Index threads finished\n",system->current_step);
            if(debug_flag>2 && system->current_step==0){
                printf("x_index = [");
                for(auto n = system->particles.begin(); n!=system->particles.end(); ++n){
                    printf("%e ",n->x[0]);
                }
                printf("]\n");
                printf("x_index_id = [");
                for(auto n = system->particles.begin(); n!=system->particles.end(); ++n){
                    printf("%i ",n->id);
                }
                printf("]\n");
            }

            // Release the worker threads to take step
            if(debug_flag) printf("[%i] Starting the Worker threads\n",system->current_step);
            pthread_barrier_wait(&begin_step_barrier);
            unsigned int nsubsteps = get_number_of_substeps();
            for(unsigned int substep=0;substep < nsubsteps; substep++){
                // Wait until worker threads are done
                pthread_barrier_wait(&end_step_barrier);
                if(debug_flag) printf("[%i] Worker threads finished substep %i/%i\n",system->current_step,substep,nsubsteps);
            }
            // Solve RDME
            if(debug_flag) printf("[%i] starting RDME simulation\n",system->current_step);
            simulate_rdme(system, system->current_step);
            if(debug_flag) printf("[%i] Finish RDME simulation\n",system->current_step);
        }
        // Record final timepoint
        if(debug_flag) printf("[%i] Starting the Output threads\n",system->current_step);
        pthread_barrier_wait(&begin_output_barrier);
        pthread_barrier_wait(&end_output_barrier);
        if(debug_flag) printf("[%i] Output threads finished\n",system->current_step);
        if(debug_flag) printf("[%i] Waiting for Async Output threads\n",system->current_step);
        pthread_barrier_wait(&begin_output_barrier); // wait for the async
        if(debug_flag) printf("[%i] Async Output threads finished\n",system->current_step);

        //clean up
        if(debug_flag) printf("Cleaning up RDME\n");
        destroy_rdme(system);

        // Kill threads and wait for them to finish
        /*
        if(debug_flag) printf("Cleaning up threads\n");
        for(i=0;i<3;i++){
            pthread_kill(sort_thread_handles[i], SIGINT);
            pthread_join(sort_thread_handles[i], NULL);
        }
        pthread_kill(output_thread_handles[0], SIGINT);
        pthread_join(output_thread_handles[0], NULL);
        for (i=0; i < num_threads; i++) {
            pthread_kill(thread_handles[i], SIGINT);
            pthread_join(thread_handles[i], NULL);
        }
        */
        // done
        if(debug_flag) printf("Simulation complete\n");;
        return;
    }
}
