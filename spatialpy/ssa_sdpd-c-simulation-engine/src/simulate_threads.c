#include "debug.h"
#include "linked_list.h"
#include "model.h"
#include "output.h"
#include "particle.h"
#include "simulate.h"
#include "simulate_rdme.h"
#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>



struct arg {
    system_t* system;
    unsigned int thread_id;
    unsigned int num_threads;
    unsigned int num_my_particles;
    //unsigned int num_my_bonds;
    node_t*my_first_particle;
    //bond*my_first_bond;
};

pthread_barrier_t begin_step_barrier;
pthread_barrier_t end_step_barrier;
pthread_barrier_t begin_sort_barrier;
pthread_barrier_t end_sort_barrier;
pthread_barrier_t begin_output_barrier;
pthread_barrier_t end_output_barrier;
unsigned int current_step;

void* output_system_thread(void* targ){
    system_t* system = (system_t*) targ;
    while(1){
        pthread_barrier_wait(&begin_output_barrier);
        //output_csv(system, current_step);
        SSA_LOG(1, "[OUT] start output_vtk__sync_step()\n");
        output_vtk__sync_step(system, current_step);
        SSA_LOG(1, "[OUT] done output_vtk__sync_step()\n");
        pthread_barrier_wait(&end_output_barrier);
        SSA_LOG(1, "[OUT] start output_vtk__async_step()\n");
        output_vtk__async_step(system);
        SSA_LOG(1, "[OUT] done output_vtk__async_step()\n");
    }
}

struct sarg {
    linked_list_t*ll;
    int sort_ndx;
};

void* sort_index_thread(void* targ_in){
    struct sarg* targ = (struct sarg*) targ_in;
    while(1){
        pthread_barrier_wait(&begin_sort_barrier);
        SSA_LOG(1, "[SORT] begin sort\n");
        linked_list_sort(targ->ll,targ->sort_ndx);
        SSA_LOG(1, "[SORT] sort complete\n");
        pthread_barrier_wait(&end_sort_barrier);
    }
}

void* run_simulation_thread(void *targ_in){
    struct arg* targ = (struct arg*)targ_in;
    system_t* system = targ->system;
    unsigned int step;
    node_t*n;
    int i;
    int count = 0;
    // each thread will take a step with each of it's particles
    for(step=0; step < system->nt; step++){
        // block on the begin barrier
        SSA_LOG(1, "[WORKER %i] waiting to begin step %i\n",targ->thread_id,step);
        pthread_barrier_wait(&begin_step_barrier);
        //---------------------------------------
        unsigned int nsubsteps = get_number_of_substeps();
        for(int substep=0;substep < nsubsteps; substep++){
            // take_step
            count = 0;
            n=targ->my_first_particle;
            for(i=0; i<targ->num_my_particles; i++){
                if(n==NULL) break;
                take_step(n->data,system,step,substep);
                count++;
                n=n->next;
            }
            // block on the end barrier
            SSA_LOG(1, "[WORKER %i] completed step %i, substep %i/%i, processed %i particles\n",targ->thread_id,step,substep,nsubsteps,count);
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
        SSA_LOG(1, "[WORKER %i] completed compute_bond_forces %i, processed %i particles\n",targ->thread_id,step,count);
        pthread_barrier_wait(&end_step_barrier);
        */
        //---------------------------------------
    } //end for(step)
    return NULL;
}

void run_simulation(int num_threads, system_t* system){
    //
    int i,j;
    int num_particles_per_thread = system->particle_list->count / num_threads;
    //int num_bonds_per_thread = system->bond_list->count / num_threads;
    // start worked threads
    struct arg*targs = (struct arg*) malloc(sizeof(struct arg)*num_threads);
    pthread_t* thread_handles = (pthread_t*) malloc(sizeof(pthread_t)*num_threads);
    // create barrier that unblocks when "num_threads+1" call wait() on it
    pthread_barrier_init(&begin_step_barrier, NULL, num_threads+1);
    pthread_barrier_init(&end_step_barrier, NULL, num_threads+1);
    // create all the worker threads
    int num_particles_left = system->particle_list->count;
    node_t*particle_list_ittr = system->particle_list->head;
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
                particle_list_ittr = particle_list_ittr->next;
                if(particle_list_ittr == NULL){
                    printf("ERROR: particle_list_ittr was unexpectly NULL\n");
                    exit(1);
                }
            }
        }
        SSA_LOG(1, "Creating thread to process %i particles\n", targs[i].num_my_particles);
        pthread_create(&thread_handles[i], NULL, run_simulation_thread, &targs[i]);
    }
    // Create threads to sort indexes
    pthread_t* sort_thread_handles = (pthread_t*) malloc(sizeof(pthread_t)*3);
    pthread_barrier_init(&begin_sort_barrier, NULL, 2);
    pthread_barrier_init(&end_sort_barrier, NULL, 2);
    struct sarg sort_args[1];
    sort_args[0].ll = system->x_index;
    sort_args[0].sort_ndx = 0;
    pthread_create(&sort_thread_handles[0], NULL, sort_index_thread, &sort_args[0]);
    SSA_LOG(1, "Creating thread to update x-position index\n");
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
    SSA_LOG(1, "Creating thread to create output files\n");

    // Start simulation, coordinate simulation
    int step,next_output_step = 0;
    for(step=0; step < system->nt; step++){
        // Release the Sort Index threads
        SSA_LOG(1, "[%i] Starting the Sort Index threads\n",step);
        pthread_barrier_wait(&begin_sort_barrier);
        // Output state
        if(step >= next_output_step){
            // Release output thread
            current_step = step;
            SSA_LOG(1, "[%i] Starting the Output threads\n",step);
            pthread_barrier_wait(&begin_output_barrier);
            next_output_step += system->output_freq;
            // Wait until output threads are done
            pthread_barrier_wait(&end_output_barrier);
            SSA_LOG(1, "[%i] Output threads finished\n",step);

        }
        // Wait until Sort Index threads are done
        pthread_barrier_wait(&end_sort_barrier);
        SSA_LOG(1, "[%i] Sort Index threads finished\n",step);
        if(debug_flag>2 && step==0){
            printf("x_index = [");
            node_t*n;
            for(n = system->x_index->head; n!=NULL; n=n->next){
                printf("%e ",n->data->x[0]);
            }
            printf("]\n");
            printf("x_index_id = [");
            for(n = system->x_index->head; n!=NULL; n=n->next){
                printf("%i ",n->data->id);
            }
            printf("]\n");
        }

        // Release the worker threads to take step
        SSA_LOG(1, "[%i] Starting the Worker threads\n",step);
        pthread_barrier_wait(&begin_step_barrier);
        unsigned int nsubsteps = get_number_of_substeps();
        for(int substep=0;substep < nsubsteps; substep++){
            // Wait until worker threads are done
            pthread_barrier_wait(&end_step_barrier);
            SSA_LOG(1, "[%i] Worker threads finished substep %i/%i\n",step,substep,nsubsteps);
        }
        // Solve RDME
        SSA_LOG(1, "[%i] starting RDME simulation\n",step);
        simulate_rdme(system, step);
        SSA_LOG(1, "[%i] Finish RDME simulation\n",step);
    }
    // Record final timepoint
    current_step = step;
    SSA_LOG(1, "[%i] Starting the Output threads\n",step);
    pthread_barrier_wait(&begin_output_barrier);
    pthread_barrier_wait(&end_output_barrier);
    SSA_LOG(1, "[%i] Output threads finished\n",step);
    SSA_LOG(1, "[%i] Waiting for Async Output threads\n",step);
    pthread_barrier_wait(&begin_output_barrier); // wait for the async
    SSA_LOG(1, "[%i] Async Output threads finished\n",step);

    //clean up
    SSA_LOG(1, "Cleaning up RDME\n");
    destroy_rdme(system);

    // Kill threads and wait for them to finish
    /*
    SSA_LOG(1, "Cleaning up threads\n");
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
    SSA_LOG(1, "Simulation complete\n");;
    return;
}



