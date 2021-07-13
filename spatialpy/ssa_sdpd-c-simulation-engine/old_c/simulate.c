/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#include "linked_list.h"
#include "simulate.h"
#include "output.h"
#include "particle.h"
#include "simulate_rdme.h"
#include "model.h"
#include <errno.h>
#include "pthread_barrier.h"
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void take_step(void*me, system_t*system, unsigned int step, unsigned int substep){
    //particle_t*me2 = (particle_t*)me;
    //printf("take_step(me->id=%i step=%i substep=%i)\n",me2->id, step, substep);
    //fflush(stdout);
    if(substep==0){
        take_step1((particle_t*)me,system,step);
    }else if(substep==1){
        compute_forces((particle_t*)me,system,step);
    }else if(substep==2){
        take_step2((particle_t*)me,system,step);
    }else{
        printf("ERROR, substep=%u\n",substep);
        exit(1);
    }
}

unsigned int get_number_of_substeps(){
    return 3;
}


// Step 1/3: First part of time step computation
void take_step1(particle_t* me, system_t* system, unsigned int step) {
    int i;
    //printf("particle id=%i Q[0]=%e\n",me->id,me->Q[0]);

    // Step 1.1: 
    if(step==0 || system->static_domain == 0){
        find_neighbors(me, system);
    }

    // Step 1.2: Predictor step
    
    // Update half-state
    if (me->solidTag == 0 && system->static_domain == 0) {
       for (i = 0; i < 3; i++) {
           // Update velocity using forces
            me->v[i] = me->v[i] + 0.5 * system->dt * me->F[i];
            // Update transport velocity using background pressure force
            me->vt[i] = me->v[i] + 0.5 * system->dt * me->Fbp[i];
            // Update position using previous velocity
            me->x[i] = me->x[i] + system->dt * me->vt[i];
        }
        // Update density using continuity equation 
        me->rho = me->rho + 0.5 * system->dt * me->Frho;
    }
    // update half-state of chem rxn
    if(step > 0){
        for(i=0; i< system->num_chem_species; i++){
            me->C[i] += me->Q[i] * system->dt * 0.5;
        }
    }

    // Apply boundary conditions
    applyBoundaryConditions(me, system);


    //  Clean forces
    for (i = 0; i < 3; i++) {
        // Clean momentum force
        me->F[i] = system->gravity[i];
        // Clean background pressure force
        me->Fbp[i] = 0.0;
    }
    // Clean mass flux term
    me->Frho = 0.0;
    // Clean chem rxn flux
    for(i=0; i< system->num_chem_species; i++){
        me->Q[i] = 0.0;
    }


}

// Step 2/3: Update the 'F' field of each partile attached to a bond
// void compute_bond_forces(bond_t*bond, system_t*system, unsigned int step){
//}


// Step 2/3: Compute forces
void compute_forces(particle_t* me, system_t* system, unsigned int step) {

    //printf("compute_forces() particle id=%i Q[0]=%e\n",me->id,me->Q[0]);
    // Step 2.2: Find nearest neighbors
    if(step>0 && system->static_domain == 0){
        find_neighbors(me, system);
    }

    // Step 2.3: Compute forces
    pairwiseForce(me, system);

}


// Step 3/3: Compute the final state
void take_step2(particle_t* me, system_t* system, unsigned int step)
{
    int i;

    // Step 3.1: Corrector step
    if (me->solidTag == 0 && system->static_domain == 0) {
        // Update velocity using forces
        for (i = 0; i < 3; i++) {
            me->v[i] = me->v[i] + 0.5 * system->dt * me->F[i];
        }

        // Update density using continuity equation and Shepard filter 
        if (step % 20 == 0) {
            filterDensity(me, system);
        }
        me->rho = me->rho + 0.5 * system->dt * me->Frho;

      // Solid (wall) particles should change density
    }else if (me->solidTag == 1 && system->static_domain == 0) {
        // Filter density field (for fixed solid particles)
        if (step % 20 == 0) {
            filterDensity(me, system);
        }
    }


    // Step 3.2: Compute boundary volume fractions (bvf)
    // Step 3.3: Apply BVF 
    if (me->solidTag == 0 && system->static_domain == 0) {
        computeBoundaryVolumeFraction(me,system);
        applyBoundaryVolumeFraction(me, system);
    }

    //  Solve deterministic/stochastic reaction-diffusion system
    // update half-state of chem rxn
    for(i=0; i< system->num_chem_species; i++){
        me->C[i] += me->Q[i] * system->dt * 0.5;
    }
    // Apply boundary conditions
    applyBoundaryConditions(me, system);

}
