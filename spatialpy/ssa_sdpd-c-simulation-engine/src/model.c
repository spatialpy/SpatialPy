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

Based on a Matlab program by Bruno Jacob (UCSB)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "particle_system.hpp"

namespace Spatialpy{

    void pairwiseForce(Particle* me, ParticleSystem* system)
    {
        // F, Frho and Fbp are output
        std::vector<NeighborNode> neighbors = me->neighbors;
        //printf("pairwiseForce(id=%i)\n",me->id);
        //fflush(stdout);

        double R, Pj, dv[3], dx[3], r, dWdr, dv_dx, fp, fv, fbp,
            transportTensor[3][3], ft[3], pressure_gradient;
        dx[0] = 0.0;
        dx[1] = 0.0;
        dx[2] = 0.0;
        double h = system->h;
        double rho0 = system->rho0;
        double P0 = system->P0;
        double c0 = system->c0;
        double Pi = P0 * (me->rho / rho0 - 1.0);
        int i, j ;
        long unsigned int s, rxn;
        Particle* pt_j;

        // Kernel function parameter
        //if (system->dimension == 3) {
        //    alpha = 105 / (16 * M_PI * h * h * h);
        //}
        //else if (system->dimension == 2) {
        //    alpha = 5 / (M_PI * h * h);
        //}
        //else {
        //    printf("Error, only 3D or 2D domains are supported\n");
        //    exit(1);
        //}

        // Zero tensors
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                transportTensor[i][j] = 0.0;
            }
        }

        //printf("pairwiseForce(id=%i) before compute force\n",me->id);
        //fflush(stdout);

        // Compute force from each neighbor
        for (auto n = neighbors.begin(); n != neighbors.end(); n = ++n) {
            pt_j = n->data;

            //r = particle_dist(me, pt_j);
            r = n->dist;
            R = r / h;
            //printf("pairwiseForce(id=%i) pt_j->id=%i r=%e R=%e\n",me->id,pt_j->id,r,R);
            //fflush(stdout);
            if (R > 1.0)
                continue; // outside kernel support
            if (r == 0.0)
                continue; // ignore sigularities
            // Compute weight function and weight function derivative
            // dWdr = (5/(pi*(h^2))) * (-12*r/(h^2)) * (1 - r/h)^2;
            //dWdr = alpha * (-12 * r / (h * h)) * pow(1 - R, 2);
            dWdr = n->dWdr;
            // Spatial deriviatives
            dv_dx = 0.0;
            for (i = 0; i < system->dimension; i++) {
                dx[i] = (me->x[i] - pt_j->x[i]);
                dv[i] = (me->v[i] - pt_j->v[i]);
                dv_dx += dv[i] * dx[i];
            }

            // Pressure of particle j
            Pj = P0 * (pt_j->rho / rho0 - 1.0);

            // Check if pressure gradient sign
            pressure_gradient = Pi / (me->rho * me->rho) + Pj / (pt_j->rho * pt_j->rho);
            if (pressure_gradient < 0) pressure_gradient = -Pi / (me->rho * me->rho) + Pj / (pt_j->rho * pt_j->rho);

            // Compute pressure force
            fp = -1.0 * pt_j->mass * pressure_gradient * dWdr / (r + 0.001 * h);

            // Compute viscous force
            fv = pt_j->mass * (2.0 * (me->nu * pt_j->nu) / (me->nu + pt_j->nu)) * (1 / (r + 0.001 * h)) * dWdr / ((me->rho * pt_j->rho));

            // Compute background pressure (bp) force
            fbp = -10.0 * P0 * (1.0 / me->mass) * (pow(me->mass / me->rho, 2) + pow(pt_j->mass / pt_j->rho, 2)) * dWdr / (r + 0.001 * h);

            // Compute transport force and transport tensor
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    transportTensor[i][j] = 0.5 * ((me->rho * me->v[i] * (me->vt[j] - me->v[j])) + (pt_j->rho * pt_j->v[i] * (pt_j->vt[j] - pt_j->v[j])));
                }
            }
            for (i = 0; i < 3; i++) {
                ft[i] = (1.0 / me->mass) * (pow(me->mass / me->rho, 2) + pow(pt_j->mass / pt_j->rho, 2)) *
                                           (transportTensor[i][0] * dx[0] + transportTensor[i][1] * dx[1] + transportTensor[i][2] * dx[2]) * dWdr / (r + 0.001 * h);
            }

            // Sum forces
            for (i = 0; i < system->dimension; i++) {
                me->F[i] += fp * dx[i] + fv * dv[i] + ft[i];
                me->Fbp[i] += fbp * dx[i];
            }
            //printf("pairwiseForce(id=%i) me->F = [%e, %e, %e]\n",me->id,me->F[0],me->F[1],me->F[2]);
            //fflush(stdout);

            // Compute density variation
            me->Frho = me->Frho + me->rho * (pt_j->mass / pt_j->rho) * dv_dx * (1 / (r + 0.001 * h)) * dWdr
                          - 0.0 * h * me->rho * c0 * pt_j->mass * 2.0 * (pt_j->rho / me->rho - 1.0) * (r * r / (r * r + 0.01 * h * h)) * (1.0 / (r + 0.001 * h)) * dWdr / pt_j->rho
                          - (pt_j->mass / pt_j->rho) * (me->rho * ((me->v[0] - me->vt[0]) * dx[0] + (me->v[1] - me->vt[1]) * dx[1] + (me->v[2] - me->vt[2]) * dx[2])
                          + pt_j->rho * ((pt_j->v[0] - pt_j->vt[0]) * dx[0] + (pt_j->v[1] - pt_j->vt[1]) * dx[1] + (pt_j->v[2] - pt_j->vt[2]) * dx[2])) * (1.0 / (r + 0.001 * h)) * dWdr;

            //printf("pairwiseForce(id=%i) me->Frho = [%e]\n",me->id,me->Frho);
            //fflush(stdout);

            // Compute Chem Rxn Flux (diffusion part)
            double wfd = (1.0 / (r + 0.001 * h)) * dWdr;
            //printf("pairwiseForce(id=%i) wfd = %e\n",me->id,wfd);
            //fflush(stdout);
            double dQc_base = 2.0* ((me->mass*pt_j->mass)/(me->mass+pt_j->mass)) * ((me->rho+pt_j->rho)/(me->rho*pt_j->rho)) * (r*r) * wfd / ((r*r) + 0.01*h*h); // (Tartakovsky et. al., 2007, JCP)
            //printf("pairwiseForce(id=%i) dQc_base = %e\n",me->id,dQc_base);
            //fflush(stdout);

            //printf("pairwiseForce(id=%i) system->num_chem_species = %i\n",me->id,system->num_chem_species);
            //fflush(stdout);
            for(s=0; s < system->num_chem_species; s++){
                // Note about below:  types start at  1
                int k = (system->num_chem_species) * (me->type - 1) + s;
                //printf("pairwiseForce(id=%i) s=%i, k=%i num_types=%i me->type=%i\n",me->id,s,k,system->num_types,me->type);
                //fflush(stdout);
                double dQc = system->subdomain_diffusion_matrix[k] * (me->C[s] - pt_j->C[s]) * dQc_base;
                //printf("pairwiseForce(id=%i) dQc = %e (me->C[s]=%e - pt_j->C[s]=%e) diffusion=%e \n",me->id,dQc,me->C[s],pt_j->C[s],system->subdomain_diffusion_matrix[k]);
                //fflush(stdout);
                me->Q[s] += dQc;
            }
            //printf("pairwiseForce(id=%i) me->Q = [%e]\n",me->id,me->Q[0]);
            //fflush(stdout);
        }
        //printf("pairwiseForce(id=%i) num_chem_rxns=%i\n",me->id,system->num_chem_rxns);
        //fflush(stdout);

        // after processing all neighbors
        // process chemical reactions

        double vol = (me->mass / me->rho);
        double cur_time = system->current_step * system->dt;
        for(rxn=0; rxn < system->num_chem_rxns; rxn++){
            //TODO: t (2nd arg) set to zero, fix
            //TODO, vol (3rd arg) set to zero, fix
            //TODO: data (4th arg) set to NULL, fix
            double flux = (*system->chem_rxn_rhs_functions[rxn])(me->C, cur_time, vol , me->data_fn, me->type);
            for(s=0; s< system->num_chem_species; s++){
                int k = system->num_chem_rxns * rxn + s;
                me->Q[s] += system->stoichiometric_matrix[k] * flux;
            }
        }

    }


    void filterDensity(Particle* me, ParticleSystem* system)
    {
      std::vector<NeighborNode> neighbors = me->neighbors;
      Particle* pt_j;
      double r, R, Wij, alpha, num, den;
      double h = system->h;
      if(system->dimension==3){
          alpha = 105/(16*M_PI*h*h*h);
      }else if(system->dimension==2){
          alpha = 5/(M_PI*h*h);
      }else{
          printf("Not 1D or 2D\n");exit(1);
      }
      num = 0.0;
      den = 0.0;

      for(auto n=neighbors.begin(); n!=neighbors.end(); n = ++n){
          pt_j = n->data;
          r = n->dist;
          R = r/h;
          if(R > 1.0) continue;
          if(r == 0.0) continue;

          // Compute weight function
          //Wij = (5/(M_PI*pow(h,2))) * (1+3*r/h) * pow((1-r/h),3) ;
          Wij = alpha*( (1+3*R) * pow( 1-R,3));

          //Compute numerator of Shepard filter
          num += pt_j->rho * Wij;

          // Compute denominator of Shepard filter
          den += Wij;
      }

      // Update density
      me->rho = num / den;

    }


    void computeBoundaryVolumeFraction(Particle* me, ParticleSystem* system)
    {
        std::vector<NeighborNode> neighbors = me->neighbors;
        Particle* pt_j;
        double r, R, Wij, dWdr, alpha, vos, vtot, nw[3], dx[3], norm_nw;
        int i;
        double h = system->h;
        if (system->dimension == 3) {
            alpha = 105 / (16 * M_PI * h * h * h);
        }
        else if (system->dimension == 2) {
            alpha = 5 / (M_PI * h * h);
        }
        else {
            printf("Not 1D or 2D\n");
            exit(1);
        }

        for (i = 0; i < 3; i++) {
            nw[i] = 0.0;
            dx[i] = 0.0;
        }
        for (i = 0; i < system->dimension; i++) {
            me->normal[i] = 0.0;
            me->vt[i] = 0.0;
        }

        me->bvf_phi = 0.0;
        vos = 0.0;
        vtot = 0.0;
        for (auto n = neighbors.begin(); n != neighbors.end(); ++n) {
            pt_j = n->data;
            r = n->dist;
            R = r / h;
            if (R > 1.0)
                continue;
            if (r == 0.0)
                continue;
            // Compute weight function and weight function derivative
            // Wij = (5/(M_PI*pow(h,2))) * (1+3*r/h) * pow((1-r/h),3) ;
            Wij = alpha * ((1 + 3 * R) * pow(1 - R, 3));
            //dWdr = alpha * (-12 * r / (h * h)) * pow(1 - R, 2);
            dWdr = n->dWdr;

            for (i = 0; i < system->dimension; i++) {
                dx[i] = (me->x[i] - pt_j->x[i]);
            }

            // Volume of solid (vos) around particle i
            if (pt_j->solidTag)
                vos += pow(pt_j->mass / pt_j->rho, 2) * Wij;

            // Total volume (vtot) around particle i
            vtot += pow(pt_j->mass / pt_j->rho, 2) * Wij;

            // Numerator of normal vectors (pointing outwards the nearby solid wall)
            if (pt_j->solidTag) {
                for (i = 0; i < 3; i++)
                    nw[i] += pow(pt_j->mass / pt_j->rho, 2) * dx[i] * dWdr / (r + 0.001 * h);
            }
        }

        // Apply denominator and normalization to construct normal vectors
        for (i = 0; i < 3; i++)
            nw[i] = nw[i] / vtot;

        // Compute norm of nw
        norm_nw = sqrt(nw[0] * nw[0] + nw[1] * nw[1] + nw[2] * nw[2]);

        // Update normals to normalized normals
        for (i = 0; i < 3; i++)
            me->normal[i] = -nw[i] / norm_nw;

        // Compute bvf_phi (boundary volume fraction) for particle i
        if (me->solidTag)
            me->bvf_phi = 0.0;
        else
            me->bvf_phi = fabs(vos / vtot);
    }


    void applyBoundaryVolumeFraction(Particle* me, ParticleSystem* system)
    {
        int i;
        double v_dot_normal;

        // Bounce-back condition for fluid particles
        v_dot_normal = me->v[0] * me->normal[0] + me->v[1] * me->normal[1] + me->v[2] * me->normal[2];

        if (me->bvf_phi >= 0.5) {
            for (i = 0; i < system->dimension; i++) {
                me->v[i] = -me->v[i] + 2.0 * fmax(0.0, v_dot_normal) * me->normal[i];
            }
        }
    }

}
