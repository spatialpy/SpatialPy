/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#include "linked_list.h"
#include "particle.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int add_to_neighbor_list(particle_t*me, particle_t*neighbor, system_t*system){
    double a = me->x[0] - neighbor->x[0];
    double b = me->x[1] - neighbor->x[1];
    double c = me->x[2] - neighbor->x[2];
    double r2 =  ( a*a + b*b + c*c);
    double r = sqrt(r2);

    if(r > system->h){ return 0; } // do not add, out side support radius

    // calculate dWdr
    double h = system->h;
    double R = r / h;
    double alpha = 105 / (16 * M_PI * h * h * h); // 3D
    //double alpha = 5 / (M_PI * h * h); // 2D
    //double dWdr = alpha * (-12 * r / (h * h)) * pow(1 - R, 2);
    double dWdr = alpha * (-12 * r / (h * h)) * ((1 - R)* (1 - R));
    // calculate D_i_j

    // Eq (13-14), Drawert et al 2019
    double ih = 1.0 / h;
    double ihsq = ih * ih;
    double dhr = h - r;
    double wfd = -25.066903536973515383e0 * dhr * dhr * ihsq * ihsq * ihsq * ih; //3D
    // Eq 28 of Drawert et al 2019, Tartakovsky et. al., 2007, JCP
    double D_i_j = -2.0*(me->mass*neighbor->mass)/(me->mass+neighbor->mass)*(me->rho+neighbor->rho)/(me->rho*neighbor->rho) * r2 * wfd / (r2+0.01*h*h);
    if(isnan(D_i_j)){
        printf("Got NaN calculating D_i_j for me=%i, neighbor=%i\n",me->id, neighbor->id);
        printf("r=%e ",r);
        printf("h=%e ",h);
        printf("alpha=%e ",alpha);
        printf("dWdr=%e ",dWdr);
        particle_t*p = me;
        printf("me->mass=%e ",p->mass);
        printf("me->rho=%e ",p->rho);
        p = neighbor;
        printf("n->mass=%e ",p->mass);
        printf("n->rho=%e ",p->rho);

        exit(1);
    }

    neighbor_node_t*n = neighbor_list_add( me->neighbors, neighbor );
    n->dist = r;
    n->dWdr = dWdr;
    n->D_i_j = D_i_j;

    return 1;
}


void find_neighbors(particle_t* me, system_t* system){
    node_t*n;
    //clean out previous neighbors
    //printf("find_neighbors.empty_linked_list\n");
    empty_neighbor_list(me->neighbors);
    // search for points forward: (assumes the list is sorted ascending)
    //printf("find_neighbors.search forward\n");
    for(n = me->x_index->next; n!=NULL; n=n->next){
        if(n->data->x[0] > (me->x[0] + system->h)) break; //stop searching
        if( (n->data->x[1] > (me->x[1] + system->h)) || (n->data->x[1] < (me->x[1] - system->h) ) ) continue;
        if( (n->data->x[2] > (me->x[2] + system->h)) || (n->data->x[2] < (me->x[2] - system->h) ) ) continue;
        add_to_neighbor_list(me, n->data, system);
        if(debug_flag>2){ 
            printf("find_neighbors(%i) forward found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
            me->id,n->data->id, particle_dist(me,n->data),
            me->x[0] - n->data->x[0],
            me->x[1] - n->data->x[1],
            me->x[2] - n->data->x[2]);
        }
    }
    // check for wrap around points x
    /*if(n==NULL && system->boundary_conditions[0] == 'p'){
        double max_x = system->h - (me->x[0] - system->xhi);
        for(n = system->x_index->head; n!=NULL; n=n->next){
            if(n->data->x[0] > max_x) break; //stop searching
            if( (n->data->x[1] > (me->x[1] + system->h)) || (n->data->x[1] < (me->x[1] - system->h) ) ) continue;
            if( (n->data->x[2] > (me->x[2] + system->h)) || (n->data->x[2] < (me->x[2] - system->h) ) ) continue;
            add_to_neighbor_list(me, n->data, system);
            if(debug_flag){ 
                printf("find_neighbors(%i) forward wrap around found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
                me->id,n->data->id, particle_dist(me,n->data),
                me->x[0] - n->data->x[0],
                me->x[1] - n->data->x[1],
                me->x[2] - n->data->x[2]);
            }
        }
    }
    // check for wrap around points y
    if(system->boundary_conditions[1] == 'p' && (me->x[1] + system->h) > system->yhi){
        //double max_y = system->ylo + ((me->x[1] + system->h) - system->yhi);
        fprintf(stderr,"periodic boundary in y not supported yet\n");
        exit(1);
        //TODO
    }
    if(system->boundary_conditions[1] == 'p' && (me->x[1] - system->h) < system->ylo){
        fprintf(stderr,"periodic boundary in y not supported yet\n");
        exit(1);
        //TODO
    }
    // check for wrap around points z
    if(system->boundary_conditions[2] == 'p' && (me->x[2] + system->h) > system->zhi){
        fprintf(stderr,"periodic boundary in z not supported yet\n");
        exit(1);
        //TODO
    }
    if(system->boundary_conditions[2] == 'p' && (me->x[2] - system->h) < system->zlo){
        fprintf(stderr,"periodic boundary in z not supported yet\n");
        exit(1);
        //TODO
    }
    */

    //search for points backward
    for(n = me->x_index->prev; n!=NULL; n=n->prev){
        //if(n->data->x[0] > (me->x[0] + system->h)) break; //stop searching forward
        /*if(me->id == 0){ 
            printf("find_neighbors(%i) backwards looking at  %i dist: %e    dx: %e   dy: %e   dz: %e\n",
            me->id,n->data->id, particle_dist(me,n->data),
            me->x[0] - n->data->x[0],
            me->x[1] - n->data->x[1],
            me->x[2] - n->data->x[2]);
        }*/
        if(n->data->x[0] < (me->x[0] - system->h)){
            //if(me->id==0){ printf("\tnode x is less than x[me]-h\n");}
            break; //stop searching backward
        }
        if( (n->data->x[1] > (me->x[1] + system->h)) || (n->data->x[1] < (me->x[1] - system->h) ) ){
            //if(me->id==0){ printf("\tnode y is outside\n");}
            continue;
        }
        if( (n->data->x[2] > (me->x[2] + system->h)) || (n->data->x[2] < (me->x[2] - system->h) ) ){ 
            //if(me->id==0){ printf("\tnode z is outside\n");}
            continue;
        }
        add_to_neighbor_list(me, n->data, system);
        if(debug_flag>2){ 
            printf("find_neighbors(%i) backwards found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
            me->id,n->data->id, particle_dist(me,n->data),
            me->x[0] - n->data->x[0],
            me->x[1] - n->data->x[1],
            me->x[2] - n->data->x[2]);
        }
    }
    /*if(n==NULL && system->boundary_conditions[0] == 'p'){
        double min_x = system->xhi - (system->h - (system->xlo - me->x[0]));
        for(n = system->x_index->tail; n!=NULL; n=n->prev){
            if(n->data->x[0] < min_x) break; //stop searching
            if( (n->data->x[1] > (me->x[1] + system->h)) || (n->data->x[1] < (me->x[1] - system->h) ) ) continue;
            if( (n->data->x[2] > (me->x[2] + system->h)) || (n->data->x[2] < (me->x[2] - system->h) ) ) continue;
            add_to_neighbor_list(me, n->data, system);
            if(debug_flag){ 
                printf("find_neighbors(%i) backwards wrap around found %i dist: %e    dx: %e   dy: %e   dz: %e\n",
                me->id,n->data->id, particle_dist(me,n->data),
                me->x[0] - n->data->x[0],
                me->x[1] - n->data->x[1],
                me->x[2] - n->data->x[2]);
            }
        }
    }*/

    // Sort the list:
    //neighbor_list_sort(me->neighbors);
}


system_t* create_system(size_t num_types, size_t num_chem_species, size_t num_chem_rxns, 
                         size_t num_stoch_species, size_t num_stoch_rxns,size_t num_data_fn){
    system_t*s = malloc(sizeof(system_t));
    s->particle_list = create_linked_list();
    s->x_index = create_linked_list();
    //s->y_index = create_linked_list();
    //s->z_index = create_linked_list();
    s->dimension = 3;
    s->boundary_conditions[0] = 'n';
    s->boundary_conditions[1] = 'n';
    s->boundary_conditions[2] = 'n';
    s->rdme = NULL;
    s->static_domain = 0;
    s->num_stoch_species = num_stoch_species;
    s->num_stoch_rxns = num_stoch_rxns;
    s->num_chem_species = num_chem_species;
    s->num_chem_rxns = num_chem_rxns;
    s->num_types = num_types;
    s->gravity = calloc(3,sizeof(double));
    s->num_data_fn = num_data_fn;
    return s;
}

particle_t* create_particle(int id){
    particle_t* me = malloc(sizeof(particle_t));
    me->id = id;
    me->nu = 0.01; 
    me->mass = 1;
    me->rho = 1;
    me->solidTag = 0;
    me->x[0] = me->x[1] = me->x[2] = 0.0;
    me->v[0] = me->v[1] = me->v[2] = 0.0;
    return me;
}

void add_particle(particle_t* me, system_t* system){
    linked_list_add(system->particle_list, me);
    me->x_index = linked_list_add(system->x_index, me);
    //me->y_index = linked_list_add(system->y_index, me);
    //me->z_index = linked_list_add(system->z_index, me);
    me->neighbors = create_neighbor_list();
    me->Q = (double*) calloc(system->num_chem_species, sizeof(double));
    me->C = (double*) calloc(system->num_chem_species, sizeof(double));
    me->data_fn = (double*) calloc(system->num_data_fn, sizeof(double));
}

double particle_dist(particle_t* p1, particle_t*p2){
    double a = p1->x[0] - p2->x[0];
    double b = p1->x[1] - p2->x[1];
    double c = p1->x[2] - p2->x[2];
    return sqrt( a*a + b*b + c*c);
}
double particle_dist_sqrd(particle_t* p1, particle_t*p2){
    double a = p1->x[0] - p2->x[0];
    double b = p1->x[1] - p2->x[1];
    double c = p1->x[2] - p2->x[2];
    return ( a*a + b*b + c*c);
}

//bond_t* create_bond(particle_t*p1, particle_t*p2, double k, double rest_distance){
//    bond_t* me = malloc(sizeof(bond_t));
//    me->p1 = p1;
//    me->p2 = p2;
//    me->param_k = k;
//    me->rest_distance = rest_distance;
//    return me;
//}




