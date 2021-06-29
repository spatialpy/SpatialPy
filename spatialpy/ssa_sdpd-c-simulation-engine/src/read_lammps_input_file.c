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
#include <stdio.h>
#include <stdlib.h>

#include "read_lammps_input_file.h"


void read_lammps_input_file(const char*filename, ParticleSystem*system){
    char* s;
    FILE* fp = fopen(filename, "r");
    if(fp == NULL){
        perror("Error opening input file");
        exit(EXIT_FAILURE);
    }
    char buffer[1024];
    //discard first two lines
    s=fgets(buffer,1024,fp);
    if(s==NULL){printf("Error reading input filen\n");exit(1);}
    s=fgets(buffer,1024,fp);
    s=fgets(buffer,1024,fp);
    int natoms;
    int r = sscanf(buffer,"%d atoms",&natoms);
    printf("r=%i natoms=%i\n",r,natoms);
    s=fgets(buffer,1024,fp);
    s=fgets(buffer,1024,fp);
    int ntypes;
    r = sscanf(buffer,"%d atom types",&ntypes);
    printf("r=%i ntypes=%i\n",r,ntypes);
    s=fgets(buffer,1024,fp);
    s=fgets(buffer,1024,fp);
    double xlow,xhigh,ylow,yhigh,zlow,zhigh;
    r = sscanf(buffer,"%lf %lf xlo xhi",&xlow,&xhigh);
    printf("r=%i xlow=%f xhi=%f\n",r,xlow,xhigh);
    s=fgets(buffer,1024,fp);
    r = sscanf(buffer,"%lf %lf ylo yhi",&ylow,&yhigh);
    printf("r=%i ylow=%f yhi=%f\n",r,ylow,yhigh);
    s=fgets(buffer,1024,fp);
    r = sscanf(buffer,"%lf %lf zlo zhi",&zlow,&zhigh);
    printf("r=%i zlow=%f zhi=%f\n",r,zlow,zhigh);

    system->xlo = xlow;
    system->xhi = xhigh;
    system->ylo = ylow;
    system->yhi = yhigh;
    system->zlo = zlow;
    system->zhi = zhigh;

    s=fgets(buffer,1024,fp);
    s=fgets(buffer,1024,fp);
    //printf("Is this 'Masses': '%s'\n",buffer);
    s=fgets(buffer,1024,fp);
    double *masses_of_types = malloc(sizeof(double)*ntypes);
    int i,j;
    double tmp_d;
    for(i=0;i<ntypes;i++){
        s=fgets(buffer,1024,fp);
        r = sscanf(buffer,"%d %lf",&j, &tmp_d);
        if(r!=2){ printf("error format mismatch on line : %s",buffer);exit(0);}
        masses_of_types[j-1] = tmp_d;
        printf("atom type %i has mass %lf\n",j,tmp_d);
    }
    s=fgets(buffer,1024,fp);
    s=fgets(buffer,1024,fp);
    printf("Is this 'Atoms': '%s'\n",buffer);
    s=fgets(buffer,1024,fp);
    int id,type,is_solid;
    double rho,x,y,z;
    Particle*me;
    for(i=0;i<natoms;i++){
        s=fgets(buffer,1024,fp);
        r = sscanf(buffer, "%d 0 %d %lf %lf %lf %lf %d",&id,&type,&rho,&x,&y,&z,&is_solid);
        if(r!=7){printf("error format mismatch on atom %i: %s",i,buffer);exit(0);}
        me = create_particle(id);
        me->x[0] = x; me->x[1] = y; me->x[2] = z;
        me->v[0] = 0.0; me->v[1] = 0.0; me->v[2] = 0.0;
        me->type = type;
        me->mass = masses_of_types[type-1];
        me->rho = rho;
        me->solidTag = is_solid;

        add_particle(me,system);

        printf("atom %i at %lf,%lf,%lf\n",id,x,y,z);
    }

}
