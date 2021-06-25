/* *****************************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3.
See the file LICENSE.txt for details.
***************************************************************************************** */
#ifndef read_lammps_input_file_h
#define read_lammps_input_file_h

#include "particle_system.hpp"

void read_lammps_input_file(const char*filename, ParticleSystem*system);

#endif //read_lammps_input_file_h
