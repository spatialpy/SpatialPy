/* *****************************************************************************
SSA-SDPD simulation engine
Copyright 2018 Brian Drawert (UNCA)

This program is distributed under the terms of the GNU General Public License.
See the file LICENSE.txt for details.
***************************************************************************** */
#ifndef simulate_hpp
#define simulate_hpp

namespace Spatialpy{
    void run_simulation(int num_threads, ParticleSystem *system);
    void take_step(Particle* me, ParticleSystem* system, unsigned int step, unsigned int substep);
    unsigned int get_number_of_substeps();
    void take_step1(Particle* me, ParticleSystem* system, unsigned int step);
    void take_step2(Particle* me, ParticleSystem* system, unsigned int step);
    void compute_forces(Particle* me, ParticleSystem* system, unsigned int step);
    //void compute_bond_forces(bond_t* this_bond, system_t*system, unsigned int step);
}

#endif // simulate_h
