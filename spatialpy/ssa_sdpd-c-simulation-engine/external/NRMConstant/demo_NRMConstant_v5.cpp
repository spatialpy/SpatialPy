// compile: g++ *.cpp
// run: a.out
#include <iostream>
#include <random>
#include "NRMConstant_v5.hpp"

int main() {
    std::cout << "testing v5...(doesn't use ParticleSystem!)" << "\n";

    /* SET UP RANDOM NUMBER GENERATOR AND "MODEL" */
    // gen is the random number generator
    std::mt19937 gen;
    //    gen.seed(123);
    std::random_device rd;
    gen.seed(rd());
    
    // create a "MODEL" (just some random propensity values)
    std::size_t modelsize=10;//total number of propensities
    std::vector<double> propensities(modelsize); //vector of all propensities  (reaction and diffusion)
    double propensitySum=0.0;
    std::size_t activeChannels=0; // initial number of nonzero propensities (typically improves initial performance)
    
    
    //let's fill the propensities with random values using a uniform continuous distribution from 0 to 100 (mean 50)
    std::uniform_real_distribution<double> uniform100(0.0,100.0);
    for (int i=0; i<modelsize; i++) {
        double next_propensity = uniform100(gen); // propensity i would normally be calculated based on state
        std::cout << "propensity " << i << " = " << next_propensity << "\n";
        propensities[i]=next_propensity;
        propensitySum+=next_propensity;
        activeChannels++;
    }
    std::cout << "propensitySum=" << propensitySum << "\n";
    std::cout << "activeChannels=" << activeChannels << "\n";
    
    /* HERE IS THE PRIORITY QUEUE "NRMConstant_v5" */
    NRMConstant_v5 reactionGenerator;
    // set up the priority queue:
    reactionGenerator.build(propensities, gen, propensitySum, activeChannels);
    
    /* RUN A "SIMULATION" */
    // selectReaction() returns a pair: the event firing time and the event index (index of propensity)
    std::pair<double,int> timeRxnPair;
    int reactionIndex;
    double currentTime=0;
    double endTime=.1;
    while (currentTime<endTime) {
        timeRxnPair=reactionGenerator.selectReaction();
        currentTime=timeRxnPair.first;
        reactionIndex=timeRxnPair.second;
        std::cout << "firing reaction " << reactionIndex << " at t = " << currentTime << "\n";
        // here is where we would update state, recalculate all affected propensities
        // suppose reaction i affects (itself and) reaction i+1 (would typically be done with a dependency graph)
        int affectedReaction = (reactionIndex+1) % modelsize;
        // update reaction that fired
        propensities[reactionIndex]=uniform100(gen); // propensity would normally be calculated based on state
        // update affected reaction
        propensities[affectedReaction]=uniform100(gen); // propensity would normally be calculated based on state
        // update the priority queue for all changed propensities
        // update takes the affected reaction index, the new propensity value, current time, and random generator
        reactionGenerator.update(reactionIndex,propensities[reactionIndex],currentTime,gen);
        reactionGenerator.update(affectedReaction,propensities[affectedReaction],currentTime,gen);
    }
    
    return 0;
}
