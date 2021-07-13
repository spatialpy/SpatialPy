Sanft added NRMConstant_v4 class (NRMConstant_v4.hpp, .cpp)
The class is a priority queue that should scale roughly O(1)
for large discrete stochastic reaction-diffusion models.

demo_NRMConstant.cpp shows how to use it. Compile with:
g++ demo_NRMConstant.cpp NRMConstant_v4.cpp

Note the key functions:
NRMConstant_v4 reactionGenerator; // constructor
reactionGenerator.build(propensities, gen, propensitySum, activeChannels); // set up priority queue from initial propensities
reactionGenerator.selectReaction(); // get the item (next reaction or diffusion event) that will occur next
reactionGenerator.update(reactionIndex,propensities[reactionIndex],currentTime,gen); // update event time

Note that the event that fires and all affected propensities should have a corresponding call to update before
calling selectReaction

