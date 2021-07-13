#include <cmath>

#include "NRMConstant_v5.hpp"

namespace Spatialpy{

    NRMConstant_v5::NRMConstant_v5() : previousFiringTime(0.0), exponential(1) {}

    void NRMConstant_v5::build(std::vector<double>& propensities,
               double propensitySum, std::size_t activeChannels,
               std::mt19937_64& rng, double timeOffset,
               double simulationEndTime) {
            //  std::cout << "in build..." << std::endl;

        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
        previousFiringTime=timeOffset; //initialize (could set to nan or something instead)

        // nextFiringTime.resize(propensities.size());
        // binIndexAndPositionInBin.resize(propensities.size());

        if (propensitySum==0.0) {
            std::cerr << "ERROR: propensity sum is zero on initial build. Terminating." << std::endl;
            exit(1);
        }

        setBinNumberAndBounds(timeOffset,propensitySum,activeChannels); //choose lowerBound, upperBound, and numberOfBins
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;

        theHashTable.resize(numberOfBins,emptyVector);

        std::size_t active_channel_count=0; //for verifying build input

        for (std::size_t i=0; i!=propensities.size(); ++i) {
            double firingTime;
            if (propensities[i]==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            }
            else {
                firingTime=exponential(rng)/propensities[i]+timeOffset;
                ++active_channel_count;
            }

            nextFiringTime[i]=firingTime; //.insert(std::make_pair<reactionIndexT,firingTimeT>(i,firingTime));
            //  std::cout << "nextFiringTime["<<i<<"]="<<nextFiringTime[i]<<"\n";
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,i)); //place this rxn at back of bin
                binIndexAndPositionInBin[i]=std::make_pair(bin,theHashTable[bin].size()-1);
            } else {
                binIndexAndPositionInBin[i]=std::make_pair<int,int>(-1,-1); //bin (and index within bin) is -1 if not in the hash table
            }
        }

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound); //keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);

        // printHashTable();

        if (active_channel_count!=activeChannelCounter) {
            std::cout << "ERROR: active channel count is inconsistent.\n";
            exit(1);
        }

    //  std::cout << "exiting build...\n";
    }

    void NRMConstant_v5::build(std::vector<Spatialpy::Particle>& particles,
           double propensitySum, std::size_t activeChannels,
           std::mt19937_64& rng, double timeOffset,
           double simulationEndTime) {
    //  std::cout << "in build..." << std::endl;
        activeChannelCounter=activeChannels;
        endTime=simulationEndTime;
        previousFiringTime=timeOffset; //initialize (could set to nan or something instead)

        // nextFiringTime.resize(propensities.size());
        // binIndexAndPositionInBin.resize(propensities.size());

        if (propensitySum==0.0) {
            std::cerr << "ERROR: propensity sum is zero on initial build. Terminating." << std::endl;
            exit(1);
        }

        setBinNumberAndBounds(timeOffset,propensitySum,activeChannels); //choose lowerBound, upperBound, and numberOfBins
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;

        theHashTable.resize(numberOfBins,emptyVector);

        std::size_t active_channel_count=0; //for verifying build input

        std::cout << "?? is index into particles vector a valid unique id?\n";

        for(auto p = particles.begin(); p!=particles.end(); p++){
            double firingTime;
            double this_propensity = p->srrate+p->sdrate;
            unsigned this_id = p->id;
            if (this_propensity==0.0) {
                firingTime=std::numeric_limits<double>::infinity();
            }
            else {
                firingTime=exponential(rng)/this_propensity+timeOffset;
                ++active_channel_count;
            }

            nextFiringTime[this_id]=firingTime; //.insert(std::make_pair<reactionIndexT,firingTimeT>(i,firingTime));
            //  std::cout << "nextFiringTime["<<i<<"]="<<nextFiringTime[i]<<"\n";
            //insert into hashTable
            int bin=computeBinIndex(firingTime);
            if (bin>=0) {
                theHashTable[bin].push_back(std::make_pair(firingTime,this_id)); //place this rxn at back of bin
                binIndexAndPositionInBin[this_id]=std::make_pair(bin,theHashTable[bin].size()-1);
            } else {
                binIndexAndPositionInBin[this_id]=std::make_pair<int,int>(-1,-1); //bin (and index within bin) is -1 if not in the hash table
            }
        }

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(0);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound); //keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);

        // printHashTable();

        if (active_channel_count!=activeChannelCounter) {
            std::cout << "ERROR: active channel count is inconsistent.\n";
            exit(1);
        }

    //  std::cout << "exiting build...\n";

    }

    // reactionIndex is particle id
    void NRMConstant_v5::update(std::size_t reactionIndex, double newPropensity, double currentTime, std::mt19937_64& rng) {

        double firingTime;
        if (newPropensity <= 0.0) {
            firingTime=std::numeric_limits<double>::infinity();
            if (nextFiringTime.at(reactionIndex)!=std::numeric_limits<double>::infinity()) {
                activeChannelCounter--;
            }
        }
        else {
            firingTime=exponential(rng)/newPropensity+currentTime;
            if (nextFiringTime.at(reactionIndex)==std::numeric_limits<double>::infinity()) {
                activeChannelCounter++;
            }
        }
    //  std::cout << "new firingTime=" << firingTime << std::endl;
        nextFiringTime[reactionIndex]=firingTime;
        int newBin=computeBinIndex(firingTime);

    //  std::cout << "newBin=" << newBin << std::endl;
        int oldBin=binIndexAndPositionInBin[reactionIndex].first;
    //  std::cout << "oldBin=" << oldBin << std::endl;
        int oldPositionInBin=binIndexAndPositionInBin[reactionIndex].second;
        if (newBin!=oldBin) {
            if (oldBin>=0) {
                //remove from old bin
                if (theHashTable[oldBin].size()>1) {

                    //take last element in old bin and place in this element's spot
                    theHashTable[oldBin][oldPositionInBin]=theHashTable[oldBin].back();
                    std::size_t movedElementIndex=theHashTable[oldBin][oldPositionInBin].second;
                    //update old last element's ...within bin index
                    binIndexAndPositionInBin[movedElementIndex].second=oldPositionInBin;
                }
                theHashTable[oldBin].pop_back();
            }
            binIndexAndPositionInBin[reactionIndex].first=newBin;
            if (newBin>=0) {
                theHashTable[newBin].push_back(std::pair<double,std::size_t>(firingTime,reactionIndex));
                binIndexAndPositionInBin[reactionIndex].second=theHashTable[newBin].size()-1;
            }
        }
        else {
            //just update firing time
            if (newBin>=0) {
                theHashTable[newBin][binIndexAndPositionInBin[reactionIndex].second].first=firingTime;
            }
        }

    }

    void NRMConstant_v5::printHashTable() {
        std::cout << "minBin=" << minBin << std::endl;
        std::cout << "lower bound="<<currentLowerBound << std::endl;
        std::cout << "upper bound="<<currentUpperBound << std::endl;
        std::cout << "approx bin width=" << (currentUpperBound-currentLowerBound)/numberOfBins << std::endl;
        for (std::size_t i=0; i!=theHashTable.size(); ++i) {
            std::cout << "[i]=" << i << ": (size=" << theHashTable[i].size() << ")\n";
            for (std::size_t j=0; j!=theHashTable[i].size(); ++j) {
                std::cout << theHashTable[i][j].second << "(" << theHashTable[i][j].first << "), ";
            }
            std::cout << "\n";
        }
    }

    NRMConstant_v5::timeIndexPair NRMConstant_v5::selectReaction() {
    //  std::cout << "in selectReaction...\n";
        //while front of queue is empty, pop it
        while (theHashTable[minBin].size()==0) {
            ++minBin;
            if (minBin==theHashTable.size()) {
            //  std::cout << "rebuilding...\n";
                if (!rebuild()) {
                    return std::make_pair(std::numeric_limits<double>::max(),-1);
                }
            }
        }
        //  std::cout << "minBin="<<minBin<<std::endl;
        //now that we have a bin with elements, find minimum
        int minTimeRxnIndex=0;
        for (std::size_t i=1;i<theHashTable[minBin].size(); ++i) {
            if (theHashTable[minBin][i].first<theHashTable[minBin][minTimeRxnIndex].first) {
                minTimeRxnIndex=i;
            }
            /*if (theHashTable[minBin][i].first < 0.0) {
                printf("Index: %i tt: %f\n",theHashTable[minBin][i].second, theHashTable[minBin][i].first);
            }*/
        }
        /*if (theHashTable[minBin][minTimeRxnIndex].first < 0.0) {
            printf("First: %f\n", theHashTable[minBin][minTimeRxnIndex].first);
            printf("Selecting: %i\n", theHashTable[minBin][minTimeRxnIndex].second);
        }*/
        previousFiringTime=theHashTable[minBin][minTimeRxnIndex].first;
        ++rxnCountThisBuildOrRebuild;
    //  std::cout << "returning time=" <<theHashTable[minBin][minTimeRxnIndex].first<<", reaction index=" << theHashTable[minBin][minTimeRxnIndex].second << " in selectReaction\n";

        return theHashTable[minBin][minTimeRxnIndex];
    }

    NRMConstant_v5::~NRMConstant_v5() {
    }

    int NRMConstant_v5::computeBinIndex(double firingTime) {
        if (firingTime>currentUpperBound) {
        //  std::cout << "computeBinIndex("<<firingTime<<") returning binIndex=" << -1 << std::endl;
            return -1;
        }
        int binIndex=(int)((firingTime-currentLowerBound)/(currentUpperBound-currentLowerBound)*numberOfBins);
        int returnVal=std::max(binIndex,(int)minBin);
        //  std::cout << "computeBinIndex("<<firingTime<<") returning binIndex=" << binIndex << std::endl;
        return returnVal;
    }

    //returns false if all propensities are 0
    bool NRMConstant_v5::rebuild() {
    //  std::cout << "in rebuild...\n";

        //estimate propensitySum based on number of steps since last build or rebuild
        double propensitySumEstimate;
        if (rxnCountThisBuildOrRebuild>0) {
        //  std::cout << "fired " << rxnCountThisBuildOrRebuild << " since last build or rebuild\n";
            propensitySumEstimate=(double)rxnCountThisBuildOrRebuild/(currentUpperBound-currentLowerBound);
            //  std::cout << "propensity sum estimate = " << propensitySumEstimate << "\n";
        }
        else {
            //this shouldn't happen...if we got here, 0 reactions fired since last rebuild!
            //will only arise in toy problems
            //but need a strategy...
            double previousBinWidth=(currentUpperBound-currentLowerBound)/(double)numberOfBins;
            propensitySumEstimate=1.0/previousBinWidth;
            std::cout << "WARNING: 0 reactions fired before rebuild.\n";
            //  std::cout << "previousBinWidth was " << previousBinWidth <<"\n";
        }

        setBinNumberAndBounds(currentUpperBound,propensitySumEstimate,activeChannelCounter);//
        minBin=0;

        std::vector<std::pair<double, std::size_t> > emptyVector;

        //theHashTable.clear();
        theHashTable.resize(numberOfBins,emptyVector);

        bool allPropensitiesZero=true;
        //std::cout << "inserting into hash table...\n";
        //for (std::size_t i=0; i!=nextFiringTime.size(); ++i) {
        for (auto it = nextFiringTime.cbegin(); it != nextFiringTime.end(); it++) {
            //insert into hashTable

            if (allPropensitiesZero && it->second!=std::numeric_limits<double>::infinity()) {
                allPropensitiesZero=false;
            }

            int bin=computeBinIndex(it->second); // it->second is firing time
            if (bin>=0) {
                //  std::cout << "inserting into bin " << bin << std::endl;
                theHashTable[bin].push_back(std::make_pair(it->second,it->first)); //place this rxn at back of bin
                binIndexAndPositionInBin[it->first]=std::make_pair(bin,theHashTable[bin].size()-1);
            }
            else {
                binIndexAndPositionInBin[it->first]=std::make_pair<int,int>(-1,-1);
            }
        }

        //record info from build
        rxnsBetweenBuildOrRebuild.push_back(rxnCountThisBuildOrRebuild);
        lowerBoundsBuildOrRebuild.push_back(currentLowerBound); //keep entry for each build/rebuild; lowerBound.back() corresponds to currentLowerBound
        upperBoundsBuildOrRebuild.push_back(currentUpperBound);

        //set rxn counter to 0
        rxnCountThisBuildOrRebuild=0;

        //  std::cout << "after rebuild, hash table:\n";
        //  printHashTable();

        return !allPropensitiesZero;
    }

    //fixed number of bins
    //currentUpperBound=
    void NRMConstant_v5::setBinNumberAndBounds(double newLowerBound, double propensitySum, int activeChannels) {
        if (activeChannels==0) {
            std::cout << "ERROR: setBinNumberAndBounds not set up to handle activeChannels=0\n";
            exit(1);
        }

        //  std::cout << "...in setBinNumberAndBounds...newLowerBound=" << newLowerBound << std::endl;
        currentLowerBound=newLowerBound;
        if (currentLowerBound>endTime) {
            std::cerr << "ERROR: calling rebuild when simulation end time exceeded. Terminating.\n";
            exit(1);
        }

        //  std::cout << "propensitySum is " << propensitySum << ", so step size is " << 1.0/propensitySum << "\n";
        double binWidth=16.0/propensitySum;
        //  std::cout << "binWidth=" << binWidth << "\n";

        // heuristic
        std::size_t maxNumberOfBins=(std::size_t)20.0*sqrt(activeChannels);

        //  std::cout << "maxNumberOfBins=" << maxNumberOfBins << "\n";

        //how many bins to get beyond end time?
        std::size_t endTimeBins=std::numeric_limits<std::size_t>::max();
        //if we know the end time:
        if (endTime<std::numeric_limits<double>::max()) {
        //  std::cout << "endTime is set, how many bins needed to exceed...\n";
            endTimeBins=(std::size_t)ceil( (endTime-currentLowerBound)/binWidth);
            //  std::cout << "endTimeBins="<<endTimeBins << "\n";
            //use number of bins that takes us past endTime, or maxNumberOfBins, whichever is smaller
            numberOfBins=std::min(endTimeBins,maxNumberOfBins);
        }
        else {
            numberOfBins=maxNumberOfBins;
        }
        currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;

        //if max bins gets us very close to the simulation end time, bump it up beyond max to avoid
        //an additional rebuild
        if (endTime<std::numeric_limits<double>::max()) {
        //  std::cout << "checking to see if we can use fewer than max bins...\n";
            if (currentUpperBound<endTime && ((double)(endTimeBins-maxNumberOfBins)/((double)numberOfBins))<.2) {
                numberOfBins=endTimeBins;
                currentUpperBound=currentLowerBound+binWidth*(double)numberOfBins;
            }
        }

    //  std::cout << "currentLowerBound=" << currentLowerBound << ", currentUpperBound=" << currentUpperBound << std::endl;
    //  std::cout << "binWidth=" << binWidth << std::endl;
    //  std::cout << "numberOfBins=" << numberOfBins << std::endl;
    }
}
