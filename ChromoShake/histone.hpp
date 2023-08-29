//
//  histone.hpp
//  ChromoShake
//
//  Created by HeYunyan on 9/26/17.
//  Updated by Ben Walker on 3/1/21

#ifndef histone_hpp
#define histone_hpp
#include <stdio.h>
#include <unordered_set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <random>

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPath3D.h"
#include "DTProgress.h"
#include "DTSeriesPath3D.h"
#include "DTRandom.h"
#include "DTUtilities.h"
#include "DTPointCollection3D.h"
#include "DTRandom.h"

#include "Model.hpp"

// The structure of each histone is modeled in the following struct. The main objects it contains are the timer and the position it bound to.
// If the histone is active, the startBeadIndex is the bead index it binds to, hence >= 0; If it is inactive, then startBeadIndex < 0.
struct HistoneUnit
{
    //constructor
    HistoneUnit() {start_bead_index = 0; remaining_time = 0; time_started = 0;};
    HistoneUnit(int i, double j) {start_bead_index = i; remaining_time = j; time_started = 0;};
    
    int start_bead_index;
    double remaining_time;
    double time_started;
};

// The class contains all the histone objects in the loopConfig variable. Besides, it stores some universal variables sharing across all histone units.
// The structure keeps track of all the beads available as a new histone attachment site candidate in availableBeads variable, in preparation for newly activated histone units to attach to the chromosome without overlapping.
class Histone : public Module
{
public:
    Histone() { };
    Histone(Model *_m, const DTDictionary &histone_params, int number_of_beads, mt19937_64 &generator, const DTDoubleArray &chains_array, double base_spring_constant);

    ~Histone() {};
    
    void updateTime(const double dt);
    void updateForce();
    void histoneDetach(const int index, const double detach_time);
    void histoneInsert(const int index, const double binding_time);
    void unbindUpdate();
    
    // loopConfig.startBeadIndex stores starting node of each loop, loopConfig.remainingTime stores loop time.
    std::vector<HistoneUnit> histone_config;
    // availableBeads stores beads that are candidate loop start.
    std::unordered_set<int> available_beads;
    // beadsNotInLoop stores beads that are not in loops.
    std::unordered_set<int> beads_not_in_loop;
    // whether a bead is theoretically available to bond just based on geometry
    std::vector<bool> bondable;
    // where a connection is to if it is made
    std::vector<int> partner;
    
    int number_of_histones, number_of_beads, beads_in_histone_loop;
    double mean_on, mean_off;
    double unattach_time, unattach_start_time;
private:
    // Model we are attached to
    Model *m;
    
    int getAvailableBead();
    void updateAfterAttachment(int index);
    void updateAfterDetachment(int index);
    
    
    bool random_placement;
    
    // for random number generator
    std::mt19937_64 *rng_p;
    std::uniform_real_distribution<double> dist_01_cont;
    std::exponential_distribution<double> dist_off;
    std::exponential_distribution<double> dist_on;
    
    bool tension_unbind;
    double unattach_threshold;
    double spring_constant;
};

#endif /* histone_hpp */
