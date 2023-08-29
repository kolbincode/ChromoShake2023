//
//  condensin.hpp
//  ChromoShake
//  Updated by Ben Walker on 3/1/21

#ifndef condensin_hpp
#define condensin_hpp
#include <stdio.h>
#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPath3D.h"
#include "DTProgress.h"
#include "DTSeriesPath3D.h"
#include "DTRandom.h"
#include "DTUtilities.h"
#include "DTPointCollection3D.h"
#include <algorithm>
#include <random>
#include "DTRandom.h"
#include <set>
#include "unitConversion.h"
#include <unordered_set>
#include "Model.hpp"



// The structure of each condensin is modeled in the following struct condensinUnit. The main objects it contains are the timers and the positions it binds to: strongSiteIndex and weakSiteIndex.
// flag determines if it is activated.
struct CondensinUnit
{
    // condensin constructor
    CondensinUnit()
    {      
        strong_site_index = -1; weak_site_index = -1; judge_time = 0; attach_time = 0; flag = 1;
    };
    
    CondensinUnit(int i, int j, double t1, double t2, double t3, bool _polarity)
    {
        strong_site_index = i;
        weak_site_index = j;
        judge_time = t1;
        attach_time = t2;
        polarity_switch_time = t3;
        polarity = _polarity;
        flag = 1;
    };
    
    int strong_site_index;
    int weak_site_index;
    double judge_time;
    double attach_time;
    double polarity_switch_time;
    bool flag;
    // true - advances to next bead. false - advances to previous bead
    bool polarity;
};

// The class contains all the condensin objects in the condensinConfig variable. Besides, it stores some universal variables sharing across all condensin units.
// condensin_array keeps a copy of a collection of binding sites of all condensins at the moment. Used for outputing condensin positions.

class Condensin : public Module
{
public:
    Condensin (){};
    Condensin(Model *_m, const DTDictionary &condensin_params, const int _number_of_beads, const double seed, const DTDoubleArray &chains_array, mt19937_64 &random, double base_spring_constant);
    
    int number_of_condensin, number_of_beads;
    int number_of_active_condensins;                    // NOT USED
    double judge_time, critical_length, initial_time;
    double increment;
    double unbind_time;
    double unbind_time_var;
    double rebind_time;
    double polarity_switch_time;
    bool unbind;
    bool dynamic_extrusion;
    double extrusion_rate_power;
    DTIntArray condensin_array;
    
    void updateTime(const double dt);
    void condensinTimeUpdate(const int condensin_index);
    void updateForce();
    // Updates value of class variable condensin_array
    void condensinArrayUpdate();
    void addCondensin();
    
    
    // condesninConfig stores condensins objects with relevant information about individual condensinUnits(strongSiteIndex = i; weakSiteIndex = j; judgeTime = t1; attachTime = t2; flag = 1; judgeTimeTracker = 0)
    std::vector<CondensinUnit> condensin_config;
    // stores a pair of available beads for condensin to bind.
    std::pair<int,int> getAvailableBeadPair();
    // availableBeads stores beads that are candidate loop start.
    std::unordered_set<int> available_beads;
    // whether a bead is theoretically available to bond based on initialization geometry
    std::vector<bool> bondable;
    
private:
    Model *m;
    
    // for randomness
    std::mt19937_64 *rng_p;
    std::uniform_int_distribution<int> dist_01;
    std::uniform_int_distribution<int> dist_20;
    std::exponential_distribution<double> dist_on;
    std::exponential_distribution<double> dist_off;
    std::exponential_distribution<double> dist_polarity;
    std::uniform_real_distribution<double> dist_judge;

    double spring_constant;
};

#endif /* condensin_hpp */
