//
//  condensin.cpp
//  ChromoShake
//
//  Created by HeYunyan on 10/2/17.
//  Updated by Ben Walker on 3/1/21

#include "condensin.hpp"
#include <exception>
typedef pair<int,int> ipair;

// Initialization.
Condensin::Condensin(Model *_m, const DTDictionary &condensin_params, const int _number_of_beads, const double seed, const DTDoubleArray &chains_array, mt19937_64 &random, double base_spring_constant) : bondable(_number_of_beads), m(_m)
{
    int _number_of_condensin = condensin_params("numberOfCondensins");
    judge_time = condensin_params("judgeTimeSeconds");     // rename to the mean of exp distribution
    increment  = condensin_params("incrementNanometers");
    critical_length = condensin_params("criticalLengthNanometers");
    unbind_time = condensin_params("unbindTime");
    unbind_time_var = condensin_params("unbindTimeVar");
    rebind_time = condensin_params("rebindTime");
    polarity_switch_time = condensin_params("polaritySwitchTime");
    unbind = condensin_params("unbind");
    dynamic_extrusion = condensin_params("dynamicExtrusionRate");     // this is a boolean flag
    extrusion_rate_power = condensin_params("extrusionRatePower");
    
    double condensin_modulus = condensin_params("condensin_modulus");
    spring_constant = condensin_modulus * base_spring_constant;
    number_of_beads = _number_of_beads;
    
    initial_time = judge_time;
    
    rng_p = &random;
    // distribution for average wait time until unbinding
    dist_on = std::exponential_distribution<double>(1.0/unbind_time);
    dist_off = std::exponential_distribution<double>(1.0/rebind_time);
    dist_polarity = std::exponential_distribution<double>(1.0/polarity_switch_time);
    dist_01 = std::uniform_int_distribution<int>(0, 1);
    dist_judge = std::uniform_real_distribution<double>(0, 1);
    
    // initialized at 0, incremented inside addCondensin loop
    number_of_condensin = 0;
    number_of_active_condensins = number_of_condensin;
    
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    /////////////     Initializing condensins in the beginning     //////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    int n = chains_array.n();
    vector<int> tmp;
    
    for (int i=0; i<number_of_beads; ++i)
    {
        // not used, but good to keep around
        bondable[i] = false;
    }
    
    for (int k=0; k<n; ++k)
    {
        // check if this chain is eligible (DNA chain) to have condensins bound to it
        // chains_array(7,k) =0 for DNA chains;  =1 not DNA, don't include beads in those chains in bondable
        if (chains_array(7,k))
        {
            continue;   // if not, no need to do anything
        }
        
        int chain_start = chains_array(0,k);
        int chain_end   = chains_array(1,k);
        
        for (int i=chain_start; i<=chain_end; ++i)
        {
            // reset bondanble to true for available beads in DNA chains
            bondable[i] = true;
            // save indices of bondable beads
            available_beads.insert(i);
            tmp.push_back(i);
        }
    }
    
    
    for (int i=0;i<_number_of_condensin;i++)
    {
        // this function increments number_of_condensins
        addCondensin();
    }
    
    condensinArrayUpdate();
}// end of constructor



// Update condensin_array by extracting indices of beads to which condensins are attached and save to output
void Condensin::condensinArrayUpdate()
{
    DTMutableIntArray tmp(2, number_of_condensin);
    
    for (int i = 0; i < number_of_condensin; ++i)
    {
        if (condensin_config[i].flag)
        {
            tmp(0,i) = condensin_config[i].weak_site_index;
            tmp(1,i) = condensin_config[i].strong_site_index;
        }
        else
        {
            // represent unbound with a value of -1
            tmp(0,i) = -1;
            tmp(1,i) = -1;
        }
    }
    
    condensin_array = tmp.Copy();
}


pair<int,int> Condensin::getAvailableBeadPair()
{
    int chosen_bead = sampleFromSet(available_beads, rng_p);
    int other_bead;
    
    if (m->nextBead(chosen_bead) == -1 && m->prevBead(chosen_bead) == -1)
    {
        DTErrorMessage("Both previous and next beads do not exist for a current index.  Error in setting available beads to bind to");
        throw std::runtime_error("Both previous and next beads do not exist for a current index.  Error in setting available beads to bind to");
    }
    // may need to change how condensins are seeded (i.e. using proximity not chain index) ::DK
    if (m->prevBead(chosen_bead) == -1)
    {
        other_bead = m->nextBead(chosen_bead);
    }
    else if (m->nextBead(chosen_bead) == -1)
    {
        other_bead = m->prevBead(chosen_bead);
    }
    else
    {
        if (dist_01(*rng_p) == 1)
        {
            other_bead = m->prevBead(chosen_bead);
        }
        else
        {
            other_bead = m->nextBead(chosen_bead);
        }
    }
    
    return ipair(chosen_bead, other_bead);
}// end of getAvailableBeadPair()


// awkward helper function to figure out which way a condensin should move so the strong
// site is going away from the weak site
bool get_polarity(int a, int b)
{
    if (a-b == 1)
    {
        return true;
    }
    else if (b-a == 1)
    {
        return false;
    }
    else
    {
        // note: this will stop working if geometries other than line/loop are used
        return b > a;
    }
}


// Add a newly active condensin. condensin_config will be updated.
void Condensin::addCondensin()
{
    //std::normal_distribution<double> norm_dist(unbind_time, unbind_time_var);
    
    auto pair = getAvailableBeadPair();
    bool polarityIs = get_polarity(pair.first, pair.second);
    
    CondensinUnit new_condensin(pair.first, pair.second, judge_time, dist_on(*rng_p), dist_polarity(*rng_p), polarityIs);
    condensin_config.push_back(new_condensin);
    
    number_of_condensin++;
    number_of_active_condensins++;
    condensinArrayUpdate();
} // add condensin if necessary


// condesninConfig stores condensins objects with relevant information about individual condensinUnits(strongSiteIndex = i; weakSiteIndex = j; judgeTime = t1; attachTime = t2; flag = 1; judgeTimeTracker = 0)

void Condensin::updateTime(const double dt)
{
    double dist, new_x, new_y, new_z;  // what did we use this for?
    int strong_bead, weak_bead;
    bool condensin_changed = false;
    //std::normal_distribution<double> norm_dist(unbind_time, unbind_time_var);
    
    for (int i=0; i<number_of_condensin; i++)
    {
        condensin_config[i].attach_time -= dt;
        condensin_config[i].polarity_switch_time -= dt;
        
        // Condition 1: condensin is active and condition change. Condensin needs to unbind and reset the timer.
        // bool unbind is set in parameters
        // all three have to be true,   When is attachTime less than ==> when condensins time is up
        if (unbind && condensin_config[i].attach_time < 0 && condensin_config[i].flag == 1)
        {
            // set condensin to be OFF for the duration of time
            condensin_config[i].flag = 0;
            // obtain time which will take condensin to rebind
            condensin_config[i].attach_time = dist_off(*rng_p);
            // indicate change has been made
            condensin_changed = true;
            // decrement condensins  (NOT USED currently)
            number_of_active_condensins--;
            continue;
        }
        // Condition 2: condensin is inactive and condition change. Condensin needs to rebind the reset the timer.
        else if (unbind && condensin_config[i].attach_time < 0 && condensin_config[i].flag == 0)
        {
            //dettached condensin rebinds on the chain.
            condensin_config[i].flag = 1;
            auto pair = getAvailableBeadPair();
            condensin_config[i].strong_site_index = pair.first;
            condensin_config[i].weak_site_index = pair.second;
            // determines how long does Condensin stays attached
            condensin_config[i].attach_time = dist_on(*rng_p);
            condensin_config[i].judge_time = judge_time;
            condensin_config[i].polarity = get_polarity(pair.first, pair.second);
            condensin_config[i].polarity_switch_time = dist_polarity(*rng_p);
            condensin_changed = true;
            // increment condensins  (NOT USED currently)
            number_of_active_condensins++;
            continue;
        }
        // polarity is checked as the timer is expired and switched from current polarity to opposite
        else if (condensin_config[i].polarity_switch_time < 0)
        {
            condensin_config[i].polarity = !condensin_config[i].polarity;
            condensin_config[i].polarity_switch_time = dist_polarity(*rng_p);
        }
        
        // Otherwise, no condition change. Work on extrusion.
        condensinTimeUpdate(i);
        
        // draw a probability of extrusion at the timestep dt
        double break_prob = 1 - exp(-dt/condensin_config[i].judge_time);
        
        // Take one step of extrusion.
        double rand_value = dist_judge(*rng_p);
        if (condensin_config[i].flag == 1 && rand_value < break_prob)
        {
            strong_bead = condensin_config[i].strong_site_index;
            weak_bead = condensin_config[i].weak_site_index;
            dist = 0;
            double diff[3];
            // get distance between strong and weak sites
            for (int k=0; k<3; ++k)
            {   // [3] next bead - weak bead  ... with the polarity on
                diff[k] = m->position(k, strong_bead) - m->position(k, weak_bead);
                dist += diff[k]*diff[k];
            }
            dist = sqrt(dist);
            
            // if this distance is larger than critical, unbind the weak site and rebind in proximity of strong site
            // direction of stepping can change, but polarity_timer should not be affected
            if (dist>critical_length)
            {
                condensin_config[i].weak_site_index = m->findClosestBead(m->position(0,strong_bead),
                                                                         m->position(1,strong_bead),
                                                                         m->position(2,strong_bead),
                                                                         weak_bead, strong_bead, bondable);
            }
            else
            {
                new_x = m->position(0,strong_bead) + increment*diff[0]/dist;
                new_y = m->position(1,strong_bead) + increment*diff[1]/dist;
                new_z = m->position(2,strong_bead) + increment*diff[2]/dist;
                // To allow for re-binding to the strong bead pass -1 as in example below...    ::D
                // condensin_config[i].strong_site_index = m->findClosestBead(new_x, new_y, new_z, -1, weak_bead, bondable);
                condensin_config[i].strong_site_index = m->findClosestBead(new_x, new_y, new_z, strong_bead, weak_bead, bondable);
                
                // Condensin is randomly reseeded when it reaches the end beads of the chain
                if (m->prevBead(condensin_config[i].strong_site_index)==-1 || m->nextBead(condensin_config[i].strong_site_index)==-1 )
                {
                                        condensin_config[i].flag = 0;
                                        condensin_config[i].attach_time = dist_off(*rng_p);
                                        number_of_active_condensins--;
                }
                
                ///////////////////////////////  INDEX BASED STEPPING WITH POLARITY ///////////////////////////////
                //                // take as indexed-based step maintaining polarity
                //                int new_bead;
                //                if (condensin_config[i].polarity)
                //                {
                //                    new_bead = m->nextBead(condensin_config[i].strong_site_index);
                //                }
                //                else
                //                {
                //                    new_bead = m->prevBead(condensin_config[i].strong_site_index);
                //                }
                //
                //                if (new_bead != -1)
                //                {
                //                    condensin_config[i].strong_site_index = new_bead;
                //                }
                //                else
                //                {
                //                    // [1] either switch polarity...
                //                    // condensin_config[i].polarity = !condensin_config[i].polarity;
                //
                //                    // [2] ...or unbind instead,
                //                    // unbind flag has to be set to 1 in ImageTank for rebinding to occur
                //                    condensin_config[i].flag = 0;
                //                    condensin_config[i].attach_time = dist_off(*rng_p);
                //                    number_of_active_condensins--;
                //                }
                //////////////////////   INDEX-BASED STEPPING ///////////////////
                
            }
            condensin_changed = true;
        }
    }
    // removed condensinArrayUpdate() outside of the for-loop
    if (condensin_changed)
    {
        condensinArrayUpdate();
    }
}// end of updateTime()



// Dynamic extrusion rate ONLY. Update the condensin timer based on temporary tension.
void Condensin::condensinTimeUpdate(const int condensin_index)
{
    // Detect greatest tension on the beads in the DNA chain
    if (dynamic_extrusion)
    {
        int strong_bead = condensin_config[condensin_index].strong_site_index;
        int weak_bead = condensin_config[condensin_index].weak_site_index;
        double max1 = max(m->norm_array(strong_bead), m->norm_array(weak_bead));
        
        double tmp1 = 0, tmp2 = 0;
        // check if prevbead ==-1  for strong and weak
        if (m->prevBead(strong_bead) != -1)
        {
            tmp1 = m->norm_array(m->prevBead(strong_bead));
        }
        
        if (m->prevBead(weak_bead) != -1)
        {
            tmp2 =  m->norm_array(m->prevBead(weak_bead));
        }
        
        double max2 = max(tmp1, tmp2);
        max1 = max(max1, max2);
        
        if (max1 < m->spring_length)
        {
            condensin_config[condensin_index].judge_time = initial_time;
        }
        else
        {
            condensin_config[condensin_index].judge_time = initial_time * pow((max1 / m->spring_length), extrusion_rate_power);
        }
    }
    else
    {
        condensin_config[condensin_index].judge_time = initial_time;
    }
}// end of condensinTimeUpdate()


// Update condensin-generated force to forceArray.
void Condensin::updateForce()
{
    double norm, scale;
    
    for (int i=0;i<number_of_condensin;i++)
    {
        if (!condensin_config[i].flag)
        {
            continue;
        }
        int strong_site_index = condensin_config[i].strong_site_index;
        int weak_site_index = condensin_config[i].weak_site_index;
        
        if (strong_site_index == weak_site_index)
        {
            continue;
        }
        
        double diff[3];
        
        norm = 0;
        for (int k=0;k<3;++k)
        {
            diff[k] = m->position(k, weak_site_index) - m->position(k, strong_site_index);
            norm += diff[k]*diff[k];
        }
        norm = sqrt(norm);
        
        double subscale = spring_constant*(1 - m->spring_length/norm);
        
        for (int k=0;k<3;++k)
        {
            scale = subscale * diff[k];
            m->force_array(k, strong_site_index) += scale;
            m->force_array(k, weak_site_index) -= scale;
        }
    }
    // assert(m->nan_check());
} //end of updateForce()
