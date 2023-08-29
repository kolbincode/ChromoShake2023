//
//  cohesin.hpp
//  ChromoShake
//
//  Created by Benjamin Walker on 3/1/21.
//

#ifndef cohesin_hpp
#define cohesin_hpp

#include <stdio.h>
#include "Model.hpp"
#include <random>
#include <unordered_set>

class Cohesin : public Module
{
public:
    Cohesin(Model *_m, const DTDoubleArray &chainsArray, mt19937_64 *_rng_p, const DTDictionary &cohesinParams, double base_spring_constant);
    void updateTime(const double dt);
    void updateForce();
private:
    Model *m;
    
    int number_of_cohesins;
    std::vector<std::pair<int, int> > ring_indices;   // start and end index of each cohesin ring
    std::vector<int> bound_to;          // which ring each ring is bound to (-1 = not bound)
    std::vector<int> active_bead;       // if bound, which bead is the one that is holding the bond
    
    // stores the chain array index of each cohesin
    std::vector<int> chain_id;
    
    // in order to tell if rings are on DNA we need to know which beads are DNA
    std::unordered_set<int> dna_beads;
    std::vector<int> dna_chains;
    
    // location of COM of each ring
    DTMutableDoubleArray com_loc;
    // pairwise distance between rings (when first index > second index)
    DTMutableDoubleArray com_dist;
    
    void updateComLoc();
    void updateComDist();
    void tailBeadForce();
    
    std::exponential_distribution<double> dist_on;
    
    // Mersenne twister type of random generator
    std::mt19937_64 *rng_p;
    
    // spring constant for cohesin bonds (make parameter)
    double spring_constant;
    
    // rings closer than this will bond
    double bond_dist;
    
    // average duration of bonds
    double mean_on;
    
    // how long bonds have remaining
    std::vector<double> time_till_break;
    
    // whether we try to put cohesins that have fallen off back on the chain
    bool restore_cohesin;
    
    // beads at end of DNA chains
    std::vector<int> tail_beads;
    
    // interval at which we check for beads falling off
    double restore_interval;       // make parameter
    double time_since_check = 0;
    
    double tail_collision_dist;
};



/*
1. Excluded volume for cohesin
2. make slipped rings disappear ...  generate new rings
 
This code in in MODEL.cpp and invokes *tmp_spring_table*
3. Cohesins rings have same springs as DNA (2Gpa) but are floppier, no hinge forces (Lp5)
 .... set as parameter in ImageTank
 .... TODO: To set cohesin stiffness separately will need to modify code and script.
 
4. Careful with DisplacementCap
5. no interpolation done.   Check tolerance in save step
 

*/

#endif /* cohesin_hpp */
