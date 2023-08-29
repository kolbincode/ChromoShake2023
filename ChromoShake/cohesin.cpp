//
//  cohesin.cpp
//  ChromoShake
//
//  Created by Benjamin Walker on 3/1/21.
//

#include "cohesin.hpp"
#include <unordered_set>


Cohesin::Cohesin(Model *_m, const DTDoubleArray &chains_array, mt19937_64* _rng_p, const DTDictionary &cohesin_params, double base_spring_constant) : rng_p(_rng_p), m(_m)
{
    // load parameters from dictionary
    mean_on = cohesin_params("meanOn");
    bond_dist = cohesin_params("bond_dist");
    double cohesin_modulus = cohesin_params("cohesin_modulus");
    spring_constant = base_spring_constant * cohesin_modulus;
    restore_cohesin = cohesin_params("restore");
    restore_interval = cohesin_params("restore_interval");
    tail_collision_dist = cohesin_params("tail_collision_dist");
    
    // go through chains array and identify cohesin rings (labeled 1)
    int n = chains_array.n();
    number_of_cohesins = 0;
    for (int k=0;k<n;++k)
    {
        if (chains_array(7,k) == 0)
        {
            // this is a DNA chain, add to relevant array
            for (int i=chains_array(0,k); i<=chains_array(1,k);++i)
            {
                dna_beads.insert(i);
            }
            
            dna_chains.push_back(k);
            
            if (chains_array(2,k) == 0)
            {
                // mark end beads for linear chains only (DK)
                tail_beads.push_back(chains_array(0, k));
                tail_beads.push_back(chains_array(1, k));
            }
        }
        else if (chains_array(7,k) == 1)
        {
            // this is a cohesin ring, process accordingly
            ++number_of_cohesins;
            bound_to.push_back(-1);    // this ring currently unbound
            ring_indices.push_back(pair<int,int>(chains_array(0,k), chains_array(1,k)));
            chain_id.push_back(k);
        }
    }
    
    active_bead.resize(number_of_cohesins);
    time_till_break.resize(number_of_cohesins);
    
    dist_on = std::exponential_distribution<double>(1.0/mean_on);
    
    // center of mass array is 3 x nCohesin
    com_loc = DTMutableDoubleArray(3, number_of_cohesins);
    com_dist = DTMutableDoubleArray(number_of_cohesins, number_of_cohesins);
}


void Cohesin::updateTime(const double dt)
{
    // rescue any rings that have fallen off the DNA
    // TODO: do we want to do this every single timestep?
    time_since_check += dt;
    if (restore_cohesin && time_since_check >= restore_interval)
    {
        time_since_check = 0;
        m->updateTree();
        std::vector<int> intersect_list;
        for (int i=0; i<number_of_cohesins;++i)
        {
            // go through the list and check for presence of DNA beads
            bool flag = true;
            for (auto v : dna_chains)
            {
                if (m->chainIntersect(chain_id[i], v))
                {
                    // there is an intersecting DNA bead so we're good
                    flag = false;
                    break;
                }
            }
            // if flag is still true there weren't any DNA beads
            if (flag)
            {
                // we fell off so do something
                
                // pick a random DNA bead and put it there
                int bead = sampleFromSet(dna_beads, rng_p);
                // resample to ensure that we get a bead with both a previous and a next
                while (!(m->prevBead(bead)!=-1 && m->nextBead(bead)!=-1))
                {
                    bead = sampleFromSet(dna_beads, rng_p);
                }
                // compute the tangent direction
                double diff[3];
                int next_bead = m->nextBead(bead), prev_bead = m->prevBead(bead);
                for (int k=0;k<3;++k)
                {
                    diff[k] = m->position(k, next_bead) - m->position(k, prev_bead);
                }
                double com[3];
                for (int k=0; k<3; ++k)
                {
                    com[k] = 0;
                }
                for (int j=ring_indices[i].first;j<=ring_indices[i].second;++j)
                {
                    for (int k=0;k<3;++k)
                    {
                        com[k] += m->position(k, j);
                    }
                }
                const int beads_in_ring = ring_indices[i].second - ring_indices[i].first + 1;
                for (int k=0;k<3;++k)
                {
                    com[k] /= beads_in_ring;
                }
                
                // TODO: actually rotate it or something
                
                for (int j=ring_indices[i].first;j<=ring_indices[i].second;++j)
                {
                    for (int k=0;k<3;++k)
                    {
                        m->position(k, j) += m->position(k, bead) - com[k];
                    }
                }
            }
        }
    }
    
    // update the COM information
    //updateComLoc();
    //updateComDist();
    
    //double break_prob = 1 - exp(-dt/mean_on);    // fix it
    
    // unbind step - check if any bonds are past their death time
    for (int i=0; i<number_of_cohesins; ++i)
    {
        int partner = bound_to[i];
        if (partner != -1)
        {
            time_till_break[i] -= dt;
            if (time_till_break[i] <= 0)
            {
                // break the bond
                bound_to[i] = -1;
                bound_to[partner] = -1;
            }
        }
    }
    
    
    // bind step - see if we form any bonds between unbound rings
    // TODO: figure out the best way of doing this
    
    const double *mPositionD = m->position.Pointer();
    
    const double bond_dist2 = bond_dist*bond_dist;
    for (int i=0;i<number_of_cohesins-1; ++i)
    {
        if (bound_to[i] != -1)
        {
            // already bound, nothing to do
            continue;
        }
        for (int j=i+1;j<number_of_cohesins;++j)
        {
            //if (m->chainIntersect(chain_id[i], chain_id[j]))
            //if (com_dist(i, j) < bond_dist)
            if (m->chainSeparation(i, j) < bond_dist)
            {
                // figure out what pair of beads is closest
                int i1=-1, i2=-1;
                // max finite number in C++
                double min_dist = __DBL_MAX__;
                for (int b1 = ring_indices[i].first; b1 <= ring_indices[i].second; ++b1)
                {
                    for (int b2 = ring_indices[j].first; b2 <= ring_indices[j].second; ++b2)
                    {
                        double dist = 0;
                        for (int k=0;k<3;++k)
                        {
                            // double tmp = m->position(k,b2) - m->position(k,b1);
                            double tmp = mPositionD[k+3*b2] - mPositionD[k+3*b1];
                            dist += tmp*tmp;
                        }
                    
                        if (dist<min_dist)
                        {
                            min_dist = dist;
                            i1 = b1;
                            i2 = b2;
                        }
                    }
                }

                if (min_dist > bond_dist2)
                {
                    continue;
                }
                
                // form a bond between i and j
                bound_to[i] = j;
                bound_to[j] = i;
                // i1 and i2 hold the indices of the minimum distance pair
                active_bead[i] = i1;
                active_bead[j] = i2;
                
                // draw the remaining time
                double bond_duration = dist_on(*rng_p);
                time_till_break[i] = bond_duration;
                time_till_break[j] = bond_duration;
            }
        }
    }
}//end of updateTime(dt,*m)


void Cohesin::updateForce()
{
    for (int i=0; i<number_of_cohesins; ++i)
    {
        // check if bound
        int partner = bound_to[i];
        
        // to prevent double-counting force, only apply force from the numerically lesser id
        if (partner != -1 && partner > i)
        {
            // find the force on this ring's bead - the opposite force is computed in that ring's iteration of this loop
            // spring force based on difference between active beads
            const int this_bead = active_bead[i];
            const int partner_bead = active_bead[partner];
            for (int k=0;k<3;++k)
            {
                //m->force_array(k, active_bead[i]) += spring_constant*(m->position(k, active_bead[partner]) - m->position(k, active_bead[i]));
                
                double norm = 0;
                // stores dx,dy,dz
                double diff[3];
                for (int k=0;k<3;++k)
                {
                    diff[k] = m->position(k, partner_bead) - m->position(k, this_bead);
                    norm += diff[k]*diff[k];
                }
                
                norm = sqrt(norm);
                
                double subscale = spring_constant*(1-m->spring_length/norm);
                for (int k=0;k<3;++k)
                {
                    double scale = subscale*diff[k];
                    m->force_array(k,this_bead) += scale;
                    m->force_array(k,partner_bead) -= scale;
                }
            }
        }
    }
    
    // force to keep rings from falling off the end of linear chains
    tailBeadForce();
}


// Compute a repulsive force between cohesin rings and DNA tail beads
void Cohesin::tailBeadForce()
{
    // right now just a simple N^2 implementation
    for (auto p : ring_indices)
    {
        for (int i=p.first; i<p.second; ++i)
        {
            for (int j : tail_beads)
            {
                double dist=0;
                double diff[3];
                for (int k=0;k<3;++k)
                {
                    diff[k] = m->position(k, j)-m->position(k, i);
                    dist += diff[k]*diff[k];
                }
                const double norm = sqrt(dist);
                // TODO: make these parameters configurable?
                // radius of tail bead
                //const double tail_collision_dist = 30;       // is now a parameter
                if (norm < tail_collision_dist)
                {
                    // cohesin spring constant cross linking
                    const double scale = spring_constant*(tail_collision_dist-norm)/norm;
                    for (int k=0;k>3;++k)
                    {
                        m->force_array(k,i) += scale*diff[k];
                        m->force_array(k,j) -= scale*diff[k];
                    }
                }
            }
        }
    }
}


void Cohesin::updateComLoc()
{
    // update the positions of center-of-mass of all cohesin rings
    com_loc = 0;
    for (int i=0;i<number_of_cohesins;++i)
    {
        int ringLength = ring_indices[i].second - ring_indices[i].first + 1;
        for (int j=ring_indices[i].first; j<=ring_indices[i].second; ++j)
        {
            for (int k=0;k<3;++k)
            {
                com_loc(k, i) += m->position(k, j);
            }
        }
        for (int k=0;k<3;++k)
        {
            com_loc(k, i) /= ringLength;
        }
    }
}


void Cohesin::updateComDist()
{
    // TODO: use collision tree info somehow here?
    com_dist = 0;
    for (int i=0;i<number_of_cohesins-1;++i)              
    {
        for (int j=i+1;j<number_of_cohesins;++j)
        {
            for (int k=0;k<3;++k)
            {
                double tmp = com_loc(k, j) - com_loc(k, i);
                com_dist(j, i) += tmp*tmp;
            }
            com_dist(j, i) = sqrt(com_dist(j, i));
        }
    }
}
