//
//  histone.cpp
//  ChromoShake
//
//  Created by HeYunyan on 9/26/17.
//  Updated by Ben Walker on 3/1/21
//

#include "histone.hpp"
#include <iterator>

// Initialization function. Histone loops are guaranteed not to overlap.
Histone::Histone(Model *_m, const DTDictionary &histone_params, int _number_of_beads, mt19937_64 &generator, const DTDoubleArray &chains_array, double base_spring_constant) : bondable(_number_of_beads), partner(_number_of_beads), rng_p(&generator), m(_m)
// periodic = plasmid
{
    random_placement = histone_params("random_placement");
    unattach_threshold = histone_params("unattachThreshold");
    tension_unbind = histone_params("unattach");
    
    number_of_histones = histone_params("numberOfHistones");
    beads_in_histone_loop = histone_params("beads_in_histone_loop");
    
    number_of_beads = _number_of_beads;
    mean_on = histone_params("halfLifeOn");
    mean_off =  histone_params("halfLifeOff");
    // Generate initial histone configuration.
    unattach_time = histone_params("unattachTime");;
    unattach_start_time = histone_params("unattachStartTime");
    
    double histone_modulus = histone_params("histone_modulus");
    spring_constant = histone_modulus * base_spring_constant;
    
    dist_off = exponential_distribution<double>(1.0l/mean_off);
    dist_on = exponential_distribution<double>(1.0l/mean_on);
    
    int n = chains_array.n();
    
    // beads not explicitly marked as bindable below can't hold bonds
    for (int i=0;i<_number_of_beads;++i)
    {
        bondable[i] = false;
    }
    
    // set beads to be available to bind
    for (int k=0;k<n;++k)
    {
        // check if this chain is eligible to have histone bonds
        if (chains_array(7,k) != 0)
        {
            // it's not so no need to do anything
            continue;
        }
        int chain_start = chains_array(0,k);
        int chain_end = chains_array(1,k);
        int chain_loop = chains_array(2, k);      // boolean but can be an int
        int beads_in_chain = chains_array(1,k) - chains_array(0,k) + 1;
        int shift = beads_in_histone_loop;
        for (int i=chain_start;i<=chain_end;++i)
        {
            // Should we change 6 to user input DT value 'beadContainedInLoop' - 1
            // Loop though j = 'beadContainedInLoop in subsequent for loops
            
            
            if (i+(shift-1) <= chain_end)
            {
                bondable[i] = true;
                available_beads.insert(i);
                partner[i] = i+(shift-1);
            }
            else if (chain_loop && beads_in_chain >= shift)
            {
                // we're at the end, but it is a loop
                bondable[i] = true;
                available_beads.insert(i);
                partner[i] = i + (shift-1) - beads_in_chain;
            }
            // what does the following line do? Collects bead indices at the end of the chain?
            beads_not_in_loop.insert(i);
        }
    }
    
    for (int i=0; i<number_of_histones; i++)
    {
        
        if (!available_beads.empty())
        {
            // randomize half life of histone binding
            double initial_portion = dist_on(*rng_p);
            int start_location = getAvailableBead();
            
            HistoneUnit tmp_unit(start_location, initial_portion);
            histone_config.push_back(tmp_unit);
            updateAfterAttachment(start_location);
        }
        else
        {
            double initial_portion = dist_off(*rng_p);
            HistoneUnit tmp_unit(-1, initial_portion);
            histone_config.push_back(tmp_unit);
        }
    }
} //end of constructor

void Histone::updateAfterAttachment(int index)
{
    int current_bead = index;
    int shift = beads_in_histone_loop;
    
    for (int j=0; j<shift; j++)
    {
        beads_not_in_loop.erase(current_bead);
        available_beads.erase(current_bead);
        current_bead = m->nextBead(current_bead);
        if (current_bead == -1)
        {
            // -1 = no connection
            break;
        }
    }
    current_bead=index;
    
    for (int j=1; j<shift; j++)
    {
        current_bead = m->prevBead(current_bead);
        if (current_bead == -1)
        {
            // -1 = no connection
            break;
        }
        available_beads.erase(current_bead);
    }
}

void Histone::updateAfterDetachment(int index)
{
    int current_bead = index;
    // get the bead after the end of this loop
    int comp_bead = m->nextBead(partner[index]);
    int shift = beads_in_histone_loop;
    
    for (int j=0; j<shift; j++)
    {
        if (comp_bead == -1)
        {
            break;
        }
        beads_not_in_loop.insert(current_bead);
        // for these ones, need to check if there's a loop further down the line blocking these
        if (bondable[current_bead] && beads_not_in_loop.count(comp_bead))
        {
            available_beads.erase(current_bead);
        }
        current_bead = m->nextBead(current_bead);
        comp_bead = m->nextBead(comp_bead);
    }
    current_bead=index;
    
    for (int j=1; j<shift; j++)
    {
        current_bead = m->prevBead(current_bead);
        if (current_bead == -1)
        {
            break;
        }
        if (bondable[current_bead])
        {
            available_beads.insert(current_bead);
        }
    }
}

// Update histone timers.
void Histone::updateTime(const double dt) {
    // if enabled, do tension-based unattach
    if (tension_unbind)
    {
        unbindUpdate();
    }
    
    // flag indicates if there are new histones attach to arm.
    bool flag = 0;
    
    for (int i=0;i<number_of_histones;i++) {
        // If unattach is turned on, histones may unattach the chromosome arm if the histone is streched long enough.
        histone_config[i].remaining_time -= dt;
        histone_config[i].time_started += dt;
        if (histone_config[i].remaining_time<0 && histone_config[i].start_bead_index<0)
        {
            flag = 1;
        }
        else if (histone_config[i].remaining_time<0 && histone_config[i].start_bead_index>=0)
        {
            double off_time = dist_off(*rng_p);
            histoneDetach(i, off_time);
        }
    }
    if (!flag) return;
    
    for (int i=0;i<number_of_histones;i++)
    {
        if (histone_config[i].remaining_time<0 && histone_config[i].start_bead_index<0)
        {
            //New loops are forming. In loopConfig, assign a new starting bead to the first entry of the histone, and assign its on time on the second entry. Update availableBeads and beadsNotInLoop.
            if (available_beads.empty())
            {
                break;
            }
            double remaining_time = dist_on(*rng_p);
            histoneInsert(i, remaining_time);
        }
    }
}

void Histone::updateForce()
{
    double norm,scale;
    for (int i=0; i<number_of_histones; i++)
    {
        int loop_start_bead = histone_config[i].start_bead_index;
        if (loop_start_bead<0)
        {
            // this indicates histone is not attached
            continue;
        }
        
        int loop_end_bead = partner[loop_start_bead];
        
        norm = 0;
        // stores dx,dy,dz
        double diff[3];
        for (int k=0;k<3;++k)
        {
            diff[k] = m->position(k, loop_end_bead) - m->position(k, loop_start_bead);
            norm += diff[k]*diff[k];
        }
        
        norm = sqrt(norm);
        
        double subscale = spring_constant*(1-m->spring_length/norm);
        for (int k=0;k<3;++k)
        {
            scale = subscale*diff[k];
            m->force_array(k,loop_start_bead) += scale;
            m->force_array(k,loop_end_bead) -= scale;
        }
    }
}



// This function helps to detach a dying histone. The timer is reset, and availableBeads & beadsNotInLoop variables are updated.
void Histone::histoneDetach(const int index, const double detach_time)
{
    // Looped but is breaking up. In loopConfig, assign a negative value to the first entry of the histone, and assign its off time on the second entry. Update availableBeads and beadsNotInLoop.
    int removed_index = histone_config[index].start_bead_index;
    histone_config[index].start_bead_index -= number_of_beads;
    histone_config[index].remaining_time = detach_time;
    // Since one histone pops off, we add beads on this histone to beadsNotInLoop, and we add beads within range of 6 from removed bead into availableBeads.
    updateAfterDetachment(removed_index);
}

// This function helps to attach a newly activated histone. The position is chosen from availableBeads to avoid overlapping. The timer is reset, and availableBeads & beadsNotInLoop variables are updated.
void Histone::histoneInsert(const int i, const double binding_time)
{
    int tmp_bead = getAvailableBead();
    
    histone_config[i].start_bead_index = tmp_bead;
    histone_config[i].time_started = 0;
    histone_config[i].remaining_time = binding_time;
    
    updateAfterAttachment(tmp_bead);
}


int Histone::getAvailableBead()
{
    if (random_placement)
    {
        return sampleFromSet(available_beads, rng_p);
    }
    else
    {
        // do we even want to keep this option?
        return *available_beads.begin();
    }
}


void Histone::unbindUpdate()
{
    for (int i=0;i<number_of_histones;++i)
    {
        if (histone_config[i].time_started >= unattach_start_time && histone_config[i].start_bead_index >= 0)
        {
            int start_bead = histone_config[i].start_bead_index;
            int start_bead_spring_index = m->prevBead(start_bead);
            
            double max1;
            if (start_bead_spring_index != -1)
            {
                // larger norm of the two springs connecting to head of histone
                max1 = max(m->norm_array(start_bead), m->norm_array(start_bead_spring_index));
            }
            else
            {
                max1 = m->norm_array(start_bead);
            }

            int end_bead = partner[start_bead];
            int end_bead_spring_index = m->prevBead(end_bead);
            double max2;
            if (end_bead_spring_index != -1)
            {
                max2 = max(m->norm_array(end_bead), m->norm_array(end_bead_spring_index));
            }
            else
            {
                max2 = m->norm_array(end_bead);
            }
            double dist = max(max1, max2);
            if (dist > unattach_threshold) {
                histoneDetach(i, unattach_time);
            }
        }
    }
}
