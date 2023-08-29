//
//
//  ChromoShake
//
//  Created by HeYunyan on 8/30/17.
//  Updated by Ben Walker on 3/1/21
//


#include <stdio.h>
#include <random>
#include <algorithm>
#include <deque>
#include <unordered_set>
#include <cassert>  // DK make compiler agreeable

#include "Model.hpp"
#include "unitConversion.h"


// TODO: check use of these functions (and where they come from)
double findMax(double *u, int start, int end);
double acos(double x);
double innerProduct(double *u, double *v);
void computeCylinderDistance(double x11,double x12, double x13, double x21, double x22, double x23, double y11, double y12, double y13, double y21, double y22, double y23,double *toReturn);
void computePointLineDistance(double x11,double x12, double x13, double y11, double y12, double y13, double y21, double y22, double y23, double *toReturn);



int sampleFromSet(const std::unordered_set<int> &s, std::mt19937_64 *rng_p)
{
    uniform_int_distribution<int> dist(0, s.size()-1);
    int bead_shift = dist(*rng_p);
    auto iter = s.begin();
    std::advance(iter, bead_shift);
    return *iter;
}

// debug
void printBeadsInTree(CollisionTreeNode *p)
{
    if (p->left)
    {
        printBeadsInTree(p->left);
    }
    if (p->right)
    {
        printBeadsInTree(p->right);
    }
    if (p->bead_id>=0)
    {
        printf("%d ", p->bead_id);
    }
}


Model::Model(double t, double timestep,
             const DTDictionary &coefficients,
             const DTDictionary &evolution,
             const DTDictionary &flags,
             const DTDictionary &cohesin,
             const double seed,
             const DTDoubleArray &chains, const DTDoubleArray &anchors,
             const DTDoubleArray &beads_to_displace,
             int _number_of_beads,
             mt19937_64 &_rng) :  rng(_rng), prev_bead(_number_of_beads), next_bead(_number_of_beads)
{
    //Read in parameters
    double node_separation_nanometers = coefficients("mass_separation_nanometers");
    double temperature_Celsius = coefficients("temperature_Celsius");
    double viscosity_centiPoise = coefficients("viscosity_centiPoise");
    double modulus_gigaPascal = coefficients("modulus_gigaPascal");
    double modulus_nanometers = gigaPascal_to_Pascal(modulus_gigaPascal) * 1e-9;
    
    double const_dna_modulus_gigaPascal = coefficients("const_dna_modulus_gigaPascal");
    double const_dna_modulus_nanometers = gigaPascal_to_Pascal(const_dna_modulus_gigaPascal) * 1e-9;
    
    // ??
    double dna_radius_nanometers = coefficients("dna_radius_nanometers");
    // factor by which the radius of the squishy sphere is computed relative to spring length
    // radius of sphere = 1/2 * factor * spring_length
    double collision_radius_factor = coefficients("collision_radius_factor");
    // factor to scale down a spring constant
    double collision_force_factor = coefficients("collision_force_factor");
    // ??
    double damping_radius_factor = coefficients("damping_radius_factor");
    // factor to scale down hinge force making cohesin rings floppier
    double cohesin_factor = cohesin("cohesin_factor");
    
    
    //Unit conversion
    //double histone_radius_meters = microns_to_meters(histone_radius_microns);
    double temperature_kelvin = Celsius_to_Kelvin(temperature_Celsius);
    double viscosity_Pascal_seconds_nanometers = centiPoise_to_Pascal_seconds(viscosity_centiPoise) * 1e-9;
    double effective_damping_radius = damping_radius_factor * node_separation_nanometers;
    
    drag_coeff = mass_damping_equivalent(viscosity_Pascal_seconds_nanometers,
                                         effective_damping_radius);
    double spring_const = dna_tensile_spring_constant(
                                                      node_separation_nanometers,
                                                      const_dna_modulus_nanometers,
                                                      dna_radius_nanometers
                                                      );
    spring_const /= drag_coeff;
    base_spring_constant = spring_const / const_dna_modulus_gigaPascal;
    
    double hinge_const = dna_bending_spring_constant(
                                                     node_separation_nanometers,
                                                     modulus_nanometers,
                                                     dna_radius_nanometers
                                                     );
    hinge_const /= drag_coeff;
    random_scale = 1e9 * Brownian_force_equivalent(
                                                   node_separation_nanometers,
                                                   temperature_kelvin,
                                                   viscosity_Pascal_seconds_nanometers
                                                   );
    random_scale /= drag_coeff;
    
    
    //  default=  0.9 x 10/2
    double mass_radius = collision_radius_factor * (node_separation_nanometers / 2.0);
    double collision_spring_constant = spring_const * collision_force_factor;
    spring_length = node_separation_nanometers;
    collision_radius = mass_radius;
    collision_scale = collision_spring_constant;
    
    // Save bead indices that are to be fixed
    anchoring = flags("anchor");
    anchors_array = anchors.Copy();
    setAnchors();
    
    // Weighted beads are always anchored beads
    weighted = flags("weighted");
    weight_factor = coefficients("weight_factor");
    
    displacement = flags ("displacement");
    displaced_beads_array = beads_to_displace.Copy();
    setDisplaceBeads();
    
    
    // set up random distributions we will use
    dist_over_beads = uniform_int_distribution<int>(0, number_of_beads-1);
    dist_01 = uniform_int_distribution<int>(0, 1);
    
    // Initialize geometry and assign mass
    number_of_beads = _number_of_beads;
    //handling of beads according to chains_array
    DTMutableDoubleArray ini_geo(3,number_of_beads),
    ini_force_array(3,number_of_beads),
    ini_dist_array(3,number_of_beads),
    ini_norm_array(number_of_beads);
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////     Initialize geometry of bead chains    //////
    ////////////////////////////////////////////////////////////////////////////
    ////////  introduced:   chains_array(first, last, loop_bool, length_nanometers, x_initial, y_initial, z_initial, DNA_bool)
    ////////                       anchored_beads_array
    
    DTMutableDoubleArray tmp_spring_table(2, number_of_beads);
    
    // Initial configuration as a loop (circle) is indicated by loop_bool which is second value of the chainsArray
    // referring to array: m first, n second, o third dimentions
    
    chains_array = chains.Copy();
    int num_chains = chains_array.n();
    for (int k=0; k<num_chains; k++)
    {
        // check that we aren't exceeding the total number of beads
        if (chains_array(1,k) >= number_of_beads)
        {
            DTErrorMessage("Chains array bead indices exceed total bead number");
            exit(1);
        }
        
        const int beads_in_cur_chain = chains_array(1,k) - chains_array(0,k) + 1;
        
        // initialize the vectors for previous and next beads
        const int chain_start = chains_array(0,k);
        const int chain_end = chains_array(1,k);
        const bool chain_loop = chains_array(2,k);
        for (int i=chain_start; i<=chain_end; ++i)
        {
            
            // set up spring constants for all beads
            // TODO: remember hinge_constant is scaled by cohesin_factor for chains that are "not DNA" regardless of cohesin flag
            // (this feature is going to affect linked chains with heterogeneous properties)
            
            // set up spring constants for all beads
            tmp_spring_table(0,i) = spring_const;
            if (chains_array(7,k)!=1)
            {
                // not cohesin
                tmp_spring_table(1,i) = hinge_const;
            }
            else if (chains_array(7,k)==1)
            {
                // cohesin - make it floppy by reducing hinge forces
                // const double cohesin_factor = 0.1;
                
                // 2022.09.28
                // IMAGE TANK, parameter. Set cohesin factor to 0 in order to make cohesins floppy (Lp5)
                // TODO: will need to rewrite for any other Lp for cohesin rings
                tmp_spring_table(1,i) = cohesin_factor * hinge_const;
            }
            
            // compute indices of previous and next beads accounting for chain geometry
            prev_bead[i] = i-1;
            next_bead[i] = i+1;
            
            if (i==chain_start && chain_loop)
            {   // if plasmid then index ends in a circle
                prev_bead[i] = chain_end;
            }
            else if (i==chain_start && !chain_loop)
            {    // -1 = first bead doesn't have prev in a line geometry
                prev_bead[i] = -1;
            }
            
            if (i==chain_end && chain_loop)
            {   // if plasmid then index ends in a circle
                next_bead[i] = chain_start;
            }
            else if (i==chain_end && !chain_loop)
            {   // -1 = last bead doesn't have next in a line geometry
                next_bead[i] = -1;
            }
        }
        
        
        // if plasmid or loop, initialize beads in a circle along the circumference which equals length_nanometers of DNA
        if (chain_loop)
        {
            double angle = 2 * M_PI / beads_in_cur_chain;
            
            for (int i=chains_array(0,k); i<=chains_array(1,k); ++i)
            {
                int current_bead = i - chains_array(0,k);
                // code to initialize a single plasmid loop ...
                if (k==0)
                {
                    // x-coordinate with offset stored in chainsArray
                    ini_geo(0,i) = chains_array(4,k) + sin(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                    // y-coordinate with offset stored in chainsArray
                    ini_geo(1,i) = chains_array(5,k) + cos(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                    // z-coordinate with offset stored in chainsArray
                    ini_geo(2,i) = chains_array(6,k);
                }
                else if (num_chains > 1 && chains_array(2,0))
                    // ... generate cohesin rings around plasmid;   generate a ring in along the x-axis then rotate accordingly ...
                {
                    double segment = 2 * M_PI / (num_chains - 1);
                    // x-plane
                    double temp_x = chains_array(4,k) + sin(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                    // y-coordinate
                    double temp_y = chains_array(5,k);
                    // z-coordinate
                    double temp_z = chains_array(6,k) + cos(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                    
                    ini_geo(0,i) = temp_x * cos(segment*k) - temp_y * sin(segment*k);
                    ini_geo(1,i) = temp_x * sin(segment*k) + temp_y * cos(segment*k);
                    ini_geo(2,i) = temp_z;
                }
                else
                    // ... generate cohesin rings around linear chain
                {
                    // x-coordinate with offset stored in chainsArray
                    ini_geo(0,i) = chains_array(4,k);
                    // y-coordinate
                    ini_geo(1,i) = chains_array(5,k) + sin(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                    // z-coordinate
                    ini_geo(2,i) = chains_array(6,k) + cos(angle*current_bead) * chains_array(3,k) / (2 * M_PI);
                }
                
                for (int j=0; j<3; j++) {
                    ini_force_array(j,i) = 0;
                    ini_dist_array(j,i) = 0;
                }
                ini_norm_array(i) = 0;
            }
        }
        else
            // initialize beads in a line corresponding to length of DNA in nanometers
        {
            for (int i=chains_array(0,k); i<=chains_array(1,k); ++i)
            {
                // x-coord beads are evenly spaced on the chain/DNA
                ini_geo(0,i) = chains_array(4,k) + (i-chains_array(0,k))*chains_array(3,k)/beads_in_cur_chain;
                // y-coord with offset
                ini_geo(1,i) = chains_array(5,k);
                // z-coord with offset
                ini_geo(2,i) = chains_array(6,k);
                
                for (int j=0; j<3; j++)
                {
                    ini_force_array(j,i) = 0;
                    ini_dist_array(j,i) = 0;
                }
                ini_norm_array(i) = 0;
            }
        }
    }
    
    position = ini_geo.Copy();
    force_array = ini_force_array.Copy();
    distance_array = ini_dist_array.Copy();
    norm_array = ini_norm_array.Copy();
    
    // TODO: code to displace beads in x-axis over time of simuation
    double time_finish = evolution("until");
    double time_step = evolution("timeStep");  // not used
    double final_bead_distance = coefficients("final_bead_distance");
    //num_of_step not use since switch to adaptive time stepping dependent on dt
    num_of_steps = time_finish/time_step;
    
    end_distance = final_bead_distance;
    bead_a = displaced_beads_array(0);
    bead_b = displaced_beads_array(1);
    double start_distance = fabs(position(0,bead_b) - position(0,bead_a));
    double path_distance = (start_distance - end_distance)/2.0;
    displacement_rate = path_distance/time_finish;  // multiplied by dt get the distance to move
    
    
    position_pointer = position.Pointer();
    force_array_pointer = force_array.Pointer();
    distance_array_pointer = distance_array.Pointer();
    norm_array_pointer = norm_array.Pointer();
    
    spring_properties = tmp_spring_table.Copy();
    spring_properties_pointer = spring_properties.Pointer();
    
    
    // debug check of collision tree
    n2_check = false;
    constructTree();
    updateTree();
    
    cr2 = 2*collision_radius;
} // end of constructor


Model::~Model()
{
    for (auto p : modules)
    {
        delete p;
    }
    free(collision_tree_storage);
}

void Model::addModule(Module *m) {
    modules.push_back(m);
}


void Model::constructTree()
{
    // pre-allocate space in collision tree storage
    // TODO: don't waste memory here
    const int cts_alloc_size = 3*number_of_beads;
    collision_tree_storage = (CollisionTreeNode*) malloc(cts_alloc_size*sizeof(CollisionTreeNode));
    // Go through each chain and construct a collision tree for that chain
    int num_chains = chains_array.n();
    tree_pointer_array.resize(num_chains);
    
    //typedef std::vector<CollisionTreeNode>::iterator it_type;
    int node_counter = cts_alloc_size-1;
    for (int k=0; k<num_chains; k++)
    {
        int beads_in_cur_chain = chains_array(1,k) - chains_array(0, k) + 1;
        // temp storage for nodes that we're putting into the tree
        std::vector<CollisionTreeNode*> tmp_a, tmp_b;
        tmp_a.reserve(beads_in_cur_chain);
        tmp_b.reserve(beads_in_cur_chain);
        
        // create a leaf node for each bead in this chain
        for (int i=chains_array(0, k); i <= chains_array(1,k); ++i)
        {
            // create the new node
            auto node = collision_tree_storage+(node_counter--);
            node->left = NULL;
            node->right = NULL;
            node->bead_id = i;
            // accepts the radius
            node->radius = collision_radius;
            tmp_a.push_back(node);
        }
        // build the first layer, getting number of beads to a power of 2
        // first work out how many beads to pair
        int num_nodes = tmp_a.size();
        
        // right now assuming there is actually a tree to make
        assert(num_nodes > 1);
        
        int pairs_to_make = (num_nodes-pow(2,floor(log2(double(num_nodes)))));
        for (int i=0;i<2*pairs_to_make;i+=2)
        {
            // create node combining the two children
            auto node = collision_tree_storage+(node_counter--);
            node->left = tmp_a[i];
            node->right = tmp_a[i+1];
            node->bead_id = -1;
            // start building up second list
            tmp_b.push_back(node);
        }
        // copy the rest of the nodes over as is
        for (int i=2*pairs_to_make;i<tmp_a.size();++i)
        {
            // TODO: more efficient copy
            tmp_b.push_back(tmp_a[i]);
        }
        tmp_a.swap(tmp_b);
        tmp_b.clear();
        
        // now iteratively build down to a single node
        while (tmp_a.size() > 1)
        {
            assert(tmp_a.size() % 2 == 0);
            for (int i=0;i<tmp_a.size(); i+=2)
            {
                // TODO: create inline helper function for this?
                auto node = collision_tree_storage+(node_counter--);
                node->left = tmp_a[i];
                node->right = tmp_a[i+1];
                node->bead_id = -1;
                // start building up second list
                tmp_b.push_back(node);
            }
            tmp_a.swap(tmp_b);
            tmp_b.clear();
        }
        // now we've got the single head node
        auto v = tmp_a.front();
        tree_pointer_array[k] = v;
    }
} // end of constructTree()


DTPath3D Model::convertPointsToPath() const
{
    int n = position.n();
    int num_chains = chains_array.n();
    DTMutableDoubleArray loop(3,n+num_chains);
    
    //MemoryCopyColumns(loop,1,points,DTRange(0,n));
    
    // current offset in path array
    int c = 0;
    for (int k=0; k<num_chains; ++k)
    {
        int beads_in_chain = chains_array(1,k)-chains_array(0,k)+1;
        loop(0,c) = 0;
        loop(1,c) = 0;
        loop(2,c) = beads_in_chain;
        ++c;
        auto r = DTRange(chains_array(0,k),beads_in_chain);
        MemoryCopyColumns(loop,c,position,r);
        c += beads_in_chain;
    }
    
    auto path = DTPath3D(loop);
    return path;
}


//Compute distance vector and its norm between two adjacent beads. Prepare to use in other functions.
// temporarily stores distance vector and its norm between adjacent beads.
//  - get distances or paired beads.
//  - keep in mind the loop (first + last distance)
void Model::computeDistance()
{
    double scale;
    int num_chains = int(chains_array.n());
    
    for (int k=0; k<num_chains; k++)
    {   //  k = chain index; i = bead index; j = coordinate index (x,y,z)
        int until = chains_array(1,k);
        for (int i = chains_array(0,k); i<until; ++i)
        {
            norm_array(i) = 0;
            for (int j=0; j<3; j++)
            {
                scale = position_pointer[(i+1)*3+j] - position_pointer[i*3+j];
                distance_array_pointer[j+3*i] = scale; // distance_array(j,i) = scale;
                norm_array(i) += scale*scale;
            }
            // norm of distances
            norm_array(i) = sqrt(norm_array(i));
        }
        
        // TODO:     YourTerminateNowCall();
        //           exit(-1); - terminates now!
        // Handle the last bead separately
        int first_bead_in_chain = chains_array(0,k);
        int last_bead_in_chain = chains_array(1,k);
        
        // computes norm for a last bead in the chain only when chain is a circle
        if (chains_array(2,0))
        {
            norm_array(last_bead_in_chain) = 0;
            for (int j=0; j<3; j++)
            {
                scale = position_pointer[first_bead_in_chain*3+j] - position_pointer[last_bead_in_chain*3+j];
                distance_array(j,last_bead_in_chain) = scale;
                norm_array(last_bead_in_chain) += scale*scale;
            }
            norm_array(last_bead_in_chain) = sqrt(norm_array(last_bead_in_chain));
        }
        else
        {
            continue;
        }
    }
} // end of computeDistance()

// Compute the spring force and update forceArray.
void Model::springForceUpdate()
{
    double scale;
    int num_chains = int(chains_array.n());
    //  k = chain index; i = bead index; j = coordinate index (x,y,z)
    for (int k=0; k<num_chains; k++)
    {
        int until = chains_array(1,k);
        for (int i=chains_array(0,k); i<until; i++)
        {
            for (int j=0; j<3; j++)
            {
                // scale = spring_properties_pointer[2*i] * distance_array(j,i) * (1 - spring_length/norm_array(i));
                scale = spring_properties_pointer[2*i] * distance_array_pointer[j+3*i] * (1 - spring_length/norm_array_pointer[i]);
                // write the force (and equal opposite force) to beads i and i+1 respectively
                force_array_pointer[j+3*i] += scale; // force_array(j,i) += scale;
                force_array_pointer[j+3*(i+1)] -= scale; // force_array(j,i+1) -= scale;
            }
        }
        
        // if in a loop, add the spring from the end back to beginning
        if (chains_array(2,k))
        {
            int first_bead_in_chain = chains_array(0,k);
            int last_bead_in_chain = chains_array(1,k);
            double dist[3];
            double norm = 0;
            for (int j=0; j<3; ++j)
            {
                dist[j] = position(j, first_bead_in_chain) - position(j, last_bead_in_chain);
                norm += dist[j]*dist[j];
            }
            norm = sqrt(norm);
            
            for (int j=0; j<3; j++)
            {
                scale = spring_properties_pointer[2*last_bead_in_chain] * dist[j] *
                (1 - spring_length/norm);
                force_array(j,last_bead_in_chain) += scale;
                force_array(j,first_bead_in_chain) -= scale;
            }
        }
    }
} // end of springForceUpdate()

// Update the hinge force by computing the pairwise angle of the |chromosome arm| and calculate the force. Update forceArray.
// Forces a linear geometry on the set of three beads in a row.
void Model::hingeForceUpdate()
{
    //cross distance, forward distance, backward distance
    double dxyz1[3],dxyz2[3];
    double norm1,norm2,scale,phi;
    int num_chains = int(chains_array.n());
    // int m_chains = chains_array.m();
    
    for (int k=0; k<num_chains; k++)
    {
        int first_bead_in_chain = chains_array(0,k);
        int last_bead_in_chain = chains_array(1,k);
        
        //  k = chain index; i = bead index; j = coordinate index (x,y,z)
        int until = chains_array(1,k);
        for (int i=chains_array(0,k); i<=until; i++) {
            // if chain is not NOT a loop AND dealing with first OR last bead
            if ((i==first_bead_in_chain || i==last_bead_in_chain) && (!chains_array(2,k)))
            {
                // without the loop, no hinge force on first or last bead
                continue;
            }
            
            for (int j=0; j<3; j++)
            {
                if (i == first_bead_in_chain)
                {
                    dxyz1[j] = distance_array_pointer[3 * last_bead_in_chain + j];
                }
                else
                {
                    dxyz1[j] = distance_array_pointer[3 * (i - 1) + j];
                }
                
                dxyz2[j] = distance_array_pointer[3 * i + j];
            }
            
            if (i == first_bead_in_chain)
            {
                norm1 = norm_array_pointer[last_bead_in_chain];
            }
            else
            {
                norm1 = norm_array_pointer[i - 1];
            }
            
            norm2 = norm_array_pointer[i];
            // dot product of two vectors divided by the product of their norms
            scale = (dxyz1[0]*dxyz2[0] + dxyz1[1]*dxyz2[1] + dxyz1[2]*dxyz2[2]) / (norm1*norm2); // cos(phi)
            phi = acos(scale);              // approximation of arccos
            
            //now compute direction of hinge force
            dxyz1[0] = -dxyz1[0]/norm1 + dxyz2[0]/norm2;
            dxyz1[1] = -dxyz1[1]/norm1 + dxyz2[1]/norm2;
            dxyz1[2] = -dxyz1[2]/norm1 + dxyz2[2]/norm2;
            
            norm1  = sqrt(dxyz1[0]*dxyz1[0] + dxyz1[1]*dxyz1[1] + dxyz1[2]*dxyz1[2]);
            
            if (norm1 > 1e-10)
            {
                for (int j=0; j<3; j++)
                {
                    scale = spring_properties_pointer[2*i+1] * phi * dxyz1[j]/norm1;
                    force_array_pointer[j+3*i] += scale; // force_array(j,i) += scale;
                    
                    if (i == first_bead_in_chain)
                    {
                        force_array(j,last_bead_in_chain) -= scale/2;
                    }
                    else
                    {
                        force_array(j,i-1) -= scale/2;
                    }
                    
                    if (i == last_bead_in_chain)
                    {
                        force_array(j,first_bead_in_chain) -= scale/2;
                    }
                    else
                    {
                        force_array(j,i+1) -= scale/2;
                    }
                }
            }
        }
    }
} // end of hingeForceUpdate()


// (Thermal noise) Update random force by generating random numbers and update forceArray.
void Model::randomForceUpdate(double dt)
{
    DTMutableDoubleArray random_noise_work_array(3,number_of_beads);
    double scale = random_scale/sqrt(dt);
    normal_distribution<double> norm_dist(0, 1);
    
    double *random_noise_work_arrayD = random_noise_work_array.Pointer();
    for (int i=0;i<3*number_of_beads;++i)
    {
        random_noise_work_arrayD[i] = norm_dist(rng); // random_noise_work_array(i) = norm_dist(rng);
    }
    random_noise_work_array *= scale;
    
    force_array += random_noise_work_array;
}


// Collision force update. Compute pairwise distance and assign collision forces to each bead. Update forceArray.
void Model::excludedVolumeUpdate()
{
    
    double norm, scale;
    double dx,dy,dz;
    int i3 = 0;
    int j3;
    double compare_dist_squared = 4 * collision_radius*collision_radius;
    
    int c=0;
    std::vector<ipair> nearby;
    nearby.reserve(2*number_of_beads);
    updateTree();
    getIntersectList(nearby);
    for (auto p : nearby)
    {
        int i = p.first, j = p.second;
        i3 = 3*i; j3 = 3*j;
        dx = position_pointer[i3] - position_pointer[j3];
        dy = position_pointer[i3+1] - position_pointer[j3+1];
        dz = position_pointer[i3+2] - position_pointer[j3+2];
        norm = dx*dx + dy*dy + dz*dz;
        // note: technically nearby can contain things that aren't actually within the distance but they mostly are
        // TODO: that might actually no longer be true
        if (norm<compare_dist_squared) {
            ++c;
            norm = sqrt(norm);
            scale = collision_scale*(cr2-norm)/norm;
            force_array(0,i) += scale*dx;
            force_array(0,j) -= scale*dx;
            force_array(1,i) += scale*dy;
            force_array(1,j) -= scale*dy;
            force_array(2,i) += scale*dz;
            force_array(2,j) -= scale*dz;
        }
    }
} // end of excludedVolumeUpdate()



void Model::setAnchors()
{
    if (anchoring)
    {
        int num_anchors = anchors_array.m();
        for (int i=0; i<num_anchors; i++)
        {
            anchored_beads.insert(anchors_array(i));
        }
    }
} // end of setAnchors()

void Model::setDisplaceBeads()
{
    if (displacement)
    {
        int num_to_store = displaced_beads_array.m();
        for (int i=0; i<num_to_store; i++)
        {
            displaced_beads.insert(displaced_beads_array(i));
        }
    }
} // end of setDisplaceBeads()


// Right now only 2 beads move toward each other along their x-coordinates
void Model::positionUpdate(double dt)
{
    
    for (int i=0; i<number_of_beads; i++)
    {
        
        if (displaced_beads.count(i))
            // keep for posterity
        {
            // change the x-coordinate only
            position(0,bead_a) +=  dt*displacement_rate/2.0;
            position(0,bead_b) -=  dt*displacement_rate/2.0;
            //  y and z are explicitely not updated and stay at their starting positions.
        }
        
        else
        {
            double scale = 1;
            if (anchored_beads.count(i))
            {
                if (!weighted) scale = 0;
                else scale = 1/weight_factor;
            }
            
            for (int j=0;j<3;j++) position(j,i) +=  dt*scale*force_array_pointer[3*i+j];
        }
    }
    
    // Implementation of boundary conditions (boundary_flag) using DTMesh3DGrid.  Space is divided into 3D grid
    // and boundary is a STL shape or DT equation embedded in the mesh.  Inside shape values are positive
    // outside shape values are negative, shape boundary is zero.
    /*
     if (1) {
     for (int i=0; i<number_of_beads; i++)
     {
     double h = 1; // grid size
     DTPoint3D origin(0,0,0);
     double x = (position_pointer[0+3*i]-origin.x)/h;
     double y = (position_pointer[1+3*i]-origin.y)/h;
     double z = (position_pointer[2+3*i]-origin.z)/h;
     
     int xi = floor(x);
     double xfract = x-xi;
     int yi = floor(y);
     double yfract = y-yi;
     int zi = floor(z);
     double zfract = z-zi;
     DTFloatArray values;
     double outputValue = (
     values(xi,yi,zi)*(1-xfract)*(1-yfract)*(1-zfract) +
     values(xi+1,yi,zi)*xfract*(1-yfract)*(1-zfract) +
     values(xi,yi+1,zi)*(1-xfract)*yfract*(1-zfract) +
     values(xi+1,yi+1,zi)*xfract*yfract*(1-zfract) +
     values(xi,yi,zi+1)*(1-xfract)*(1-yfract)*zfract +
     values(xi+1,yi,zi+1)*xfract*(1-yfract)*zfract +
     values(xi,yi+1,zi+1)*(1-xfract)*yfract*zfract +
     values(xi+1,yi+1,zi+1)*xfract*yfract*zfract
     );
     if (outputValue>0) {
     // Outside, need to move this point back.
     // previous = position_pointer[j+3*i] -  dt*force_array_pointer[3*i+j];
     }
     
     }
     
     }
     */
} // end of positionUpdate()



double Model::findMaxForce() const
{
    return *std::max_element(force_array_pointer, force_array_pointer+3*number_of_beads);
}

void Model::moduleTimeUpdate(const double dt)
{
    for (auto p : modules)
    {
        p->updateTime(dt);
    }
}

void Model::moduleForceUpdate()
{
    for (auto p : modules)
    {
        p->updateForce();
    }
}


// This function is triggered at the beginning of every time step. Clear forceArray, distanceArray and normArray.
void Model::initializeForce()
{
    force_array = 0;
    distance_array = 0;
    norm_array = 0;
}


double findMax(double *u, int start, int end)
{
    if (end-start<2) return fabs(u[end])>fabs(u[start])?fabs(u[end]):fabs(u[start]);
    double m1 = findMax(u,start,floor(end/2));
    double m2 = findMax(u,floor(end/2)+1,end);
    return m1>m2?m1:m2;
}


// alt version may use bb_tree
int Model::findClosestBead(double x,double y, double z, int excluded_bead, int itself, std::vector<bool> bondable) const
{
    double min_dist = __DBL_MAX__, dx, dy, dz, min_dist_temp;
    int i3;
    int ret = -1;
    
    for (int i=0; i<number_of_beads; i++)
    {
        i3 = 3*i;
        if (i == excluded_bead || i == itself || !bondable[i])
        {
            continue;
        }
        else
        {
            // check if correct
            dx = x-position_pointer[i3];
            dx = dx*dx;
            if (min_dist<dx) continue;
            dy = y-position_pointer[i3+1];
            dy = dy*dy;
            if (min_dist<dy) continue;
            dz = z-position_pointer[i3+2];
            dz = dz*dz;
            if (min_dist<dz) continue;
            min_dist_temp = dx+dy+dz;
            if (min_dist_temp<min_dist)
            {
                min_dist = min_dist_temp;
                ret = i;
            }
        }
    }
    return ret;
} // end of findClosestBead()

// Update all collision trees based on current positions
void Model::updateTree()
{
    for (auto p : tree_pointer_array)
    {
        updateTreeSub(p);
    }
}

inline void boundingSphereUpdate(CollisionTreeNode *node);
// Update the collision tree given by the specified head pointer
void Model::updateTreeSub(CollisionTreeNode* head)
{
    // TODO: iterative implementation?
    
    // check if we're a leaf
    if (head->bead_id >= 0)
    {
        // copy position from data
        for (int k=0;k<3;++k)
        {
            head->pos[k] = position_pointer[k + 3*(head->bead_id)]; // head->pos[k] = position(k, head->bead_id);
        }
    }
    else
    {
        // all non-leaves should have two children
        assert(head->left && head->right);
        
        // recursively update children
        updateTreeSub(head->left);
        updateTreeSub(head->right);
        boundingSphereUpdate(head);
    }
}

// helper function to do the math for updating bounding sphere assuming children are updated
inline void boundingSphereUpdate(CollisionTreeNode *node)
{
    CollisionTreeNode *c1 = node->left, *c2 = node->right;
    // two inner nodes - spheres have arbitrary radius
    double r1 = c1->radius;
    double r2 = c2->radius;
    
    // swap if necessary so larger sphere comes first
    if (r2 > r1)
    {
        std::swap(r1, r2);
        std::swap(c1, c2);
    }
    double *a = c1->pos, *b = c2->pos;
    
    double sep=0;
    double diff[3];
    for (int i=0;i<3;++i)
    {
        double tmp = b[i] - a[i];
        diff[i] = tmp;
        sep += tmp*tmp;
    }
    sep = sqrt(sep);
    // accounts for whether sphere 1 entirely encompasses sphere 2
    double v = std::max(r1, r2+sep);
    
    // new center
    double scale = (v-r1)/(2*sep);
    for (int i=0;i<3;++i)
    {
        node->pos[i] = a[i]+scale*diff[i];
    }
    // new radius
    node->radius = (r1 + v)/2;
}

bool Model::sphereIntersect(CollisionTreeNode* n1, CollisionTreeNode *n2) const
{
    const double *a = n1->pos;
    const double *b = n2->pos;
    
    double comp = n1->radius + n2->radius;
    double comp2 = comp*comp;
    
    double dx = b[0]-a[0];
    double dy = b[1]-a[1];
    double dz = b[2]-a[2];
    double sep = dx*dx+dy*dy+dz*dz;
    
    return (sep<=comp2);
}

// Returns the separation (amount of empty space) between the spheres - this is a lower bound on the distance between points in each bounding sphere
double Model::sphereSeparation(CollisionTreeNode* n1, CollisionTreeNode *n2) const
{
    double *a = n1->pos;
    double *b = n2->pos;
    double sep = 0;
    
    double tmp;
    double comp = n1->radius + n2->radius;
    //double comp2 = comp*comp;
    for (int i=0;i<3;++i)
    {
        tmp = b[i]-a[i];
        sep += tmp*tmp;
    }
    sep = sqrt(sep);
    
    sep -= comp;
    
    // if sep < 0 the spheres intersect
    if (sep < 0)
    {
        sep = 0;
    }
    return sep;
}

inline void Model::conditional_push(CollisionTreeNode* n1, CollisionTreeNode *n2, std::deque<std::pair<CollisionTreeNode*,CollisionTreeNode*> > &s) const
{
    if (sphereIntersect(n1, n2))
    {
        s.push_back(std::make_pair(n1, n2));
    }
}

static inline void push_single(CollisionTreeNode* node, std::deque<std::pair<CollisionTreeNode*,CollisionTreeNode*> > &s)
{
    // don't push leaf nodes onto queue
    if (node->bead_id==-1)
    {
        s.push_back(std::make_pair(node,(CollisionTreeNode*) NULL));
    }
}

void Model::getIntersectList(vector<ipair> &neighbors)
{
    neighbors.clear();
    
    NodePairDeque queue;
    
    // first figure out which pairs of trees intersect overall
    for (int i=0;i<tree_pointer_array.size();++i)
    {
        auto n1 = tree_pointer_array[i];
        push_single(n1, queue);
        for (int j=i+1;j<tree_pointer_array.size();++j)
        {
            auto n2 = tree_pointer_array[j];
            conditional_push(n1, n2, queue);
        }
    }
    
    std::sort(queue.begin(), queue.end());
    
    // process queue
    while (queue.size() > 0)
    {
        auto cur_pair = queue.front();
        queue.pop_front();
        if (!cur_pair.second)
        {
            // second is null -> just a single node
            auto cur = cur_pair.first;
            
            // nothing to do for leaves on single stack
            if (cur->bead_id >= 0)
            {
                // single beads shouldn't have gotten pushed
                assert(false);
            }
            
            // both children go on the single stack
            push_single(cur->left, queue);
            push_single(cur->right, queue);
            
            // if they intersect, they go on the pair
            conditional_push(cur->left, cur->right, queue);
        }
        else
        {
            auto n1 = cur_pair.first, n2 = cur_pair.second;
            // check if both are leaves
            if (n1->bead_id >= 0 && n2->bead_id >= 0)
            {
                if (sphereIntersect(n1, n2))
                {
                    neighbors.push_back(std::make_pair(n1->bead_id, n2->bead_id));
                    continue;
                }
            }
            if (n1->left && n1->right && n2->left && n2->right)
            {
                conditional_push(n1->left, n2->left, queue);
                conditional_push(n1->left, n2->right, queue);
                conditional_push(n1->right, n2->left, queue);
                conditional_push(n1->right, n2->right, queue);
            }
            else
            {
                if (n1->left)
                {
                    conditional_push(n1->left, n2, queue);
                }
                if (n1->right)
                {
                    conditional_push(n1->right, n2, queue);
                }
                if (n2->left)
                {
                    conditional_push(n1, n2->left, queue);
                }
                if (n2->right)
                {
                    conditional_push(n1, n2->right, queue);
                }
            }
        }
    }
    
    if (n2_check)
    {
        double dx, dy, dz, norm;
        int i3, j3;
        double compare_dist_squared = collision_radius*collision_radius*4;
        double compare_dist = 2 * collision_radius;
        int c = 0;
        for (int i=0;i<number_of_beads - 1;i++) {
            i3 = 3*i;
            
            for (int j=i+1;j<number_of_beads;j++) {
                j3 = j*3;
                dx = position_pointer[i3] - position_pointer[j3];
                // checking these conditions as we go allows us to save some time if it's clearly too far even without checking all coordinates
                if (dx > compare_dist || dx < -compare_dist)
                {
                    continue;
                }
                dy = position_pointer[i3+1] - position_pointer[j3+1];
                if (dy > compare_dist || dy < -compare_dist)
                {
                    continue;
                }
                dz = position_pointer[i3+2] - position_pointer[j3+2];
                if (dz > compare_dist || dz < -compare_dist)
                {
                    continue;
                }
                norm = dx*dx + dy*dy + dz*dz;
                if (norm<compare_dist_squared) {
                    ++c;
                }
            }
        }
        
        assert(c == neighbors.size());
    }
}


// returns a list of all beads intersecting with a given node (excluding its descendents)
void Model::getIntersectList(CollisionTreeNode* node, std::vector<int> &ret) const
{
    ret.clear();
    
    // TODO: check if there is any way to make this have common code with full collision listing
    // go through all active tree pointers and start looking for matches
    
    // right now going to use existing pair-based code but have the second element always be node
    NodePairDeque queue;
    for (auto p : tree_pointer_array)
    {
        // exclude node itself
        if (p == node)
        {
            continue;
        }
        conditional_push(p, node, queue);
    }
    // process queue
    while (queue.size() > 0)
    {
        // node in question is always first element of pair (second is node from argument)
        auto cur = queue.front().first;
        queue.pop_front();
        
        // this definitely shouldn't be a leaf node
        assert(cur->left && cur->right);
        
        // check each child
        if (sphereIntersect(cur->left, node))
        {
            if (cur->left->bead_id >= 0)
            {
                // leaf node intersects, push on
                ret.push_back(cur->left->bead_id);
                continue;
            }
            else if (cur->left != node)
            {
                // not leaf, push onto queue
                queue.push_back(std::make_pair(cur->left, node));
            }
        }
        
        if (sphereIntersect(cur->right, node))
        {
            if (cur->right->bead_id >= 0)
            {
                // leaf node intersects, push on
                ret.push_back(cur->right->bead_id);
                continue;
            }
            else if (cur->right != node)
            {
                queue.push_back(std::make_pair(cur->right, node));
            }
        }
    }
}


// Check for intersection between two chains by looking for whether any beads of chain b are inside the bounding sphere of chain a
bool Model::chainIntersect(int a, int b) const
{
    auto node = tree_pointer_array[a], other = tree_pointer_array[b];
    
    // recursively look for intersections
    NodePairDeque queue;
    conditional_push(other, node, queue);
    
    // process queue
    while (queue.size() > 0)
    {
        // node in question is always first element of pair
        auto cur = queue.front().first;
        queue.pop_front();
        
        // this definitely shouldn't be a leaf node
        assert(cur->left && cur->right);
        
        // check each child
        if (sphereIntersect(cur->left, node))
        {
            if (cur->left->bead_id >= 0)
            {
                // leaf node intersects
                return true;
            }
            else if (cur->left != node)
            {
                // not leaf, push onto queue
                queue.push_back(std::make_pair(cur->left, node));
            }
        }
        
        if (sphereIntersect(cur->right, node))
        {
            if (cur->right->bead_id >= 0)
            {
                // leaf node intersects
                return true;
            }
            else if (cur->right != node)
            {
                queue.push_back(std::make_pair(cur->right, node));
            }
        }
    }
    
    return false;
}

double Model::bead_dist(int i, int j) const
{
    double dist = 0;
    for (int k=0;k>3;++k)
    {
        double tmp = position(k, i) - position(k, j);
        dist += tmp*tmp;
    }
    dist = sqrt(dist);
    return dist;
}


// returns true if there aren't any nan's
//bool Model::nan_check() const
//{
//    for (int i=0;i<3*number_of_beads;++i)
//    {
//        if (std::isnan(position_pointer[i]))
//        {
//            return false;
//        }
//    }
//    return true;
//}
