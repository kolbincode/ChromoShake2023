//
//  Model_h
//  ChromoShake
//
//  Created by HeYunyan on 8/29/17.
//  Updated by Ben Walker on 3/1/21
//

#ifndef Model_h
#define Model_h
#include <stdio.h>
#include <vector>
#include <random>
#include <deque>
#include <unordered_set>
#include <cassert>   //added to make compiler agreeable

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPath3D.h"
#include "DTProgress.h"
#include "DTSeriesPath3D.h"
#include "DTRandom.h"
#include "DTUtilities.h"
#include "DTPointCollection3D.h"

#include "unitConversion.h"


typedef std::pair<int, int> ipair;  //ipair


// TODO: figure out if there's a better file to put this function
int sampleFromSet(const std::unordered_set<int> &s, std::mt19937_64 *rng_p);


class Model;

class Module
//ask again about Module-Model relationship
{
public:
    virtual ~Module() {};
    virtual void updateTime(const double dt) {};
    virtual void updateForce() {};
};


// collision bounding-sphere tree is made of these
struct CollisionTreeNode
{
    //CollisionTreeNode *left, *right, *parent;
    CollisionTreeNode *left, *right;
    double pos[3];
    double radius;
    // if this is a leaf, bead_id is the bead it represents
    int bead_id;
    
    CollisionTreeNode(int id) { left = NULL; right = NULL; bead_id = id; };
    CollisionTreeNode() { left = NULL; right = NULL; bead_id = -1; };
};
typedef std::deque<std::pair<CollisionTreeNode*, CollisionTreeNode*> > NodePairDeque;
void printBeadsInTree(CollisionTreeNode *p);


class Model {
public:
    Model(double t, double timestep,
          const DTDictionary &coefficients,
          const DTDictionary &evolution,
          const DTDictionary &flags,
          const DTDictionary &cohesin,
          const double seed,
          const DTDoubleArray &chains, const DTDoubleArray &anchors,
          const DTDoubleArray &beads_to_displace,
          int _number_of_beads,
          mt19937_64 &rng);
    
    ~Model();
    
    // member functions
    int findClosestBead(double x, double y, double z, int excluded_bead, int itself, std::vector<bool> bondable) const;
    void initializeForce();
    void computeDistance();
    void springForceUpdate();
    void hingeForceUpdate();
    void randomForceUpdate(double dt);
    void dragForceUpdate();
    void excludedVolumeUpdate();
    void condensinForceUpdate();
    void positionUpdate(double dt);
    double findMaxForce() const;
    // handle the modules (extra things going on)
    void moduleTimeUpdate(const double dt);
    void moduleForceUpdate();

    double get_base_spring_constant() const { return base_spring_constant; };
    double cohesin_factor;
    
    DTPath3D convertPointsToPath() const;

// private:
    DTMutableDoubleArray chains_array;  // indices making up the DNA chains, and other bead-spring structures
    DTMutableDoubleArray anchors_array;       // indices of fixed beads   <== TODO: do we need this?
    DTMutableDoubleArray displaced_beads_array;       // indices of fixed beads   <== TODO: do we need this?
    DTMutableDoubleArray position;      // array contains 3d positions. (3,numberOfBeads)
    DTMutableDoubleArray force_array;   // array contains instant 3d force applied on each bead. (3,numberOfBeads)
    DTMutableDoubleArray distance_array; // array contains the pairwise 3D distance vector of beads. For speeding up. (3,numberOfBeads). [dist(bead1,bead0), dist(bead2,bead1),...,dist(bead0,beadLast)].
    DTMutableDoubleArray norm_array; // Norm array of distanceArray. (1,numberOfBeads).
    
    double spring_length;
    
// public:
    void addModule(Module *m);
    
    void getIntersectList(CollisionTreeNode* node, std::vector<int> &ret) const;
    CollisionTreeNode* getTreeNode(int i) const { return tree_pointer_array[i]; };
    // update the bounding spheres based on current positions
    void updateTree();
    // return whether chains a and b intersect based on bounding spheres
    bool chainIntersect(int a, int b) const; // { return sphereIntersect(tree_pointer_array[a], tree_pointer_array[b]);};
    double chainSeparation(int a, int b) const { return sphereSeparation(tree_pointer_array[a], tree_pointer_array[b]);};
    
    
    // Helper functions for geometry
    int nextBead(int bead) const { return next_bead[bead]; };
    int prevBead(int bead) const { return prev_bead[bead]; };
    
    // debug function to find issues that make nan's in the position
    bool nan_check() const;
    
    // compute distance between two beads
    double bead_dist(int i, int j) const;
    
    int bead_a;
    int bead_b;
    double end_distance;
    double num_of_steps;
    double displacement_rate;
    
private:
    // Total number of beads
    int number_of_beads;
    
    // Set up Anchoring 
    bool anchoring;
    // anchored_beads stores beads that are fixed.
    std::unordered_set<int> anchored_beads;
    void setAnchors();
    
    // Set up weighted beads with factor
    bool weighted;
    double weight_factor;
    
    // Set up Displacement of beads

    bool displacement;
    std::unordered_set<int> displaced_beads;
    void setDisplaceBeads();
    
    std::vector<Module*> modules;
    
    // Parameters
    double collision_radius;
    double collision_scale;
    double random_scale;
    double drag_coeff;
    double base_spring_constant;
    
    // These pointer arrays are for speeding up.
    double* position_pointer;
    double* distance_array_pointer;
    double* norm_array_pointer;
    double* force_array_pointer;
    const double* spring_properties_pointer;
    
    std::vector<int> next_bead, prev_bead;
    
    
    // random engine we use
    std::mt19937_64 &rng;
    // distributions used in random number generation
    std::uniform_int_distribution<int> dist_over_beads; // 0 to numberofbeads-1
    std::uniform_int_distribution<int> dist_01; // 0 or 1
    
    double cr2;
    
    //DTMutableDoubleArray bb_tree;
    // This is where the collision tree nodes are stored (think about memory efficiency)
    CollisionTreeNode* collision_tree_storage;
    // Points to head of tree for each chain
    std::vector<CollisionTreeNode*> tree_pointer_array;
    
    
    // Spring const and length table which is a 2*num_of_beads table. (0,i) stores spring const, (1,i) stores hinge const.
    DTDoubleArray spring_properties;
    
    void constructTree();
    void updateTreeSub(CollisionTreeNode* head);
    bool sphereIntersect(CollisionTreeNode* a, CollisionTreeNode *b) const;
    // compute the separation between two spheres (0 if they intersect)
    double sphereSeparation(CollisionTreeNode* n1, CollisionTreeNode *n2) const;
    void conditional_push(CollisionTreeNode* n1, CollisionTreeNode *n2, std::deque<std::pair<CollisionTreeNode*,CollisionTreeNode*> > &s) const;
    
    void getIntersectList(vector<ipair> &neighbors);
    
    bool n2_check;
};
#endif /* ChromosomeArm_h */
