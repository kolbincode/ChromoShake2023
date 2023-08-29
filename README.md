# ChromoShake2023
 
Simulation program for small scale polymer bead-spring model of Chromatin


## Description of parameters
The Chromoshake is configured with a set of parameters with default values. Parameters are grouped in dictionaries. Parameters and their values used for the experiment are highlighted with ***

### Parameters not in dictionaries
1. seed: Sets integer used to seed the random number generator used in the simulation
2. chains: Provides 8-by-k matrix that describes the geometry of the bead-spring model (todo: describe in more detail)
            _chainArrays_ (first_bead_index, last_bead_index, LOOP_flag, length_nanometers, X_pos_first, Y_pos_first, Z_pos_first, DNA_flag) 
3. anchor list:    Provides a list of bead indices to be fixed in position (remember bead indices start at 0 end on n-1).
            
            
## Flags
These parameters are boolean flags toggling certain parts of the model
1. histone = 0                                                       // [0/1], includes histone module in the model  (cross liker with histone spring)
2. condensin = 1                                                     // [0/1], includes condensin module in the model  (cross liker with condensin spring)
3. cohesin = 0                                                       //  [0/1], includes cohesin module in the model  (bead spring rings around main chain)          
4. collision = 1                                                          //  [0/1],  Volume exclusion force
5. spring = 1                                                             // [0/1],  Spring forces 
6. thermo_motion = 1                                               // [0/1],  Browinian forces
7. hinge = 1                                                          // [0/1],  Force that keeps DNA linear related to Lp
8. anchor = 0  or  1***                                   // [0/1],  allows to fix positions of indicated beads, which become tethers if 1.
9. displacement = 0                                       // [0/1], are the end beads (tethers) being displaced?
10. weighted = 1  or  0***                                  // [0/1], Do the end beads (anchors, tethers) provide resitance/drag when force is applied (act as weights)


## Coefficients 
These parameters contain physical constants describing the model.
1. mass_separation_nanometers = 10                  // resting length of the spring between beads = 10nm
2. temperature_Celsius = 25                               // Room temperature. C is converted to K for physical calculations
3. viscosity_centiPoise = 1                               // viscosity 10 Poise = 1 cPoise = 1 Pa * sec
4. modulus_gigaPascal = 0, 3.25, 16***                      // (hinge factor) determines stiffness of the chain by scaling hinge force
5. const_dna_modulus_gigaPascal = 2                        // Young's modulus for DNA:  2 GPa
6. dna_radius_nanometers = 0.6                             //  (0.6)   to calculate tensile stiffness of DNA spring when it's diameter when used in formula with Youngs modulus
7. damping_radius_factor = 0.8                                //  (0.8)  Does this scale the drag force on DNA 
8. collision_radius_factor = 0.9                                //  (0.9) sets the diameter of the vol exclusion spheres around the beads    SET >1 and check for rings slipping.
9. collision_force_factor = 0.75                                //  (0.25) scales collision force
10. final_bead_distance = 50                                   // final distance between beads if moving beads are displaced
11. weight_factor = 10e0... 10e9***                                 // factor to scale force down o end beads/anchors/tethers whe weighted_flag is on.                                         


## Evolution
These parameters govern aspects of the numerical simulation of the model.
1. timeStep = 2e-10                            // (2e-9)    time step in seconds    ~ default: 1 nanosecond
2. until = 35e-3                               //  35 milliseconds,    duration of simulation    ~ default:   until 0.008 seconds
3. save_stride = 2e-6                         // (2e-6)    interval in second for saving position data     ~ default: 1 microsecond
4. numberOfNodes = 101                   //  1um chain, number of nodes is number of bead determines length of DNA segment modeled
5. displacementCap = 3                   //  (5) Max distance nanometers a bead is allowed to move in dt.     


### histone NOT USED
1. random_placement = 1                       // [0/1]  flag that turns on/off random placement or first available histone
2. numberOfHistones = 80                     // default: 80    (for 80 full occupancy is on 11.7 kb ~ 390 Nodes) 
3. beads_in_histone_loop = 7                // default: 7
4. histone_modulus = 5                          // default: 5   specifies a factor by which _spring_constant_ for histone crosslink is scaled 
5. halfLifeOn = 5e-3                                  // (5e-3)  mean of exponential dist of how long histones stays on
6. halfLifeOff = 0                                       // (0)  mean of exponential dist of how long histones stays on
7. unattach = 1                                         // [0/1]  flag turns on/off histone unattachment
8. unattachTime = 0                                 //
9. unattachStartTime = 1e-7                  //
10. unattachThreshold = 30                     //


### condensin
1. numberOfCondensins = 0...6***           // default: 6     3 condensins for each 10kb (arms),  6 condensins for 10kb (pericentric)
2. condesin_modulus = 0.02/ 0.2/ 2.0***    // default: 5   specifies a factor by which _spring_constant_ for condensin crosslink is scaled 
3. dynamicExtrusionRate = 1                 // [0/1]  flag that turns on/off dynamicExtrusionRate for condensins
4. unbind = 1                                   // [0/1]  flag for condensins to be allowed to unbind/rebind or stay fixed
5. unbindTime =      0.01135             //  (160/14100)      plan to make it an *exponential random* variable
6. unbindTimeVar = 0.00425          //   (60/14100)    plan to make it an *exponential random* variable
7. rebindTime = 1e-8                               // (1e-8)     plan to make it an *exponential random* variable
8. judgeTimeSeconds = 2e-5           //
9. criticalLengthNanometers = 30         // max of condensins crosslinking spring, beyond which it breaks at weakSite
10. incrementNanometers = 10          // 10nm goes into direction of StrongSite, moves StrongSite to closest bead
11. extrusionRatePower = 1             //  increased will reduce extrusion, decrease will increase extrusion affecting judge time btwn extrusions.
12. polaritySwitchTime = 2e-2          // Frequence of condensin switching directions (weak site and strong site interchanged)
 

### cohesin  NOT USED
1. meanOn = 5e-3                   // (5e-3)    Time in sec two cohesins are bonded
2. bond_dist = 20                   // (15)  Distance in nanometers below which cohesins are allowed to bond
3. cohesin_modulus = 1         // (1)  specifies the strength of spring (factor of _spring_constant_) that crosslinks cohesins
4. restore = 0                          // (1)  allows slipped cohesin rings to pop back on the DNA chain

