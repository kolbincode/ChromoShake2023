#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTProgress.h"
#include "unitConversion.h"
#include "Model.hpp"
#include "histone.hpp"
#include "condensin.hpp"
#include "cohesin.hpp"

#include "DTSeriesPath3D.h"
#include "DTSeriesTable.h"

void Computation(const DTDictionary &flags,const DTDictionary &coefficients,
                 double seed,const DTDictionary &evolution,
                 const DTDictionary &histone,const DTDictionary &condensin,
                 const DTDictionary &cohesin,const DTTable &chainsTable,
                 const DTTable &beads_anchored,const DTTable &beads_to_displace,
                 DTDataFile &output) // Write all output to this file
{
    //Read in parameters and flags.
    DTTableColumnNumber values = beads_anchored("bead index");
    DTDoubleArray anchors = values.DoubleVersion();

    values = beads_to_displace("bead index");
    DTDoubleArray displace = values.DoubleVersion();

    DTProgress progress;
    
    DTSeriesPath3D outputBeads(output,"Chromosomes");
    DTSeriesTable condensinOutput(output,"Condensin");
    
    bool histone_flag = flags("histone");
    bool condensin_flag = flags("condensin");
    bool cohesin_flag = flags("cohesin");
    bool spring_flag = flags("spring");
    bool thermomotion_flag = flags("thermo_motion");
    bool collision_flag = flags("collision");
    bool hinge_flag = flags("hinge");
    //  NO Anchor flag here??
    
    double until = evolution("until");
    double save_stride = evolution("save_stride");
    double dt_base = evolution("timeStep");
    int number_of_nodes = evolution("numberOfNodes");
    double displacement_cap = evolution("displacementCap");
    

    // start time at 0
    double t = 0;
    // set current time to 0
    double time_since_save = t;
    // time step actually taken to calculate position updates (less than dt)
    double dt_current_step;
    // tolerance needed to account for small numerical addition errors
    double tolerance = dt_base * 1e-4;
    double save_stride_judge = save_stride - tolerance;
    double current_max_force;
    // displacementCap set to 5 nm by default
    double force_cap = displacement_cap/dt_base;
    
    mt19937_64 generator(seed);
    
    // Construct a chains array from the chainsTable object. Should really hand in the DTTable object
    DTMutableDoubleArray chains(chainsTable.NumberOfColumns(),chainsTable.NumberOfRows());
    chains = NAN;
    for (int colN=0;colN<chainsTable.NumberOfColumns();colN++) {
        DTTableColumnNumber singleColumn = chainsTable(colN);
        DTDoubleArray values = singleColumn.Values();
        for (int rowN=0;rowN<values.Length();rowN++) {
            chains(colN,rowN) = values(rowN);
        }
    }
    /*
    Model(double t, double timestep,
          const DTDictionary &coefficients,
          const DTDictionary &evolution,
          const DTDictionary &flags,
          const DTDictionary &cohesin,
          const double seed,
          const DTDoubleArray &chains, const DTDoubleArray &anchors,
          int _number_of_beads,
          mt19937_64 &rng);
     */
    
    // Initialize data structures including the chromosome arm, histones and condensins.
    Model model = Model(t, dt_base, coefficients, evolution, flags, cohesin, seed, chains, anchors, displace,number_of_nodes, generator);
    
    double base_spring_constant = model.get_base_spring_constant();
    
    if (histone_flag)
    {
        Histone *histone_temp = new Histone(&model, histone, number_of_nodes, generator, chains, base_spring_constant);
        model.addModule(histone_temp);
    }
    
    
    Condensin *condensin_info = NULL;
    if (condensin_flag)
    {
        condensin_info = new Condensin(&model, condensin, number_of_nodes, seed, chains, generator, base_spring_constant);
        model.addModule(condensin_info);
    }

    
    if (cohesin_flag)
    {
        Cohesin *cohesin_temp = new Cohesin(&model, chains, &generator, cohesin, base_spring_constant);
        model.addModule(cohesin_temp);
    }
    
    //Add first time step to output.
    DTMutableDoubleArray to_update = model.position.Copy();
    to_update *= 1e-9;

    auto tmp = model.convertPointsToPath();
    
    DTIntArray to_update_condensin;
    if (condensin_flag)
    {
        to_update_condensin = condensin_info->condensin_array;
    }
    
    // this gets set to true to force a save after the timestep is reduced to match the save stride
    bool force_save = false;
    
    outputBeads.Add(tmp,t);
    //to_update_group.Var = tmp;
    //to_update_group.condensin = to_update_condensin.Copy();
    //computed.Add(to_update_group,t);
     
///////////////////////////////////////////////////////
/////////             Main loop           //////
///////////////////////////////////////////////////////
    ///
    ///////////// use lines 140 and 144 for debugging behavior at a certain time ///////
    ///
    // double stopAt = 0.000202;
    while (t < until)
    { //compute dynamic time step to avoid instablity lead by all forces.
        
        // if (t>stopAt) {  std::cerr << "stopped at...\n"; }
        model.initializeForce();
        model.computeDistance();
        
        if (spring_flag) model.springForceUpdate();
        if (hinge_flag) model.hingeForceUpdate();
        if (collision_flag) model.excludedVolumeUpdate();

        model.moduleForceUpdate();
        
        
        //Determine dynamic timestep for this run
        current_max_force = model.findMaxForce();
        dt_current_step = dt_base;
        // Enforce displacement cap condition
        if (current_max_force > force_cap) {
            dt_current_step = displacement_cap/current_max_force;
        }
        // Avoid overstepping save time
        if (time_since_save + dt_current_step >= save_stride_judge)
        {
            dt_current_step = save_stride_judge - time_since_save;
            // to avoid issues with floating point numbers not quite lining up
            force_save = true;
        }
        
        //Determine random motion parameters based on time step.
        if (thermomotion_flag) model.randomForceUpdate(dt_current_step);
        
        // Now do the update. First update the positions of the chromosome and cohesin bead-springs. Then update histones and condensins, including their built-in timers, attaching sites and conditions. The update is done through member functions.
        model.positionUpdate(dt_current_step);
        
        t += dt_current_step;
        time_since_save += dt_current_step;            // to store time period for output.
        
        model.moduleTimeUpdate(dt_current_step);

        // Add new instances to the output file.
        if (force_save) {
            force_save = false;
            progress.UpdatePercentage(t/until);
            
            if (condensin_flag && condensin_info)
            {
                to_update_condensin = condensin_info->condensin_array.Copy();
                // Array
                
                DTMutableList<DTTableColumn> condensinTableColumns(2);
                ssize_t condensinCount = to_update_condensin.n();
                DTMutableIntArray fromList(condensinCount);
                DTMutableIntArray toList(condensinCount);
                for (int ptN=0;ptN<condensinCount;ptN++) {
                    fromList(ptN) = to_update_condensin(0,ptN);
                    toList(ptN) = to_update_condensin(1,ptN);
                }

                condensinTableColumns(0) = CreateTableColumn("from",fromList);
                condensinTableColumns(1) = CreateTableColumn("to",toList);

                DTTable condensinTable(condensinTableColumns);
                condensinOutput.Add(condensinTable,t);
            }
            
            outputBeads.Add(model.convertPointsToPath(),t);

            //to_update_group.condensin = to_update_condensin.Copy();
            //computed.Add(to_update_group,t);
            time_since_save -= save_stride;
            // do we need cohesins crosslinking updated here?   Apparently not... they are handled in Model
        }
    } // end of main
}
