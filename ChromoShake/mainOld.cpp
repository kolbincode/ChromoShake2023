// #include "DTSource.h"
#include "DTSaveError.h"
#include "DTSeriesGroup.h"
#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTDoubleArray.h"
#include "DTMatlabDataFile.h"
#include "DTPath3D.h"
#include "DTProgress.h"
#include "DTSeriesPath3D.h"
#include "DTRandom.h"
#include "unitConversion.h"
#include "Model.hpp"
#include "histone.hpp"
#include "condensin.hpp"
#include "cohesin.hpp"


//////////////////////////////////////////////////////////////////////////////
//    DT_RetGroup
//////////////////////////////////////////////////////////////////////////////

struct DT_RetGroup {
    DTPath3D Var;
    DTDoubleArray centromereForce;
    DTIntArray condensin;
    
    void pinfo(void) const;
    void pinfoIndent(string) const;
    
    static void WriteStructure(DTDataStorage &,string);
};

// pinfo??. position information
void DT_RetGroup::pinfo(void) const
{
    pinfoIndent("");
}

// what is pinfoIndent?.. position information
void DT_RetGroup::pinfoIndent(string pad) const
{
    cerr << pad << "Var = "; Var.pinfo();
    cerr << pad << "centromereForce = "; centromereForce.pinfo();
    cerr << pad << "condensin = "; condensin.pinfo();
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save("Var",name+"_1N");
    output.Save("Path3D",name+"_1T");
    output.Save("centromereForce",name+"_2N");
    output.Save("Array",name+"_2T");
    output.Save("condensin",name+"_3N");
    output.Save("Array",name+"_3T");
    output.Save(3,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
    Write(output,name+"_Var",var.Var);
    Write(output,name+"_centromereForce",var.centromereForce);
    Write(output,name+"_condensin",var.condensin);
    Write(output,name,DTDoubleArray()); // So that DataTank can see the variable.
}

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTDictionary &Coefficients,const DTDictionary &flags,
                 const DTDictionary &evolution,
                 double seed, const DTDictionary &histone,
                 const DTDictionary &condensin, const DTDictionary &cohesin,
                 const DTDoubleArray &chains, const DTDoubleArray &anchors,
                 const DTDoubleArray &beads_to_displace,
                 DTSeriesGroup<DT_RetGroup> &computed);

int main(int argc, const char *argv[])
{
    DTSetArguments(argc, argv);
    
    DTMatlabDataFile inputFile("Input.mat",DTFile::ReadOnly);
    DTMatlabDataFile outputFile("Output.mat",DTFile::NewReadWrite);
    
    // Input variables.
    DTDoubleArray chains = inputFile.ReadDoubleArray("chains");
    DTDoubleArray anchors = inputFile.ReadDoubleArray("beads_anchored");
    DTDoubleArray beads_to_displace = inputFile.ReadDoubleArray("beads_to_displace");
    double seed = inputFile.ReadNumber("seed");
    
    DTDictionary Coefficients;
    Read(inputFile,"Coefficients", Coefficients);
    DTDictionary flags;
    Read(inputFile,"flags", flags);
    DTDictionary evolution;
    Read(inputFile,"evolution", evolution);
    DTDictionary histone;
    Read(inputFile,"histone", histone);
    DTDictionary condensin;
    Read(inputFile,"condensin",condensin);
    DTDictionary cohesin;
    Read(inputFile,"cohesin", cohesin);
    
    
    // Output series.
    DTSeriesGroup<DT_RetGroup> computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput"))
    { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"Coefficients",Coefficients);
        WriteOne(outputFile,"flags",flags);
        WriteOne(outputFile,"evolution",evolution);
        WriteOne(outputFile,"seed",seed);
        WriteOne(outputFile,"histone",histone);
        WriteOne(outputFile,"condensin",condensin);
    }
    
///////////////////////////////////////////////////////////
////////////    The computation   /////////////
//////////////////////////////////////////////////////////
    clock_t t_before = clock();
    Computation(Coefficients,flags,evolution,seed,histone,condensin,cohesin,chains,anchors,beads_to_displace,computed);

    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);  //double casts to double
    printf("Execution completed in %f seconds\n", exec_time);
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    return 0;
}


//////////////////////////////////////////////////////////////////////////////
///        Computational routine
//////////////////////////////////////////////////////////////////////////////

DTPath3D ConvertPointsToPath(const DTMutableDoubleArray &points);

DTPath3D ConvertPointsToPath(const DTMutableDoubleArray &points)
{
    int n = points.n();
    DTMutableDoubleArray loop(3,n+1);
    // loop = NAN;
    loop(0,0) = 0;
    loop(1,0) = 0;
    loop(2,0) = n;
    
    MemoryCopyColumns(loop,1,points,DTRange(0,n));
    
    return DTPath3D(loop);
}

void Computation(const DTDictionary &coefficients,const DTDictionary &flags,
                 const DTDictionary &evolution,
                 double seed, const DTDictionary &histone,
                 const DTDictionary &condensin, const DTDictionary &cohesin,
                 const DTDoubleArray &chains, const DTDoubleArray &anchors,
                 const DTDoubleArray &beads_to_displace,
                 DTSeriesGroup<DT_RetGroup> &computed)
{
    //Read in parameters and flags.
    
    DTProgress progress;
    DT_RetGroup to_update_group;
    
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
    
    
    // Initialize data structures including the chromosome arm, histones and condensins.
    Model model = Model(t, dt_base, coefficients, evolution, flags, cohesin, seed, chains, anchors, beads_to_displace, number_of_nodes, generator);
    // model.SetCoeffients(coefficients);
    // model.SetImplicitBoundary(distance);
    
    double base_spring_constant = model.get_base_spring_constant();
    
    if (histone_flag)
    {
        Histone *histone_temp = new Histone(&model, histone, number_of_nodes, generator, chains, base_spring_constant);
        model.addModule(histone_temp);
    }
    
    
    Condensin *condensin_info;
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
    
    DTMutableIntArray to_update_condensin;
    if (condensin_flag)
    {
        to_update_condensin = condensin_info->condensin_array.Copy();
    }
    
    // this gets set to true to force a save after the timestep is reduced to match the save stride
    bool force_save = false;
    
    
    to_update_group.Var = tmp;
    to_update_group.condensin = to_update_condensin.Copy();
    computed.Add(to_update_group,t);
    
///////////////////////////////////////////////////////
/////////             Main loop           //////
///////////////////////////////////////////////////////
    while (t < until)
    { //compute dynamic time step to avoid instablity lead by all forces.
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
        
        /* DK: use dt_current_step, t/until, or progress.UpdatePercentage to accurately compute
         rate displacement for beadsInMotion */

        // Add new instances to the output file.
        // TODO: look example of force applied to centromeres from old YY code
        if (force_save) {
            force_save = false;
            progress.UpdatePercentage(t/until);

//            // old code used to apply force to beads representing centromeres
//            if (anchoring)
//            {
//                for (int i = 0; i < 3; ++i) {
//                    to_update_centromere_force(i,0) = model.force_array(i,start_bead);
//                    to_update_centromere_force(i,1) = model.force_array(i,end_bead);
//                }
//                to_update_centromere_force *= 1e-9;
//            }
            
            if (condensin_flag)
            {
                to_update_condensin = condensin_info->condensin_array.Copy();
            }
            
            to_update_group.Var = model.convertPointsToPath();
//            to_update_group.centromereForce = to_update_centromere_force.Copy();
            to_update_group.condensin = to_update_condensin.Copy();
            computed.Add(to_update_group,t);
            time_since_save -= save_stride;
            // do we need cohesins crosslinking updated here?   Apparently not... they are handled in Model
        }
    } // end of main
}// end of Computation

