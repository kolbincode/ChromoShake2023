// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDictionary flags, coefficients, evolution, histone, condensin, cohesin;
    DTTable chains, beads_anchored, beads_to_displace;
    double seed;

    {
        // Inside a scope so that the input data file will be closed before the computation is called.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        DTDataFile variableDataFile;

        Read(inputDataFile,"flags",flags);
        Read(inputDataFile,"coefficients",coefficients);
        seed = inputDataFile.ReadNumber("seed");
        Read(inputDataFile,"evolution",evolution);
        Read(inputDataFile,"histone",histone);
        Read(inputDataFile,"condensin",condensin);
        Read(inputDataFile,"cohesin",cohesin);
        variableDataFile = DTDataFile("chains.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"chains",chains);
        variableDataFile = DTDataFile("beads_anchored.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"beads_anchored",beads_anchored);
        variableDataFile = DTDataFile("beads_to_displace.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"beads_to_displace",beads_to_displace);
    }

    // Create the output file.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    //DTTimer timer;
    //timer.Start();
    Computation(flags,coefficients,seed,evolution,histone,condensin,cohesin,chains,beads_anchored,
                beads_to_displace,outputFile); // Write all output to this file

    //timer.Stop(); // Use timer.Time() to get the elapsed time

    // To speed up reading.
    outputFile.SaveIndex();

    return 0;
}
