// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"
#include "DTError.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTTable EndtoEnd_forAvg;
    int NumSim;

    {
        // Inside a scope so that the input data file will be closed before the computation is called.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        if (inputDataFile.IsOpen()==false) {
            std::cerr << "No input file found. Might have to save input for debugging." << std::endl;
        }
        DTDataFile variableDataFile;

        variableDataFile = DTDataFile("EndtoEnd_forAvg.dtbin",DTFile::ReadOnly);
        Read(variableDataFile,"EndtoEnd_forAvg",EndtoEnd_forAvg);
        NumSim = inputDataFile.ReadNumber("NumSim");
    }

    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    //DTTimer timer;
    //timer.Start();
    DTTable output = Computation(EndtoEnd_forAvg,NumSim);

    //timer.Stop(); // Use timer.Time() to get the elapsed time
    if (DTHowManyErrors()>0) outputFile.Save(DTHowManyErrors(),"ErrorCount"); // For error logging

    WriteOne(outputFile,"Var",output);
    // The structure, to make it easy to open the output file
    output.WriteStructure(outputFile,"SeqInfo_Var");
    outputFile.Save("Table","Seq_Var");

    // To speed up reading.
    outputFile.SaveIndex();

    return 0;
}