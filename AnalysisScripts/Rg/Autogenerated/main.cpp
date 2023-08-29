// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#include "computation.h"

#include "DTArguments.h"
#include "DTTimer.h"
#include "DTDataFile.h"
#include "DTError.h"

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTPath3D BeadsPath;
    double stride;

    {
        // Inside a scope so that the data files will be closed before the computation starts.
        DTDataFile inputDataFile("Input.dtbin",DTFile::ReadOnly);
        if (inputDataFile.IsOpen()==false) {
            std::cerr << "No input file found. Might have to save input for debugging." << std::endl;
        }

        stride = inputDataFile.ReadNumber("stride");
        Read(inputDataFile,"BeadsPath",BeadsPath);
    }

    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    //DTTimer timer;
    //timer.Start();
    DTTable output = Computation(stride,BeadsPath);

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
