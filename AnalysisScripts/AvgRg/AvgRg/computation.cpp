#include "computation.h"
#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &BeadPos,int NumSim)
{
    if (BeadPos.IsEmpty()) { DTErrorMessage("Bead Position is empty");   return DTTable();}
    // time
    DTTableColumnNumber timeValues = BeadPos("time");
    DTDoubleArray time = timeValues.DoubleVersion();
    // Avg
    DTTableColumnNumber AvgValues = BeadPos("Rg");
    DTDoubleArray Avg = AvgValues.DoubleVersion();
    
    double size = Avg.m();
    int sizeOut = size / NumSim;
    DTMutableDoubleArray timeOut(sizeOut, 1);
    DTMutableDoubleArray AverageOut(sizeOut,1);
    DTMutableDoubleArray SigmaOut(sizeOut,1);
    
    double Var;
    
    double sum;
    for (int j = 0 ; j < sizeOut ; ++j)
    {
        sum = 0;
        int counter = 0;
        
        for (int i = 0; i< NumSim ; ++i)
        {
            sum += Avg(j+counter);
            counter += sizeOut;
        }
        
        // Get variance and SD
        AverageOut(j) = sum/NumSim;
        Var = 0;
        counter = 0;
    
        for (int i = 0; i< NumSim ; ++i)
        {
            Var += (Avg(j+counter)-AverageOut(j))*(Avg(j+counter)-AverageOut(j));
            counter += sizeOut;
        }
        
        Var = Var/(NumSim-1);
        SigmaOut(j) = sqrt(Var);
        timeOut(j) = timeValues(j);
    }
    
    // Create the table
    DTMutableList<DTTableColumn> columns(3);
    columns(0) = CreateTableColumn("time",timeOut);
    columns(1) = CreateTableColumn("avgRg",AverageOut);
    columns(2) = CreateTableColumn("sigmaRg",SigmaOut);
    
    return DTTable(columns);
} // end of computation
