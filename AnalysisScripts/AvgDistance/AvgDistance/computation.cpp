#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &EndtoEnd_forAvg,int NumSim)
{
    if (EndtoEnd_forAvg.IsEmpty())
    {
        return DTTable();
    }
    
    DTTableColumnNumber timeValues = EndtoEnd_forAvg("time");
    DTDoubleArray time = timeValues.DoubleVersion();
    // distance
    DTTableColumnNumber avegrageValues = EndtoEnd_forAvg("distance");
    DTDoubleArray Avg = avegrageValues.DoubleVersion();
    
    double size = Avg.m();
    int sizeOut = size / NumSim;
    DTMutableDoubleArray timeOut(sizeOut, 1);
    DTMutableDoubleArray AverageOut(sizeOut,1);
    
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
        AverageOut(j) = sum/NumSim;
        timeOut(j) = timeValues(j);
    }
    
    // Create the table
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("time",timeOut);
    columns(1) = CreateTableColumn("AvgDistance",AverageOut);
    return DTTable(columns);
}
