#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &RG_combined_table)
{
    
    // t
    //   DTTableColumnNumber values = RG_combined_table("t");
    //   DTDoubleArray t = values.DoubleVersion();
    // timeStep
    DTTableColumnNumber timeColumn = RG_combined_table("timeStep");
    DTDoubleArray timeStepValues = timeColumn.DoubleVersion();
    // aveRg
    DTTableColumnNumber RgColumn = RG_combined_table("aveRg");
    DTDoubleArray RgValues = RgColumn.DoubleVersion();
    
    DTMutableDoubleArray output(4);
    DTMutableDoubleArray timeOut(4);
    ssize_t totalLength = RgValues.Length();
    ssize_t  pos = 0;
    int groupNumber = 0;
    while (pos<totalLength)
    {
        double groupDT = timeStepValues(pos);
        double sum = 0;
        int count = 0;
        while (pos<totalLength && timeStepValues(pos)==groupDT)
        {
            sum += RgValues(pos);
            pos++;
            count++;
        }
        if (output.Length()==groupNumber)
        {
            output = IncreaseSize(output,output.Length());
            timeOut = IncreaseSize(timeOut,output.Length());
        }
        
        output(groupNumber) = sum/count;
        timeOut(groupNumber) = groupDT;
        groupNumber++;
    }
    output = TruncateSize(output,groupNumber);
    timeOut = TruncateSize(timeOut,groupNumber);

    // Create the output table
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("timeStep",timeOut);
    columns(1) = CreateTableColumn("aveRg",output);
    return DTTable(columns);
    
    
    //    for (i=0;i<len;i++) {
    //        block
    //    }
    //
    //    i = 0;
    //    while (i<length) {
    //        block
    //        i++;
    //    }
    
}
