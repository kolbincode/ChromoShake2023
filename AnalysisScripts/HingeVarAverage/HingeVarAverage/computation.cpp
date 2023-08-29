#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

Group Computation(const DTTable &Rgyration)
{
    
    DTTableColumnNumber RgColumn = Rgyration("Rg");
    DTDoubleArray RgValues = RgColumn.Values();
    
    ssize_t num_rows = RgValues.Length();
    
    double sum = 0;
    for (int i=0; i<num_rows; i++ )
    {
        sum += RgColumn(i);
    }
   
    
    Group toReturn;

    toReturn.averageRG = sum/num_rows;

    return toReturn;
}

