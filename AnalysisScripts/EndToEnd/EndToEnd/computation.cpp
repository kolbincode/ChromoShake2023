#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"
#include "DTSet.h"

DTTable Computation(const DTSet<DTPath3D> &BeadPos,int Stride)
{
    // BeadPos is a set
    if (BeadPos.IsEmpty()) {  DTErrorMessage("No data in beads");   return DTTable();  }
    
    //Parameters
    DTTable timetable = BeadPos.Parameters();
    DTTableColumnNumber timecolumn = timetable("time");
    DTDoubleArray timevalues = timecolumn.Values();
    
    ssize_t num_rows = timevalues.Length()/Stride;
    
    DTMutableDoubleArray timeValuesOut(num_rows);
    DTMutableDoubleArray distOut(num_rows);
    
    for (int i=0; i<num_rows; i++ )
    {
        DTPath3D path = BeadPos(i*Stride);
        DTMutableDoubleArray points = Points(path);
        // end to end distance computation
        ssize_t numPoints = points.n();
        int lastPtIndex = (int)numPoints - 1;
        int firstPtIndex = 0;
        
        double dx = points(0,firstPtIndex) - points(0,lastPtIndex);
        double dy = points(1,firstPtIndex) - points(1,lastPtIndex);
        double dz = points(2,firstPtIndex) - points(2,lastPtIndex);
        
        distOut(i) = sqrt (dx*dx + dy*dy + dz*dz);
        timeValuesOut(i) = timevalues(i*Stride);
    }
    
    // Create the table
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("time",timeValuesOut);
    columns(1) = CreateTableColumn("distance",distOut);
    return DTTable(columns);
}
