#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTSet<DTPath3D> &BeadPositions,int Stride)
{
    // BeadPositions is a set
    
    // Error check
    if (BeadPositions.IsEmpty()) {  DTErrorMessage("No data in beads");   return DTTable();  }
    
    //Parameters
    DTTable timetable = BeadPositions.Parameters();
    DTTableColumnNumber timecolumn = timetable("time");
    DTDoubleArray timevalues = timecolumn.Values();
    
    ssize_t num_rows = timevalues.Length()/Stride;
    
    DTMutableDoubleArray timeValuesOut(num_rows);
    DTMutableDoubleArray comOut(3,num_rows);
    
    for (int i=0; i<num_rows; i++ )
    {
        DTPath3D path = BeadPositions(i*Stride);
        DTMutableDoubleArray points = Points(path);
        // COM computation
        double sumx = 0;
        double sumy = 0;
        double sumz = 0;
        ssize_t numPoints = points.n();
        for (int ptN=0;ptN<numPoints;ptN++) {
            sumx += points(0,ptN);
            sumy += points(1,ptN);
            sumz += points(2,ptN);
        }
        comOut(0,i) = sumx/numPoints;
        comOut(1,i) = sumy/numPoints;
        comOut(2,i) = sumz/numPoints;
        timeValuesOut(i) = timevalues(i*Stride);
    }
    
    // Create the table
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("time",timeValuesOut);
    columns(1) = CreateTableColumn("COM",DTPointCollection3D(comOut));
    return DTTable(columns);
}
