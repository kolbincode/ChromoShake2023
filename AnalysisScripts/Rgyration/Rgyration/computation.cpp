#include "computation.h"
#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTSet<DTPath3D> &PositionBeads,int Stride)
{
    // PositionBeads is a set
    //Error check
    if (PositionBeads.IsEmpty()) {  DTErrorMessage("No data in beads");   return DTTable();  }
    
    //Parameters
    DTTable timetable = PositionBeads.Parameters();
    DTTableColumnNumber timecolumn = timetable("time");
    DTDoubleArray timevalues = timecolumn.Values();
    
    ssize_t num_rows = timevalues.Length()/Stride;
    
    DTMutableDoubleArray timeValuesOut(num_rows);
    DTMutableDoubleArray comOut(3,num_rows);
    DTMutableDoubleArray RgOut(num_rows);
    
    for (int i=0; i<num_rows; i++ )
    {
        DTPath3D path = PositionBeads(i*Stride);
        DTMutableDoubleArray points = Points(path);
        // int totalPoints = path.NumberOfPoints();
        
        // COM computation
        double sumx = 0;
        double sumy = 0;
        double sumz = 0;
        ssize_t numPoints = points.n();
        for (int ptN=0; ptN<numPoints; ptN++)
        {
            sumx += points(0,ptN);
            sumy += points(1,ptN);
            sumz += points(2,ptN);
        }
        comOut(0,i) = sumx/numPoints;
        comOut(1,i) = sumy/numPoints;
        comOut(2,i) = sumz/numPoints;
        
        double Rgx = 0;
        double Rgy = 0;
        double Rgz = 0;
        
        for (int ptN=0; ptN<numPoints; ptN++)
        {
            Rgx += (points(0,ptN) - comOut(0, i))*(points(0,ptN) - comOut(0, i));
            Rgy += (points(1,ptN) - comOut(1, i))*(points(1,ptN) - comOut(1, i));
            Rgz += (points(2,ptN) - comOut(2, i))*(points(2,ptN) - comOut(2, i));
        }
        Rgx = Rgx/numPoints;
        Rgy = Rgy/numPoints;
        Rgz = Rgz/numPoints;
        // get the norm of radius of gyration vector
        RgOut(i) = sqrt(Rgx + Rgy + Rgz);
        timeValuesOut(i) = timevalues(i*Stride);
    }
    
    // Create the table
    DTMutableList<DTTableColumn> columns(2);
    columns(0) = CreateTableColumn("time",timeValuesOut);
    columns(1) = CreateTableColumn("Rg",RgOut);
    return DTTable(columns);
}
