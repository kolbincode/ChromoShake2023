#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(double stride,const DTPath3D &BeadsPath)
{
    //Error check
    if (BeadsPath.IsEmpty()) {  DTErrorMessage("No data in beads");   return DTTable();  }
    
    //only one row for output
    int num_rows = 1;
    DTMutableDoubleArray comOut(3,num_rows);
    DTMutableDoubleArray RgOut(num_rows);
    
    // loop runs once :)
    for (int i=0; i<num_rows; i++ )
    {
        DTPath3D path = BeadsPath;
        DTMutableDoubleArray points = Points(path);
        ssize_t numPoints = points.n();
        //int totalPoints = path.NumberOfPoints();
        
        // COM computation
        double sumx = 0;
        double sumy = 0;
        double sumz = 0;
        
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
    }
    
    // Create the table for output
    DTMutableList<DTTableColumn> columns(1);
    columns(0) = CreateTableColumn("Rg",RgOut);
    return DTTable(columns);
}
