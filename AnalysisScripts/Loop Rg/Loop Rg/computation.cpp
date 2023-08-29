#include "computation.h"

#include <math.h>
#include "DTDoubleArrayOperators.h"
#include "DTUtilities.h"

DTTable Computation(const DTTable &table,const DTSet<DTPath3D> &chromosomes)
{
    ssize_t rowCount = table.NumberOfRows();
    DTTableColumnNumber run = table("Run");
    DTTableColumnNumber time = table("time");
    
    DTTableColumnNumber from = table("from");
    DTTableColumnNumber to = table("to");
    
    DTTable chromosomesMetaData = chromosomes.Parameters();
    DTTableColumnNumber runFromMeta = chromosomesMetaData("Run");
    DTDoubleArray runValues = runFromMeta.Values();
    DTTableColumnNumber timeFromMeta = chromosomesMetaData("time");
    DTDoubleArray timeValues = timeFromMeta.Values();
    
    ssize_t rowMeta,rowCountInMeta = chromosomesMetaData.NumberOfRows();
    
    DTMutableDoubleArray RgOutColumn(rowCount);
    DTMutableDoubleArray comOut(3,rowCount);
    // double dt = timeValues(1)-timeValues(0);
    // int howManyTimes = find the first index >0 s.that timeValues(index)=0
    //                  = find where runCountInMeta(index)==1
    
    
    for (ssize_t row=0;row<rowCount;row++) {
        double runNumber = run(row);
        double timeValue = time(row);
        
        // timeIndex = round(timeValue/dt)
        // find the row of (runNumber,timeValue) from the run/time FromMeta
        
        for (rowMeta=0;rowMeta<rowCountInMeta;rowMeta++) {
            if (runValues(rowMeta)==runNumber && timeValues(rowMeta)==timeValue) {
                break;
            }
        }
        
        if (rowMeta==rowCountInMeta) {
            // Failed to find it
            DTErrorMessage("Failed to find");
            continue;
        }
        
        int fromIndex = from(row);
        int toIndex = to(row);
        if (fromIndex>toIndex) std::swap(fromIndex,toIndex);
        
        DTPath3D chromosome = chromosomes(rowMeta);
        DTMutableDoubleArray points = Points(chromosome);
        
        // COM computation
        double sumx = 0;
        double sumy = 0;
        double sumz = 0;
        
        ssize_t numPoints = toIndex-fromIndex+1;
        for (int ptN=fromIndex; ptN<toIndex+1; ptN++)
        {
            sumx += points(0,ptN);
            sumy += points(1,ptN);
            sumz += points(2,ptN);
        }
        comOut(0,row) = sumx/numPoints;
        comOut(1,row) = sumy/numPoints;
        comOut(2,row) = sumz/numPoints;
        
        double Rgx = 0;
        double Rgy = 0;
        double Rgz = 0;
        
        for (int ptN=fromIndex; ptN<toIndex+1; ptN++)
        {
            Rgx += (points(0,ptN) - comOut(0, row))*(points(0,ptN) - comOut(0,row));
            Rgy += (points(1,ptN) - comOut(1, row))*(points(1,ptN) - comOut(1,row));
            Rgz += (points(2,ptN) - comOut(2, row))*(points(2,ptN) - comOut(2,row));
        }
        Rgx = Rgx/numPoints;
        Rgy = Rgy/numPoints;
        Rgz = Rgz/numPoints;
        // get the norm of radius of gyration vector
        RgOutColumn(row) = sqrt(Rgx + Rgy + Rgz);  // the actual value
    }
    
    
    // Table is a list of columns
    return DTTable({
        CreateTableColumn("Run",run),
        CreateTableColumn("time",time),
        CreateTableColumn("Rg",RgOutColumn)
    });
}
