// Auto-generated by ImageTank, can overwrite using the C++ template gear menu.

#ifndef IT_computation_h
#define IT_computation_h

#include "structures.h"

#include "DTProgress.h"

#include "DTDataFile.h"

void Computation(const DTDictionary &flags,const DTDictionary &coefficients,double seed,
                 const DTDictionary &evolution,const DTDictionary &histone,
                 const DTDictionary &condensin,const DTDictionary &cohesin,const DTTable &chains,
                 const DTTable &beads_anchored,const DTTable &beads_to_displace,
                 DTDataFile &output); // Write all output to this file

#endif /* IT_computation_h */ 
