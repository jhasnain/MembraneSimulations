#ifndef MESHBIAS_H
#define MESHBIAS_H

#include "membranesys.h"

void (*ComputeBiasPotential)();
double (*MeshBiasChange)();
double (*MeshBiasGlobalChange)();
void (*StoreBiasChange)();
double (*PrintBiasPotential)();

void ComputeBiasPotentialMinDist(sys *Sys);
double MeshBiasChangeMinDist(sys *Sys);
void StoreBiasChangeMinDist(sys *Sys, double DeltaEBias);

void ComputeBiasPotentialAverage(sys *Sys);
double MeshBiasChangeAverage(sys *Sys);
void StoreBiasChangeAverage(sys *Sys, double DeltaEBias);

void ComputeBiasPotentialHarmonic(sys *Sys);
double MeshBiasChangeHarmonic(sys *Sys);
double MeshBiasGlobalHarmonic(sys *Sys, double dzSuggested);
void StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias);

double PrintBiasPotentialAverage(sys *Sys);
double PrintBiasPotentialHarmonic(sys *Sys);
#endif
