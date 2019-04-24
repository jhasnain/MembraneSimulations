#ifndef MESHBIAS_H
#define MESHBIAS_H

#include "membranesys.h"
#include "effectiveparticlepotential.h"

void (*ComputeBiasPotential)();
double (*MeshBiasChange)();
double (*MeshBiasGlobalChange)();
void (*StoreBiasChange)();
double (*PrintBiasPotential)();
int (*ErrorCheck_Bias)();

void ComputeBiasPotentialMinDist(sys *Sys);
double MeshBiasChangeMinDist(sys *Sys);
void StoreBiasChangeMinDist(sys *Sys, double DeltaEBias);

void ComputeBiasPotentialAverage(sys *Sys);
double MeshBiasChangeAverage(sys *Sys);
void StoreBiasChangeAverage(sys *Sys, double DeltaEBias);

void ComputeBiasPotentialHarmonic(sys *Sys);
double MeshBiasChangeHarmonic(sys *Sys);
double MeshBiasGlobalHarmonic(sys *Sys);
void StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias);
int  ErrorCheckHarmonicBias(sys *Sys);

void ComputeBiasPotentialHarmonicPlusVeff(sys *Sys);
double MeshBiasChangeHarmonicPlusVeff(sys *Sys);
double MeshBiasGlobalHarmonicPlusVeff(sys *Sys);
void StoreBiasChangeHarmonicPlusVeff(sys *Sys, double DeltaEBias);
int  ErrorCheckHarmonicPlusVeffBias(sys *Sys);

double PrintBiasPotentialAverage(sys *Sys);
double PrintBiasPotentialHarmonic(sys *Sys);
double PrintBiasPotentialHarmonicPlusVeff(sys *Sys);
#endif
