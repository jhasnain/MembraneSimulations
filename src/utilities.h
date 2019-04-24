#ifndef UTILITIES_H
#define UTILITIES_H

#include "membranesys.h"
#include <string.h>

void PrintConfiguration(FILE *fp, sys *Sys);
void ReadConfiguration(FILE *fp, sys *Sys);
void PrintData(FILE *fp, sys *Sys, long currentMCstep);
void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep);
void PrintVec(vec target, char *note);
void ConfigurationRestart(FILE *fp, sys *Sys);

int computeindexdistsqr(int di, int dj, int Nx, int Ny);

void CopyVec(vec *target, vec source);
double ReturnHeightStd(sys *Sys);
void ComputeHHist(sys *Sys);
void ComputeCurrentCorrelations(sys *Sys);
void UpdateCorrelations(sys *Sys);
void UpdateCompHist(sys *Sys);

#endif
