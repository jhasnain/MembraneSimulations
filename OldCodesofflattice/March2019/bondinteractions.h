#ifndef BONDS_H
#define BONDS_H

#include "membranesys.h"
#include "utilities.h"
#include <math.h>

void InitializeBondFunctions(sys *Sys);

int (*ErrorCheck_Bonds)(sys *Sys);
void (*ComputeAllParticleBonds)(sys *Sys);
double (*BondParticleMove)(sys *Sys);
double (*BondMeshMove)(sys *Sys);
double (*BondMeshShift)(sys *Sys, double suggestedz);
double (*BondGCMC)(sys *Sys);
void (*BondStoreGCMC)(sys *Sys);
void (*BondMove)(sys *Sys);

void (*UpdateBondLists)(sys *Sys, int opt);
void RearrangeBonds(sys *Sys);
void PrintBonds(sys *Sys);
void PrintBondList(particlebonds Bonds, int p);
#endif
