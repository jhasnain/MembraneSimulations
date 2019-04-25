#ifndef PARTIC_H
#define PARTIC_H

#include "membranesys.h"
#include "utilities.h"
#include "math.h"

void InitializeParticleEnergies(sys *Sys);
void AssignCells(sys *Sys);
int PerformParticleMove(sys *Sys);

double CheckOverlap(particle Moved, particle Partner, double Lx, double Ly, int *loop);
double CheckForBond(particle MovedParticle, particle  PartnerParticle, double Lx, double Ly, sys *Sys);
double ComputeDeltaMPEnergy(sys *Sys, particle *MyParticle, int MyParticleIndex);
double MeshParticleInteraction(double dRsqr);

double CollisionTestAndParticleBonds(sys *Sys);
void Compute_Z_nhat_PhaseFactor(sys *Sys, particle *My_Particle);
double ParticleMeshEnergy(sys *Sys, particle *MyParticle);
void ComputeMeshParticleEnergy(sys *Sys);

double PrintAllParticleEnergies(sys *Sys);
double PrintParticleHeight(sys *Sys, int p);

double ParticleEnergy(sys *Sys, particle TestParticle);
double ReturnTotalParticleEnergies(sys *Sys);

int ErrorCheck_CellLists(sys *Sys);
int ErrorCheck_Z_nhat_PhaseFactor(sys *Sys);
int ErrorCheck_ParticleEnergies(sys *Sys);

#endif
