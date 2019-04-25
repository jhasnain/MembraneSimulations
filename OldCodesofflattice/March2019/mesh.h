#ifndef MESH_H
#define MESH_H

#include "membranesys.h"
#include "particles.h"
#include "meshbiasfuncs.h"
#include <math.h>

double ReturnVanillaMeshEnergy(sys *Sys);
double ReturnVanillaBiasEnergy(sys *Sys);
void ReturnNodeNormal(double *Norm, vec *Target, MeshPoint **MyMesh, int i0, int j0, sys *Sys, int Orient);

void InitializeMeshEnergies(sys *Sys);
int PerformMeshMove(sys *Sys);

int ErrorCheck_Mesh(sys *Sys);
#endif
