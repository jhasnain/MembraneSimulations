#ifndef MESH_H
#define MESH_H

#include "membranesys.h"
#include "utilities.h"
#include "meshbiasfuncs.h"
#include <math.h>

double ReturnVanillaMeshEnergy(sys *Sys);
double ReturnVanillaBiasEnergy(sys *Sys);

void ReturnNodeNormal(double *Norm, vec *Target, MeshPoint **MyMesh, int i0, int j0, sys *Sys, int Orient);
void MeshPBC(int *left, int *right, int *out, int *in, sys *Sys, int i, int j);
void InitializeMeshEnergies(sys *Sys);

int PerformMeshMove(sys *Sys);
int ErrorCheck_Mesh(sys *Sys);


#endif
