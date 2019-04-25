#ifndef MESH_H
#define MESH_H

#include "membranesys.h"
#include <math.h>


double ReturnVanillaMeshEnergy(sys *Sys);
double ReturnVanillaBiasEnergy(sys *Sys);

void InitializeMeshEnergies(sys *Sys);
int PerformMeshMove(sys *Sys);

int ErrorCheck_MeshEnergies_BiasEnergy(sys *Sys);

#endif
