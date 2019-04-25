#ifndef ERROR_H
#define ERROR_H

#include "mesh.h"
#include "membranesys.h"
#include "effectiveparticlepotential.h"
#include <math.h>

void (*ErrorCheckIntermittent)();
void (*ErrorCheckEvery)();

void ErrorCheckConfiguration(sys *Sys, int mciteration, int acc);
void DontErrorCheck(sys *Sys, int mciteration, double movechoice);



#endif
