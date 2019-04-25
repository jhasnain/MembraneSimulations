#ifndef ERROR_H
#define ERROR_H

#include "mesh.h"
#include "particles.h"
#include "membranesys.h"
#include <math.h>

void (*ErrorCheck)();
void ErrorCheckConfiguration(sys *Sys, int mciteration, double movechoice, double mrange, double prange, double gcrange);
void DontErrorCheck(sys *Sys, int mciteration, double movechoice, double mrange, double prange, double gcrange);

#endif
