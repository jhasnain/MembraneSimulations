#ifndef ERROR_H
#define ERROR_H

#include "mesh.h"
#include "particles.h"
#include "membranesys.h"
#include "bondinteractions.h"
#include <math.h>

void (*ErrorCheckIntermittent)();
void (*ErrorCheckEvery)();

void ErrorCheckConfiguration(sys *Sys, int mciteration, double movechoice, 
                             double mrange, double prange, double gcrange,
                             double accM, double accP, double accG);
void DontErrorCheck(sys *Sys, int mciteration, double movechoice, 
                             double mrange, double prange, double gcrange,
                             double accM, double accP, double accG);

#endif
