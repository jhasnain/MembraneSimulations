#ifndef INIT_H
#define INIT_H

#include "parseparameters.h"
#include "membranesys.h"
#include "meshbiasfuncs.h"
#include "errorcheck.h"
#include "bondinteractions.h"
#include "utilities.h"

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <dirent.h>

void InitializeStructs(sys *Sys, struct arguments *arguments);
void WriteRestart(sys *Sys, struct arguments *arguments);
#endif
