#ifndef SYSDEFS_H
#define SYSDEFS_H

#include <stdlib.h>
#include <argp.h>
#include "string.h"

struct arguments{
/* Struct which stores all input arguments. */

// MC properties
  int MCSweeps;
  double MDisp;
  double BiasStr, z0, theta0;
  char BiasType[1024];
  long int RandSeed;
  
// MeshProperties
  double kc, Lx, Ly;
  int    Nx, Ny;
  
// ParticleProperties
  int NumPType;
  char PTypePop[1024], LPType[1024], DPType[1024], BindPType[1024], RanPType[1024];
  
// OutPutProperties
  char OutputDir[1024], Tag[1024]; 
  int framenums, datnums, Histbnums;
  
// GCMCProperties
  char ChemPots[2024];
  
// Restart/Equilibration option
  char Inputfile[1024], Equil[3];
  
// ErrorCheckRun
  char ErrorCheck[1024];
};

static struct argp_option options[] = {
/* The options we understand.
First is the parameter name, then is _a single letter with single quotes_
for the code that is to be used in parse_opt, the variable type description,
and the help message. End it with an empty item {0}.
*/

  {"Lx"        , 'a', "DOUBLE", 0,  "Mesh length in x direction [nm]" },
  {"Ly"        , 'b', "DOUBLE", 0,  "Mesh length in y direction [nm]" },
  {"Nx"        , 'c', "INT"   , 0,  "Meshpoint nums in x direction" },
  {"Ny"        , 'd', "INT"   , 0,  "Meshpoint nums in y direction" },
  {"z0"        , 'e', "DOUBLE", 0,  "Initial height of top membrane [nm]" },
  {"MeshDisp"  , 'f', "DOUBLE", 0,  "Size of Mesh height suggestions [nm]" },
  {"kc"        , 'g', "DOUBLE", 0,  "Bending modulus [kBT]" },
  
  {"NumPType"  , 'h', "INT"   , 0, "Number of protein flavors" },
  {"LPType"    , 'j', "LIST"  , 0, "Length of each protein type" },
  {"DPType"    , 'k', "LIST"  , 0, "Diameter of each protein type" },
  {"BindPType" , 'l', "LIST"  , 0, "Binding strength of each protein type" },
  {"RanPType"  , 'm', "LIST"  , 0, "Range of binding of each protein type" },
  
  {"ChemPots"  , 'n', "LIST"  , 0, "Chemical potential of each species on the membrane, the first [NumPType] entries are assigned to the bottom mesh, the second set to the top" },
  
  {"OutDir"    , 'o', "DIR"   , 0, "Output directory" },
  {"Tag"       , 'p', "STR"   , 0, "Cosmetic tag to add to filenames" },
  
  {"framenums" , 'q', "INT"   , 0, "Number of frames in video" },
  {"datnums"   , 'r', "INT"   , 0, "Number of datapoints in datafile"},
  {"histbnums" , 's', "INT"   , 0, "Number of bins in dz histogram"},
  
  {"Input"     , 't', "FILE"  , 0, "Input filename" },
  {"Equil"     , 'u', "OPT"   , 0, "Equilibration option [YES/NO]"},
  
  {"MCSweeps"  , 'v', "INT"   , 0, "Number of MC sweeps to perform" },
  
  {"BiasStr"   , 'w', "DOUBLE", 0, "Strength of bias applied on mesh distance [kBT]"},
  {"BiasType"  , 'x', "STR"   , 0, "Type of Bias to apply on membranes: Harmonic, Point, or Linear"},
  {"BiasCenter", 'y', "DOUBLE", 0, "Shift parameter of bias potential"},
  
  {"ErrorCheck", 'z', "STR"   , 0, "Check all configurations"},
  {"SEED"      , '1', "LONG INT", 0, "Random seed for rng"},
  {0}
};

void InterpretInputCommands( struct arguments *arguments, int argc, char **argv);
void PrintInputCommands(FILE *fp, struct arguments *arguments);
void ReadInputCommands(FILE *fp, struct arguments *arguments);

#endif
