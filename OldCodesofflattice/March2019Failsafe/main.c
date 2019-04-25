#include <stdio.h>	
#include <stdlib.h>
#include <argp.h>
#include <math.h>

#include "parseparameters.h"
#include "initialize.h"
#include "membranesys.h"
#include "utilities.h"
#include "mesh.h"
#include "particles.h"
#include "gcmc.h"
#include "meshbiasfuncs.h"
#include "errorcheck.h"


void ConductMC(sys *Sys, struct arguments *args){
//   Extract MC vars for shorter code
  long MCStepNums    = Sys->MCSteps; 
  long *ptr_SweepIndex=&(Sys->SweepIndex);
  long datarate       = Sys->datarate;
  long vidrate        = Sys->framerate;
  long printfrate     = max(1, Sys->SweepNum/500);
   
  gsl_rng *RNG = Sys->rng;
  double mrange =Sys->MRate;
  double prange =Sys->PRate  + mrange;
  double gcrange=Sys->GCRate + prange;
  
//  Helper vars
  long step;
  double movechoice, accM, accP, accGC;
  
//  Initialize all energies and lists
  InitializeMeshEnergies(Sys);
  AssignCells(Sys);
  InitializeParticleEnergies(Sys);
  Sys->TotEnergy=Sys->TotParticleEnergy + Sys->TotPhaseFactorEnergy + Sys->MeshEnergy + Sys->BiasEnergy + Sys->TotMeshParticleEnergy;
  
//  Perform MC simulation
  for (step=0;step<MCStepNums;step++){
    movechoice=gsl_rng_uniform(RNG);
    accM=0.0;accP=0.0;accGC=0.0;
    if (movechoice < mrange){
      accM=PerformMeshMove(Sys);
      Sys->MeshAttempt++;
    }
    else if ( movechoice < prange ){
      accP=PerformParticleMove(Sys);
      Sys->ParticleAttempt++;
    }
    else if ( movechoice < gcrange ){
      accGC=PerformGCMCMove(Sys);
      Sys->GCAttempt++;
    }

    BondMove(Sys);
    ErrorCheckEvery(Sys, step, movechoice, mrange, prange, gcrange, accM, accP, accGC);
    
    Sys->AcceptMesh     += accM;
    Sys->AcceptPart     += accP;
    Sys->AcceptGC       += accGC;
    Sys->AcceptanceRate += accM + accP + accGC;
    
    if(step % Sys->DoF == 0){
      ErrorCheckIntermittent(Sys, step, movechoice, mrange, prange, gcrange, accM, accP, accGC);
      if(*ptr_SweepIndex%printfrate==0) {PrintDataToScreen(stdout, Sys, step);WriteRestart(Sys, args);}
      if(*ptr_SweepIndex%datarate==0)   PrintData(Sys->fileptr_data, Sys, step);
      if(*ptr_SweepIndex%vidrate==0)    PrintConfiguration(Sys->fileptr_traj, Sys, 1);
      (*ptr_SweepIndex)++;
    }
  }
  
  PrintDataToScreen(stdout, Sys, step);
  PrintData(Sys->fileptr_data, Sys, step);
  PrintConfiguration(Sys->fileptr_traj, Sys, 1);
  
  Sys->AcceptanceRate /= MCStepNums;
  Sys->AcceptMesh     /= Sys->MeshAttempt;
  Sys->AcceptPart     /= Sys->ParticleAttempt;
  Sys->AcceptGC       /= Sys->GCAttempt;
}

int main(int argc, char **argv){
  struct arguments arguments;
  sys Sys;
  InterpretInputCommands(&arguments, argc, argv);
  InitializeStructs(&Sys, &arguments);
  PrintInputCommands(Sys.fileptr_data, &arguments);
  ConductMC(&Sys, &arguments);
  PrintInputCommands(stdout, &arguments);
}
    
   


