#include <stdio.h>	
#include <stdlib.h>
#include <argp.h>
#include <math.h>

#include "parseparameters.h"
#include "initialize.h"
#include "membranesys.h"
#include "utilities.h"
#include "mesh.h"
#include "meshbiasfuncs.h"
#include "errorcheck.h"
#include "effectiveparticlepotential.h"

void ConductMC(sys *Sys, struct arguments *args){
//   Extract MC vars for shorter code
  long MCStepNums      = Sys->MCSteps; 
  long *ptr_SweepIndex = &(Sys->SweepIndex);
  long datarate        = Sys->datarate;
  long vidrate         = Sys->framerate;
  long printfrate      = max(1, Sys->SweepNum/500);
  long corrate         = Sys->SweepNum/100;
   
//  Helper vars
  long step;
  double acc;
  
//  Initialize all energies and lists
  InitializeMeshEnergies(Sys);
  Sys->TotEnergy=Sys->MeshEnergy + Sys->BiasEnergy + Sys->TotVeff;
  
//  Perform MC simulation
  for (step=0;step<MCStepNums;step++){
    acc=PerformMeshMove(Sys);
    
    ErrorCheckEvery(Sys, step, acc);
    Sys->AcceptanceRate += acc;
    if(Sys->MeshShift == 1) Sys->ShiftAccept+=acc;
    else Sys->NodeAccept+=acc;
    Sys->ShiftAttempt+=Sys->MeshShift;
      
    if(step % Sys->DoF == 0){
      if(*ptr_SweepIndex%printfrate==0 || *ptr_SweepIndex%datarate==0 ){ComputeParticlePops(Sys);}
      ErrorCheckIntermittent(Sys, step, acc);
      if(*ptr_SweepIndex%printfrate==0) {PrintDataToScreen(stdout, Sys, step);WriteRestart(Sys, args);}
      if(*ptr_SweepIndex%datarate==0)    PrintData(Sys->fileptr_data, Sys, step);
      if(*ptr_SweepIndex%vidrate==0)     PrintConfiguration(Sys->fileptr_traj, Sys);
      if(corrate!=0 && *ptr_SweepIndex%corrate==0) {UpdateCompHist(Sys);UpdateCorrelations(Sys);}
      (*ptr_SweepIndex)++;
    }
  }
  
  ComputeParticlePops(Sys);
  PrintDataToScreen(stdout, Sys, step);
  PrintData(Sys->fileptr_data, Sys, step);
  PrintConfiguration(Sys->fileptr_traj, Sys);
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
    
   


