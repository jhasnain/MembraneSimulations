#include "errorcheck.h"


void DontErrorCheck(sys *Sys, int mciteration, double movechoice, double mrange, double prange, double gcrange){}

int ParseError(int Error){
  char routine[256];
  
       if (Error==0) return 0;
  else if (Error==1) sprintf(routine, "CellLists\n");
  else if (Error==2) sprintf(routine, "PhaseFactor, Height, and nhat\n");
  else if (Error==3) sprintf(routine, "Particle Energies\n");
  else if (Error==4) sprintf(routine, "MeshEnergy\n");
//   else if (Error==5) sprintf(routine, "\n");
  
  printf("ERROR:\nMismatch in computed energies!\nFlag thrown by %s\n", routine);
  return 1;
}

void ErrorCheckConfiguration(sys *Sys, int mciteration, double movechoice, double mrange, double prange, double gcrange){
  char lastmove[256];
  int TotErrs=0;
  if(mciteration>10000){ 
    TotErrs+=ParseError(ErrorCheck_CellLists(Sys));
    TotErrs+=ParseError(ErrorCheck_Z_nhat_PhaseFactor(Sys));
    TotErrs+=ParseError(ErrorCheck_ParticleEnergies(Sys));
    TotErrs+=ParseError(ErrorCheck_MeshEnergies_BiasEnergy(Sys));
    
    if(TotErrs!=0){
      if      ( movechoice < mrange )  sprintf(lastmove, "MeshMove");
      else if ( movechoice < prange )  sprintf(lastmove, "ParticleMove");
      else if ( movechoice < gcrange ) sprintf(lastmove, "GCMCMove");
      printf("Last applied move %s\nMC iteration: %d\nExiting...\n", lastmove, mciteration);
      exit(1);
    }
  }
}
