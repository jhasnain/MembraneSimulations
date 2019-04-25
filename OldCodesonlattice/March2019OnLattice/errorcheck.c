#include "errorcheck.h"


void DontErrorCheck(sys *Sys, int mciteration, double movechoice){}

void MoveInfo(char *lastmove, sys *Sys, int acc){
  if(strcmp(lastmove, "MeshMove Meshpoint")==0){
    MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
    printf("MeshPoint Info:\n");
    printf("  %s\n", (acc > 1e-6) ? "Acc" : "Rej");
    printf("  MeshChoice %s\n", (Sys->choice == TOP) ? "Top" : "Bot");
    printf("  Indices    %d %d\n", Sys->i0, Sys->j0);
    printf("  Position   %lf %lf %lf\n", MeshToModify[Sys->i0][Sys->j0].Pos.x, MeshToModify[Sys->i0][Sys->j0].Pos.y, MeshToModify[Sys->i0][Sys->j0].Pos.z);
    printf("  dz         %lf\n", Sys->suggestedz);
  }
  
  if(strcmp(lastmove, "MeshMove MeshShift")==0){
    printf("Move Info:\n");
    printf("  %s\n", (acc > 1e-6) ? "Acc" : "Rej");
    printf("  MeshChoice %s\n", (Sys->choice == TOP) ? "Top" : "Bot");
    printf("  dz         %lf\n", Sys->suggestedz);
    printf("\n");
  }
}

int ParseError(int Error){
  char routine[256];
       if (Error==0) return 0;
  else if (Error==1) sprintf(routine, "MeshEnergy\n");
  else if (Error==2) sprintf(routine, "BiasEnergy\n");
  else if (Error==3) sprintf(routine, "Veff\n");
  printf("\nERROR: Flag thrown by %s", routine);
  return 1;
}

void ErrorCheckConfiguration(sys *Sys, int mciteration, int acc){
  char lastmove[256];
  int TotErrs=0;
  if(mciteration>0){ 
    TotErrs+=ParseError(ErrorCheck_Mesh(Sys));
    TotErrs+=ParseError(ErrorCheck_Bias(Sys));
    TotErrs+=ParseError(ErrorCheck_Veff(Sys));
    
    if(TotErrs>0){  
      if(Sys->MeshShift==0) sprintf(lastmove, "MeshMove Meshpoint");
      else                  sprintf(lastmove, "MeshMove MeshShift");
      printf("Last applied move %s\n", lastmove);
      MoveInfo(lastmove, Sys, acc);
      printf("MC iteration: %d\nExiting...\n", mciteration);
      exit(1);
    }
  }
}
