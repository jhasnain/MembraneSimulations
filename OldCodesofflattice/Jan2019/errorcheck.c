#include "errorcheck.h"


void DontErrorCheck(sys *Sys, int mciteration, double movechoice, 
                             double mrange, double prange, double gcrange,
                             double accM, double accP, double accG){}

void MoveInfo(char *lastmove, sys *Sys, double accM, double accP, double accG){
  if(strcmp(lastmove, "MeshMove Meshpoint")==0){
    MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
    printf("MeshPoint Info:\n");
    printf("  %s\n", (accM > 1e-6) ? "Acc" : "Rej");
    printf("  MeshChoice %s\n", (Sys->choice == TOP) ? "Top" : "Bot");
    printf("  Indices    %d %d\n", Sys->i0, Sys->j0);
    printf("  Position   %lf %lf %lf\n", MeshToModify[Sys->i0][Sys->j0].Pos.x, MeshToModify[Sys->i0][Sys->j0].Pos.y, MeshToModify[Sys->i0][Sys->j0].Pos.z);
    printf("  dz         %lf\n", Sys->suggestedz);
    
    printf("  Affected Num %d\n", Sys->AffectedParticleNum);
    for(int i=0;i<Sys->AffectedParticleNum;i++){
      printf("    Id           %d\n", Sys->Unique_dummy_index[Sys->AffectedParticleNum]);
    }
    printf("\n");
  }
  
  if(strcmp(lastmove, "MeshMove MeshShift")==0){
    printf("MeshPoint Info:\n");
    printf("  %s\n", (accM > 1e-6) ? "Acc" : "Rej");
    printf("  MeshChoice %s\n", (Sys->choice == TOP) ? "Top" : "Bot");
    printf("  dz         %lf\n", Sys->suggestedz);
    printf("\n");
  }
  
  if(strcmp(lastmove, "ParticleMove")==0){
    printf("PMove Info:\n");
    printf("  %s\n", (accP > 1e-6) ? "Acc" : "Rej");
    printf("  Index    %d\n", Sys->p0);
    printf("  Base     %lf %lf %lf\n", Sys->Particles[Sys->p0].Pos.x, Sys->Particles[Sys->p0].Pos.y, Sys->Particles[Sys->p0].Pos.z);
    printf("  MidPoint %lf %lf %lf\n", Sys->Particles[Sys->p0].MidPoint.x, Sys->Particles[Sys->p0].MidPoint.y, Sys->Particles[Sys->p0].MidPoint.z);
    printf("\n");
  }
  
  if(strcmp(lastmove, "GCMCMove AddParticle")==0){
    printf("AddParticle Info:\n");
    printf("  %s\n", (accG > 1e-6) ? "Acc" : "Rej");
    printf("  Index    %d\n", Sys->p0);
    printf("  Base     %lf %lf %lf\n", Sys->Particles[Sys->p0].Pos.x, Sys->Particles[Sys->p0].Pos.y, Sys->Particles[Sys->p0].Pos.z);
    printf("  MidPoint %lf %lf %lf\n", Sys->Particles[Sys->p0].MidPoint.x, Sys->Particles[Sys->p0].MidPoint.y, Sys->Particles[Sys->p0].MidPoint.z);
    printf("\n");
  }
  
  if(strcmp(lastmove, "GCMCMove RemoveParticle")==0){
    printf("RemoveParticle Info:\n");
    printf("  %s\n", (accG > 1e-6) ? "Acc" : "Rej");
    printf("  Index    %d\n", Sys->p0);
    printf("  BondBuddy %d\n", Sys->BondList.ParticleBonds[Sys->p0].currbond);
    printf("  MidPoint %lf %lf %lf\n", Sys->Particles[Sys->p0].MidPoint.x, Sys->Particles[Sys->p0].MidPoint.y, Sys->Particles[Sys->p0].MidPoint.z);
    printf("\n");
  }
  
  
}

int ParseError(int Error){
  char routine[256];
       if (Error==0) return 0;
  else if (Error==1) sprintf(routine, "CellLists\n");
  else if (Error==2) sprintf(routine, "Particle phaseFactor, Height, and nhat\n");
  else if (Error==3) sprintf(routine, "Particle Overlap\n");
  else if (Error==4) sprintf(routine, "MeshEnergy\n");
  else if (Error==5) sprintf(routine, "Bonds List\n");
  else if (Error==6) sprintf(routine, "BiasEnergy\n");
  else if (Error==7) sprintf(routine, "Particle properties\n");
  else if (Error==8) sprintf(routine, "Mesh Particle\n");
  printf("\nERROR: Flag thrown by %s", routine);
  return 1;
}

void ErrorCheckConfiguration(sys *Sys, int mciteration, double movechoice, 
                             double mrange, double prange, double gcrange,
                             double accM, double accP, double accG){
  char lastmove[256];
  int TotErrs=0;
  if(mciteration>0){ 
    TotErrs+=ParseError(ErrorCheck_CellLists(Sys));
    TotErrs+=ParseError(ErrorCheck_Z_nhat_PhaseFactor(Sys));
    TotErrs+=ParseError(ErrorCheck_ParticleOverlap(Sys));
    TotErrs+=ParseError(ErrorCheck_Mesh(Sys));
    TotErrs+=ParseError(ErrorCheck_Bias(Sys));
    TotErrs+=ParseError(ErrorCheck_Bonds(Sys));
    TotErrs+=ParseError(ErrorCheck_ParticleProps(Sys));
    TotErrs+=ParseError(ErrorCheck_ParticleMesh(Sys));
    
    if(TotErrs!=0){
      if      ( movechoice < mrange ){
        if(Sys->MeshShift==0) sprintf(lastmove, "MeshMove Meshpoint");
        else                  sprintf(lastmove, "MeshMove MeshShift");
      }
      else if ( movechoice < prange ){
                              sprintf(lastmove, "ParticleMove");
      }
      else if ( movechoice < gcrange ){
        if(Sys->GC_AddRemove==0) sprintf(lastmove, "GCMCMove AddParticle");
        else                     sprintf(lastmove, "GCMCMove RemoveParticle");
      }
      printf("Number of errors %d\n\nLast applied move %s\n", TotErrs, lastmove);
      MoveInfo(lastmove, Sys, accM, accP, accG);
      printf("MC iteration: %d\nExiting...\n", mciteration);
      exit(1);
    }
  }
}
