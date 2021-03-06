#include "meshbiasfuncs.h"

void ComputeBiasPotentialMinDist(sys *Sys){
 double minz0=10000000000, dummy;

  for(int i=0; i<Sys->Nx;i++){
    for(int j=0; j<Sys->Ny;j++){
      dummy = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if (dummy<minz0){
        Sys->MeshBias.imin=i;Sys->MeshBias.jmin=j;
        minz0=dummy;
        Sys->MeshBias.thetaTop=minz0;
      }
    }
  }
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(Sys->MeshBias.thetaTop - Sys->MeshBias.theta0);
}

double MeshBiasChangeMinDist(sys *Sys){
  double dummy, deltaEBias;
  int *idummy = &(Sys->MeshBias.idummy), i0=Sys->i0;
  int *jdummy = &(Sys->MeshBias.jdummy), j0=Sys->j0;
  double *thetadummy = &(Sys->MeshBias.thetadummy);
  
  double dzSuggested = (Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->BotMesh[i0][j0].Pos.z: 
                       Sys->TopMesh[i0][j0].Pos.z - Sys->suggestedz;
  
  if(i0 == Sys->MeshBias.imin && j0 == Sys->MeshBias.jmin){
    *thetadummy=dzSuggested;*idummy=i0;*jdummy=j0;
    for(int i=0; i<Sys->Nx;i++){
      for(int j=0; j<Sys->Ny;j++){
        dummy = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
        if (dummy<*thetadummy && i!=i0 && j!=j0){
          *idummy=i;
          *jdummy=j;
          *thetadummy=dummy;
        }
      }
    }
    deltaEBias=Sys->MeshBias.BiasStr*(*thetadummy - Sys->MeshBias.thetaTop);
  }
  
  else{
    dummy=dzSuggested;
    if( dummy < Sys->MeshBias.thetaTop ){
      *idummy=i0;
      *jdummy=j0;
      *thetadummy=dummy;
      deltaEBias=Sys->MeshBias.BiasStr*( *thetadummy - Sys->MeshBias.thetaTop );
    }
    else{
      deltaEBias=0.0;
      *idummy=Sys->MeshBias.imin;
      *jdummy=Sys->MeshBias.jmin;
      *thetadummy=Sys->MeshBias.thetaTop;
    }
  }
  return deltaEBias;
}

void StoreBiasChangeMinDist(sys *Sys, double DeltaEBias){
  Sys->MeshBias.imin = Sys->MeshBias.idummy;
  Sys->MeshBias.jmin = Sys->MeshBias.jdummy;
  Sys->MeshBias.thetaTop= Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}


double PrintBiasPotentialAverage(sys *Sys){
  double zave=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  return Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0);
}

void ComputeBiasPotentialAverage(sys *Sys){
  double zave=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  Sys->MeshBias.thetaTop = zave;
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0);
}

double MeshBiasChangeAverage(sys *Sys){
  double dzSuggested = (Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z: 
                       Sys->BotMesh[Sys->i0][Sys->j0].Pos.z - Sys->suggestedz;
  
  Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
  return Sys->MeshBias.BiasStr*dzSuggested;      
}

void StoreBiasChangeAverage(sys *Sys, double DeltaEBias){
  Sys->MeshBias.thetaTop= Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}


int  ErrorCheckHarmonicBias(sys *Sys){
  int Err=0;
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  if(fabs(ztop-Sys->MeshBias.thetaTop)>1e-6){
    printf("displacement of top is %lf stored %lf\ntheta0: %lf\n", ztop, Sys->MeshBias.thetaTop, Sys->MeshBias.theta0);
    Err=6; 
  }
  if(fabs(zbot-Sys->MeshBias.thetaBot)>1e-6){
    printf("displacement of bot is %lf stored %lf\ntheta0: %lf\n", zbot, Sys->MeshBias.thetaBot,-Sys->MeshBias.theta0);
    Err=6; 
  }
  double En = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
  if(fabs(En-Sys->BiasEnergy)>1e-6){
    printf("Bias energy off computed %lf stored %lf\n", En, Sys->BiasEnergy);
    Err=6;
  }
  return Err;
}

double PrintBiasPotentialHarmonic(sys *Sys){
  double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  return Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

void ComputeBiasPotentialHarmonic(sys *Sys){
 double ztop=0.0, zbot=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      ztop += Sys->TopMesh[i][j].Pos.z/Sys->MTot;
      zbot += Sys->BotMesh[i][j].Pos.z/Sys->MTot;
    }
  }
  ztop-=Sys->MeshBias.theta0;
  zbot+=Sys->MeshBias.theta0;
  
  Sys->MeshBias.thetaTop = ztop;
  Sys->MeshBias.thetaBot = zbot;
  Sys->BiasEnergy        = Sys->MeshBias.BiasStr*(zbot*zbot + ztop*ztop);
}

double MeshBiasChangeHarmonic(sys *Sys){
  double dzSuggested;
  if(Sys->choice == TOP){
    dzSuggested = (Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaTop);
  }
  else{
    dzSuggested = (Sys->suggestedz - Sys->BotMesh[Sys->i0][Sys->j0].Pos.z)/Sys->MTot;
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot);
  }
}

double MeshBiasGlobalHarmonic(sys *Sys, double dzSuggested){
  if(Sys->choice == 0){
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaBot + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaBot);
  }
  else{
    Sys->MeshBias.thetadummy = Sys->MeshBias.thetaTop + dzSuggested;
    return Sys->MeshBias.BiasStr*dzSuggested*(dzSuggested + 2.0*Sys->MeshBias.thetaTop);
  }
}

void StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias, int choice){
  if (choice == 0) Sys->MeshBias.thetaBot = Sys->MeshBias.thetadummy;
  else             Sys->MeshBias.thetaTop = Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}
