#include "meshbiasfuncs.h"

void ComputeBiasPotentialMinDist(sys *Sys){
 double minz0=10000000000, dummy;

  for(int i=0; i<Sys->Nx;i++){
    for(int j=0; j<Sys->Ny;j++){
      dummy = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if (dummy<minz0){
        Sys->MeshBias.imin=i;Sys->MeshBias.jmin=j;
        minz0=dummy;
        Sys->MeshBias.theta=minz0;
      }
    }
  }
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(Sys->MeshBias.theta - Sys->MeshBias.theta0);
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
    deltaEBias=Sys->MeshBias.BiasStr*(*thetadummy - Sys->MeshBias.theta);
  }
  
  else{
    dummy=dzSuggested;
    if( dummy < Sys->MeshBias.theta ){
      *idummy=i0;
      *jdummy=j0;
      *thetadummy=dummy;
      deltaEBias=Sys->MeshBias.BiasStr*( *thetadummy - Sys->MeshBias.theta );
    }
    else{
      deltaEBias=0.0;
      *idummy=Sys->MeshBias.imin;
      *jdummy=Sys->MeshBias.jmin;
      *thetadummy=Sys->MeshBias.theta;
    }
  }

  return deltaEBias;
}

void StoreBiasChangeMinDist(sys *Sys, double DeltaEBias){
  Sys->MeshBias.imin = Sys->MeshBias.idummy;
  Sys->MeshBias.jmin = Sys->MeshBias.jdummy;
  Sys->MeshBias.theta= Sys->MeshBias.thetadummy;
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
  Sys->MeshBias.theta = zave;
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0);
}

double MeshBiasChangeAverage(sys *Sys){
  double dzSuggested = (Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z: 
                       Sys->BotMesh[Sys->i0][Sys->j0].Pos.z - Sys->suggestedz;
  
  Sys->MeshBias.thetadummy = Sys->MeshBias.theta + dzSuggested;
  return Sys->MeshBias.BiasStr*dzSuggested;      
}

void StoreBiasChangeAverage(sys *Sys, double DeltaEBias){
  Sys->MeshBias.theta= Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}


double PrintBiasPotentialHarmonic(sys *Sys){
  double zave=0.0;
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  zave/=Sys->MTot;
  return Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0)*(zave - Sys->MeshBias.theta0);
}

void ComputeBiasPotentialHarmonic(sys *Sys){
  double zave=0.0;
  
  for(int i = 0;i<Sys->Nx;i++){
    for(int j = 0;j<Sys->Ny;j++){
      zave+=(Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z);
    }
  }
  zave/=Sys->MTot;
  Sys->MeshBias.theta = zave - Sys->MeshBias.theta0;
  Sys->BiasEnergy=Sys->MeshBias.BiasStr*(zave - Sys->MeshBias.theta0)*(zave - Sys->MeshBias.theta0);
}

double MeshBiasChangeHarmonic(sys *Sys){
  double dzSuggested =(Sys->choice == TOP) ? 
                       Sys->suggestedz - Sys->TopMesh[Sys->i0][Sys->j0].Pos.z: 
                       Sys->BotMesh[Sys->i0][Sys->j0].Pos.z - Sys->suggestedz;
  dzSuggested/=Sys->MTot;        
  Sys->MeshBias.thetadummy = Sys->MeshBias.theta + dzSuggested;
  return Sys->MeshBias.BiasStr*(dzSuggested + 2.0*Sys->MeshBias.theta)*dzSuggested;
}

double MeshBiasGlobalHarmonic(sys *Sys, double dzSuggested){
  Sys->MeshBias.thetadummy = Sys->MeshBias.theta + dzSuggested;
  return Sys->MeshBias.BiasStr*(dzSuggested + 2.0*Sys->MeshBias.theta)*dzSuggested;
}

void StoreBiasChangeHarmonic(sys *Sys, double DeltaEBias){
  Sys->MeshBias.theta= Sys->MeshBias.thetadummy;
  Sys->BiasEnergy+=DeltaEBias;
}
