#include "utilities.h"
#include "meshbiasfuncs.h"
#include "mesh.h"

void CopyVec(vec *target, vec source){
 target->x=source.x; 
 target->y=source.y;
 target->z=source.z; 
}

double ReturnHeightStd(sys *Sys){
  double sigmasqr=0.0, kx, ky;
  int i,j;
  for(i=0;i<Sys->Nx;i++){
    kx=2.0*M_PI*(i-Sys->Nx/2)/Sys->Lx;
    for(j=0;j<Sys->Ny;j++){
      ky=2.0*M_PI*(j-Sys->Ny/2)/Sys->Ly;
    if(kx*kx + ky*ky>0.0)sigmasqr += 1.0/(Sys->A*Sys->kc*(kx*kx + ky*ky)*(kx*kx + ky*ky));
    }
  }
  return sqrt(sigmasqr);
}

double ReturnTotalArea(sys *Sys, int choice){
  MeshPoint **MyMesh = (choice == 0 ? Sys->BotMesh: Sys->TopMesh);
  double Area=0.0;
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      Area+=MyMesh[i][j].nNorm;
    }
  }
  return Area*Sys->dx*Sys->dy;
}

void GetHHist(sys *Sys){
  hist *Hhist=&(Sys->Hhist);
  int i, j, entry;
  double deltaH;
  for(i=0;i<Hhist->binnums;i++)Hhist->Hhist[1][i]=0.0;
  if(Sys->SweepIndex>0){
    for(i=0;i<Sys->Nx;i++){
      for(j=0;j<Sys->Ny;j++){
        deltaH= Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
        entry = (int)round( (deltaH-Hhist->zmin)/Hhist->dz);
        if(entry<0)entry=0;
        if(entry>Hhist->binnums-1)entry=Hhist->binnums-1;
        Hhist->Hhist[1][entry]+=1.0/Hhist->histnorm;
      }
    }
  }
}

void UpdateCompHist(sys *Sys){
  hist *Hhist=&(Sys->Hhist);
  int samplenums = Sys->SweepIndex/Sys->datarate;
  int i;
  if(Sys->SweepIndex>0){
    for(i=0;i<Hhist->binnums;i++)Hhist->CompHhist[1][i]+=Hhist->Hhist[1][i];
    FILE *fp=fopen(Sys->Histfile, "w");
    fprintf(fp, "#%d\n", samplenums);
    for(i=0;i<Hhist->binnums;i++)fprintf(fp, "%lf %lf\n", Hhist->CompHhist[0][i], Hhist->CompHhist[1][i]/samplenums);
    fclose(fp);
  }
}

void PrintVec(vec target, char *note){printf("%s: %lf %lf %lf\n", note, target.x, target.y, target.z);}

void PrintConfiguration(FILE *fp, sys *Sys){
  int i, j;
  
  fprintf(fp, "%ld\nSweep no: %ld out of %ld\n", Sys->DoF, Sys->SweepIndex, Sys->SweepNum);
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      Sys->TopMesh[i][j].Pos.x,  Sys->TopMesh[i][j].Pos.y,  Sys->TopMesh[i][j].Pos.z,
      Sys->TopMesh[i][j].nhat.x, Sys->TopMesh[i][j].nhat.y, Sys->TopMesh[i][j].nhat.z,
      Sys->TopMesh[i][j].lc_x, Sys->TopMesh[i][j].lc_y
    );
   }
  }
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      Sys->BotMesh[i][j].Pos.x,  Sys->BotMesh[i][j].Pos.y,  Sys->BotMesh[i][j].Pos.z,
      Sys->BotMesh[i][j].nhat.x, Sys->BotMesh[i][j].nhat.y, Sys->BotMesh[i][j].nhat.z,
      Sys->BotMesh[i][j].lc_x, Sys->BotMesh[i][j].lc_y
    );
   }
  }
}

void PrintData(FILE *fp, sys *Sys, long currentMCstep){
  int i, shift = 13;
  
  if (currentMCstep==0) {
    fprintf(fp, "#Sweep(1) TAcc(2) TotE(3) MeshE(4) BiasE(5) Veff(6) ThetaValBot(7) ThetaValTop(8) TotAreaBot(9) TotAreaTop(10) TotPop(11) TotPopBot(12) TotPopTop(13) ");
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "NBot%d(%d) ", i+1, shift+1+i);
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "NTop%d(%d) ", i+1,  shift + 2 +   Sys->NumPTypes + i);
    for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "Bonds%d(%d) ", i+1, shift + 2 + 2*Sys->NumPTypes + i);
    fprintf(fp, "\n");
  }
  
  fprintf( fp, "%ld %.2lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
    Sys->SweepIndex, Sys->AcceptanceRate/(currentMCstep+1), 
    Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotVeff,
    Sys->MeshBias.thetaBot, Sys->MeshBias.thetaTop, ReturnTotalArea(Sys, 0), ReturnTotalArea(Sys, 1),
    Sys->TotCurPop[0] + Sys->TotCurPop[1], Sys->TotCurPop[0], Sys->TotCurPop[1]
  );
  
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->CurParPop[0][i]);
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->CurParPop[1][i]);
  for(i = 0; i<Sys->NumPTypes; i++) fprintf(fp, "%lf ", Sys->BondNums[i]);
  fprintf(fp,"\n");
  
  UpdateCompHist(Sys);
}

void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep){
  if (currentMCstep==0) fprintf(fp, "Sweep(1) TAcc(2) TotE(3) MeshE(4) BiasE(5) Veff(6) TotPopBot(7) TotPopTop(8)\n");
  
  fprintf( fp,
    "%ld %.2lf %lf %lf %lf %lf %.0lf %.0lf\n", 
    Sys->SweepIndex, 
    Sys->AcceptanceRate/(currentMCstep+1), Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotVeff, Sys->TotCurPop[0], Sys->TotCurPop[1]
  );
}








