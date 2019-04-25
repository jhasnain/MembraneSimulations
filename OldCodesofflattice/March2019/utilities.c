#include "utilities.h"
#include "particles.h"
#include "meshbiasfuncs.h"
#include "mesh.h"

void twodvecdiff(vec *target, vec v2, vec v1){ target->x = v2.x-v1.x;target->y = v2.y-v1.y;target->z=0.0;}
void twodvecsum(vec *target, vec v2, vec v1){ target->x = v2.x+v1.x;target->y = v2.y+v1.y;target->z=0.0;}
double twod_NormSqr(vec v){return v.x*v.x + v.y*v.y;}

void threedvecdiff(vec *target, vec v2, vec v1){ target->x = v2.x-v1.x;target->y = v2.y-v1.y;target->z=v2.z-v1.z;}
void threedvecsum(vec *target, vec v2, vec v1){ target->x = v2.x+v1.x;target->y = v2.y+v1.y;target->z=v2.z+v1.z;}
void threedvecsumscalar(vec *target, vec v2, vec v1, double a){ target->x = v2.x+a*v1.x;target->y = v2.y+a*v1.y;target->z=v2.z+a*v1.z;}
void rescalevec(vec *target, double a){ target->x*=a;target->y*=a;target->z*=a;}
double threed_NormSqr(vec v){return v.x*v.x + v.y*v.y + v.z*v.z;}

void Position_To_Index(int *i, double x, double dx){*i=(int)round(x/dx);}
void Index_To_Position(double *x, int i, double dx){*x=i*dx;}
void Distance_To_IndexRan(int *IR, double R, double dr){*IR=(int)ceil(R/dr);}

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

void ShiftVec(vec *target, vec *source, double x, double y, double z){
  target->x = source->x + x;
  target->y = source->y + y;
  target->z = source->z + z;
}

void TranslateParticle(particle *Target, particle *Source, double x, double y, double z){
  ShiftVec(&(Target->Pos),      &(Source->Pos),      x, y ,z);
  ShiftVec(&(Target->MidPoint), &(Source->MidPoint), x, y ,z);
  ShiftVec(&(Target->EndPoint), &(Source->EndPoint), x, y ,z);
}

void Position_To_Index_PBC(int *i, double x, double dx, int maxx){
  Position_To_Index(i, x, dx);*i = IndexPBC(*i, maxx);
}

void threed_NormalizeVec(vec *target){
  double tmp = sqrt(threed_NormSqr(*target));
  target->x/=tmp;
  target->y/=tmp;
  target->z/=tmp;
}
void crossprod(vec *target, vec v1, vec v2){ 
  target->x = v1.y*v2.z - v1.z*v2.y;
  target->y = v1.z*v2.x - v1.x*v2.z;
  target->z = v1.x*v2.y - v1.y*v2.x;
}
double dotprod(vec v1, vec v2){return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;}

void CopyVec(vec *target, vec source){
 target->x=source.x; 
 target->y=source.y;
 target->z=source.z; 
}

void NearestImageConv(vec *myvec, double Lx, double Ly ){
  
  if(myvec->x<-Lx/2)myvec->x+=Lx;
  if(myvec->y<-Ly/2)myvec->y+=Ly;
  
  if(myvec->x> Lx/2)myvec->x-=Lx;
  if(myvec->y> Ly/2)myvec->y-=Ly;
}

void ScalarNearestImageConv(double *val, double L ){
  if(*val<-L/2)*val+=L;
  if(*val> L/2)*val-=L;
}

int IndexPBC(int index, int Nmax){
  if(index<0    )index+=Nmax;
  if(index>=Nmax)index-=Nmax;
  return index;
}

void PosPBC(vec *myvec, double Lx, double Ly){
  
  if(myvec->x<0.0)myvec->x+=Lx;
  if(myvec->y<0.0)myvec->y+=Ly;
  
  if(myvec->x>=Lx)myvec->x-=Lx;
  if(myvec->y>=Ly)myvec->y-=Ly;
}

double deltaRsqr(vec v1, vec v2, double Lx, double Ly){
  vec diff;
  threedvecdiff(&diff, v1, v2);
  NearestImageConv(&diff, Lx, Ly); 
  return threed_NormSqr(diff);
}

void InterpolateProperty(double *Target, vec diffvec, double xy, double xpy, double xyp, double xpyp){
 
//   printf("dx %lf dy %lf\n %lf %lf %lf %lf\n\n", diffvec.x, diffvec.y, xy, xpy, xyp, xpyp);
  
  *Target = xy 
        - diffvec.x*(xy - xpy)
        - diffvec.x*(xy - xyp)
        + diffvec.x*diffvec.y*(xy - xpy - xyp + xpyp);
}

double ReturnAverageCurvature(sys *Sys, const char *Tag){
  
  MeshPoint **MyMesh = (strcmp(Tag, "Top") ? Sys->TopMesh :  Sys->BotMesh);
  double AveCurvature=0.0;
  
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      AveCurvature+=MyMesh[i][j].lc_x + MyMesh[i][j].lc_y;
    }
  }

  return AveCurvature/Sys->MTot;
}

void PrintVec(vec target, char *note){printf("%s: %lf %lf %lf\n", note, target.x, target.y, target.z);}

const char * AtomNames[] = {"N", "O", "C", "Na", "S"};

void PrintConfiguration(FILE *fp, sys *Sys, const int opt){
  int i, j, mylength, mytype;
  
  if(Sys->NumPTypes>4){
    printf("ERROR: number of atom types exceeds number of labels for atom types in AtomNames struct.\n"
    "Define some new ones and you will be fine\n"
    "Exiting\n");
    exit(1);
  }
  
  fprintf(fp, "%d\nSweep no: %ld out of %ld\n", Sys->vmdDoF, Sys->SweepIndex, Sys->SweepNum);
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
  
  for(mytype = 0;mytype<Sys->NumPTypes;mytype++){
    mylength = opt*(Sys->vmdLength[0][mytype] - 1) + 1;
//     printf("%d %d\n", Sys->TotCurPop[0], Sys->TotCurPop[1]);
    for(j = 0;j<Sys->TotCurPop[0];j++){
      if(Sys->Particles[j].type == mytype ){
        for(int k=0;k<mylength;k++){
          fprintf(
            fp, "%s %lf %lf %lf %lf %lf %lf\n", 
            AtomNames[mytype+1],
            Sys->Particles[j].Pos.x + k*Sys->DPType[0][mytype]*Sys->Particles[j].nhat.x,  
            Sys->Particles[j].Pos.y + k*Sys->DPType[0][mytype]*Sys->Particles[j].nhat.y,  
            Sys->Particles[j].Pos.z + k*Sys->DPType[0][mytype]*Sys->Particles[j].nhat.z,  
            Sys->Particles[j].nhat.x,
            Sys->Particles[j].nhat.y,
            Sys->Particles[j].nhat.z
          );
        }  
      } 
    }
    
    for(j=0;j<Sys->maxpop - Sys->CurParPop[0][mytype];j++){
      for(int k=0;k<mylength;k++){
        fprintf(
          fp, "%s %lf %lf %lf %lf %lf %lf\n", 
          AtomNames[mytype+1], -100.0, -100.0, -100.0, -100.0, -100.0, -100.0
        );
      }  
    }
    
    mylength = opt*(Sys->vmdLength[1][mytype] - 1) + 1;
    for(j = Sys->maxpop;j<Sys->TotCurPop[1] + Sys->maxpop;j++){
      if(Sys->Particles[j].type == mytype ){
        for(int k=0;k<mylength;k++){
          fprintf(
            fp, "%s %lf %lf %lf %lf %lf %lf\n", 
            AtomNames[mytype+1],
            Sys->Particles[j].Pos.x + k*Sys->DPType[1][mytype]*Sys->Particles[j].nhat.x,  
            Sys->Particles[j].Pos.y + k*Sys->DPType[1][mytype]*Sys->Particles[j].nhat.y,  
            Sys->Particles[j].Pos.z + k*Sys->DPType[1][mytype]*Sys->Particles[j].nhat.z,  
            Sys->Particles[j].nhat.x,
            Sys->Particles[j].nhat.y,
            Sys->Particles[j].nhat.z
          );
        }  
      } 
    }
    
    for(j=0;j<Sys->maxpop - Sys->CurParPop[1][mytype];j++){
      for(int k=0;k<mylength;k++){
        fprintf(
          fp, "%s %lf %lf %lf %lf %lf %lf\n", 
          AtomNames[mytype+1], -100.0, -100.0, -100.0, -100.0, -100.0, -100.0
        );
      }  
    }
  }
}

void PrintData(FILE *fp, sys *Sys, long currentMCstep){
  int i, shift = 17;
  
  if (currentMCstep==0) {
    fprintf(fp, "#Sweep(1) TAcc(2) PAcc(3) MAcc(4) GCAcc(5) TotE(6) MeshE(7) BiasE(8) PartE(9) PhaseE(10) PartMeshE(11) TotAreaBot(12) TotAreaTop(13) ThetaValTop(14) ThetaValBot(15) TotPop(16) TotBotPop(17) ");
    for(i = 0; i<Sys->NumPTypes; i++){fprintf(fp, "NBot%d(%d) ", i+1, shift+1+i);}
    fprintf(fp, "TotTopPop(%d) ", shift + 1 + Sys->NumPTypes);
    for(i = 0; i<Sys->NumPTypes; i++){fprintf(fp, "NTop%d(%d) ", i+1,  shift + 2 +   Sys->NumPTypes + i);}
    for(i = 0; i<Sys->NumPTypes; i++){fprintf(fp, "Bonds%d(%d) ", i+1, shift + 2 + 2*Sys->NumPTypes + i);}
    fprintf(fp, "\n");
  }
  
  fprintf( fp, "%ld %.2lf %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
    Sys->SweepIndex, 
    Sys->AcceptanceRate/(currentMCstep+1), Sys->AcceptPart/Sys->ParticleAttempt, Sys->AcceptMesh/Sys->MeshAttempt,
    Sys->AcceptGC/Sys->GCAttempt,
    Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotParticleEnergy, Sys->TotPhaseFactorEnergy, Sys->TotMeshParticleEnergy, 
    ReturnTotalArea(Sys, 0), ReturnTotalArea(Sys, 1),
    Sys->MeshBias.thetaTop, Sys->MeshBias.thetaBot
  );
  
  fprintf( fp, "%d " , Sys->TotCurPop[0] + Sys->TotCurPop[1]);
  fprintf( fp, "%d " , Sys->TotCurPop[0]);
  for(i = 0; i<Sys->NumPTypes; i++){fprintf(fp, "%d ", Sys->CurParPop[0][i]);}
  fprintf( fp, "%d " , Sys->TotCurPop[1]);
  for(i = 0; i<Sys->NumPTypes; i++){fprintf(fp, "%d ", Sys->CurParPop[1][i]);}
  
  int k=0;
  for(i=0;i<Sys->NumPTypes;i++){ 
    if (Sys->BindPType[i]>1e-10){
      fprintf( fp, "%.0lf " , Sys->BondList.Curbondnums[k]);
      k++;
    }
    else{fprintf( fp, "%d " , 0);}
  }
  fprintf(fp,"\n");
}

void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep){
  if (currentMCstep==0) fprintf(fp, "Sweep(1) TAcc(2) PAcc(3) MAcc(4) GCAcc(5) TotE(6) MeshE(7) BiasE(8) PartE(9) PhaseE(10) PartMeshE(11) TotBotPop(12) TotTopPop(13)\n");
  
  fprintf( fp,
    "%ld %.2lf %.2lf %.2lf %.2lf %lf %lf %lf %lf %lf %lf %d %d\n", 
    Sys->SweepIndex, 
    Sys->AcceptanceRate/(currentMCstep+1), 
    Sys->AcceptPart/Sys->ParticleAttempt, Sys->AcceptMesh/Sys->MeshAttempt, Sys->AcceptGC/Sys->GCAttempt,
    Sys->TotEnergy, Sys->MeshEnergy, Sys->BiasEnergy, Sys->TotParticleEnergy, Sys->TotPhaseFactorEnergy, Sys->TotMeshParticleEnergy, Sys->TotCurPop[0], Sys->TotCurPop[1]
  );
}

void PrintJustMesh(FILE *fp, FILE *gp, sys *Sys){
  int i,j;
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "%lf %lf %lf %lf %lf %lf\n", 
      Sys->BotMesh[i][j].Pos.x,  Sys->BotMesh[i][j].Pos.y,  Sys->BotMesh[i][j].Pos.z,
      Sys->BotMesh[i][j].nhat.x, Sys->BotMesh[i][j].nhat.y, Sys->BotMesh[i][j].nhat.z
    );
   }
   fprintf(fp, "\n");
  }
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      gp, "%lf %lf %lf %lf %lf %lf\n", 
      Sys->TopMesh[i][j].Pos.x,  Sys->TopMesh[i][j].Pos.y,  Sys->TopMesh[i][j].Pos.z,
      Sys->TopMesh[i][j].nhat.x, Sys->TopMesh[i][j].nhat.y, Sys->TopMesh[i][j].nhat.z
    );
   }
   fprintf(gp, "\n");
  }
}

void WriteGPScript(char *outname, char *label, particle P1, particle P2){
  FILE *gp = fopen(outname, "w");
  
  vec a, b, c;
  a.x=1.0;a.y=0.0;a.z=0.0;
  crossprod(&a, P1.nhat, a);
  double tmp = sqrt(threed_NormSqr(a) );
  if(tmp<1e-5){
    a.x=0.0;a.y=1.0;a.z=0.0;
    b.x=0.0;b.y=0.0;b.z=1.0;
  }
  else{
    crossprod(&b, P1.nhat, a);
    threed_NormalizeVec(&b);
    threed_NormalizeVec(&a);
  }
    
  CopyVec(&c, P1.MidPoint);
  double MyRad=P1.MyRad;
  double Mylengthon2=P1.MyLength/2;
  
  fprintf(gp,"clear;set terminal qt font 'Helvetica, 24';set parametric;\nset isosample 100;\n");
  fprintf(gp, "set vrange [-10:10];\n");
  fprintf(gp, "set xrange [-10:10];");
  fprintf(gp, "set yrange [-10:10];");
  fprintf(gp, "set zrange [-10:10];\n");
  
  
  fprintf(gp, "sp ");
  fprintf(gp, "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN), "
              "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN), "
              "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN) title 'P1', ",
                c.x, a.x*MyRad, b.x*MyRad, P1.nhat.x, Mylengthon2,
                c.y, a.y*MyRad, b.y*MyRad, P1.nhat.y, Mylengthon2,
                c.z, a.z*MyRad, b.z*MyRad, P1.nhat.z, Mylengthon2
  );
  
  a.x=1.0;a.y=0.0;a.z=0.0;
  crossprod(&a, P2.nhat, a);
  tmp = sqrt(threed_NormSqr(a) );
  if(tmp<1e-5){
    a.x=0.0;a.y=1.0;a.z=0.0;
    b.x=0.0;b.y=0.0;b.z=1.0;
  }
  else{
    crossprod(&b, P1.nhat, a);
    threed_NormalizeVec(&b);
    threed_NormalizeVec(&a);
  }
  crossprod(&b, P2.nhat, a);
  CopyVec(&c, P2.MidPoint);
  MyRad=P2.MyRad;
  Mylengthon2=P2.MyLength/2;
  
  fprintf(gp, "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN), "
              "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN), "
              "%lf + %lf*cos(u) + %lf*sin(u) + %lf*(abs(v) < %lf ? v:NaN) title 'P2' ",
                c.x, a.x*MyRad, b.x*MyRad, P2.nhat.x, Mylengthon2,
                c.y, a.y*MyRad, b.y*MyRad, P2.nhat.y, Mylengthon2,
                c.z, a.z*MyRad, b.z*MyRad, P2.nhat.z, Mylengthon2
  );
  
  fprintf(gp, ";set title '%s'\npause -1\n", label);
  fclose(gp);
}

void TestCollisionDetection(sys *Sys){
  int loop;
  
  particle P1, P2;
  
  P1.MidPoint.x=1.0;P1.MidPoint.y=-4.0;P1.MidPoint.z=0.0;
  P1.nhat.x=0.0;P1.nhat.y=0.0;P1.nhat.z=3.0;
  threed_NormalizeVec(&P1.nhat);
  P1.MyRad = 2.0;
  P1.MyLength= 10.0;
  
  P2.MidPoint.x=0.0;P2.MidPoint.y=-10.0;P2.MidPoint.z=2.0;
  P2.nhat.x=2.0;P2.nhat.y=2.0;P2.nhat.z=1.0;
  threed_NormalizeVec(&P2.nhat);
  P2.MyRad = 2.0;
  P2.MyLength= 10.0;

  double overlap=CheckOverlap(P2, P1, 10000, 10000, &loop);
    
  char outname[1024], result[1024];
  strcpy(outname, "mygpscript.gp");
  if (overlap == ALOT)sprintf(result, "overlap %d", loop);
  else sprintf(result, "notoverlap %d", loop);
  WriteGPScript(outname, result, P1, P2);
}









