#include "particles.h"

int ErrorCheck_CellLists(sys *Sys){
  celllist *CellList = &(Sys->CellList);
  
  int nx,ny, p, tmpcx, tmpcy, index, Err=0;
  int botpoptest=0, toppoptest=0;
  
  for(nx=0;nx<CellList->NumCellsx;nx++){
    for(ny=0;ny<CellList->NumCellsy;ny++){
      for(p=0;p<CellList->BotPopulations[nx][ny];p++){
        index=CellList->BotParticleIds[nx][ny][p];
        
        Position_To_Index_PBC(&tmpcx, Sys->Particles[index].Pos.x, CellList->dCx, CellList->NumCellsx);
        Position_To_Index_PBC(&tmpcy, Sys->Particles[index].Pos.y, CellList->dCy, CellList->NumCellsy);
        
        if(nx!=tmpcx || ny!=tmpcy){
          printf("\nCellList ERROR!\nParticle %d missassigned! Bot particle Pos:%f %f\nAssigned to %d %d entroNo %d popsize: %d\nActually in %d %d\nTotalVals %d %d Lx %lf Ly %lf dcx %lf dcy %lf\n",
          index, Sys->Particles[index].Pos.x, Sys->Particles[index].Pos.y, nx, ny, p, CellList->BotPopulations[nx][ny],tmpcx, tmpcy, CellList->NumCellsx, CellList->NumCellsy, Sys->Lx, Sys->Ly, CellList->dCx, CellList->dCy
          );
          Err=4;
        }
      }
      for(p=0;p<CellList->TopPopulations[nx][ny];p++){
        index=CellList->TopParticleIds[nx][ny][p];
        
        Position_To_Index_PBC(&tmpcx, Sys->Particles[index].Pos.x, CellList->dCx, CellList->NumCellsx);
        Position_To_Index_PBC(&tmpcy, Sys->Particles[index].Pos.y, CellList->dCy, CellList->NumCellsy);
        
        if(nx!=tmpcx || ny!=tmpcy){
          printf("\nCellList ERROR!\nParticle %d Top missassigned! particle Pos:%f %f\nAssigned to %d %d popsize: %d\nActually in %d %d\nTotalVals %d %d Lx %lf Ly %lf dcx %lf dcy %lf\n",
            index, Sys->Particles[index].Pos.x, Sys->Particles[index].Pos.y, nx, ny, CellList->TopPopulations[nx][ny], tmpcx, tmpcy, CellList->NumCellsx, CellList->NumCellsy, Sys->Lx, Sys->Ly, CellList->dCx, CellList->dCy
          );
          Err=1;
        }
      }
    }
  }
  
  for(nx=0;nx<CellList->NumCellsx;nx++){
    for(ny=0;ny<CellList->NumCellsy;ny++){
      botpoptest+=CellList->BotPopulations[nx][ny];
      toppoptest+=CellList->TopPopulations[nx][ny];
    }
  }
  if(botpoptest!=Sys->TotCurPop[0]){
    printf("\nCellList ERROR!\nMismatch in bot population: clist pop: %d Actual pop:%d\n", botpoptest, Sys->TotCurPop[0]);Err=1;
  }
  
  if(toppoptest!=Sys->TotCurPop[1]){
    printf("\nCellList ERROR!\nMismatch in top population: clist pop: %d Actual pop:%d\n", toppoptest, Sys->TotCurPop[1]);Err=1;
  }
  
  return Err;
}

int ErrorCheck_Z_nhat_PhaseFactor(sys *Sys){
  int i, p, inew, jnew, xNeigh, yNeigh, Err=0;
  MeshPoint **MyMesh;
  vec testnhat, diffvec, v1, v2 , v3, v4;
  double Zval, TotPhaseEnergy=0.0, PhaseEnergy;
  
  for(i=0;i<2;i++){
    for(p=i*Sys->maxpop;p<Sys->TotCurPop[i]+i*(Sys->maxpop);p++){
      
      if(i==0) MyMesh = Sys->BotMesh;
      if(i==1) MyMesh = Sys->TopMesh;
      
      Position_To_Index_PBC(&inew, Sys->Particles[p].Pos.x, Sys->dx, Sys->Nx);
      Position_To_Index_PBC(&jnew, Sys->Particles[p].Pos.y, Sys->dy, Sys->Ny);
      
      twodvecdiff(&(diffvec), Sys->Particles[p].Pos, MyMesh[inew][jnew].Pos);
      NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
      diffvec.x/=Sys->dx;
      diffvec.y/=Sys->dy;
      
      xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
      yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
      
      InterpolateProperty(&Zval, diffvec, 
                                 MyMesh[inew][jnew].Pos.z, 
                                 MyMesh[xNeigh][jnew].Pos.z, 
                                 MyMesh[inew][yNeigh].Pos.z, 
                                 MyMesh[xNeigh][yNeigh].Pos.z 
                         );
    
      if(fabs(Zval - Sys->Particles[p].Pos.z)>1e-6){
        printf("\nZval Error!\n");Err=2;
        printf("Id: %d\n", p);
        printf("computed height %lf\n", Zval);
        printf("particle height %lf\n", Sys->Particles[p].Pos.z);
      }
      
      
      CopyVec(&v1, MyMesh[inew][jnew].nhat);     rescalevec(&v1, MyMesh[inew][jnew].nNorm    );
      CopyVec(&v2, MyMesh[xNeigh][jnew].nhat);   rescalevec(&v2, MyMesh[xNeigh][jnew].nNorm  );
      CopyVec(&v3, MyMesh[inew][yNeigh].nhat);   rescalevec(&v3, MyMesh[inew][yNeigh].nNorm  );
      CopyVec(&v4, MyMesh[xNeigh][yNeigh].nhat); rescalevec(&v4, MyMesh[xNeigh][yNeigh].nNorm);
      
      InterpolateProperty(&(testnhat.x), diffvec, v1.x, v2.x, v3.x, v4.x );
      InterpolateProperty(&(testnhat.y), diffvec, v1.y, v2.y, v3.y, v4.y );
      InterpolateProperty(&(testnhat.z), diffvec, v1.z, v2.z, v3.z, v4.z );
      
      PhaseEnergy = -0.5*log(1.0 
                            + testnhat.x*testnhat.x
                            + testnhat.y*testnhat.y);
      TotPhaseEnergy+=PhaseEnergy;
      
      if(fabs(PhaseEnergy - Sys->Particles[p].PhaseFactorEnergy)>1e-6){
        printf("\nPhaseError!\n");Err=2;
        printf("EComp: %lf EStored: %lf EDiff: %lf\n", PhaseEnergy, Sys->Particles[p].PhaseFactorEnergy, PhaseEnergy - Sys->Particles[p].PhaseFactorEnergy);
        printf("Id: %d\n", p);
        PrintVec(Sys->Particles[p].Pos,"Position");
      }
      
      threed_NormalizeVec( &testnhat );
      if(fabs(testnhat.x -Sys->Particles[p].nhat.x)>1e-6 || fabs(testnhat.y -Sys->Particles[p].nhat.y)>1e-6 || fabs(testnhat.z -Sys->Particles[p].nhat.z)>1e-6){
        printf("\nNhat Error!\n");Err=2;
        printf("Id: %d\n", p);
        PrintVec(testnhat,"computed orientation");
        PrintVec(Sys->Particles[p].nhat,"particle orientation");
      }
      threedvecsumscalar( &(testnhat), Sys->Particles[p].Pos, Sys->Particles[p].nhat, Sys->Particles[p].MyLength/2);
      if(fabs(testnhat.x -Sys->Particles[p].MidPoint.x)>1e-6 || fabs(testnhat.y -Sys->Particles[p].MidPoint.y)>1e-6 || fabs(testnhat.z -Sys->Particles[p].MidPoint.z)>1e-6){
        printf("\nMidpoint Error!\n");Err=2;
        printf("Id: %d\n", p);
        PrintVec(testnhat,"computed midpoint");
        PrintVec(Sys->Particles[p].MidPoint,"particle midpoint");
      }
      threedvecsumscalar( &testnhat, Sys->Particles[p].Pos, Sys->Particles[p].nhat, Sys->Particles[p].MyLength);
      if(fabs(testnhat.x -Sys->Particles[p].EndPoint.x)>1e-6 || fabs(testnhat.y -Sys->Particles[p].EndPoint.y)>1e-6 || fabs(testnhat.z -Sys->Particles[p].EndPoint.z)>1e-6){
        printf("\nEndpoint Error!\n");Err=2;
        printf("Id: %d %lf %lf\n", p, Zval, Sys->Particles[p].Pos.z);
        PrintVec(testnhat,"computed Endpoint");
        PrintVec(Sys->Particles[p].EndPoint, "particle Endpoint");
        break;
      }
    }
  }
  
  if(fabs(TotPhaseEnergy - Sys->TotPhaseFactorEnergy)>1e-6){printf("\nError Mismatch in PhaseE %lf TotalE %lf\n", TotPhaseEnergy, Sys->TotPhaseFactorEnergy);printf("p1: %d p2: %d\n", Sys->TotCurPop[0], Sys->TotCurPop[1]);Err=2;}
  
  return Err;
}

int ErrorCheck_ParticleEnergies(sys *Sys){
  int Err = 0;
  
  if ( fabs(Sys->TotParticleEnergy - ReturnTotalParticleEnergies(Sys))>1e-6 ){ printf("\nStored Particle Energy: %lf Computed: %lf\n", Sys->TotParticleEnergy, ReturnTotalParticleEnergies(Sys));Err=3;}
  return Err;
}

double ReturnTotalParticleEnergies(sys *Sys){
  int tmp;
  particle CurrentParticle, PartnerParticle;
  double Lx=Sys->Lx, Ly=Sys->Ly, l_current_on2, l_partner_on2, cutoffsqr;
  
  int MyMeshID, MyTypeID;
  double DeltaEPart=0.0, MyintRangesqr, MyIntStr;
  
  for(int p0=0;p0<Sys->twomaxpop;p0++){
    if(p0 == Sys->TotCurPop[0]) p0=Sys->maxpop;
    if(p0 == Sys->TotCurPop[1] + Sys->maxpop) break;
    
    CurrentParticle = Sys->Particles[p0];
    l_current_on2 = CurrentParticle.MyLength/2;
    MyMeshID = CurrentParticle.MeshId;
    MyTypeID=CurrentParticle.type;
    MyintRangesqr=CurrentParticle.MyIntCutoffsqr;
    MyIntStr=CurrentParticle.MyIntStr;
    
    for(int index = p0+1; index<Sys->twomaxpop;index++){
      if(index == Sys->TotCurPop[0]) index=Sys->maxpop;
      if(index == Sys->TotCurPop[1] + Sys->maxpop) break;
      
      PartnerParticle = Sys->Particles[index];
      l_partner_on2 = PartnerParticle.MyLength/2;
      cutoffsqr= (l_partner_on2 + l_current_on2)*(l_partner_on2 + l_current_on2);
      
      if( deltaRsqr(CurrentParticle.MidPoint, PartnerParticle.MidPoint, Lx, Ly) <  cutoffsqr){
        if ( CheckOverlap(PartnerParticle, CurrentParticle, Lx, Ly, &tmp) == ALOT) {
          PrintVec(CurrentParticle.MidPoint, "r1");
          PrintVec(PartnerParticle.MidPoint, "r2");
          printf("%d %d\n", p0, index);
          printf("%lf %lf\n", l_partner_on2, l_current_on2);
          printf("%d %lf\n", tmp, sqrt(deltaRsqr(CurrentParticle.MidPoint, PartnerParticle.MidPoint, Lx, Ly)));
          return ALOT;
        }
      }
      
      if( MyIntStr>1e-10                     &&
          PartnerParticle.MeshId != MyMeshID && 
          PartnerParticle.type == MyTypeID   &&
          deltaRsqr(CurrentParticle.EndPoint, PartnerParticle.EndPoint, Lx, Ly)< MyintRangesqr
          ){
//           printf("%d %d\n");
          DeltaEPart-=MyIntStr;
      }
    }
  }
  return DeltaEPart;
}

double SimpleParticleInteractions(sys *Sys){
  vec diffnew, diffold;
  if      (Sys->p0 == 0 ){ 
    threedvecdiff(&diffnew, Sys->dummy_particle.Pos, Sys->Particles[1].Pos);
    threedvecdiff(&diffold, Sys->Particles[0].Pos, Sys->Particles[1].Pos);
  }
  else if (Sys->p0 == 1 ) {
    threedvecdiff(&diffnew, Sys->dummy_particle.Pos, Sys->Particles[0].Pos);  
    threedvecdiff(&diffold, Sys->Particles[0].Pos, Sys->Particles[1].Pos);
  }
  else if (Sys->p0 == 2 ) {
    threedvecdiff(&diffnew, Sys->dummy_particle.Pos, Sys->Particles[3].Pos);  
    threedvecdiff(&diffold, Sys->Particles[3].Pos, Sys->Particles[2].Pos);
  }
  else if (Sys->p0 == 3 ) {
    threedvecdiff(&diffnew, Sys->dummy_particle.Pos, Sys->Particles[2].Pos);  
    threedvecdiff(&diffold, Sys->Particles[3].Pos, Sys->Particles[2].Pos); 
  }
    
  NearestImageConv(&diffnew, Sys->Lx, Sys->Ly);
  NearestImageConv(&diffold, Sys->Lx, Sys->Ly);
  
  return ( threed_NormSqr(diffnew) - threed_NormSqr(diffold) );
}

void Compute_Z_nhat_PhaseFactor(sys *Sys, particle *My_Particle){
  double length = My_Particle->MyLength;
  MeshPoint **MyMesh = (My_Particle->MeshId == TOP) ? Sys->TopMesh : Sys->BotMesh;
  
  vec diffvec, v1, v2 , v3, v4;
  int inew, jnew, xNeigh, yNeigh;
  
  Position_To_Index_PBC(&inew, My_Particle->Pos.x, Sys->dx, Sys->Nx);
  Position_To_Index_PBC(&jnew, My_Particle->Pos.y, Sys->dy, Sys->Ny);
  
  twodvecdiff(&(diffvec), My_Particle->Pos, MyMesh[inew][jnew].Pos);
  NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
  diffvec.x/=Sys->dx;
  diffvec.y/=Sys->dy;
  
  xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
  yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
  
  InterpolateProperty(&(My_Particle->Pos.z), diffvec, 
                                             MyMesh[inew][jnew].Pos.z, 
                                             MyMesh[xNeigh][jnew].Pos.z, 
                                             MyMesh[inew][yNeigh].Pos.z, 
                                             MyMesh[xNeigh][yNeigh].Pos.z 
                     );
  
  CopyVec(&v1, MyMesh[inew][jnew].nhat);     rescalevec(&v1, MyMesh[inew][jnew].nNorm    );
  CopyVec(&v2, MyMesh[xNeigh][jnew].nhat);   rescalevec(&v2, MyMesh[xNeigh][jnew].nNorm  );
  CopyVec(&v3, MyMesh[inew][yNeigh].nhat);   rescalevec(&v3, MyMesh[inew][yNeigh].nNorm  );
  CopyVec(&v4, MyMesh[xNeigh][yNeigh].nhat); rescalevec(&v4, MyMesh[xNeigh][yNeigh].nNorm);
  
  InterpolateProperty(&(My_Particle->nhat.x), diffvec, v1.x, v2.x, v3.x, v4.x );
  InterpolateProperty(&(My_Particle->nhat.y), diffvec, v1.y, v2.y, v3.y, v4.y );
  InterpolateProperty(&(My_Particle->nhat.z), diffvec, v1.z, v2.z, v3.z, v4.z );
  My_Particle->PhaseFactorEnergy = -0.5*log(1.0 
                                            + My_Particle->nhat.x*My_Particle->nhat.x
                                            + My_Particle->nhat.y*My_Particle->nhat.y);//Use unNormalized orientation vector!
  
  threed_NormalizeVec( &My_Particle->nhat );
  threedvecsumscalar(  &(My_Particle->MidPoint), My_Particle->Pos, My_Particle->nhat, length/2);
  threedvecsumscalar(  &(My_Particle->EndPoint), My_Particle->Pos, My_Particle->nhat, length);
}

void ComputeAllParticleBonds(sys *Sys){
  particle pi, pj;
  double cutoffsqr, mystrength, Lx = Sys->Lx, Ly = Sys->Ly;
  int myparticletype, myMeshId;

  for(int i = 0;i<Sys->twomaxpop;i++){
    if(i == Sys->TotCurPop[0]) i=Sys->maxpop;
    if(i == Sys->TotCurPop[1] + Sys->maxpop) break;
    
    pi = Sys->Particles[i];
    mystrength= pi.MyIntStr;
    if(mystrength>1e-10){
      myparticletype = pi.type;
      cutoffsqr = pi.MyIntCutoffsqr;
      myMeshId=pi.MeshId;
      
      for(int j = i+1; j<Sys->twomaxpop;j++){
        if(j == Sys->TotCurPop[0]) j=Sys->maxpop;
        if(j == Sys->TotCurPop[1] + Sys->maxpop) break;
        
        pj = Sys->Particles[j];
        if( pj.MeshId != myMeshId &&  
            pj.type == myparticletype && 
            deltaRsqr(pi.EndPoint, pj.EndPoint, Lx, Ly)< cutoffsqr 
          )
        {Sys->TotParticleEnergy -= mystrength;}
      }
    }
  }
}

double MeshParticleInteraction(double dRsqr){
  return ALOT;
}

double ParticleMeshEnergy(sys *Sys, particle *MyParticle){
  int MyMeshId = MyParticle->MeshId;
  double MPCutoffsqr=MyParticle->MeshParticleCutoffsqr;
  vec MyPos; 
  CopyVec( &(MyPos), MyParticle->EndPoint);
  MeshPoint **OtherMesh = (MyMeshId == TOP) ? Sys->BotMesh : Sys->TopMesh;
  
  vec diffvec, v;
  int inew, jnew, xNeigh, yNeigh;
  
  Position_To_Index_PBC(&inew, MyPos.x, Sys->dx, Sys->Nx);
  Position_To_Index_PBC(&jnew, MyPos.y, Sys->dy, Sys->Ny);
  
  twodvecdiff(&(diffvec), MyPos, OtherMesh[inew][jnew].Pos);
  NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
  v.x=diffvec.x/Sys->dx;
  v.y=diffvec.y/Sys->dy;
  
  xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
  yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
  
  InterpolateProperty( &(diffvec.z), 
                       v, 
                       OtherMesh[inew][jnew].Pos.z, 
                       OtherMesh[xNeigh][jnew].Pos.z, 
                       OtherMesh[inew][yNeigh].Pos.z, 
                       OtherMesh[xNeigh][yNeigh].Pos.z 
                     );
  
  diffvec.z -= MyPos.z;
  double dRsqr = threed_NormSqr(diffvec), MeshParticleEnergy = 0.0;
  
  if( MyMeshId*diffvec.z < 0.0)  MeshParticleEnergy+= ALOT;
  else if( dRsqr < MPCutoffsqr ) MeshParticleEnergy+=MeshParticleInteraction(dRsqr);
  
  return MeshParticleEnergy;
}

void ComputeMeshParticleEnergy(sys *Sys){
  Sys->TotMeshParticleEnergy = 0.0;
  double MeshParticleEnergy;
  for(int p = 0; p<Sys->twomaxpop;p++){
    if(p == Sys->TotCurPop[0]) p=Sys->maxpop;
    if(p == Sys->TotCurPop[1] + Sys->maxpop) break;
    
    MeshParticleEnergy = ParticleMeshEnergy(Sys, &(Sys->Particles[p]));
    Sys->Particles[p].MeshParticleEnergy = MeshParticleEnergy;
    Sys->TotMeshParticleEnergy  += MeshParticleEnergy;
  }
}

double ComputeDeltaMPEnergy(sys *Sys, particle *MyParticle, int p0){
  double MeshParticleEnergy = ParticleMeshEnergy(Sys, MyParticle);
  MyParticle->MeshParticleEnergy = MeshParticleEnergy;
  return  MeshParticleEnergy - Sys->Particles[p0].MeshParticleEnergy;
}

void InitializeParticleEnergies(sys *Sys){
  Sys->TotParticleEnergy=0.0;Sys->TotPhaseFactorEnergy=0.0;
  for(int p=0;p<Sys->twomaxpop;p++){
    if(p == Sys->TotCurPop[0]) p=Sys->maxpop;
    if(p == Sys->TotCurPop[1] + Sys->maxpop) break;
    Sys->Particles[p].MeshParticleEnergy = 0.0;
    Compute_Z_nhat_PhaseFactor(Sys, &(Sys->Particles[p]));
    Sys->TotPhaseFactorEnergy += Sys->Particles[p].PhaseFactorEnergy;
  }
  ComputeMeshParticleEnergy(Sys);
  ComputeAllParticleBonds(Sys);
}

void ParticleToDummyParticle(sys *Sys, int p0){
  Sys->dummy_particle.type        = Sys->Particles[p0].type;
  Sys->dummy_particle.MeshId      = Sys->Particles[p0].MeshId;
  Sys->dummy_particle.MyDiam      = Sys->Particles[p0].MyDiam;
  Sys->dummy_particle.MyRad       = Sys->Particles[p0].MyRad;
  Sys->dummy_particle.MyRadSqr    = Sys->Particles[p0].MyRadSqr;
  Sys->dummy_particle.MyIntStr    = Sys->Particles[p0].MyIntStr;
  Sys->dummy_particle.MyIntCutoff = Sys->Particles[p0].MyIntCutoff;
  Sys->dummy_particle.MyIntCutoffsqr = Sys->Particles[p0].MyIntCutoffsqr;
  Sys->dummy_particle.MeshParticleCutoffsqr = Sys->Particles[p0].MeshParticleCutoffsqr;
  Sys->dummy_particle.MyLength = Sys->Particles[p0].MyLength;
}

double SuggestParticleMove(sys *Sys){
  
  if(Sys->TotCurPop[0] + Sys->TotCurPop[1]>0){
    int p0 = gsl_rng_uniform_int(Sys->rng, Sys->TotCurPop[0] + Sys->TotCurPop[1] );
    p0 = p0 + (p0 < Sys->TotCurPop[0] ? 0 : Sys->maxpop - Sys->TotCurPop[0]);
    Sys->dummy_particle.Pos.x  = Sys->Particles[p0].Pos.x + (gsl_rng_uniform(Sys->rng) - 0.5)*Sys->PDisp;
    Sys->dummy_particle.Pos.y  = Sys->Particles[p0].Pos.y + (gsl_rng_uniform(Sys->rng) - 0.5)*Sys->PDisp;
    PosPBC(&(Sys->dummy_particle.Pos), Sys->Lx, Sys->Ly);
    ParticleToDummyParticle(Sys, p0);
    Sys->p0 = p0;
    Compute_Z_nhat_PhaseFactor(Sys, &(Sys->dummy_particle) );
    return Sys->dummy_particle.PhaseFactorEnergy - Sys->Particles[p0].PhaseFactorEnergy;
  }
  else{return ALOT;}
}

double PerformLineTest(double fpx, double fpy, double l_moved_on2, double r_moved, vec MySlope){
  double s, t;
  
  t = (-r_moved - fpx)/MySlope.x;
  if(fabs(t)<1.0){
    s = (MySlope.y * t + fpy)/l_moved_on2;
    if(fabs(s)<1.0){return ALOT;}
  }
  t = (r_moved - fpx)/MySlope.x;
  if(fabs(t)<1.0){
    s = (MySlope.y * t + fpy)/l_moved_on2;
    if(fabs(s)<1.0){return ALOT;}
  }
  t = (-l_moved_on2 - fpy)/MySlope.y;
  if(fabs(t)<1.0){
    s = (MySlope.x * t + fpx)/r_moved;
    if(fabs(s)<1.0){return ALOT;}
  }
  t = ( l_moved_on2 - fpy)/MySlope.y;
  if(fabs(t)<1.0){
    s = (MySlope.x * t + fpx)/r_moved;
    if(fabs(s)<1.0){return ALOT;}
  }
  return 0.0;
}

double CheckOverlap(particle Moved, particle Partner, double Lx, double Ly, int *loop){
  double tmp;
  double l_moved_on2 = Moved.MyLength/2, l_partner_on2 = Partner.MyLength/2;
  double r_moved = Moved.MyRad, r_partner = Partner.MyRad;
  
  vec axperp_both_nhat, axperp_moved_nhat, r12, nhat_partner;

  CopyVec(&nhat_partner, Partner.nhat);
  threedvecdiff(&r12, Moved.MidPoint, Partner.MidPoint);
  NearestImageConv(&r12, Lx, Ly );
  
  crossprod(&axperp_both_nhat, Moved.nhat, nhat_partner);
  tmp = sqrt(threed_NormSqr(axperp_both_nhat) );
 
  //If cylinders are parrallel, use simple test
  if (tmp<1e-5){
    tmp = dotprod(r12, nhat_partner);
    axperp_both_nhat.x = r12.x - tmp*nhat_partner.x;
    axperp_both_nhat.y = r12.y - tmp*nhat_partner.y;
    axperp_both_nhat.z = r12.z - tmp*nhat_partner.z;
    if( fabs(tmp) > l_partner_on2 + l_moved_on2 ||
        sqrt(threed_NormSqr(axperp_both_nhat)) > r_partner + r_moved){*loop=-1;return 0.0;}
    else{*loop=1;return ALOT;}
  }
  
  else{
    
    axperp_both_nhat.x/=tmp;
    axperp_both_nhat.y/=tmp;
    axperp_both_nhat.z/=tmp;
    //Check if cylinders are too high to intersect
    if ( fabs(dotprod(r12, axperp_both_nhat))> r_partner + r_moved ){*loop=-2;return 0.0;}
    
    //Transform to coordinate system with moved particle parallel to y axis and centered at origin
    crossprod(&axperp_moved_nhat, Moved.nhat, axperp_both_nhat);
    threed_NormalizeVec(&axperp_moved_nhat);
  
    vec trans_r12, trans_partner_n, trans_partner_n_perp;
    
    trans_r12.x = dotprod(r12, axperp_moved_nhat);
    trans_r12.y = dotprod(r12, Moved.nhat);
    
    trans_partner_n.x = l_partner_on2*dotprod(nhat_partner, axperp_moved_nhat);
    trans_partner_n.y = l_partner_on2*dotprod(nhat_partner, Moved.nhat);
    
    trans_partner_n_perp.x = -trans_partner_n.y*(r_partner/l_partner_on2);
    trans_partner_n_perp.y =  trans_partner_n.x*(r_partner/l_partner_on2);
    
    /***********If necessary, add corner tests to analysis************/
    
    //Line test for rotated rectangles
    
    double fpx, fpy;
    
    fpx = trans_r12.x + trans_partner_n_perp.x;
    fpy = trans_r12.y + trans_partner_n_perp.y;
    if(PerformLineTest(fpx, fpy, l_moved_on2, r_moved, trans_partner_n) == ALOT){*loop=2;return ALOT;}

    fpx = trans_r12.x - trans_partner_n_perp.x;
    fpy = trans_r12.y - trans_partner_n_perp.y;
    if(PerformLineTest(fpx, fpy, l_moved_on2, r_moved, trans_partner_n) == ALOT){*loop=3;return ALOT;}

    fpx = trans_r12.x + trans_partner_n.x;
    fpy = trans_r12.y + trans_partner_n.y;
    if(PerformLineTest(fpx, fpy, l_moved_on2, r_moved, trans_partner_n_perp) == ALOT){*loop=4;return ALOT;}

    fpx = trans_r12.x - trans_partner_n.x;
    fpy = trans_r12.y - trans_partner_n.y;
    if(PerformLineTest(fpx, fpy, l_moved_on2, r_moved, trans_partner_n_perp) == ALOT){*loop=5;return ALOT;}

    *loop=-3;
    return 0.0;
  }
}

double CheckForBond(particle TestParticle, particle  PartnerParticle, double Lx, double Ly, sys *Sys){
  int MyMeshID = TestParticle.MeshId, MyTypeID=TestParticle.type;
  double MyintRangesqr=TestParticle.MyIntCutoffsqr, MyIntStr=TestParticle.MyIntStr, E=0.0;
    
  if( PartnerParticle.MeshId != MyMeshID && PartnerParticle.type == MyTypeID){
    if(deltaRsqr(TestParticle.EndPoint, PartnerParticle.EndPoint, Lx, Ly)< MyintRangesqr){
      E=-MyIntStr;
    }
  }
  return E;
}

double CollisionTestAndParticleBonds(sys *Sys){
  int p0 = Sys->p0;
  particle MovedParticle = Sys->dummy_particle;
  celllist *CellList = &(Sys->CellList);
  double Lx=Sys->Lx, Ly=Sys->Ly, l_moved_on2 = MovedParticle.MyLength/2;
  double MyIntStr=MovedParticle.MyIntStr;
  
  particle PartnerParticle;
  int tmp, MyCx, MyCy, Cx, Cy, partner, k;
  double l_partner_on2, cutoffsqr, deltaE = 0.0;
  
  Position_To_Index_PBC(&MyCx, MovedParticle.Pos.x, CellList->dCx, CellList->NumCellsx);
  Position_To_Index_PBC(&MyCy, MovedParticle.Pos.y, CellList->dCy, CellList->NumCellsy);
  
  for(int i=0;i<9;i++){
    Cx = CellList->Neighbors[MyCx][MyCy][2*i];
    Cy = CellList->Neighbors[MyCx][MyCy][2*i + 1];
    
    for(k=0;k<CellList->BotPopulations[Cx][Cy];k++){
      partner = CellList->BotParticleIds[Cx][Cy][k];
      if(partner!=p0){
        PartnerParticle = Sys->Particles[partner];
        l_partner_on2 = PartnerParticle.MyLength/2;
        cutoffsqr= (l_partner_on2 + l_moved_on2)*(l_partner_on2 + l_moved_on2);
        if( deltaRsqr(MovedParticle.MidPoint, PartnerParticle.MidPoint, Lx, Ly) <  cutoffsqr){
          if ( CheckOverlap(PartnerParticle, MovedParticle, Lx, Ly, &tmp) == ALOT) return ALOT;
        }
        if( MyIntStr>1e-10 ){
          deltaE += CheckForBond(MovedParticle     , PartnerParticle, Lx, Ly, Sys) 
                   -CheckForBond(Sys->Particles[p0], PartnerParticle, Lx, Ly, Sys);
        }
      }
    }
    
    for(k=0;k<CellList->TopPopulations[Cx][Cy];k++){
      partner = CellList->TopParticleIds[Cx][Cy][k];
      if(partner!=p0){
        PartnerParticle = Sys->Particles[partner];
        l_partner_on2 = PartnerParticle.MyLength/2;
        cutoffsqr= (l_partner_on2 + l_moved_on2)*(l_partner_on2 + l_moved_on2);
        if( deltaRsqr(MovedParticle.MidPoint, PartnerParticle.MidPoint, Lx, Ly) <  cutoffsqr){
          if ( CheckOverlap(PartnerParticle, MovedParticle, Lx, Ly, &tmp) == ALOT) return ALOT;
        }
        if( MyIntStr>1e-10 ){
          deltaE += CheckForBond(MovedParticle     , PartnerParticle, Lx, Ly, Sys) 
                   -CheckForBond(Sys->Particles[p0], PartnerParticle, Lx, Ly, Sys);
        }
      }
    }
  }
  return deltaE;
}

double ParticleEnergy(sys *Sys, particle TestParticle){
  double E = 0.0;
  if( TestParticle.MyIntStr>1e-10 ){
    for(int i = 0; i<Sys->twomaxpop;i++){
      if(i == Sys->TotCurPop[0]) i=Sys->maxpop;
      if(i == Sys->TotCurPop[1] + Sys->maxpop) break;
      if(i!=Sys->p0){E += CheckForBond( TestParticle, Sys->Particles[i], Sys->Lx, Sys->Ly, Sys);}
    }
  }
  return E;
}

void AssignCells(sys *Sys){
  celllist *CellList = &(Sys->CellList);
  int i, j;
  for(int p=0;p<Sys->maxpop;p++){
    if(p == Sys->TotCurPop[0])break;
    
    Position_To_Index_PBC(&i, Sys->Particles[p].Pos.x, CellList->dCx, CellList->NumCellsx);
    Position_To_Index_PBC(&j, Sys->Particles[p].Pos.y, CellList->dCy, CellList->NumCellsy);
    
    CellList->BotParticleIds[i][j][CellList->BotPopulations[i][j]] = p;
    CellList->BotPopulations[i][j]++;
  }
  
  for(int p=Sys->maxpop;p<Sys->twomaxpop;p++){
    if(p == Sys->TotCurPop[1] + Sys->maxpop)break;
    
    Position_To_Index_PBC(&i, Sys->Particles[p].Pos.x, CellList->dCx, CellList->NumCellsx);
    Position_To_Index_PBC(&j, Sys->Particles[p].Pos.y, CellList->dCy, CellList->NumCellsy);
    
    CellList->TopParticleIds[i][j][CellList->TopPopulations[i][j]] = p;
    CellList->TopPopulations[i][j]++;
  }
}

void UpdateCellList(sys *Sys){
  celllist *CellList = &(Sys->CellList);
  int p0 = Sys->p0;
  
  int NewCx, NewCy, OldCx, OldCy, k;
  
  Position_To_Index_PBC(&NewCx, Sys->dummy_particle.Pos.x, CellList->dCx, CellList->NumCellsx);
  Position_To_Index_PBC(&NewCy, Sys->dummy_particle.Pos.y, CellList->dCy, CellList->NumCellsy);
  
  Position_To_Index_PBC(&OldCx, Sys->Particles[Sys->p0].Pos.x, CellList->dCx, CellList->NumCellsx);
  Position_To_Index_PBC(&OldCy, Sys->Particles[Sys->p0].Pos.y, CellList->dCy, CellList->NumCellsy);
  
  if(NewCx != OldCx || NewCy != OldCy){
    
    int *NewPopulation = (Sys->Particles[p0].MeshId == TOP) ? 
                         &(CellList->TopPopulations[NewCx][NewCy]): 
                         &(CellList->BotPopulations[NewCx][NewCy]);
    
    int *NewIds        = (Sys->Particles[p0].MeshId == TOP) ? 
                         CellList->TopParticleIds[NewCx][NewCy]: 
                         CellList->BotParticleIds[NewCx][NewCy];
    
    if(*NewPopulation==CellList->MaxPopCell){
      printf("Update CellList ERROR! (particles.c)\nTrying to jam too many particles into Cell\n");
      printf("%d %d\n%f %f\n%d %d\n", NewCx, NewCy, Sys->dummy_particle.Pos.x, Sys->dummy_particle.Pos.y, *NewPopulation, CellList->MaxPopCell);
      exit(1);
    }
                         
    NewIds[*NewPopulation] = Sys->p0;
    (*NewPopulation)++;
    
    int *OldPopulation = (Sys->Particles[p0].MeshId == TOP) ? 
                         &(CellList->TopPopulations[OldCx][OldCy]): 
                         &(CellList->BotPopulations[OldCx][OldCy]);
    
    int *OldIds        = (Sys->Particles[p0].MeshId == TOP) ? 
                         CellList->TopParticleIds[OldCx][OldCy]: 
                         CellList->BotParticleIds[OldCx][OldCy];
    
    for(k = 0;k<*OldPopulation;k++){
      if(OldIds[k] == Sys->p0){
        (*OldPopulation)--;
        OldIds[k] = OldIds[*OldPopulation];
        OldIds[*OldPopulation]=EMPTY;
        break; 
      }
    }
  }
  
}

void StoreParticleChange(sys *Sys, double DeltaEPhase, double DeltaEPart, double DeltaMPEnergy){
  int p0=Sys->p0;
  UpdateCellList(Sys);
  
  CopyVec(&(Sys->Particles[p0].Pos)     , Sys->dummy_particle.Pos );
  CopyVec(&(Sys->Particles[p0].nhat)    , Sys->dummy_particle.nhat);
  CopyVec(&(Sys->Particles[p0].MidPoint), Sys->dummy_particle.MidPoint );
  CopyVec(&(Sys->Particles[p0].EndPoint), Sys->dummy_particle.EndPoint );

  Sys->Particles[p0].MeshParticleEnergy = Sys->dummy_particle.MeshParticleEnergy; 
  Sys->Particles[p0].PhaseFactorEnergy  = Sys->dummy_particle.PhaseFactorEnergy; 
  
  Sys->TotPhaseFactorEnergy  += DeltaEPhase;
  Sys->TotParticleEnergy     += DeltaEPart;
  Sys->TotMeshParticleEnergy += DeltaMPEnergy;
}

int PerformParticleMove(sys *Sys){
  int acc=0;
  double DeltaEPhase   = SuggestParticleMove(Sys);
  if(DeltaEPhase<ALOT){
    double DeltaMPEnergy = ComputeDeltaMPEnergy(Sys, &(Sys->dummy_particle), Sys->p0);
    if(DeltaMPEnergy<ALOT){
      double DeltaEPart= CollisionTestAndParticleBonds(Sys);
//       double DeltaEPart = 0.0;
      double DeltaETot = DeltaEPhase + DeltaEPart + DeltaMPEnergy;
      if ( gsl_rng_uniform(Sys->rng)< exp(-DeltaETot) ){
        StoreParticleChange(Sys, DeltaEPhase, DeltaEPart, DeltaMPEnergy );
        Sys->TotEnergy+=DeltaETot;
        acc=1;
      }
    }
  }
//   ErrorCheckCellLists(Sys, 0);
  return acc;
}


// void ComputeMyNodes(sys *Sys, int p){
//   MeshPoint **MyMesh = (Sys->Particles[p].MeshId == TOP) ? Sys->TopMesh : Sys->BotMesh;
//   double MyRad = Sys->Particles[p].MyRad, MyRadSqr=Sys->Particles[p].MyRadSqr;
//   vec MyPos, diffvec;CopyVec(&(MyPos), Sys->Particles[p].Pos);
//   double Lx=Sys->Lx, Ly=Sys->Ly;
//   double dx=Sys->dx, dy=Sys->dy;
//   int Nx=Sys->Nx, Ny=Sys->Ny;
//   int i0, j0, iran, jran, i, j;
//   
//   Position_To_Index_PBC(&i0, MyPos.x, dx);Position_To_Index_PBC(&j0, MyPos.y, dy);
//   Distance_To_IndexRan(&iran, MyRad, dx);Distance_To_IndexRan(&jran, MyRad, dy);
//   int *ptr_nodelist= Sys->NodesInParticle[p];
//   Sys->NodeNums[p]=0;
//   
//   for(int idum=i0-iran; idum<i0+iran;idum++){
//     i=IndexPBC(idum, Nx);
//     for(int jdum=j0-jran; jdum<j0+jran;jdum++){
//       j=IndexPBC(jdum, Ny);
//       twodvecdiff( &diffvec, MyMesh[i][j].Pos, MyPos);
//       NearestImageConv(&diffvec, Lx, Ly);
//       
//       if( twod_NormSqr(diffvec)<MyRadSqr ){
//         ptr_nodelist[2*Sys->NodeNums[p]  ]=i;
//         ptr_nodelist[2*Sys->NodeNums[p]+1]=j;
//         MyMesh[i][j].ParticleOnNode=p;
//         Sys->NodeNums[p]++;
//       }
//     }
//   }
// }

/// Determine if particle Pokes out of membrane    
//   int i0, j0, iran, jran, i, j;
//   vec r12;
//   
//   Position_To_Index_PBC(&i0, MyPos.x, dx, Nx);Position_To_Index_PBC(&j0, MyPos.y, dy, Ny);
//   Distance_To_IndexRan(&iran, MPCutoff, dx);Distance_To_IndexRan(&jran, MPCutoff, dy);
//   
// //   printf("%lf %lf %d\n", MPCutoffsqr, MPCutoff, iran);exit(1);
//   
//   for(int idum=i0-iran; idum<i0+iran;idum++){
//     i=IndexPBC(idum, Nx);
//     for(int jdum=j0-jran; jdum<j0+jran;jdum++){
//       j=IndexPBC(jdum, Ny);
//       threedvecdiff(&r12, OtherMesh[i][j].Pos, MyPos);
//       NearestImageConv(&r12, Lx, Ly );
//       
//       if (twod_NormSqr(r12)< MPCutoffsqr){
//         dRsqr = threed_NormSqr(r12);
//         if( MyMeshId*r12.z < 0.0){MeshParticleEnergy+= ALOT;}
//         else if( dRsqr < MPCutoffsqr ){MeshParticleEnergy+=MeshParticleInteraction(dRsqr);}
//       }
//     }
//   }

// int ErrorCheckPhaseFactorEnergies(sys *Sys){
//   int i, p, inew, jnew, xNeigh, yNeigh, Err=0;
//   MeshPoint **MyMesh;
//   vec testnhat, diffvec, v1, v2 , v3, v4;
//   double PhaseEnergy, TotPhaseEnergy=0.0;
//   
//   for(i=0;i<2;i++){
//     for(p=i*Sys->maxpop;p<Sys->TotCurPop[i]+i*(Sys->maxpop);p++){
//       
//       if(i==0) MyMesh = Sys->BotMesh;
//       if(i==1) MyMesh = Sys->TopMesh;
//       
//       Position_To_Index_PBC(&inew, Sys->Particles[p].Pos.x, Sys->dx, Sys->Nx);
//       Position_To_Index_PBC(&jnew, Sys->Particles[p].Pos.y, Sys->dy, Sys->Ny);
//       
//       twodvecdiff(&(diffvec), Sys->Particles[p].Pos, MyMesh[inew][jnew].Pos);
//       NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
//       diffvec.x/=Sys->dx;
//       diffvec.y/=Sys->dy;
//       
//       xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
//       yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
//       
//       CopyVec(&v1, MyMesh[inew][jnew].nhat);     rescalevec(&v1, MyMesh[inew][jnew].nNorm    );
//       CopyVec(&v2, MyMesh[xNeigh][jnew].nhat);   rescalevec(&v2, MyMesh[xNeigh][jnew].nNorm  );
//       CopyVec(&v3, MyMesh[inew][yNeigh].nhat);   rescalevec(&v3, MyMesh[inew][yNeigh].nNorm  );
//       CopyVec(&v4, MyMesh[xNeigh][yNeigh].nhat); rescalevec(&v4, MyMesh[xNeigh][yNeigh].nNorm);
//       
//       InterpolateProperty(&(testnhat.x), diffvec, v1.x, v2.x, v3.x, v4.x );
//       InterpolateProperty(&(testnhat.y), diffvec, v1.y, v2.y, v3.y, v4.y );
//       InterpolateProperty(&(testnhat.z), diffvec, v1.z, v2.z, v3.z, v4.z );
//       PhaseEnergy = -0.5*log(1.0 
//                             + testnhat.x*testnhat.x
//                             + testnhat.y*testnhat.y);
//       TotPhaseEnergy+=PhaseEnergy;
//       
//       if(fabs(PhaseEnergy - Sys->Particles[p].PhaseFactorEnergy)>1e-4){
//         printf("EComp: %lf EStored: %lf EDiff: %lf\n", PhaseEnergy, Sys->Particles[p].PhaseFactorEnergy, PhaseEnergy - Sys->Particles[p].PhaseFactorEnergy);
//         printf("Id: %d\n", p);
//         PrintVec(Sys->Particles[p].Pos,"Position");
//         printf("PhaseError!\n");Err=5;
//       }
//     }
//   }
//   
//   if(fabs(TotPhaseEnergy - Sys->TotPhaseFactorEnergy)>1e-4){printf("Error Mis total %lf %lf\n", TotPhaseEnergy, Sys->TotPhaseFactorEnergy);printf("p1: %d p2: %d\n", Sys->TotCurPop[0], Sys->TotCurPop[1]);Err=5;}
//   
//   return Err;
// }

// double ReturnParticleHeight(sys *Sys, int p){
//   MeshPoint **MyMesh = (Sys->Particles[p].MeshId == TOP) ? Sys->TopMesh : Sys->BotMesh;
//   particle *My_Particle = &(Sys->Particles[p]);
//   
//   vec diffvec;
//   int inew, jnew, xNeigh, yNeigh;
//   double znew;
//   
//   Position_To_Index_PBC(&inew, My_Particle->Pos.x, Sys->dx, Sys->Nx);
//   Position_To_Index_PBC(&jnew, My_Particle->Pos.y, Sys->dy, Sys->Ny);
//   
//   twodvecdiff(&(diffvec), My_Particle->Pos, MyMesh[inew][jnew].Pos);
//   NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
//   diffvec.x/=Sys->dx;
//   diffvec.y/=Sys->dy;
//   
//   xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
//   yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
//   
//   InterpolateProperty(&(znew), diffvec,
//                                MyMesh[inew][jnew].Pos.z, 
//                                MyMesh[xNeigh][jnew].Pos.z, 
//                                MyMesh[inew][yNeigh].Pos.z, 
//                                MyMesh[xNeigh][yNeigh].Pos.z 
//                      );
//   return znew;
// }

/*
void Dummy_z_nhat_PhaseFactor(sys *Sys){
  double length = Sys->dummy_particle.MyLength;
  particle *dummy_particle = &(Sys->dummy_particle);
  MeshPoint **MyMesh = (dummy_particle->MeshId == TOP) ? Sys->TopMesh : Sys->BotMesh;
  
  
  vec diffvec, v1, v2 , v3, v4;
  int inew, jnew, xNeigh, yNeigh;
  
  Position_To_Index_PBC(&inew, dummy_particle->Pos.x, Sys->dx, Sys->Nx);
  Position_To_Index_PBC(&jnew, dummy_particle->Pos.y, Sys->dy, Sys->Ny);
  
  twodvecdiff(&(diffvec), dummy_particle->Pos, MyMesh[inew][jnew].Pos);
  PosPBC(&(diffvec), Sys->Lx, Sys->Ly);
//   diffvec.x/=Sys->dx;
//   diffvec.y/=Sys->dy;
    
  xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
  yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
  
  dummy_particle->Pos.z = MyMesh[inew][jnew].Pos.z 
                        - diffvec.x*(MyMesh[inew][jnew].Pos.z - MyMesh[xNeigh][jnew].Pos.z)
                        - diffvec.y*(MyMesh[inew][jnew].Pos.z - MyMesh[inew][yNeigh].Pos.z)
                        + diffvec.x*diffvec.y*(
                          MyMesh[inew][jnew].Pos.z 
                          - MyMesh[xNeigh][jnew].Pos.z - MyMesh[inew][yNeigh].Pos.z
                          + MyMesh[xNeigh][yNeigh].Pos.z
                         );
  
  CopyVec(&v1, MyMesh[inew][jnew].nhat);     rescalevec(&v1, MyMesh[inew][jnew].nNorm    );
  CopyVec(&v2, MyMesh[xNeigh][jnew].nhat);   rescalevec(&v2, MyMesh[xNeigh][jnew].nNorm  );
  CopyVec(&v3, MyMesh[inew][yNeigh].nhat);   rescalevec(&v3, MyMesh[inew][yNeigh].nNorm  );
  CopyVec(&v4, MyMesh[xNeigh][yNeigh].nhat); rescalevec(&v4, MyMesh[xNeigh][yNeigh].nNorm);
  
  dummy_particle->nhat.x =  v1.x
                          - diffvec.x*(v1.x - v2.x)
                          - diffvec.y*(v1.x - v3.x)
                          + diffvec.x*diffvec.y*(v1.x - v2.x - v3.x + v4.x);
                        
  dummy_particle->nhat.y =  v1.y
                          - diffvec.x*(v1.y - v2.y)
                          - diffvec.y*(v1.y - v3.y)
                          + diffvec.x*diffvec.y*(v1.y - v2.y - v3.y + v4.y);
  
  dummy_particle->nhat.z =  v1.z
                          - diffvec.x*(v1.z - v2.z)
                          - diffvec.y*(v1.z - v3.z)
                          + diffvec.x*diffvec.y*(v1.z - v2.z - v3.z + v4.z);                   
  
  threed_NormalizeVec(&dummy_particle->nhat);
  threedvecsumscalar( &(dummy_particle->MidPoint), dummy_particle->Pos, dummy_particle->nhat, length/2);
  threedvecsumscalar( &(dummy_particle->EndPoint), dummy_particle->Pos, dummy_particle->nhat, length  );
  
  Sys->dummy_particle.PhaseFactorEnergy = -0.5*log(1.0 
                                                 + dummy_particle->nhat.x*dummy_particle->nhat.x
                                                 + dummy_particle->nhat.y*dummy_particle->nhat.y);
}
  
    //Check if points of partner are inside of moved particle
  if( Partner_Verts[k].x < r_moved     && Partner_Verts[k].x > -r_moved && 
      Partner_Verts[k].y < l_moved_on2 && Partner_Verts[k].x > -l_moved_on2 ) return ALOT;
    
  minx=min(minx, Partner_Verts[k].x);maxx=max(maxx, Partner_Verts[k].x);
  miny=min(miny, Partner_Verts[k].y);maxy=max(maxy, Partner_Verts[k].y);
  
 
  //Check if rectangles are too far away to intersect
  if(maxx < -r_moved || minx > r_moved || maxy < -l_moved_on2 || miny > l_moved_on2) return 0.0;*/
