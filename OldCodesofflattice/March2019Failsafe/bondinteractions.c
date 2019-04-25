#include "bondinteractions.h"

double BondInteractions(particle TestParticle, particle PartnerParticle, double Lx, double Ly, double IntRanSqr, double IntStr){
  double E=0.0;
  if(deltaRsqr(TestParticle.EndPoint, PartnerParticle.EndPoint, Lx, Ly)< IntRanSqr)E=-IntStr;
  return E;
}

void PrintBonds(sys *Sys){
  int p0, p1;
  double Lx = Sys->Lx, Ly = Sys->Ly, IntStr, IntCutoffsqr;
  
  printf("\nParticleBonds:\n");
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    IntStr = Sys->Particles[p0].MyIntStr;IntCutoffsqr = Sys->Particles[p0].MyIntCutoffsqr;
    printf("Particle %d --> Particle %d (%d) E: %lf rsqrC: %lf C: %lf\n", p0, p1, p1-Sys->maxpop,
      p1!=-1 ? BondInteractions(Sys->Particles[p0], Sys->Particles[p1], Lx, Ly, IntCutoffsqr, IntStr):0.0,
      p1!=-1 ? deltaRsqr(Sys->Particles[p0].EndPoint, Sys->Particles[p1].EndPoint, Lx, Ly):-1.0,
      IntCutoffsqr
    ); 
  }
  
  for(p0=Sys->maxpop;p0<Sys->maxpop+Sys->TotCurPop[0];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    IntStr = Sys->Particles[p0].MyIntStr;IntCutoffsqr = Sys->Particles[p0].MyIntCutoffsqr;
    printf("Particle %d (%d) --> Particle %d E: %lf rsqrC: %lf C: %lf\n", p0, p0-Sys->maxpop, p1,
      p1!=-1 ? BondInteractions(Sys->Particles[p0], Sys->Particles[p1], Lx, Ly, IntCutoffsqr, IntStr):0.0,
      p1!=-1 ? deltaRsqr(Sys->Particles[p0].EndPoint, Sys->Particles[p1].EndPoint, Lx, Ly):-1.0,
      IntCutoffsqr
    ); 
  }
}

void PrintBondList(particlebonds Bonds, int p){
 printf("particle %d nums %d currbond %d\n", p, Bonds.num, Bonds.currbond);
 for(int i=0;i<BONDNUMS;i++) printf(" %d %lf\n", Bonds.Ids[i], Bonds.rsqr[i]);
}

void SortBondList(particlebonds *BondList, int n){
    // Base case
    if (n == 1)
        return;
    // One pass of bubble sort. After
    // this pass, the largest element
    // is moved (or bubbled) to end.
    for (int i=0; i<n-1; i++)
        if (BondList->rsqr[i] > BondList->rsqr[i+1]){
          double tmprsqr      = BondList->rsqr[i];
          int tmpindex        = BondList->Ids[i];
          BondList->rsqr[i]   = BondList->rsqr[i+1];
          BondList->Ids[i]    = BondList->Ids[i+1];
          BondList->rsqr[i+1] = tmprsqr;
          BondList->Ids[i+1]  = tmpindex;
        }
    // Largest element is fixed,
    // recur for remaining array
    SortBondList(BondList, n-1);
}

void CopyBondList(particlebonds *Target, particlebonds *Source){
  Target->num=Source->num;Target->currbond = Source->currbond;
  for(int i=0;i<BONDNUMS;i++){
    Target->rsqr[i]= Source->rsqr[i];
    Target->Ids[i ]= Source->Ids[i];
  }
}

void InitBondList(particlebonds *BondList){
  BondList->num=0; //BondList->currbond=-1;
  for(int i=0;i<BONDNUMS;i++){
    BondList->rsqr[i]=ALOT;
    BondList->Ids[i]=-1;
  }
}

void PopulateTmpBondList(particlebonds *tmpBondList, particle *MovedParticle, int currbond, sys *Sys){
  InitBondList(tmpBondList);
  tmpBondList->currbond = currbond;
  
  particle *PartnerParticle;
  celllist *CellList = &(Sys->CellList);
  int MyCx, MyCy, Cx, Cy, Populations, *Ids, k;
  int mytype = MovedParticle->type, tmpnum;
  double rsqr, Lx=Sys->Lx, Ly=Sys->Ly, cutoffsqr = MovedParticle->MyIntCutoffsqr;
  
  Position_To_Index_PBC(&MyCx, MovedParticle->Pos.x, CellList->dCx, CellList->NumCellsx);
  Position_To_Index_PBC(&MyCy, MovedParticle->Pos.y, CellList->dCy, CellList->NumCellsy);
  
  for(int i=0;i<9;i++){
    Cx = CellList->Neighbors[MyCx][MyCy][2*i];
    Cy = CellList->Neighbors[MyCx][MyCy][2*i + 1];
    Populations = (MovedParticle->MeshId == TOP ? CellList->BotPopulations[Cx][Cy] : CellList->TopPopulations[Cx][Cy]);
    Ids         = (MovedParticle->MeshId == TOP ? CellList->BotParticleIds[Cx][Cy] : CellList->TopParticleIds[Cx][Cy]);
    
    for(k=0;k<Populations;k++){
      PartnerParticle = &(Sys->Particles[Ids[k]]);
      if(PartnerParticle->type == mytype){
        rsqr = deltaRsqr(PartnerParticle->EndPoint, MovedParticle->EndPoint, Lx, Ly);
        if(rsqr<cutoffsqr){
          tmpnum = tmpBondList->num;
          tmpBondList->Ids[tmpnum] = Ids[k];
          tmpBondList->rsqr[tmpnum] = rsqr;
          tmpBondList->num++;
        }
      }
    }
  }
  SortBondList(tmpBondList, BONDNUMS);
}

void ComputeAllParticleBondLists(sys *Sys){
  particlebonds *ParticleBonds = Sys->BondList.ParticleBonds;
  particle *Particles = Sys->Particles, PartnerParticle;
  celllist *CellList = &(Sys->CellList);
  int MyCx, MyCy, Cx, Cy, Populations, *Ids, k, num;
  int mytype, meshId, p0, p1;
  double IntRanSqr, rsqr, Lx=Sys->Lx, Ly=Sys->Ly;
  
  for(p0=0;p0<Sys->maxpop + Sys->TotCurPop[1];p0++){
    if(p0==Sys->TotCurPop[0]) p0=Sys->maxpop;
    if(Particles[p0].MyIntStr>1e-6)InitBondList(&(ParticleBonds[p0]));
  }
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    if(Particles[p0].MyIntStr>1e-6){
      mytype = Particles[p0].type;IntRanSqr = Particles[p0].MyIntCutoffsqr;
      Position_To_Index_PBC(&MyCx, Particles[p0].Pos.x, CellList->dCx, CellList->NumCellsx);
      Position_To_Index_PBC(&MyCy, Particles[p0].Pos.y, CellList->dCy, CellList->NumCellsy);
      meshId = Particles[p0].MeshId;
      
      for(int i=0;i<9;i++){
        Cx = CellList->Neighbors[MyCx][MyCy][2*i];
        Cy = CellList->Neighbors[MyCx][MyCy][2*i + 1];
        Populations = (meshId == TOP ? CellList->BotPopulations[Cx][Cy] : CellList->TopPopulations[Cx][Cy]);
        Ids         = (meshId == TOP ? CellList->BotParticleIds[Cx][Cy] : CellList->TopParticleIds[Cx][Cy]);
          
        for(k=0;k<Populations;k++){
          p1 = Ids[k];
          PartnerParticle = Sys->Particles[p1];
          if(PartnerParticle.type == mytype){
            rsqr = deltaRsqr(PartnerParticle.EndPoint, Particles[p0].EndPoint, Lx, Ly);
            if(rsqr<IntRanSqr){
              num=ParticleBonds[p0].num;
              ParticleBonds[p0].Ids[num]  = p1;
              ParticleBonds[p0].rsqr[num] = rsqr;
              ParticleBonds[p0].num++;
              
              num=ParticleBonds[p1].num;
              ParticleBonds[p1].Ids[num]  = p0;
              ParticleBonds[p1].rsqr[num] = rsqr;
              ParticleBonds[p1].num++;
            }
          }
        }
      }
      SortBondList(&(ParticleBonds[p0]), BONDNUMS);
    }
  }
  
  for(p1=Sys->maxpop;p1<Sys->maxpop+Sys->TotCurPop[1];p1++){
    if(Particles[p1].MyIntStr>1e-6)SortBondList(&(ParticleBonds[p1]), BONDNUMS);
  }
}

int CreateBondList(sys *Sys){
  int interactions=0;
  for(int i=0;i<Sys->NumPTypes;i++) interactions+=(Sys->BindPType[i]>1e-10 ? 1:0);
  
  Sys->BondList.Curbondnums=(double *)malloc(interactions*sizeof(double));
  for(int p=0;p<Sys->NumPTypes;p++) Sys->BondList.Curbondnums[p] = 0.0;
  
  Sys->BondList.ParticleBonds = (particlebonds *)malloc(Sys->twomaxpop*sizeof(particlebonds));
  particlebonds *BondList = Sys->BondList.ParticleBonds;
  for(int p=0;p<Sys->twomaxpop;p++){
    InitBondList(&(BondList[p]));BondList[p].currbond=-1;
  }
  
  Sys->BondList.BondableTypeNum = interactions;
  Sys->BondList.BondProbs       = (double**)malloc(2*sizeof(double *));
  for(int i=0;i<2;i++) Sys->BondList.BondProbs[i] = (double *)malloc(interactions*sizeof(double));
  
  Sys->BondList.BondableTypes = (double *)malloc(interactions*sizeof(double));
  int tmp=0;
  for(int i=0;i<Sys->NumPTypes;i++){
   if(Sys->BindPType[i]>1e-10) Sys->BondList.BondableTypes[tmp] = i;
   tmp++;
  }
  
  return interactions;
}

void OverwriteOldBondLists(int paffected, int pOld, int i_affected, sys *Sys){
  particlebonds *BondList = Sys->BondList.ParticleBonds;
  particlebonds *partnerbondlist;
  particlebonds *newbondlist = &(Sys->tmpBondList[i_affected]);
  int i, j;
  
  particlebonds *oldbondlist = &(BondList[pOld]);
  for(i=0;i<oldbondlist->num;i++){
    partnerbondlist = &(BondList[oldbondlist->Ids[i]]);
    for(j=0;j<partnerbondlist->num;j++){
      if(partnerbondlist->Ids[j] == pOld){
        partnerbondlist->Ids[j] = -1;
        partnerbondlist->rsqr[j] = ALOT;
        partnerbondlist->num--;
        break;
      }
    }
    SortBondList(partnerbondlist, BONDNUMS);
  }
  
  if(paffected!=pOld){
    oldbondlist = &(BondList[paffected]);
    for(i=0;i<oldbondlist->num;i++){
      partnerbondlist = &(BondList[oldbondlist->Ids[i]]);
      for(j=0;j<partnerbondlist->num;j++){
        if(partnerbondlist->Ids[j] == paffected){
          partnerbondlist->Ids[j] = -1;
          partnerbondlist->rsqr[j] = ALOT;
          partnerbondlist->num--;
          break;
        }
      }
      SortBondList(partnerbondlist, BONDNUMS);
    }
  }
  
  for(i=0;i<newbondlist->num;i++){
    partnerbondlist = &(BondList[newbondlist->Ids[i]]);
    partnerbondlist->Ids[partnerbondlist->num]  = paffected;
    partnerbondlist->rsqr[partnerbondlist->num] = newbondlist->rsqr[i];
    partnerbondlist->num++;
    SortBondList(partnerbondlist, BONDNUMS);
  }
  CopyBondList(&BondList[paffected], newbondlist);
}

void UpdateBondLists_NotEmpty(sys *Sys, int opt){
//   printf("%d b4: %d ", opt, Sys->BondList.ParticleBonds[1036].currbond);
  if(opt==0 && Sys->dummy_particle.MyIntStr>1e-6)
    OverwriteOldBondLists(Sys->p0, Sys->p0, 0, Sys);
  
  if(opt==1){
    int current_p;
    for(int i_affect=0;i_affect<Sys->AffectedParticleNum;i_affect++){
      current_p=Sys->Unique_dummy_index[i_affect];
      OverwriteOldBondLists(current_p, current_p, i_affect, Sys);
    }
  }
  
  if(opt==2) ComputeAllParticleBondLists(Sys);
  
  if(opt==3){
    if(Sys->GC_AddRemove==0){
      if(Sys->dummy_particle.MyIntStr>1e-6){
        InitBondList(&(Sys->tmpBondList[0]));Sys->tmpBondList[0].currbond=-1;
        PopulateTmpBondList(&(Sys->tmpBondList[0]), &(Sys->dummy_particle), -1, Sys);
        OverwriteOldBondLists(Sys->p0, Sys->p0, 0, Sys);
      }
    }
    if(Sys->GC_AddRemove==1){
      int meshindex = ((Sys->GC_TopOrBot == BOT) ? 0 : 1);
      int replaceindex = Sys->maxpop*meshindex + Sys->TotCurPop[meshindex]-1;
      particlebonds *ReplaceList = &(Sys->BondList.ParticleBonds[replaceindex]);
      
      int currbond = ReplaceList->currbond;
      CopyBondList(&(Sys->tmpBondList[0]), ReplaceList );
      if(Sys->p0==replaceindex) Sys->tmpBondList[0].num=0;

      OverwriteOldBondLists(Sys->p0, replaceindex, 0, Sys);
      if(currbond!=-1) Sys->BondList.ParticleBonds[currbond].currbond = Sys->p0;
      
      InitBondList(ReplaceList);ReplaceList->currbond=-1;
    }
  }
}

void UpdateBondLists_Empty(sys *Sys, int opt){}

int ErrorCheck_Bonds_Empty(sys *Sys){return 0;}

void ComputeParticleBonds_Empty(sys *Sys){}

double BondParticle_Empty(sys *Sys){return 0.0;}

double BondMeshMove_Empty(sys *Sys){return 0.0;}

double BondMeshShift_Empty(sys *Sys, double suggestedz){return 0.0;}

double BondGCMC_Empty(sys *Sys){return 0.0;}

void BondMove_Empty(sys *Sys){}


int ErrorCheck_Bonds_RanNetLR(sys *Sys){
  particle *Particles=Sys->Particles;
  double TmpTotalEnergy=0.0, Lx=Sys->Lx, Ly=Sys->Ly;
  double DeltaE, *tmpbonds, IntStr, IntCutoffsqr;
  int p0, p1, Err=0;
  tmpbonds = (double *)malloc(Sys->NumPTypes*sizeof(double));
  for(p0=0;p0<Sys->NumPTypes;p0++)tmpbonds[p0]=0.0;
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    if(p1!=-1){
      IntStr = Particles[p0].MyIntStr;
      IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
      
      DeltaE = BondInteractions(Particles[p0], Particles[p1], Lx, Ly, IntCutoffsqr, IntStr);
      if (DeltaE>1e-6)tmpbonds[Particles[p0].type]++;
      TmpTotalEnergy += DeltaE;
      if( (Sys->BondList.ParticleBonds[p0].currbond != p1) || 
          (Sys->BondList.ParticleBonds[p1].currbond != p0)
      ){
        printf("\nError! List integrity not given, p0: %d connected to %d, but\np1: %d connected to %d\n",
              p0, Sys->BondList.ParticleBonds[p0].currbond, p1, Sys->BondList.ParticleBonds[p1].currbond
        );
        Err=5;
      }
    }
  }
  
  for(p0=Sys->maxpop;p0<Sys->maxpop+Sys->TotCurPop[1];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    if(p1!=-1){
      if( (Sys->BondList.ParticleBonds[p0].currbond != p1) || 
          (Sys->BondList.ParticleBonds[p1].currbond != p0)
      ){
        printf("\nError! List integrity not given, p0: %d connected to %d, but\np1: %d connected to %d\n",
              p0, Sys->BondList.ParticleBonds[p0].currbond, p1, Sys->BondList.ParticleBonds[p1].currbond
        );
        Err=5;
      }
    }
  }
  
  
  if(fabs(Sys->TotParticleEnergy-TmpTotalEnergy)>1e-6){
    printf("\nLastMoved: %d\nError in bondEnergies computed %lf stored %lf\n", Sys->p0, TmpTotalEnergy, Sys->TotParticleEnergy);
    Err=5;
  }
//   for(p0=0;p0<Sys->NumPTypes;p0++){
//     if(fabs(tmpbonds[p0]-Sys->BondList.Curbondnums[p0])>1e-6){
//       printf("\nError in bond numbers computed %lf stored %lf particle type %d\n", tmpbonds[p0], Sys->BondList.Curbondnums[p0], p0);
//       Err=5;
//     }
//   }
  return Err;
}

void ComputeParticleBonds_RanNetLR(sys *Sys){
  int p0, p1;
  double Lx=Sys->Lx, Ly=Sys->Ly;
  double IntCutoffsqr, IntStr;
  particle *Particles = Sys->Particles;
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    if(p1!=-1){
      IntStr = Particles[p0].MyIntStr;
      IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
      if (deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly)<IntCutoffsqr){
        Sys->BondList.Curbondnums[Particles[p0].type]++;
        Sys->TotParticleEnergy+=IntStr;
      }
    }
  }
}

void BondMove_RanNetLR(sys *Sys){
  int inttype, type, choice, NTop, NBot;
  int p0, p1, p0p, p1p;
  double IntStr, IntCutoffsqr, deltaE;
  double Lx = Sys->Lx, Ly = Sys->Ly;
  particle *Particles=Sys->Particles;
  particlebonds *ParticleBonds = Sys->BondList.ParticleBonds;
  
  for(int counter = 0;counter<50;counter++){
    inttype = gsl_rng_uniform_int(Sys->rng, Sys->BondList.BondableTypeNum);choice=0;
    for(type=0;type<Sys->NumPTypes;type++){if(choice==inttype){break;}choice+=(Sys->BindPType[type]>1e-10 ?1:0);}
    NBot = Sys->CurParPop[0][type];
    NTop = Sys->CurParPop[1][type];
    
    if(NTop>0 && NBot>0){
      IntStr=Sys->BindPType[type];
      IntCutoffsqr = Sys->RanPType[type]*Sys->RanPType[type];
      
      choice=gsl_rng_uniform_int(Sys->rng, NBot);
      for(p0=0;p0<Sys->TotCurPop[0];p0++){
        if( Particles[p0].type==type){
          if(inttype==choice)break;
          else inttype++;
        }
      }
      
      choice=gsl_rng_uniform_int(Sys->rng, NTop);
      for(p1=Sys->maxpop;p1<Sys->TotCurPop[1];p1++){
        if( Particles[p1].type==type){
          if(inttype==choice)break;
          else inttype++;
        }
      }
      
//       printf("%d %d %d %d\n", p0, p1, ParticleBonds[p0].currbond, ParticleBonds[p1].currbond);
      if(ParticleBonds[p0].currbond==-1 && ParticleBonds[p1].currbond==-1){
        deltaE = BondInteractions(Sys->Particles[p0], Sys->Particles[p1], Lx, Ly, IntCutoffsqr, IntStr);
        if(gsl_rng_uniform(Sys->rng)< exp(-deltaE)){
          ParticleBonds[p1].currbond = p0;
          ParticleBonds[p0].currbond = p1;
//           Sys->BondList.Curbondnums[type]++;
          Sys->TotParticleEnergy+=deltaE;
          Sys->TotEnergy+=deltaE;
        }
      }
      
      else if(ParticleBonds[p1].currbond==ParticleBonds[p0].currbond){
        deltaE =-BondInteractions(Sys->Particles[p0], Sys->Particles[p1], Lx, Ly, IntCutoffsqr, IntStr);
        if(gsl_rng_uniform(Sys->rng)< exp(-deltaE)){
          ParticleBonds[p1].currbond = -1;
          ParticleBonds[p0].currbond = -1;
//           Sys->BondList.Curbondnums[type]--;
          Sys->TotParticleEnergy+=deltaE;
          Sys->TotEnergy+=deltaE;
        }
      }
      
      else {
        p0p = ParticleBonds[p0].currbond;
        p1p = ParticleBonds[p1].currbond;
        deltaE =   BondInteractions(Particles[p0], Particles[p1], Lx, Ly, IntCutoffsqr, IntStr)
                 +(p1p ==-1 || p0p ==-1 ? 0.0:BondInteractions(Particles[p0p], Particles[p1p], Lx, Ly, IntCutoffsqr, IntStr))
                 -(p0p ==-1 ? 0.0:BondInteractions(Particles[p0], Particles[p0p], Lx, Ly, IntCutoffsqr, IntStr))
                 -(p1p ==-1 ? 0.0:BondInteractions(Particles[p1], Particles[p1p], Lx, Ly, IntCutoffsqr, IntStr));
//         printf("Switch p0 %d p1 %d and p0p %d p1p %d\n", p0, p1, p0p, p1p);
//         printf("deltaE %lf test %lf\n", deltaE, Sys->TotParticleEnergy);
        if(gsl_rng_uniform(Sys->rng)< exp(-deltaE)){
          ParticleBonds[p0].currbond   = p1;
          ParticleBonds[p1].currbond   = p0;
          if(p0p!=-1) ParticleBonds[p0p].currbond = p1p;
          if(p1p!=-1) ParticleBonds[p1p].currbond = p0p;
          Sys->TotParticleEnergy+=deltaE;
          Sys->TotEnergy+=deltaE;
        }
      }
    }
  }
//   printf("%lf %lf\n", Sys->TotParticleEnergy, Sys->TotEnergy);
}

double BondParticle_RanNetLR(sys *Sys){
  particle MovedParticle = Sys->dummy_particle;
  double deltaE=0.0;
  int p1;
  
  if(MovedParticle.MyIntStr>1e-10 && Sys->BondList.ParticleBonds[Sys->p0].currbond!=-1){
    p1 = Sys->BondList.ParticleBonds[Sys->p0].currbond;
    deltaE = BondInteractions(Sys->Particles[p1], MovedParticle, Sys->Lx, Sys->Ly, MovedParticle.MyIntCutoffsqr, MovedParticle.MyIntStr) 
            -BondInteractions(Sys->Particles[p1], Sys->Particles[Sys->p0], Sys->Lx, Sys->Ly, MovedParticle.MyIntCutoffsqr, MovedParticle.MyIntStr);
  }
  return deltaE;
}

double BondMeshMove_RanNetLR(sys *Sys){
  double Lx = Sys->Lx, Ly = Sys->Ly;
  int i_affect, current_p;
  double DeltaE=0.0, IntStr, IntRanSqr;
  particle current_Particle, partner_Particle;
  
  //This routine is fine, all affected particles (at the top or bottom), only interact with unperturbed particles below.
  //So there is no need to check whether the particles are virtual or real.
  
  for(i_affect=0;i_affect<Sys->AffectedParticleNum;i_affect++){
    current_p=Sys->Unique_dummy_index[i_affect];
    current_Particle = Sys->mesh_dummy_particles[i_affect];
    if(current_Particle.MyIntStr>1e-10 && Sys->BondList.ParticleBonds[current_p].currbond!=-1){
      IntStr = current_Particle.MyIntStr;
      IntRanSqr = current_Particle.MyIntCutoffsqr;
      partner_Particle = Sys->Particles[Sys->BondList.ParticleBonds[current_p].currbond];
      DeltaE+= BondInteractions(current_Particle,          partner_Particle, Lx, Ly, IntRanSqr, IntStr)
              -BondInteractions(Sys->Particles[current_p], partner_Particle, Lx, Ly, IntRanSqr, IntStr);
    }
  }
  return DeltaE;
}

double BondMeshShift_RanNetLR(sys *Sys, double suggestedz){
  double DeltaE=0.0, EOld, IntStr, IntRanSqr, Lx=Sys->Lx, Ly=Sys->Ly;
  int pTop, pBot;
  particle CurrentParticle, PartnerParticle;
  
  for(pTop=Sys->maxpop;pTop<Sys->maxpop + Sys->TotCurPop[1];pTop++){
    if(Sys->BondList.ParticleBonds[pTop].currbond!=-1){
      pBot = Sys->BondList.ParticleBonds[pTop].currbond;
      PartnerParticle = Sys->Particles[pBot];
      CurrentParticle = Sys->Particles[pTop];
      IntStr = CurrentParticle .MyIntStr;
      IntRanSqr = CurrentParticle .MyIntCutoffsqr;
      
      EOld = BondInteractions(CurrentParticle, PartnerParticle, Lx, Ly, IntRanSqr, IntStr);
      CurrentParticle.Pos.z +=suggestedz;
      DeltaE+= BondInteractions(CurrentParticle, PartnerParticle, Lx, Ly, IntRanSqr, IntStr)
              -EOld;
      CurrentParticle.Pos.z -=suggestedz;
    }
  }
  return DeltaE;
}

double BondGCMC_RanNetLR(sys *Sys){
  double DeltaE=0.0;
  if(Sys->GC_AddRemove==1 && Sys->BondList.ParticleBonds[Sys->p0].currbond!=-1){
    int p0=Sys->p0, p1=Sys->BondList.ParticleBonds[Sys->p0].currbond;
    DeltaE = BondInteractions(
                Sys->Particles[p0], Sys->Particles[p1], 
                Sys->Lx, Sys->Ly, Sys->Particles[p0].MyIntCutoffsqr, 
                Sys->Particles[p0].MyIntStr
              );
//     printf("Removing bonded particle %d %lf\n", );
  }
  return DeltaE;
}



int ErrorCheck_Bonds_RanNetSR(sys *Sys){
  particle *Particles=Sys->Particles;
  particlebonds *ParticleBonds = Sys->BondList.ParticleBonds;
  double TmpTotalEnergy=0.0, Lx=Sys->Lx, Ly=Sys->Ly;
  double DeltaE, *tmpbonds, IntStr, IntCutoffsqr;
  int p0, p1, Err=0, k, l, tmp;
  tmpbonds = (double *)malloc(Sys->NumPTypes*sizeof(double));
  for(p0=0;p0<Sys->NumPTypes;p0++)tmpbonds[p0]=0.0;
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1=ParticleBonds[p0].currbond;
    if(p1!=-1 && (p1<Sys->maxpop || p1>Sys->maxpop + Sys->TotCurPop[1])){
      printf("\nError! Particle %d is bonded weirdly to %d maxtoppop %d\n", p0, p1, Sys->maxpop + Sys->TotCurPop[1]);
      return 5;
    }
  }
  
  for(p0=Sys->maxpop;p0<Sys->maxpop+Sys->TotCurPop[1];p0++){
    p1=ParticleBonds[p0].currbond;
    if(p1!=-1 && p1>Sys->TotCurPop[0]){
      printf("\nError! Particle %d is bonded weirdly to %d maxbotpop %d\n", p0, p1, Sys->TotCurPop[0]);
      return 5;
    }
  }
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1 = ParticleBonds[p0].currbond;
    if(p1!=-1){
      IntStr = Particles[p0].MyIntStr;
      IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
      tmpbonds[Particles[p0].type]++;
      
      if(deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly)>IntCutoffsqr){
        printf("\nError! Paricles %d %d have a distance of\n %lf vs %lf\nshould not be in each others' list\n", p0, p1, deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly), IntCutoffsqr );
        Err=5;
      }
      
      DeltaE = BondInteractions(Particles[p0], Particles[p1], Lx, Ly, IntCutoffsqr, IntStr);
      TmpTotalEnergy += DeltaE;
      if( (ParticleBonds[p0].currbond != p1) || 
          (ParticleBonds[p1].currbond != p0)
      ){
        printf("\nError! List integrity not given, p0: %d bonded to %d, but\np1: %d bonded to %d\n",
              p0, ParticleBonds[p0].currbond, p1, ParticleBonds[p1].currbond
        );
        Err=5;
      }
    }
  }
  
  if(fabs(Sys->TotParticleEnergy-TmpTotalEnergy)>1e-6){
    printf("\nLastMoved: %d\nError in bondEnergies computed %lf stored %lf\n", Sys->p0, TmpTotalEnergy, Sys->TotParticleEnergy);
    Err=5;
  }
  //Not incrementing curbondnums
//   for(p0=0;p0<Sys->NumPTypes;p0++){
//     if(fabs(tmpbonds[p0]-Sys->BondList.Curbondnums[p0])>1e-6){
//       printf("\nError in bond numbers computed %lf stored %lf particle type %d\n", tmpbonds[p0], Sys->BondList.Curbondnums[p0], p0);
//       Err=5;
//     }
//   }

  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    IntCutoffsqr = Particles[p0].MyIntCutoffsqr;tmp=0;
    for(p1=Sys->maxpop;p1<Sys->maxpop+Sys->TotCurPop[1];p1++){
      if(deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly)<IntCutoffsqr){
        l=0;
        for(k=0;k<ParticleBonds[p0].num;k++){
          if(ParticleBonds[p0].Ids[k]==p1){
            l=1;break;
          }
        }
        if(l==0){
          printf("Error! for particle %d, particle %d is in range %lf cutoff %lf\n", p0, p1, deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly), IntCutoffsqr);
          printf("but the Id is not in bondlist:\n");
          PrintBondList(ParticleBonds[p0], p0);
          Err=5;
        }
        tmp++;
        if(tmp==BONDNUMS){
          printf("Error! too many potential bond partners for particle %d partner %d %d\nAdjust BONDNUMS variable!\n", p0, p1, tmp);
          Err=5;
        }
      }
    }
  }
    

  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
    for(k=0;k<ParticleBonds[p0].num;k++){
      p1 = ParticleBonds[p0].Ids[k];
      if(deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly)>IntCutoffsqr){
        PrintBonds(Sys);
        printf("\nBondlists erronous\np0 %d p1 %d dist %lf cutoff %lf\n", p0, p1, deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly), IntCutoffsqr);
        Err=5;
      }
      
      if(Particles[p0].type != Particles[p1].type){
        PrintBonds(Sys);
        printf("\nBondlists erronous\np0 %d has type %d wheras p1 %d has type %d\n", p0, Particles[p0].type, p1, Particles[p1].type);
        Err=5;
      }
      
      tmp = 0;
      for(l=0;l<ParticleBonds[p1].num;l++){
        if(ParticleBonds[p1].Ids[l]==p0) tmp=1;
      }
      if(tmp==0){
        printf("\nBondList reflexivity broken! p0 %d p1 %d\n", p0, p1);
        PrintBondList(ParticleBonds[p0], p0);
        PrintBondList(ParticleBonds[p1], p1);
        Err=5;
      }
    }
  }
  
  for(p0=Sys->maxpop;p0<Sys->maxpop+Sys->TotCurPop[1];p0++){
    IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
    for(k=0;k<ParticleBonds[p0].num;k++){
      p1 = ParticleBonds[p0].Ids[k];
      if(deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly)>IntCutoffsqr){
        PrintBonds(Sys);
        printf("\nBondlists erronous\np0 %d p1 %d dist %lf cutoff %lf\n", p0, p1, deltaRsqr(Particles[p0].EndPoint, Particles[p1].EndPoint, Lx, Ly), IntCutoffsqr);
        Err=5;
      }
      
      if(Particles[p0].type != Particles[p1].type){
        PrintBonds(Sys);
        printf("\nBondlists erronous\np0 %d has type %d wheras p1 %d has type %d\n", p0, Particles[p0].type, p1, Particles[p1].type);
        Err=5;
      }
      
      tmp = 0;
      for(l=0;l<ParticleBonds[p1].num;l++){
        if(ParticleBonds[p1].Ids[l]==p0) tmp=1;
      }
      if(tmp==0 ){
        printf("\nBondList reflexivity broken! p0 %d p1 %d\n", p0, p1);
        PrintBondList(ParticleBonds[p0], p0);
        PrintBondList(ParticleBonds[p1], p1);
        Err=5;
      }
    }
  }
  return Err;
}

void ComputeParticleBonds_RanNetSR(sys *Sys){
  int p0, p1;
  double Lx=Sys->Lx, Ly=Sys->Ly;
  double IntCutoffsqr, IntStr, deltaE;
  particle *Particles = Sys->Particles;
  
  for(p0=0;p0<Sys->TotCurPop[0];p0++){
    p1 = Sys->BondList.ParticleBonds[p0].currbond;
    if(p1!=-1){
      IntStr = Particles[p0].MyIntStr;
      IntCutoffsqr = Particles[p0].MyIntCutoffsqr;
      Sys->BondList.Curbondnums[Particles[p0].type]++;
      deltaE = BondInteractions(Particles[p0], Particles[p1], Lx, Ly, IntCutoffsqr, IntStr);
      Sys->TotParticleEnergy+=deltaE;
    }
  }
  ComputeAllParticleBondLists(Sys);
}

void BondMove_RanNetSR(sys *Sys){
  int **CurParPop = Sys->CurParPop, BondableTypeNum = Sys->BondList.BondableTypeNum;
  double **BondProbs = Sys->BondList.BondProbs;
  
  int type, p0, p1, choice;
  double Acc, IntStr;
  particle *Particles=Sys->Particles;
  particlebonds *ParticleBonds = Sys->BondList.ParticleBonds;
  
  int i, toporbot;
  int bondableTop=0, bondableBot=0;
  
  for(i=0;i<BondableTypeNum;i++){
    type = Sys->BondList.BondableTypes[i];
    bondableBot+=CurParPop[0][type];
    bondableTop+=CurParPop[1][type];
    BondProbs[0][i]=CurParPop[0][type];
    BondProbs[1][i]=CurParPop[1][type];
  }
  
  if(bondableBot + bondableTop == 0) return;
  double probBot= ((double)bondableBot)/(bondableTop + bondableBot);
  BondProbs[0][0]/=bondableBot;
  BondProbs[1][0]/=bondableTop;
  for(i=1;i<BondableTypeNum;i++){
    BondProbs[0][i] = BondProbs[0][i]/bondableBot + BondProbs[0][i-1];
    BondProbs[1][i] = BondProbs[1][i]/bondableTop + BondProbs[1][i-1];
  }
  
  for(int counter = 0;counter<1;counter++){
    toporbot = (gsl_rng_uniform(Sys->rng)<probBot ? 0:1);
    Acc = gsl_rng_uniform(Sys->rng);type=-1;
    for(i=0;i<BondableTypeNum;i++){
      if(Acc < BondProbs[toporbot][i]){type = Sys->BondList.BondableTypes[i];break;}
    }
    choice = gsl_rng_uniform_int(Sys->rng, CurParPop[toporbot][type] );
    i=0;
    for(p0=toporbot*Sys->maxpop;p0<toporbot*Sys->maxpop+Sys->TotCurPop[toporbot];p0++){
      if( Particles[p0].type==type ){
        if(i==choice)break;
        else i++;
      }
    }
    
    if(ParticleBonds[p0].num>0){
      choice = gsl_rng_uniform_int(Sys->rng, ParticleBonds[p0].num);
      IntStr = Particles[p0].MyIntStr;
      p1 = ParticleBonds[p0].Ids[choice];
//       printf("bonding info\np0 %d bonded to %d\np1 %d bonded to %d\n", p0, ParticleBonds[p0].currbond, p1, ParticleBonds[p1].currbond);
      //Delete Bond
      if(ParticleBonds[p1].currbond==p0){
        if(gsl_rng_uniform(Sys->rng)< exp(-IntStr)){
//           printf("del complete! %d %d %lf %lf\n", p0, p1, exp(-IntStr), Sys->TotParticleEnergy);
          Sys->TotParticleEnergy += IntStr;
          Sys->TotEnergy         += IntStr;
          ParticleBonds[p1].currbond = -1;
          ParticleBonds[p0].currbond = -1;
          Sys->BondList.Curbondnums[type]--;
        }
      }
      //If partner is unbonded
      else if(ParticleBonds[p1].currbond==-1){
        //Add Bond
        if(ParticleBonds[p0].currbond==-1){
//        if(1==1)
            Sys->TotParticleEnergy +=-IntStr;
            Sys->TotEnergy         +=-IntStr;
            ParticleBonds[p1].currbond = p0;
            ParticleBonds[p0].currbond = p1;
            Sys->BondList.Curbondnums[type]++;
        }
        //Transfer Bond IFF p1 is unsaturated
        else{
//        if(1==1)
            ParticleBonds[p1].currbond = p0;
            ParticleBonds[ParticleBonds[p0].currbond].currbond = -1;
            ParticleBonds[p0].currbond = p1; 
        }
      }
    }
  }
}

double BondParticle_RanNetSR(sys *Sys){
  particle *MovedParticle = &(Sys->dummy_particle);
  if(MovedParticle->MyIntStr>1e-10){
    particle *PartnerParticle;
    particlebonds *ParticleBonds = Sys->BondList.ParticleBonds;
    double cutoffsqr = MovedParticle->MyIntCutoffsqr, Lx=Sys->Lx, Ly=Sys->Ly;
    
    if(ParticleBonds[Sys->p0].currbond!=-1){
      PartnerParticle = &(Sys->Particles[ParticleBonds[Sys->p0].currbond]);
      if(deltaRsqr(PartnerParticle->EndPoint, MovedParticle->EndPoint, Lx, Ly)>cutoffsqr){
        return ALOT;
      }
    }
    PopulateTmpBondList( &(Sys->tmpBondList[0]), MovedParticle, ParticleBonds[Sys->p0].currbond, Sys);
  }
  return 0.0;
}

double BondMeshMove_RanNetSR(sys *Sys){
  double Lx = Sys->Lx, Ly = Sys->Ly;
  int i_affect, current_p;
  double IntRanSqr;
  particle MovedParticle, PartnerParticle;
  
  //This routine is fine, all affected particles (at the top or bottom), only interact with unperturbed particles below/above.
  //So there is no need to check whether the particles are virtual or real.
  for(i_affect=0;i_affect<Sys->AffectedParticleNum;i_affect++){
    current_p=Sys->Unique_dummy_index[i_affect];
    MovedParticle = Sys->mesh_dummy_particles[i_affect];
    if(MovedParticle.MyIntStr>1e-10 && Sys->BondList.ParticleBonds[current_p].currbond!=-1){
      IntRanSqr = MovedParticle.MyIntCutoffsqr;
      PartnerParticle = Sys->Particles[Sys->BondList.ParticleBonds[current_p].currbond];
      if(deltaRsqr(PartnerParticle.EndPoint, MovedParticle.EndPoint, Lx, Ly)>IntRanSqr){return ALOT;}
    }
  }
  
  for(i_affect=0;i_affect<Sys->AffectedParticleNum;i_affect++){
    current_p=Sys->Unique_dummy_index[i_affect];
    MovedParticle = Sys->mesh_dummy_particles[i_affect];
    PopulateTmpBondList( &(Sys->tmpBondList[i_affect]), &MovedParticle, Sys->BondList.ParticleBonds[current_p].currbond, Sys);
  }
  return 0.0;
}

double BondMeshShift_RanNetSR(sys *Sys, double suggestedz){
  double IntRanSqr, Lx=Sys->Lx, Ly=Sys->Ly;
  int pTop, pBot;
  particle CurrentParticle, PartnerParticle;
  for(pBot=0;pBot<Sys->TotCurPop[0];pBot++){
    CurrentParticle = Sys->Particles[pBot];
    CurrentParticle.EndPoint.z +=suggestedz;
    if(Sys->BondList.ParticleBonds[pBot].currbond!=-1){
      pTop = Sys->BondList.ParticleBonds[pBot].currbond;
      PartnerParticle = Sys->Particles[pTop];
      IntRanSqr = CurrentParticle.MyIntCutoffsqr;
      if(deltaRsqr(CurrentParticle.EndPoint, PartnerParticle.EndPoint, Lx, Ly)>IntRanSqr){
        CurrentParticle.EndPoint.z -=suggestedz;
        return ALOT;
      }
    }
    CurrentParticle.EndPoint.z -=suggestedz;
  }
  return 0.0;
}

double BondGCMC_RanNetSR(sys *Sys){
  if(Sys->GC_AddRemove==1 && Sys->BondList.ParticleBonds[Sys->p0].currbond!=-1)
    return -ALOT;
  else return 0.0;
}

void InitializeBondFunctions(sys *Sys){
  int interactions = CreateBondList(Sys);  
  
  if(interactions==0){
    ErrorCheck_Bonds        = ErrorCheck_Bonds_Empty;
    ComputeAllParticleBonds = ComputeParticleBonds_Empty;
    BondParticleMove        = BondParticle_Empty;
    BondMeshMove            = BondMeshMove_Empty;
    BondMeshShift           = BondMeshShift_Empty;
    BondGCMC                = BondGCMC_Empty;
    BondMove                = BondMove_Empty;
    UpdateBondLists         = UpdateBondLists_Empty;
  }
  
  else if(1==2){
    ErrorCheck_Bonds        = ErrorCheck_Bonds_RanNetLR;
    ComputeAllParticleBonds = ComputeParticleBonds_RanNetLR;
    BondParticleMove        = BondParticle_RanNetLR;
    BondMeshMove            = BondMeshMove_RanNetLR;
    BondMeshShift           = BondMeshShift_RanNetLR;
    BondGCMC                = BondGCMC_RanNetLR;
    BondMove                = BondMove_RanNetLR;
    UpdateBondLists         = UpdateBondLists_NotEmpty;
  }
  
  else if(1==1){
    ErrorCheck_Bonds        = ErrorCheck_Bonds_RanNetSR;
    ComputeAllParticleBonds = ComputeParticleBonds_RanNetSR;
    BondParticleMove        = BondParticle_RanNetSR;
    BondMeshMove            = BondMeshMove_RanNetSR;
    BondMeshShift           = BondMeshShift_RanNetSR;
    BondGCMC                = BondGCMC_RanNetSR;
    BondMove                = BondMove_RanNetSR;
    UpdateBondLists         = UpdateBondLists_NotEmpty;
  }
}


//Old Switch Move
//           p0pr = BondList[p0].currbond;
//           p1pr = BondList[p1].currbond;
//           
//           DeltaE =-(p0pr==-1 ? 0.0: BondInteractions(Particles[p0], Particles[p0pr], Lx, Ly, IntRanSqr, IntStr));
//           DeltaE+=-(p1pr==-1 ? 0.0: BondInteractions(Particles[p1], Particles[p1pr], Lx, Ly, IntRanSqr, IntStr));
//           DeltaE+= (p0pr==-1 ? 0.0: BondInteractions(Particles[p1], Particles[p0pr], Lx, Ly, IntRanSqr, IntStr));
//           DeltaE+= (p1pr==-1 ? 0.0: BondInteractions(Particles[p0], Particles[p1pr], Lx, Ly, IntRanSqr, IntStr));
          
//           if(gsl_rng_uniform(Sys->rng)< exp(-DeltaE)){
//             inttype = BondList[p0].currbond;
//             BondList[p0].currbond   = p1pr;
//             BondList[p1].currbond   = p0pr;
//             if(p0pr!=-1) BondList[p0pr].currbond = p1;
//             if(p1pr!=-1) BondList[p1pr].currbond = p0;
//             Sys->TotParticleEnergy+=DeltaE;
//             Sys->TotEnergy+=DeltaE;
//           }



//  else if(1==0){
//     ErrorCheck_Bonds        = ErrorCheck_Bonds_NN;
//     ComputeAllParticleBonds = ComputeParticleBonds_NN;
//     BondParticleMove        = BondParticle_NN;
//     BondMeshMove            = BondMeshMove_NN;
//     BondMeshShift           = BondMeshShift_NN;
//     BondGCMC                = BondGCMC_NN;
//     BondMove                = BondMove_NN;
//     UpdateBondLists         = UpdateBondLists_NotEmpty;
//   }


// int CompareBondLists_NN(particlebonds *tmpBondList, particlebonds *ParticleBonds){
//   int Err=0;
//   
//   if (tmpBondList->num!=ParticleBonds->num){
//     printf("\nError candidate mismatch computed %d stored %d\n", tmpBondList->num, ParticleBonds->num);
//     Err=5;
//   }
//   
//   if (tmpBondList->currbond!=ParticleBonds->currbond){
//     printf("\nError assigned bond mismatch computed %d stored %d\n", tmpBondList->currbond, ParticleBonds->currbond);
//     Err=5;
//   }
//   
//   for(int i=0;i<BONDNUMS;i++){
//     if (tmpBondList->Ids[i]!=ParticleBonds->Ids[i]){
//       printf("\nError local environment computed %d stored %d\n", tmpBondList->Ids[i], ParticleBonds->Ids[i]);
//       Err=5;
//     }
//     if (fabs(tmpBondList->rsqr[i]-ParticleBonds->rsqr[i])>1e-6){
//       printf("\nError mutual distances computed %lf stored %lf\n", tmpBondList->rsqr[i], ParticleBonds->rsqr[i]);
//       Err=5;
//     }
//   }
//   
//   return Err;
// }
// 
// int ErrorCheck_Bonds_NN(sys *Sys){
//   int Err = 0;
//   double Lx = Sys->Lx, Ly = Sys->Ly;
//   
//   particlebonds tmpBondList;
//   particle MyParticle, PartnerParticle;
//   int p, pp, ppbot, pptop, mytype, tmpnum, meshId;
//   double cutoffsqr, mystrength, rsqr, TotParticleEnergy=0.0;
//   
//   for(p = 0;p <Sys->twomaxpop;p++){
//     if(p == Sys->TotCurPop[0]) p=Sys->maxpop;
//     if(p == Sys->TotCurPop[1] + Sys->maxpop) break;
//     MyParticle = Sys->Particles[p];
//     mystrength = MyParticle.MyIntStr;
//     
//     if(mystrength>1e-10){
//       InitBondList(&tmpBondList);
//       mytype    = MyParticle.type;
//       cutoffsqr = MyParticle.MyIntCutoffsqr;
//       meshId    = MyParticle.MeshId;
//       ppbot = Sys->maxpop*(meshId==TOP ? 0 : 1);
//       pptop = Sys->maxpop*(meshId==TOP ? 0 : 1) + Sys->TotCurPop[(meshId==TOP ? 0 : 1)];
//       
//       for(pp = ppbot;pp<pptop;pp++){
//         PartnerParticle = Sys->Particles[pp];
//         if(PartnerParticle.type == mytype){
//           rsqr = deltaRsqr(PartnerParticle.EndPoint, MyParticle.EndPoint, Lx, Ly);
//           if(rsqr<cutoffsqr){
//             tmpnum = tmpBondList.num;
//             tmpBondList.Ids[tmpnum] = pp;
//             tmpBondList.rsqr[tmpnum] = rsqr;
//             tmpBondList.num++;
//             if (tmpBondList.num==BONDNUMS){printf("Too many potential bonds computed %d permitted %d\n", tmpBondList.num, BONDNUMS);Err=5;return Err;}
//           }
//         }
//       }
//     }
//     SortBondList(&(tmpBondList), BONDNUMS);
//     Err = CompareBondLists_NN( &(tmpBondList), &(Sys->BondList.ParticleBonds[p]));
//   }
//   
//   if(fabs(TotParticleEnergy - Sys->TotParticleEnergy)>1e-6) 
//     printf("Particle Energy mismatch\n computed %lf stored %lf\n", TotParticleEnergy, Sys->TotParticleEnergy);
//   return Err;
// }
// 
// int prefPoverCurr(particlebonds Bondlist, int proposed, int existing){
//   for(int j=0;j<BONDNUMS;j++){
//     if(proposed == Bondlist.Ids[j]) return 1;
//     if(existing == Bondlist.Ids[j]) return 0;
//   }
//   return 0;
// }
// 
// void FormBonds_NN(sys *Sys){
//   int p, spot, suggestId;
//   particlebonds *BondList=Sys->BondList.ParticleBonds;
//   
//   for(p=0;p<Sys->TotCurPop[0];p++){ 
//     if(BondList[p].num>0) BondList[p].spotonlist=0;
//   }
//   
//   while(1){
//     //Determine bachelor
//     for(p=0;p<Sys->TotCurPop[0];p++){
//       if(BondList[p].num>0 && BondList[p].spotonlist<BONDNUMS && BondList[p].currbond==-1) break;
//     }
//     //Escape if none apply
//     if(p==Sys->TotCurPop[0]) break;
//     
//     for(spot=BondList[p].spotonlist;spot<BONDNUMS && BondList[p].currbond==-1;spot++){
//       suggestId=BondList[p].Ids[spot];
//       if(BondList[suggestId].currbond==-1){
//         BondList[p].currbond=suggestId;
//         BondList[suggestId].currbond=p;
//       }
//       
//       else{
//         if(prefPoverCurr(BondList[suggestId], p, BondList[suggestId].currbond)== 1){
//           BondList[BondList[suggestId].currbond].currbond = -1;
//           BondList[suggestId].currbond = p;
//           BondList[p].currbond = suggestId;
//         }
//       }
//     }
//     BondList[p].spotonlist++;
//   }
//   
//   for(p=0;p<Sys->TotCurPop[0];p++){
//     if(BondList[p].currbond!=-1){
//       Sys->TotParticleEnergy+=-Sys->Particles[p].MyIntStr;
//       Sys->BondList.Curbondnums[Sys->Particles[p].type]++;
//     }
//   }
// //   PrintBonds(Sys);
// }
// 
// void BondTestDat(sys *Sys){
//   particlebonds *BondList = Sys->BondList.ParticleBonds;
//     
//   BondList[0].num = 1;
//   BondList[0].Ids[0] = Sys->maxpop + 3;
//   BondList[0].Ids[1] = Sys->maxpop + 0;
//   BondList[0].Ids[2] = Sys->maxpop + 2;
//   BondList[0].Ids[3] = Sys->maxpop + 1;
//   
//   BondList[1].num = 1;
//   BondList[1].Ids[0] = Sys->maxpop + 3;
//   BondList[1].Ids[1] = -1;
//   BondList[1].Ids[2] = -1;
//   BondList[1].Ids[3] = -1;
//   
//   BondList[2].num = 1;
//   BondList[2].Ids[0] = Sys->maxpop + 0;
//   BondList[2].Ids[1] = -1;
//   BondList[2].Ids[2] = -1;
//   BondList[2].Ids[3] = -1;
//   
//   BondList[3].num = 1;
//   BondList[3].Ids[0] = Sys->maxpop + 2;
//   BondList[3].Ids[1] = -1;
//   BondList[3].Ids[2] = -1;
//   BondList[3].Ids[3] = -1;
//   
//   BondList[Sys->maxpop + 0].num = 1;
//   BondList[Sys->maxpop + 0].Ids[0] = 2;
//   BondList[Sys->maxpop + 0].Ids[1] = -1;
//   BondList[Sys->maxpop + 0].Ids[2] = -1;
//   BondList[Sys->maxpop + 0].Ids[3] = -1;
//   
//   BondList[Sys->maxpop + 1].num = 1;
//   BondList[Sys->maxpop + 1].Ids[0] = 0;
//   BondList[Sys->maxpop + 1].Ids[1] = -1;
//   BondList[Sys->maxpop + 1].Ids[2] = -1;
//   BondList[Sys->maxpop + 1].Ids[3] = -1;
//   
//   BondList[Sys->maxpop + 2].num = 1;
//   BondList[Sys->maxpop + 2].Ids[0] = 3;
//   BondList[Sys->maxpop + 2].Ids[1] = -1;
//   BondList[Sys->maxpop + 2].Ids[2] = -1;
//   BondList[Sys->maxpop + 2].Ids[3] = -1;
//   
//   BondList[Sys->maxpop + 3].num = 1;
//   BondList[Sys->maxpop + 3].Ids[0] = 1;
//   BondList[Sys->maxpop + 3].Ids[1] = -1;
//   BondList[Sys->maxpop + 3].Ids[2] = -1;
//   BondList[Sys->maxpop + 3].Ids[3] = -1;
//   
// }
// 
// void ComputeParticleBonds_NN(sys *Sys){
//   double Lx = Sys->Lx, Ly = Sys->Ly;
//   particlebonds *BondList = Sys->BondList.ParticleBonds;
//   
//   particle MyParticle, PartnerParticle;
//   int p, pp, ppbot, pptop, mytype, tmpnum, meshId;
//   double cutoffsqr, mystrength, rsqr;
//   
//   for(p = 0;p <Sys->twomaxpop;p++){
//     if(p == Sys->TotCurPop[0]) p=Sys->maxpop;
//     if(p == Sys->TotCurPop[1] + Sys->maxpop) break;
//     MyParticle = Sys->Particles[p];
//     mystrength = MyParticle.MyIntStr;
//     
//     if(mystrength>1e-10){
//       InitBondList(&BondList[p]);
//       mytype    = MyParticle.type;
//       cutoffsqr = MyParticle.MyIntCutoffsqr;
//       meshId    = MyParticle.MeshId;
//       ppbot = Sys->maxpop*(meshId==TOP ? 0 : 1);
//       pptop = Sys->maxpop*(meshId==TOP ? 0 : 1) + Sys->TotCurPop[(meshId==TOP ? 0 : 1)];
//       
//       for(pp = ppbot;pp<pptop;pp++){
//         PartnerParticle = Sys->Particles[pp];
//         if(PartnerParticle.type == mytype){
//           rsqr = deltaRsqr(PartnerParticle.EndPoint, MyParticle.EndPoint, Lx, Ly);
//           if(rsqr<cutoffsqr){
//             tmpnum = BondList[p].num;
//             BondList[p].Ids[tmpnum] = pp;
//             BondList[p].rsqr[tmpnum] = rsqr;
//             BondList[p].num++;
//           }
//         }
//       }
//     }
//     SortBondList(&(BondList[p]), BONDNUMS);
//   }
//   
//   for(p=0;p<Sys->TotCurPop[0];p++){PrintBondList(BondList[p], p);}
//   
// //   BondTestDat(Sys);
//   FormBonds_NN(Sys);
// }
// 
// double BondParticle_NN(sys *Sys){
//   double deltaE=0.0;
//   particle *MovedParticle   = &(Sys->dummy_particle);
//   double MyIntStr = MovedParticle->MyIntStr;
//   particle *PartnerParticle;
//   particlebonds *ParticleBonds = Sys->BondList.ParticleBonds; 
//   
//   if(MyIntStr>1e-10){
//     particlebonds *tmpBondList = &(Sys->tmpBondList[0]);
//     InitBondList(tmpBondList);
//     
//     celllist *CellList = &(Sys->CellList);
//     int MyCx, MyCy, Cx, Cy, Populations, *Ids, k;
//     int mytype = MovedParticle->type, tmpnum;
//     double cutoffsqr = MovedParticle->MyIntCutoffsqr, Lx=Sys->Lx, Ly=Sys->Ly;
//     double rsqr;
//     
//     Position_To_Index_PBC(&MyCx, MovedParticle->Pos.x, CellList->dCx, CellList->NumCellsx);
//     Position_To_Index_PBC(&MyCy, MovedParticle->Pos.y, CellList->dCy, CellList->NumCellsy);
//     
//     for(int i=0;i<9;i++){
//       Cx = CellList->Neighbors[MyCx][MyCy][2*i];
//       Cy = CellList->Neighbors[MyCx][MyCy][2*i + 1];
//       Populations = (Sys->dummy_particle.MeshId == TOP ? CellList->BotPopulations[Cx][Cy] : CellList->TopPopulations[Cx][Cy]);
//       Ids         = (Sys->dummy_particle.MeshId == TOP ? CellList->BotParticleIds[Cx][Cy] : CellList->TopParticleIds[Cx][Cy]);
//       
//       for(k=0;k<Populations;k++){
//         PartnerParticle = &(Sys->Particles[Ids[k]]);
//         if(PartnerParticle->type == mytype){
//           rsqr = deltaRsqr(PartnerParticle->EndPoint, MovedParticle->EndPoint, Lx, Ly);
//           if(rsqr<cutoffsqr){
//             tmpnum = tmpBondList->num;
//             tmpBondList->Ids[tmpnum] = Ids[k];
//             tmpBondList->rsqr[tmpnum] = rsqr;
//             tmpBondList->num++;
//           }
//         }
//       }
//     }
//     deltaE = (tmpBondList->num > 0 ? -MyIntStr:0.0) - (ParticleBonds[Sys->p0].currbond!=-1 ? -MyIntStr:0.0);
//   }
//   return deltaE;
// }
// 
// double BondMeshMove_NN(sys *Sys){return 0.0;}
// 
// double BondMeshShift_NN(sys *Sys, double suggestedz){return 0.0;}
// 
// double BondGCMC_NN(sys *Sys){return 0.0;}
// 
// void BondMove_NN(sys *Sys){}
// 
