#include "gcmc.h"

void ErrorCheckGCMC(sys *Sys){
  
 if(Sys->TotCurPop[0] == Sys->maxpop || Sys->TotCurPop[1] == Sys->maxpop ){
  printf("ERROR: You are exceeding max allowable popnums:\nbot %d top %d max %d sum %d maxsum %d",
         Sys->TotCurPop[0], Sys->TotCurPop[1], Sys->maxpop,
         Sys->TotCurPop[0]+Sys->TotCurPop[1],  Sys->twomaxpop
        );
    exit(1);
  }
}

double SuggestParticleInsertion(sys *Sys){
  int PType    = gsl_rng_uniform_int(Sys->rng, Sys->NumPTypes);
  int TopOrBot = 2*gsl_rng_uniform_int(Sys->rng, 2) - 1;
  
  Sys->GC_TopOrBot = TopOrBot;
  Sys->GC_PType = PType;
    
  Sys->dummy_particle.Pos.x = gsl_rng_uniform(Sys->rng)*Sys->Lx;
  Sys->dummy_particle.Pos.y = gsl_rng_uniform(Sys->rng)*Sys->Ly;
  
  Sys->dummy_particle.type     = PType;
  Sys->dummy_particle.MeshId   = TopOrBot;
  Sys->dummy_particle.MyDiam   = Sys->DPType[PType];
  Sys->dummy_particle.MyRad    = Sys->dummy_particle.MyDiam/2;
  Sys->dummy_particle.MyRadSqr = Sys->dummy_particle.MyRad*Sys->dummy_particle.MyRad;
  Sys->dummy_particle.MyLength = Sys->LPType[PType];
  
  Sys->dummy_particle.MyIntStr              = Sys->BindPType[PType];
  Sys->dummy_particle.MyIntCutoff           = Sys->RanPType[PType];
  Sys->dummy_particle.MyIntCutoffsqr        = Sys->dummy_particle.MyIntCutoff*Sys->dummy_particle.MyIntCutoff;
  Sys->dummy_particle.MeshParticleCutoffsqr = (Sys->dummy_particle.MyRad + Sys->dx/2)*(Sys->dummy_particle.MyRad + Sys->dx/2);
  
  Sys->p0 = Sys->TotCurPop[(TopOrBot == BOT) ? 0 : 1] + ((TopOrBot == BOT) ? 0 : Sys->maxpop);
  
  Compute_Z_nhat_PhaseFactor(Sys, &(Sys->dummy_particle) );
  return Sys->dummy_particle.PhaseFactorEnergy;
}

double SuggestParticleRemoval(sys *Sys){
  int PType    = gsl_rng_uniform_int(Sys->rng, Sys->NumPTypes);
  int TopOrBot = 2*gsl_rng_uniform_int(Sys->rng, 2) - 1;
  int MeshIndex = ((TopOrBot == BOT) ? 0 : 1);
  
  if (Sys->CurParPop[MeshIndex][PType] == 0 ) return ALOT;
  
  else{
    int particleflavornum = gsl_rng_uniform_int(Sys->rng, Sys->CurParPop[MeshIndex][PType] );
    int imin = Sys->maxpop*MeshIndex, imax = imin + Sys->TotCurPop[MeshIndex];
    int p0 = -1, i;
    
    for(i=imin; i<imax; i++){
      if(Sys->Particles[i].type == PType) p0++;
      if(p0 == particleflavornum){p0 = i;break;}
    }
    
    Sys->p0 = p0;
    Sys->GC_TopOrBot = TopOrBot;
    Sys->GC_PType = PType;
    return -Sys->Particles[p0].PhaseFactorEnergy;
  }
}

double GCMCAdditionProb(sys *Sys, double DeltaETot){
 int MeshIndex = (Sys->GC_TopOrBot == BOT) ? 0 : 1;
 int TypeIndex = Sys->GC_PType;
 return Sys->AccessVol[TypeIndex]/(Sys->CurParPop[MeshIndex][TypeIndex] + 1)*exp(Sys->ChemPots[MeshIndex][TypeIndex] - DeltaETot);
}

double GCMCRemovalProb(sys *Sys, double DeltaETot){
 int MeshIndex = (Sys->GC_TopOrBot == BOT) ? 0 : 1;
 int TypeIndex = Sys->GC_PType; 
 return (Sys->CurParPop[MeshIndex][TypeIndex]/(float)Sys->AccessVol[TypeIndex])*exp(-Sys->ChemPots[MeshIndex][TypeIndex] - DeltaETot);
}

void CreateNewParticle(sys *Sys, double DeltaEPhase, double DeltaEPart, double DeltaMPEnergy){
  int p0=Sys->p0;
  
  CopyVec(&(Sys->Particles[p0].Pos)      , Sys->dummy_particle.Pos );
  CopyVec(&(Sys->Particles[p0].nhat)     , Sys->dummy_particle.nhat);
  CopyVec(&(Sys->Particles[p0].MidPoint) , Sys->dummy_particle.MidPoint );
  CopyVec(&(Sys->Particles[p0].EndPoint) , Sys->dummy_particle.EndPoint );

  Sys->Particles[p0].MeshParticleEnergy = Sys->dummy_particle.MeshParticleEnergy;
  Sys->Particles[p0].PhaseFactorEnergy  = Sys->dummy_particle.PhaseFactorEnergy;
    
  Sys->Particles[p0].type     = Sys->dummy_particle.type;
  Sys->Particles[p0].MeshId   = Sys->dummy_particle.MeshId;
  Sys->Particles[p0].MyDiam   = Sys->dummy_particle.MyDiam;
  Sys->Particles[p0].MyRad    = Sys->dummy_particle.MyRad;
  Sys->Particles[p0].MyRadSqr = Sys->dummy_particle.MyRadSqr;
  Sys->Particles[p0].MyLength = Sys->dummy_particle.MyLength;
  
  Sys->Particles[p0].MyIntStr              = Sys->dummy_particle.MyIntStr;
  Sys->Particles[p0].MyIntCutoff           = Sys->dummy_particle.MyIntCutoff;
  Sys->Particles[p0].MyIntCutoffsqr        = Sys->dummy_particle.MyIntCutoffsqr;
  Sys->Particles[p0].MeshParticleCutoffsqr = Sys->dummy_particle.MeshParticleCutoffsqr;

  Sys->TotCurPop[(Sys->GC_TopOrBot == BOT) ? 0 : 1]++;
  Sys->CurParPop[(Sys->GC_TopOrBot == BOT) ? 0 : 1][Sys->GC_PType]++;
  
  int MyCx, MyCy;
  Position_To_Index_PBC(&MyCx, Sys->Particles[p0].Pos.x, Sys->CellList.dCx, Sys->CellList.NumCellsx);
  Position_To_Index_PBC(&MyCy, Sys->Particles[p0].Pos.y, Sys->CellList.dCy, Sys->CellList.NumCellsy);
  int  **Populations = (Sys->GC_TopOrBot == BOT) ? Sys->CellList.BotPopulations : Sys->CellList.TopPopulations;
  int ***ParticleIds = (Sys->GC_TopOrBot == BOT) ? Sys->CellList.BotParticleIds : Sys->CellList.TopParticleIds;
  ParticleIds[MyCx][MyCy][Populations[MyCx][MyCy]] = p0;
  Populations[MyCx][MyCy]++;
  
  Sys->TotPhaseFactorEnergy  += DeltaEPhase;
  Sys->TotParticleEnergy     += DeltaEPart;
  Sys->TotMeshParticleEnergy += DeltaMPEnergy;
}

void ReplaceParticle(sys *Sys, double DeltaEPhase, double DeltaEPart){
  int meshindex    = (Sys->GC_TopOrBot == BOT) ? 0 : 1;
  int p0 = Sys->p0;
  
  Sys->TotCurPop[meshindex]--;
  Sys->CurParPop[meshindex][Sys->GC_PType]--;
  Sys->TotPhaseFactorEnergy  += DeltaEPhase;
  Sys->TotParticleEnergy     -= DeltaEPart;
  Sys->TotMeshParticleEnergy -= Sys->Particles[p0].MeshParticleEnergy;
  
  int MyCx, MyCy;
  Position_To_Index_PBC(&MyCx, Sys->Particles[p0].Pos.x, Sys->CellList.dCx, Sys->CellList.NumCellsx);
  Position_To_Index_PBC(&MyCy, Sys->Particles[p0].Pos.y, Sys->CellList.dCy, Sys->CellList.NumCellsy);
  int  **Populations = (Sys->GC_TopOrBot == BOT) ? Sys->CellList.BotPopulations : Sys->CellList.TopPopulations;
  int ***ParticleIds = (Sys->GC_TopOrBot == BOT) ? Sys->CellList.BotParticleIds : Sys->CellList.TopParticleIds;
  
  Populations[MyCx][MyCy]--;
  for(int k=0;k<Populations[MyCx][MyCy] + 1;k++){
    if(ParticleIds[MyCx][MyCy][k]==p0){
      ParticleIds[MyCx][MyCy][k] = ParticleIds[MyCx][MyCy][Populations[MyCx][MyCy]];
      break;
    }
  }
  ParticleIds[MyCx][MyCy][Populations[MyCx][MyCy]] = EMPTY;
  
  int replaceindex = Sys->maxpop*meshindex + Sys->TotCurPop[meshindex];
  
  if(replaceindex != p0){
    Position_To_Index_PBC(&MyCx, Sys->Particles[replaceindex].Pos.x, Sys->CellList.dCx, Sys->CellList.NumCellsx);
    Position_To_Index_PBC(&MyCy, Sys->Particles[replaceindex].Pos.y, Sys->CellList.dCy, Sys->CellList.NumCellsy);
    for(int k=0;k<Populations[MyCx][MyCy];k++){
      if(ParticleIds[MyCx][MyCy][k]==replaceindex){
        ParticleIds[MyCx][MyCy][k] = p0;
        break;
      }
    }
  
    CopyVec(&(Sys->Particles[p0].Pos)      , Sys->Particles[replaceindex].Pos );
    CopyVec(&(Sys->Particles[p0].nhat)     , Sys->Particles[replaceindex].nhat);
    CopyVec(&(Sys->Particles[p0].MidPoint) , Sys->Particles[replaceindex].MidPoint );
    CopyVec(&(Sys->Particles[p0].EndPoint) , Sys->Particles[replaceindex].EndPoint );

    Sys->Particles[p0].MeshParticleEnergy = Sys->Particles[replaceindex].MeshParticleEnergy; 
    Sys->Particles[p0].PhaseFactorEnergy  = Sys->Particles[replaceindex].PhaseFactorEnergy; 
    
    Sys->Particles[p0].type     = Sys->Particles[replaceindex].type;
    Sys->Particles[p0].MeshId   = Sys->Particles[replaceindex].MeshId;
    Sys->Particles[p0].MyDiam   = Sys->Particles[replaceindex].MyDiam;
    Sys->Particles[p0].MyRad    = Sys->Particles[replaceindex].MyRad;
    Sys->Particles[p0].MyRadSqr = Sys->Particles[replaceindex].MyRadSqr;
    Sys->Particles[p0].MyLength = Sys->Particles[replaceindex].MyLength;
    
    Sys->Particles[p0].MyIntStr              = Sys->Particles[replaceindex].MyIntStr;
    Sys->Particles[p0].MyIntCutoff           = Sys->Particles[replaceindex].MyIntCutoff;
    Sys->Particles[p0].MyIntCutoffsqr        = Sys->Particles[replaceindex].MyIntCutoffsqr;
    Sys->Particles[p0].MeshParticleCutoffsqr = Sys->Particles[replaceindex].MeshParticleCutoffsqr;
  }
  
  Sys->Particles[replaceindex].Pos.x = EMPTY;
  Sys->Particles[replaceindex].Pos.y = EMPTY;
  Sys->Particles[replaceindex].Pos.z = EMPTY;
  
  Sys->Particles[replaceindex].MidPoint.x = EMPTY;
  Sys->Particles[replaceindex].MidPoint.y = EMPTY;
  Sys->Particles[replaceindex].MidPoint.z = EMPTY;
  
  Sys->Particles[replaceindex].EndPoint.x = EMPTY;
  Sys->Particles[replaceindex].EndPoint.y = EMPTY;
  Sys->Particles[replaceindex].EndPoint.z = EMPTY;
  Sys->Particles[replaceindex].PhaseFactorEnergy = 0.0; 
  Sys->Particles[replaceindex].MeshParticleCutoffsqr = 0.0;
}

double PerformGCMCMove(sys *Sys){
  int acc=0;
  Sys->GC_AddRemove = gsl_rng_uniform_int(Sys->rng, 2);
  
  if (Sys->GC_AddRemove == 0){
    double DeltaEPhase   = SuggestParticleInsertion(Sys);
    double DeltaMPEnergy = ParticleMeshEnergy(Sys, &(Sys->dummy_particle));
    Sys->dummy_particle.MeshParticleEnergy = DeltaMPEnergy;
    if(DeltaMPEnergy<ALOT){
      double DeltaEPart = CollisionTestAndParticleBonds(Sys);
//       double DeltaEPart = 0.0;
      if(DeltaEPart<ALOT){
        double DeltaETot = DeltaEPhase + DeltaEPart + DeltaMPEnergy;
        if ( gsl_rng_uniform(Sys->rng)< GCMCAdditionProb(Sys, DeltaETot) ){
//           ErrorCheckGCMC(Sys);
          CreateNewParticle(Sys, DeltaEPhase, DeltaEPart, DeltaMPEnergy);
          Sys->TotEnergy+=DeltaETot;
          acc=1;
        }
      }
    }
  }
  
  else{
    double DeltaEPhase = SuggestParticleRemoval(Sys);
    if(DeltaEPhase < ALOT){
      double DeltaEPart = ParticleEnergy(Sys, Sys->Particles[Sys->p0]);
//       double DeltaEPart = 0.0;
      double DeltaETot = -(DeltaEPhase + DeltaEPart + Sys->Particles[Sys->p0].MeshParticleEnergy);
      if ( gsl_rng_uniform(Sys->rng)< GCMCRemovalProb(Sys, DeltaETot) ){
//         ErrorCheckGCMC(Sys);
        ReplaceParticle(Sys, DeltaEPhase, DeltaEPart);
        Sys->TotEnergy-=DeltaETot;
        acc=1;
      }
    }
  }
  
//   ErrorCheckCellLists(Sys, 1);
  return acc;
}
  
