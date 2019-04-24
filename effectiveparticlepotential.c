#include "effectiveparticlepotential.h"


int ErrorCheck_Veff(sys *Sys){
  int Err=0;
  if (fabs(Sys->TotVeff - ReturnTotalParticlePotential(Sys))>1e-6){
    printf("\nVeffError! Stored: %lf Computed: %lf\n", Sys->TotVeff, ReturnTotalParticlePotential(Sys));
    Err=3;
  }
  return Err;
}

double ComputeVeffZero(sys *Sys){
  int p1, p2;
  double expveff=1.0;
  for(p1=0;p1<Sys->NumPTypes;p1++){
    expveff+=Sys->Fugacities[0][p1] + Sys->Fugacities[1][p1];
    for(p2=0;p2<Sys->NumPTypes;p2++){
      expveff+=Sys->Fugacities[0][p1]*Sys->Fugacities[1][p2];
    }
  }
  return -log(expveff); 
}

double ReturnParticlePotentialLatticeSite(sys *Sys, double dz, double sqrtnormbot, double sqrtnormtop){
  double **Fugacities = Sys->Fugacities;
  double  *BindPType  = Sys->BindPType;
  double **Lengths    = Sys->LPType;
  
  //Enumerate "spin" states
  double expveff=1.0; //no spins
  double interaction, sumlengths;
  int p1, p2;
  
  for(p1=0;p1<Sys->NumPTypes;p1++){
    if(dz > Lengths[0][p1]) expveff+=Fugacities[0][p1]*sqrtnormbot; //1 spin bot
    if(dz > Lengths[1][p1]) expveff+=Fugacities[1][p1]*sqrtnormtop; //1 spin top
    
    //Two spins
    for(p2=0;p2<Sys->NumPTypes;p2++){
      interaction=0.0;
      sumlengths = Lengths[0][p1] + Lengths[1][p2];
      
      if(p1==p2 && BindPType[p1]>0.000001){
        if(dz > sumlengths){        
          interaction = Fugacities[0][p1]*Fugacities[1][p2]*sqrtnormbot*sqrtnormtop;
          if( fabs(dz - sumlengths)<Sys->RanPType[p1] ){
            interaction *= exp(Sys->BindPType[p1]);
          }
        }
      }
      
      else{
        if(dz > sumlengths){
          interaction = Fugacities[0][p1]*Fugacities[1][p2]*sqrtnormbot*sqrtnormtop;
        }
      }
      expveff+=interaction;
    }
  }
  
  return -log(expveff)-Sys->VeffZero;
}

double ReturnTotalParticlePotential(sys *Sys){
  double totVeff=0.0;
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      totVeff+= ReturnParticlePotentialLatticeSite(
                  Sys, Sys->TopMesh[i][j].Pos.z-Sys->BotMesh[i][j].Pos.z, 
                  Sys->BotMesh[i][j].nNorm, Sys->TopMesh[i][j].nNorm
                );
    }
  }
  return totVeff;
}

void SetVeffLatticeSite(sys *Sys, int i0, int j0){
  double dz=Sys->TopMesh[i0][j0].Pos.z-Sys->BotMesh[i0][j0].Pos.z;
  double sqrtnormbot = Sys->BotMesh[i0][j0].nNorm, sqrtnormtop = Sys->TopMesh[i0][j0].nNorm;
  Sys->TopMesh[i0][j0].veff = Sys->BotMesh[i0][j0].veff = 0.5*ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop); 
}

void InitializeVeff(sys *Sys){
  Sys->TotVeff=0.0;Sys->VeffZero = ComputeVeffZero(Sys);
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      SetVeffLatticeSite(Sys, i, j);
      Sys->TotVeff+=Sys->TopMesh[i][j].veff + Sys->BotMesh[i][j].veff;
    }
  }
}

double VeffChange(sys *Sys){
  int i0=Sys->i0, j0=Sys->j0;
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  
  double veffold=0.0, dz, sqrtnormbot, sqrtnormtop;
  
  dz = (Sys->choice == TOP ? Sys->suggestedz - Sys->BotMesh[i0][j0].Pos.z : Sys->TopMesh[i0][j0].Pos.z - Sys->suggestedz);
  sqrtnormbot = Sys->BotMesh[i0][j0].nNorm;sqrtnormtop = Sys->TopMesh[i0][j0].nNorm;
  Sys->veffC=ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop);
  veffold+=(Sys->TopMesh[i0][j0].veff + Sys->BotMesh[i0][j0].veff);
  
  dz = Sys->TopMesh[left][j0].Pos.z - Sys->BotMesh[left][j0].Pos.z;
  sqrtnormbot = (Sys->choice == TOP ? Sys->BotMesh[left][j0].nNorm : Sys->nNormLeft);
  sqrtnormtop = (Sys->choice == TOP ? Sys->nNormLeft : Sys->TopMesh[left][j0].nNorm);
  Sys->veffLeft=ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop);
  veffold+=(Sys->TopMesh[left][j0].veff + Sys->BotMesh[left][j0].veff);
  
  dz = Sys->TopMesh[right][j0].Pos.z - Sys->BotMesh[right][j0].Pos.z;
  sqrtnormbot = (Sys->choice == TOP ? Sys->BotMesh[right][j0].nNorm : Sys->nNormRight);
  sqrtnormtop = (Sys->choice == TOP ? Sys->nNormRight : Sys->TopMesh[right][j0].nNorm);
  Sys->veffRight=ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop);
  veffold+=(Sys->TopMesh[right][j0].veff + Sys->BotMesh[right][j0].veff);
  
  dz = Sys->TopMesh[i0][in].Pos.z - Sys->BotMesh[i0][in].Pos.z;
  sqrtnormbot = (Sys->choice == TOP ? Sys->BotMesh[i0][in].nNorm : Sys->nNormIn);
  sqrtnormtop = (Sys->choice == TOP ? Sys->nNormIn : Sys->TopMesh[i0][in].nNorm);
  Sys->veffIn=ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop);
  veffold+=(Sys->TopMesh[i0][in].veff + Sys->BotMesh[i0][in].veff);
  
  dz = Sys->TopMesh[i0][out].Pos.z - Sys->BotMesh[i0][out].Pos.z;
  sqrtnormbot = (Sys->choice == TOP ? Sys->BotMesh[i0][out].nNorm : Sys->nNormOut);
  sqrtnormtop = (Sys->choice == TOP ? Sys->nNormOut : Sys->TopMesh[i0][out].nNorm);
  Sys->veffOut=ReturnParticlePotentialLatticeSite(Sys, dz, sqrtnormbot, sqrtnormtop);
  veffold+=(Sys->TopMesh[i0][out].veff + Sys->BotMesh[i0][out].veff);
  
  return Sys->veffLeft + Sys->veffRight + Sys->veffIn + Sys->veffOut + Sys->veffC - veffold;
}

void StoreVeffChange(sys *Sys, double DeltaEVeff){
  int i0=Sys->i0, j0=Sys->j0;
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0); 
  Sys->BotMesh[i0][j0].veff    =0.5*Sys->veffC;
  Sys->TopMesh[i0][j0].veff    =0.5*Sys->veffC;
  Sys->BotMesh[left][j0].veff  =0.5*Sys->veffLeft;
  Sys->TopMesh[left][j0].veff  =0.5*Sys->veffLeft;
  Sys->BotMesh[right][j0].veff =0.5*Sys->veffRight;
  Sys->TopMesh[right][j0].veff =0.5*Sys->veffRight;
  Sys->BotMesh[i0][in].veff    =0.5*Sys->veffIn;
  Sys->TopMesh[i0][in].veff    =0.5*Sys->veffIn;
  Sys->BotMesh[i0][out].veff   =0.5*Sys->veffOut;
  Sys->TopMesh[i0][out].veff   =0.5*Sys->veffOut;
  Sys->TotVeff+=DeltaEVeff;
}

void ComputeProbsOnSite(sys *Sys, int i0, int j0){
  double dz = Sys->TopMesh[i0][j0].Pos.z - Sys->BotMesh[i0][j0].Pos.z;
  double sqrtnormbot = Sys->BotMesh[i0][j0].nNorm;
  double sqrtnormtop = Sys->TopMesh[i0][j0].nNorm;
  double expveff = exp((Sys->VeffZero + Sys->TopMesh[i0][j0].veff + Sys->BotMesh[i0][j0].veff));
  
  double **Fugacities = Sys->Fugacities,  *BindPType = Sys->BindPType;
  double **Lengths    = Sys->LPType;
  
  int p1, p2, NumPTypes = Sys->NumPTypes;
  double sumlengths, ztheta;
  
  Sys->BotMesh[i0][j0].ParticleProbs[NumPTypes][NumPTypes] = 
  Sys->TopMesh[i0][j0].ParticleProbs[NumPTypes][NumPTypes] = expveff;
  
  for(p1=0;p1<NumPTypes;p1++){
    ztheta=0.0;
    if(dz > Lengths[0][p1]) ztheta = sqrtnormbot*Fugacities[0][p1]*expveff;
    Sys->BotMesh[i0][j0].ParticleProbs[p1][NumPTypes] = 
    Sys->TopMesh[i0][j0].ParticleProbs[p1][NumPTypes] = ztheta;
  
    ztheta=0.0;
    if(dz > Lengths[1][p1]) ztheta = sqrtnormtop*Fugacities[1][p1]*expveff;
    Sys->BotMesh[i0][j0].ParticleProbs[NumPTypes][p1] = 
    Sys->TopMesh[i0][j0].ParticleProbs[NumPTypes][p1] = ztheta;
    
    for(p2=0;p2<NumPTypes;p2++){
      ztheta=0.0;sumlengths=Lengths[0][p1] + Lengths[1][p2];
      
      if(p1==p2 && BindPType[p1]>0.000001){
        if(dz > sumlengths){        
          ztheta=sqrtnormtop*sqrtnormbot*Fugacities[0][p1]*Fugacities[1][p2]*expveff;
          if(fabs(dz-sumlengths)<Sys->RanPType[p1]){
            ztheta*=exp(Sys->BindPType[p1]);
            Sys->BotMesh[i0][j0].BondNums[p1] = Sys->TopMesh[i0][j0].BondNums[p1] = ztheta;
          }
          else Sys->BotMesh[i0][j0].BondNums[p1] = Sys->TopMesh[i0][j0].BondNums[p1] = 0.0;
        }
      }
      
      else if (dz> sumlengths) ztheta=sqrtnormbot*sqrtnormtop*Fugacities[0][p1]*Fugacities[1][p2]*expveff;
      
      
      Sys->BotMesh[i0][j0].ParticleProbs[p1][p2] = 
      Sys->TopMesh[i0][j0].ParticleProbs[p1][p2] = ztheta;
    }
  }
}

void ComputeParticlePops(sys *Sys){
  int i, j, p1, p2, Nx = Sys->Nx, Ny = Sys->Ny, statenums = Sys->CorrFuncs.statenums, compnums = Sys->CorrFuncs.compnums;
  double **CurParPop = Sys->CurParPop, *MeanProbs = Sys->CorrFuncs.Cur_MeanProbs, *BondNums = Sys->BondNums;
         
  for(p1=0;p1<Sys->NumPTypes;p1++) BondNums[p1] = CurParPop[0][p1] = CurParPop[1][p1] = 0.0;
  CurParPop[0][p1] = CurParPop[1][p1] = 0.0;
  for(p1=0;p1<statenums;p1++) MeanProbs[p1] = 0.0;
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      ComputeProbsOnSite(Sys, i, j);
      for(p1=0;p1<compnums;p1++){
        if(p1<Sys->NumPTypes) BondNums[p1]+= Sys->BotMesh[i][j].BondNums[p1];
        for(p2=0;p2<compnums;p2++){
          CurParPop[0][p1]  += Sys->BotMesh[i][j].ParticleProbs[p1][p2];
          CurParPop[1][p1]  += Sys->TopMesh[i][j].ParticleProbs[p2][p1];
          MeanProbs[p1*compnums + p2] += Sys->TopMesh[i][j].ParticleProbs[p1][p2];
        }
      }
    }
  }
  
  Sys->TotCurPop[0] = Sys->TotCurPop[1] = 0.0;
  for(p1=0;p1<Sys->NumPTypes;p1++){
    Sys->TotCurPop[0]   += Sys->CurParPop[0][p1];
    Sys->TotCurPop[1]   += Sys->CurParPop[1][p1];
  }
}

double VeffGlobalShiftChange(sys *Sys){
  double totVeff=0.0, dummy, znew, dz=(Sys->choice == BOT ? -1.0:1.0)*Sys->suggestedz;
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      znew = Sys->TopMesh[i][j].Pos.z-Sys->BotMesh[i][j].Pos.z + dz;
      if(znew<0.0) return ALOT;
      dummy= ReturnParticlePotentialLatticeSite(
                Sys, znew, Sys->BotMesh[i][j].nNorm, Sys->TopMesh[i][j].nNorm
             );
      Sys->TopMesh[i][j].veffdummy=0.5*dummy;
      Sys->BotMesh[i][j].veffdummy=0.5*dummy;
      totVeff+=dummy;
    }
  }
  return totVeff - Sys->TotVeff;
}

void StoreVeffGlobalShiftChange(sys *Sys, double DeltaEVeff){
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      Sys->TopMesh[i][j].veff = Sys->TopMesh[i][j].veffdummy;
      Sys->BotMesh[i][j].veff = Sys->BotMesh[i][j].veffdummy;
    }
  }
  Sys->TotVeff+=DeltaEVeff;
}









// void OldComputeParticlePops(sys *Sys){
//   int p1, p2, i, j;
//   double expveff, ntop, nbot, inttop, intbot, bond, sqrtnormbot, sqrtnormtop, dz;
//   double sumlengths, sumrads;
//   
//   double **CurParPop  = Sys->CurParPop, *BondNums   = Sys->BondNums;
//   double **Fugacities = Sys->Fugacities, *BindPType = Sys->BindPType;
//   double **Lengths    = Sys->LPType, **Diams        = Sys->DPType;
//   
//   for(p1=0;p1<Sys->NumPTypes;p1++){BondNums[p1] = CurParPop[0][p1] = CurParPop[1][p1] = 0.0;}
//   
//   for(i=0;i<Sys->Nx;i++){
//     for(j=0;j<Sys->Ny;j++){
//       dz = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
//       sqrtnormbot = Sys->BotMesh[i][j].nNorm;
//       sqrtnormtop = Sys->TopMesh[i][j].nNorm;
//       expveff = exp((Sys->VeffZero + Sys->TopMesh[i][j].veff + Sys->BotMesh[i][j].veff));
//       
//       for(p1=0;p1<Sys->NumPTypes;p1++){
//         ntop = nbot = bond = 0.0;
//         if(dz>(Lengths[0][p1] + Diams[0][p1]/2))nbot+=1.0;
//         if(dz>(Lengths[1][p1] + Diams[1][p1]/2))ntop+=1.0;
//         
//         for(p2=0;p2<Sys->NumPTypes;p2++){
//           sumlengths=Lengths[0][p1] + Lengths[1][p2];
//           
//           if(p1==p2 && BindPType[p1]>0.000001){
//             if(dz > sumlengths){        
//               intbot=sqrtnormtop*Fugacities[1][p2];
//               inttop=sqrtnormbot*Fugacities[0][p2];
//               
//               if(fabs(dz-sumlengths)<Sys->RanPType[p1]){
//                 intbot*=exp(Sys->BindPType[p1]);
//                 inttop*=exp(Sys->BindPType[p1]);
//                 bond+=sqrtnormbot*sqrtnormtop*Fugacities[0][p1]*Fugacities[1][p2]*exp(Sys->BindPType[p1]);
//               }
//               nbot+=intbot;
//               ntop+=inttop;
//             }
//           }
//           
//           else{
//             sumrads=(Diams[0][p1] + Diams[1][p2])/2.0;
//             if(dz>(sumlengths+sumrads)){
//               nbot+=sqrtnormtop*Fugacities[1][p2];
//               ntop+=sqrtnormbot*Fugacities[0][p2];
//             }
//           }
//         }
//         Sys->CurParPop[0][p1]+=nbot*sqrtnormbot*Fugacities[0][p1]*expveff;
//         Sys->CurParPop[1][p1]+=ntop*sqrtnormtop*Fugacities[1][p1]*expveff;
//         BondNums[p1] += bond*expveff;
//       }
//     }
//   }
//   
//   Sys->TotCurPop[0] = Sys->TotCurPop[1] = 0.0;
//   for(p1=0;p1<Sys->NumPTypes;p1++){
//     Sys->TotCurPop[0] += Sys->CurParPop[0][p1];
//     Sys->TotCurPop[1] += Sys->CurParPop[1][p1];
//   }
// }


// for(p1=0;p1<NumPTypes;p1++){
//     ntop = nbot = bond = 0.0;
//     if(dz>(Lengths[0][p1] + Diams[0][p1]/2))nbot+=1.0;
//     if(dz>(Lengths[1][p1] + Diams[1][p1]/2))ntop+=1.0;
// 
//     for(p2=0;p2<NumPTypes;p2++){
//       sumlengths=Lengths[0][p1] + Lengths[1][p2];
//       
//       if(p1==p2 && BindPType[p1]>0.000001){
//         if(dz > sumlengths){        
//           intbot=sqrtnormtop*Fugacities[1][p2];
//           inttop=sqrtnormbot*Fugacities[0][p2];
//           
//           if(fabs(dz-sumlengths)<Sys->RanPType[p1]){
//             intbot*=exp(Sys->BindPType[p1]);
//             inttop*=exp(Sys->BindPType[p1]);
//             bond+=sqrtnormbot*sqrtnormtop*Fugacities[0][p1]*Fugacities[1][p2]*exp(Sys->BindPType[p1]);
//           }
//           nbot+=intbot;
//           ntop+=inttop;
//         }
//       }      
//     
//       else{
//         sumrads=(Diams[0][p1] + Diams[1][p2])/2.0;
//         if(dz>(sumlengths+sumrads)){
//           nbot+=sqrtnormtop*Fugacities[1][p2];
//           ntop+=sqrtnormbot*Fugacities[0][p2];
//         }
//       }
//     }
//     
//     Sys->BotMesh[i0][j0].ParticleProbs[p1*compnums + p2] = 
//     Sys->TopMesh[i0][j0].ParticleProbs[p1*compnums + p2] = ntop*sqrtnormtop*Fugacities[1][p1]*expveff;
//     Sys->BotMesh[i0][j0].BondNums[p1]     = Sys->TopMesh[i0][j0].BondNums[p1] = bond*expveff;
//   }
//   
//   Sys->BotMesh[i0][j0].ParticlePops[p1] = Sys->TopMesh[i0][j0].ParticlePops[p1] = 1.0;
//   for(p1=0;p1<Sys->NumPTypes;p1++){
//     Sys->BotMesh[i0][j0].ParticlePops[Sys->NumPTypes] -= Sys->BotMesh[i0][j0].ParticlePops[p1];
//     Sys->TopMesh[i0][j0].ParticlePops[Sys->NumPTypes] -= Sys->TopMesh[i0][j0].ParticlePops[p1];
//   }
