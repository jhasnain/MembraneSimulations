#include "mesh.h"

int ErrorCheck_MeshEnergies_BiasEnergy(sys *Sys){
  int Err=0;
  if( fabs(Sys->BiasEnergy - PrintBiasPotential(Sys)  )>1e-6 ){printf("\nMeshBias error Stored: %lf Computed: %lf\n", Sys->BiasEnergy, PrintBiasPotential(Sys));Err=4;}
  if( fabs(Sys->MeshEnergy - ReturnVanillaMeshEnergy(Sys)  )>1e-6 ){printf("\nMeshEnergy off Stored: %lf Computed: %lf\n", Sys->MeshEnergy, ReturnVanillaMeshEnergy(Sys));Err=4;}
  return Err;
}

int ErrorCheck_MeshNodes(sys *Sys){
  int Err=0;
  MeshPoint **TheMesh;
  vec testvec;
  double Orient, norm;
  
  for(int m=0;m<2;m++){
    if(m==0){TheMesh=Sys->BotMesh;Orient= 1.0;}
    if(m==1){TheMesh=Sys->TopMesh;Orient=-1.0;}
    
    for(int i=0;i<Sys->Nx;i++){
      for(int j=0;j<Sys->Ny;j++){
        ReturnNodeNormal(&norm, &testvec, TheMesh, i, j, Sys, Orient);
        if( fabs(TheMesh[i][j].nhat.x - testvec.x)>1e-6 ||
            fabs(TheMesh[i][j].nhat.y - testvec.y)>1e-6 ||
            fabs(TheMesh[i][j].nhat.z - testvec.z)>1e-6 ||
            fabs(TheMesh[i][j].nNorm - norm)>1e-6
        ){
          printf("Error in meshpoint calculation, %s i %d j %d,\n", (m==0?"Bot":"Top"), i, j );
          PrintVec(testvec, "Computed");
          PrintVec(TheMesh[i][j].nhat, "Stored\n");Err=2;
        }
      }  
    }
  }
  
  for(int i=0;i<Sys->Nx;i++){
    for(int j=0;j<Sys->Ny;j++){
      norm = Sys->TopMesh[i][j].Pos.z - Sys->BotMesh[i][j].Pos.z;
      if(norm<0.0)
      {
        printf("Error node points have crossed!\n");
        PrintVec(Sys->TopMesh[i][j].Pos, "TopPos");
        PrintVec(Sys->BotMesh[i][j].Pos, "BotPos");
        printf("nodevals: %d %d\n\n", i, j);
        Err=2;
      }
    }  
  }
  
 return Err; 
}

int CheckAffectedPNum(sys *Sys){
  int Err=0, l;
  for(int k=0;k<Sys->AffectedParticleNum;k++){
    for(l=0;l<Sys->AffectedParticleNum;l++){
      if(Sys->Unique_dummy_index[k]==Sys->Unique_dummy_index[l] && l!=k){
        printf("\nERROR found duplicate particles affected by meshmove:\n");
        for(int m=0;m<Sys->AffectedParticleNum;m++) printf("%d %d %d\n", m, Sys->Unique_dummy_index[m], Sys->AffectedParticleNum);
        Err=4;k=2*Sys->AffectedParticleNum;break;
      }
    }
  }
  return Err;
}

int ErrorCheck_Mesh(sys *Sys){
  int Err;
  Err = ErrorCheck_MeshEnergies_BiasEnergy(Sys);
  if(Err!=0)return Err;
  Err = ErrorCheck_MeshNodes(Sys);
  Err = CheckAffectedPNum(Sys);
  return Err;
}

void MeshPBC(int *left, int *right, int *out, int *in, sys *Sys, int i, int j){
  *left  = (i-1) < 0 ? Sys->Nx -1 : i-1;
  *right = (i+1) % Sys->Nx;
  *in    = (j-1) < 0 ? Sys->Ny -1 : j-1;
  *out   = (j+1) % Sys->Ny;
}

double ReturnVanillaMeshEnergy(sys *Sys){
  int i, j, in, out, left, right;
  double dxsqr=Sys->dxsqr, dysqr=Sys->dysqr;
  double partialenergy, totenergy=0.0;
  
  MeshPoint **MeshConfig = Sys->TopMesh;
  for(i=0; i<Sys->Nx;i++){
    for(j=0; j<Sys->Ny;j++){
      MeshPBC(&left, &right, &out, &in, Sys, i, j);
      partialenergy= (  (MeshConfig[left][j].Pos.z - 2.0*MeshConfig[i][j].Pos.z + MeshConfig[right][j].Pos.z)/dxsqr
                      + (MeshConfig[i][in].Pos.z   - 2.0*MeshConfig[i][j].Pos.z + MeshConfig[i][out].Pos.z  )/dysqr
                     );
      totenergy+= partialenergy*partialenergy;
    }
  }
  
  MeshConfig = Sys->BotMesh;
  for(i=0; i<Sys->Nx;i++){
    for(j=0; j<Sys->Ny;j++){
      MeshPBC(&left, &right, &out, &in, Sys, i, j);
      partialenergy= (  (MeshConfig[left][j].Pos.z - 2.0*MeshConfig[i][j].Pos.z + MeshConfig[right][j].Pos.z)/dxsqr
                      + (MeshConfig[i][in].Pos.z   - 2.0*MeshConfig[i][j].Pos.z + MeshConfig[i][out].Pos.z  )/dysqr
                     );
      totenergy+= partialenergy*partialenergy;
    }
  }
  return Sys->MeshEnergyPrefac*totenergy;
}

void ReturnNodeNormal(double *Norm, vec *Target, MeshPoint **MyMesh, int i0, int j0, sys *Sys, int Orient){
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  
  double nx = (MyMesh[left][j0].Pos.z - MyMesh[right][j0].Pos.z)/(2.0*Sys->dx);
  double ny = (MyMesh[i0][in].Pos.z - MyMesh[i0][out].Pos.z)/(2.0*Sys->dy);
  double norm = sqrt(nx*nx + ny*ny + 1.0);
  Target->x = Orient*nx/norm;
  Target->y = Orient*ny/norm;
  Target->z = Orient/norm;
  *Norm =  norm;
}

void ComputeNodeNormal(MeshPoint **MyMesh, int i0, int j0, sys *Sys, int Orient){
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  
  double nx = (MyMesh[left][j0].Pos.z - MyMesh[right][j0].Pos.z)/(2.0*Sys->dx);
  double ny = (MyMesh[i0][in].Pos.z - MyMesh[i0][out].Pos.z)/(2.0*Sys->dy);
  double norm = sqrt(nx*nx + ny*ny + 1.0);
  MyMesh[i0][j0].nhat.x = Orient*nx/norm;
  MyMesh[i0][j0].nhat.y = Orient*ny/norm;
  MyMesh[i0][j0].nhat.z = Orient/norm;
  MyMesh[i0][j0].nNorm  = norm;
}

void ComputeAllNodeNormals(sys *Sys){
  for(int i=0; i<Sys->Nx;i++){
    for(int j=0; j<Sys->Ny;j++){
      ComputeNodeNormal(Sys->TopMesh, i, j, Sys, TOP);
      ComputeNodeNormal(Sys->BotMesh, i, j, Sys, BOT);
    }
  }
}

void ComputeLocalCurvatures(MeshPoint **MyMesh, int i0, int j0, sys *Sys ){
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  
  MyMesh[i0][j0].lc_x = (MyMesh[left][j0].Pos.z - 2.0*MyMesh[i0][j0].Pos.z + MyMesh[right][j0].Pos.z)/Sys->dxsqr;
  MyMesh[i0][j0].lc_y = (MyMesh[i0][in].Pos.z   - 2.0*MyMesh[i0][j0].Pos.z + MyMesh[i0][out].Pos.z)/Sys->dysqr;
  MyMesh[i0][j0].CurvatureEnergy = Sys->MeshEnergyPrefac
                                  *(MyMesh[i0][j0].lc_x + MyMesh[i0][j0].lc_y)
                                  *(MyMesh[i0][j0].lc_x + MyMesh[i0][j0].lc_y);
}

void ComputeAllCurvatureEnergies(sys *Sys){
  Sys->MeshEnergy=0.0;
  for(int i=0; i<Sys->Nx;i++){
    for(int j=0; j<Sys->Ny;j++){
      ComputeLocalCurvatures(Sys->TopMesh, i, j, Sys);
      ComputeLocalCurvatures(Sys->BotMesh, i, j, Sys);
      Sys->MeshEnergy+=Sys->TopMesh[i][j].CurvatureEnergy 
                       + Sys->BotMesh[i][j].CurvatureEnergy;
    }
  }
}

void InitializeMeshEnergies(sys *Sys){
  ComputeAllCurvatureEnergies(Sys);
  ComputeAllNodeNormals(Sys);
  (*ComputeBiasPotential)(Sys);
  Sys->TotEnergy+=Sys->MeshEnergy + Sys->BiasEnergy;
}

int GlobalMeshShift(sys *Sys){
  celllist *CellList = &(Sys->CellList);
  int ***Neighbors = CellList->Neighbors;  
  int NumCx = CellList->NumCellsx, NumCy = CellList->NumCellsy;
  double suggestedz=gsl_ran_flat(Sys->rng, -Sys->MDisp, Sys->MDisp);
  double Lx=Sys->Lx, Ly=Sys->Ly;
  
  int CellNeighx, CellNeighy, CxBot, CyBot, ibot, itop, pBot, pTop;
  int acc=0, popind, tmp;
  particle CurrentParticle, PartnerParticle;
  double deltaMeshParticle=0.0;
  
  for(pBot=0;pBot<Sys->TotCurPop[0];pBot++){
    CurrentParticle = Sys->Particles[pBot];
    CurrentParticle.EndPoint.z +=suggestedz;
    deltaMeshParticle+=ComputeDeltaMPEnergy(Sys, &CurrentParticle, pBot);
    CurrentParticle.EndPoint.z -=suggestedz;
    if(deltaMeshParticle>=1.0){return 0;}
  }
  
  for(pTop=Sys->maxpop;pTop<Sys->maxpop + Sys->TotCurPop[1];pTop++){
    CurrentParticle = Sys->Particles[pTop];
    CurrentParticle.EndPoint.z -=suggestedz;
    deltaMeshParticle+=ComputeDeltaMPEnergy(Sys, &CurrentParticle, pTop);
    CurrentParticle.EndPoint.z +=suggestedz;
    if(deltaMeshParticle>=1.0){return 0;}
  }
  
  for(CxBot=0; CxBot<NumCx; CxBot++){
    for(CyBot=0; CyBot<NumCy; CyBot++){
      for(ibot = 0;ibot<CellList->BotPopulations[CxBot][CyBot];ibot++){
        pBot = CellList->BotParticleIds[CxBot][CyBot][ibot];
        CurrentParticle = Sys->Particles[pBot];
        
        TranslateParticle(&CurrentParticle, &CurrentParticle, 0.0, 0.0, suggestedz);
        for(popind=0;popind<9;popind++){
          CellNeighx = Neighbors[CxBot][CyBot][2*popind];
          CellNeighy = Neighbors[CxBot][CyBot][2*popind + 1];
          
          for(itop = 0;itop<CellList->TopPopulations[CellNeighx][CellNeighy];itop++){
            pTop = CellList->TopParticleIds[CellNeighx][CellNeighy][itop];
            PartnerParticle = Sys->Particles[pTop];
            if(   deltaRsqr(CurrentParticle.MidPoint, PartnerParticle.MidPoint, Lx, Ly) 
                < MaxParticleCutoffSqr(CurrentParticle, PartnerParticle)
              ){
              if ( CheckOverlap(CurrentParticle, PartnerParticle, Lx, Ly, &tmp) == ALOT){
                TranslateParticle(&CurrentParticle, &CurrentParticle, 0.0, 0.0, -suggestedz);
                return 0;
              }
            }
          }
        }
        TranslateParticle(&CurrentParticle, &CurrentParticle, 0.0, 0.0, -suggestedz);
      }
    } 
  }
  
  for(pBot=0;pBot<Sys->Nx;pBot++){
    for(pTop=0;pTop<Sys->Ny;pTop++){
      if(Sys->TopMesh[pBot][pTop].Pos.z-Sys->BotMesh[pBot][pTop].Pos.z - suggestedz<0.0){
        return 0;
      }
    }
  }
  
  double deltaEParticle = BondMeshShift(Sys, suggestedz);
  
  int choice  = Sys->choice = (gsl_rng_uniform(Sys->rng)>0.5 ? 0:1);
  if (choice == 1) suggestedz = -suggestedz;
  double DeltaEBias= MeshBiasGlobalChange(Sys, suggestedz);
  
  if ( gsl_rng_uniform(Sys->rng)< exp(-(deltaEParticle+DeltaEBias)) ){
    MeshPoint **MyMesh = (choice==0 ? Sys->BotMesh : Sys->TopMesh);
    for(int nx=0; nx<Sys->Nx; nx++){
      for(int ny=0; ny<Sys->Ny; ny++){
        MyMesh[nx][ny].Pos.z+=suggestedz;
      }
    }
    for(int p=choice*Sys->maxpop ;p<choice*Sys->maxpop + Sys->TotCurPop[choice];p++){
      Sys->Particles[p].Pos.z+=suggestedz;
      Sys->Particles[p].MidPoint.z+=suggestedz;
      Sys->Particles[p].EndPoint.z+=suggestedz;
    }
    UpdateBondLists(Sys, 2);
    Sys->TotParticleEnergy +=deltaEParticle;
    StoreBiasChange(Sys, DeltaEBias, choice);
    Sys->TotEnergy+=deltaEParticle+DeltaEBias;
    acc=1;
  }
  return acc;
}

void SuggestMeshMove(sys *Sys){
  Sys->p0 = EMPTY;
  Sys->i0=gsl_rng_uniform_int(Sys->rng, Sys->Nx);
  Sys->j0=gsl_rng_uniform_int(Sys->rng, Sys->Ny);
  Sys->choice = 2*gsl_rng_uniform_int(Sys->rng, 2) - 1 ;
  
  MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  Sys->suggestedz=MeshToModify[Sys->i0][Sys->j0].Pos.z + gsl_ran_flat(Sys->rng, -Sys->MDisp, Sys->MDisp);
}  

void ParticleToMeshDummyParticle(sys *Sys, int dummyindex, int p0){
  Sys->mesh_dummy_particles[dummyindex].type        = Sys->Particles[p0].type;
  Sys->mesh_dummy_particles[dummyindex].MeshId      = Sys->Particles[p0].MeshId;
  Sys->mesh_dummy_particles[dummyindex].MyDiam      = Sys->Particles[p0].MyDiam;
  Sys->mesh_dummy_particles[dummyindex].MyRad       = Sys->Particles[p0].MyRad;
  Sys->mesh_dummy_particles[dummyindex].MyRadSqr    = Sys->Particles[p0].MyRadSqr;
  Sys->mesh_dummy_particles[dummyindex].MyIntStr    = Sys->Particles[p0].MyIntStr;
  Sys->mesh_dummy_particles[dummyindex].MyIntCutoff = Sys->Particles[p0].MyIntCutoff;
  Sys->mesh_dummy_particles[dummyindex].MyIntCutoffsqr = Sys->Particles[p0].MyIntCutoffsqr;
  Sys->mesh_dummy_particles[dummyindex].MeshParticleCutoffsqr = Sys->Particles[p0].MeshParticleCutoffsqr;
}

void AssignNormVecsToDummy(vec *v, int ix, int jy, sys *Sys, MeshPoint **MyMesh, const int i0, const int j0,
                           const int left, const int right, const int out, const int in){
  
  CopyVec(v, MyMesh[ix][jy].nhat);rescalevec(v, MyMesh[ix][jy].nNorm);
  
       if(   ix==left  && jy == j0 ){CopyVec(v, Sys->nhatLeft); rescalevec(v, Sys->nNormLeft);}
  else if(   ix==right && jy == j0 ){CopyVec(v, Sys->nhatRight);rescalevec(v, Sys->nNormRight);}
  else if(   ix==i0    && jy == in ){CopyVec(v, Sys->nhatIn);   rescalevec(v, Sys->nNormIn);}
  else if(   ix==i0    && jy == out){CopyVec(v, Sys->nhatOut);  rescalevec(v, Sys->nNormOut);}
  
}

void PopulateDummyMeshParticles(sys *Sys, int dummyindex, int pcurrent, int i0, int j0, int left, int right, int out, int in){
  MeshPoint **MyMesh = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  double length = Sys->Particles[pcurrent].MyLength;
  particle *MyParticle = &(Sys->mesh_dummy_particles[dummyindex]);
  
  vec diffvec, v1, v2 , v3, v4, dummy_nhat;
  int ip, jp, xNeigh, yNeigh;
  
  Position_To_Index_PBC(&ip, MyParticle->Pos.x, Sys->dx, Sys->Nx);
  Position_To_Index_PBC(&jp, MyParticle->Pos.y, Sys->dy, Sys->Ny);
  
  twodvecdiff(&(diffvec), MyParticle->Pos, MyMesh[ip][jp].Pos);
  NearestImageConv(&(diffvec), Sys->Lx, Sys->Ly);
  diffvec.x/=Sys->dx;
  diffvec.y/=Sys->dy;
  
  xNeigh = IndexPBC(ip + ((diffvec.x > 0) ? 1 : -1), Sys->Nx);
  yNeigh = IndexPBC(jp + ((diffvec.y > 0) ? 1 : -1), Sys->Ny);
  
  v1.z = (ip==i0     &&     jp == j0) ? Sys->suggestedz : MyMesh[ip][jp].Pos.z;
  v2.z = (xNeigh==i0 &&     jp == j0) ? Sys->suggestedz : MyMesh[xNeigh][jp].Pos.z;
  v3.z = (ip==i0     && yNeigh == j0) ? Sys->suggestedz : MyMesh[ip][yNeigh].Pos.z;
  v4.z = (xNeigh==i0 && yNeigh == j0) ? Sys->suggestedz : MyMesh[xNeigh][yNeigh].Pos.z;
  
  InterpolateProperty(&(MyParticle->Pos.z), diffvec, v1.z, v2.z, v3.z, v4.z);
                      
  AssignNormVecsToDummy(&v1,     ip, jp,     Sys, MyMesh, i0, j0, left, right, out, in);
  AssignNormVecsToDummy(&v2, xNeigh, jp,     Sys, MyMesh, i0, j0, left, right, out, in);
  AssignNormVecsToDummy(&v3,     ip, yNeigh, Sys, MyMesh, i0, j0, left, right, out, in);
  AssignNormVecsToDummy(&v4, xNeigh, yNeigh, Sys, MyMesh, i0, j0, left, right, out, in);
  
  InterpolateProperty(&(dummy_nhat.x), diffvec, v1.x, v2.x, v3.x, v4.x);
  InterpolateProperty(&(dummy_nhat.y), diffvec, v1.y, v2.y, v3.y, v4.y);
  InterpolateProperty(&(dummy_nhat.z), diffvec, v1.z, v2.z, v3.z, v4.z);
  MyParticle->PhaseFactorEnergy = -0.5*log(1.0 + dummy_nhat.x*dummy_nhat.x + dummy_nhat.y*dummy_nhat.y);
  
  threed_NormalizeVec(&dummy_nhat);  
  threedvecsumscalar(&(MyParticle->MidPoint), MyParticle->Pos, dummy_nhat, length/2);
  threedvecsumscalar(&(MyParticle->EndPoint), MyParticle->Pos, dummy_nhat, length);
  
  MyParticle->MyLength = length;
  CopyVec(&(MyParticle->nhat), dummy_nhat);
  ParticleToMeshDummyParticle(Sys, dummyindex, pcurrent);
}

void PopulateAffectedParticles(sys *Sys){
  MeshPoint **MyMesh = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  celllist *CellList = &(Sys->CellList);
  int i0 = Sys->i0, j0 = Sys->j0;
  int **Populations = (Sys->choice == TOP) ? CellList->TopPopulations : CellList->BotPopulations;
  int ***ParticleIds= (Sys->choice == TOP) ? CellList->TopParticleIds : CellList->BotParticleIds;
  double twodx = 2.0*Sys->dx, twody = 2.0*Sys->dy;
  
  int MyCx, MyCy, Cx, Cy, left, right, out, in;
  int i, k, p;
  vec diffvec;
  
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  Position_To_Index_PBC(&(MyCx), MyMesh[i0][j0].Pos.x, CellList->dCx, CellList->NumCellsx);
  Position_To_Index_PBC(&(MyCy), MyMesh[i0][j0].Pos.y, CellList->dCy, CellList->NumCellsy);
  
  Sys->AffectedParticleNum = 0;
  
  for(i=0;i<9;i++){
    Cx = CellList->Neighbors[MyCx][MyCy][2*i];
    Cy = CellList->Neighbors[MyCx][MyCy][2*i + 1];
    for(k=0;k<Populations[Cx][Cy];k++){
      p = ParticleIds[Cx][Cy][k];
      twodvecdiff(&diffvec, Sys->Particles[p].Pos, MyMesh[i0][j0].Pos );
      NearestImageConv(&diffvec, Sys->Lx, Sys->Ly);
      if( fabs(diffvec.x)<twodx && fabs(diffvec.y)<twody ){
        Sys->Unique_dummy_index[Sys->AffectedParticleNum]=p;
        CopyVec(&(Sys->mesh_dummy_particles[Sys->AffectedParticleNum].Pos), Sys->Particles[p].Pos);
        PopulateDummyMeshParticles(Sys, Sys->AffectedParticleNum, p, i0, j0, left, right, out, in);
        Sys->AffectedParticleNum++;
      }
    }
  }
}

void MeshParticleChange(sys *Sys, double *Energies){
  double Lx=Sys->Lx, Ly=Sys->Ly;
  int *MovedParticleIndices = Sys->Unique_dummy_index;
  celllist *CellList = &(Sys->CellList);
  int ***Neighbors = CellList->Neighbors;
  
  particle current_Particle, PartnerParticleNew, PartnerParticleOld;
  int tmp, flag, current_p, MyCx, MyCy, Cx, Cy;
  int i_affect, j_affect, i_cellNeigh, i_cellpop, i_partner;
  
  Energies[0]=0.0;Energies[1]=0.0;
  
  for(i_affect=0;i_affect<Sys->AffectedParticleNum;i_affect++){
    current_p = MovedParticleIndices[i_affect];
    current_Particle = Sys->mesh_dummy_particles[i_affect];
    
    Position_To_Index_PBC(&MyCx, current_Particle.Pos.x, CellList->dCx, CellList->NumCellsx);
    Position_To_Index_PBC(&MyCy, current_Particle.Pos.y, CellList->dCy, CellList->NumCellsy);
    
    Energies[0]+=current_Particle.PhaseFactorEnergy - Sys->Particles[current_p].PhaseFactorEnergy;
    
    for(i_cellNeigh=0;i_cellNeigh<9;i_cellNeigh++){
      Cx = Neighbors[MyCx][MyCy][2*i_cellNeigh];
      Cy = Neighbors[MyCx][MyCy][2*i_cellNeigh + 1];
      for(i_cellpop = 0;i_cellpop<CellList->TopPopulations[Cx][Cy];i_cellpop++){
        i_partner = CellList->TopParticleIds[Cx][Cy][i_cellpop];
        if(current_p != i_partner){
          flag=2;
          for(j_affect = 0;j_affect < Sys->AffectedParticleNum; j_affect++){
            if (MovedParticleIndices[j_affect] == i_partner){
              if (j_affect < i_affect) {flag=0;break;}
              else{
                PartnerParticleNew = Sys->mesh_dummy_particles[j_affect];
                PartnerParticleOld = Sys->Particles[i_partner];
                flag=1;break;
              }
            }
          }
          if (flag == 2){PartnerParticleNew = PartnerParticleOld = Sys->Particles[i_partner];}
          if(flag>0){
            if(   deltaRsqr(current_Particle.MidPoint, PartnerParticleNew.MidPoint, Lx, Ly) 
                < MaxParticleCutoffSqr(current_Particle, PartnerParticleNew)
              ){
              if ( CheckOverlap(current_Particle, PartnerParticleNew, Lx, Ly, &tmp) == ALOT){Energies[1]=ALOT;return;}
            }
          } 
        }
      }
      
      for(i_cellpop = 0;i_cellpop<CellList->BotPopulations[Cx][Cy];i_cellpop++){
        i_partner = CellList->BotParticleIds[Cx][Cy][i_cellpop];
        if(current_p != i_partner){
          flag=2;
          for(j_affect = 0;j_affect < Sys->AffectedParticleNum; j_affect++){
            if (MovedParticleIndices[j_affect] == i_partner){
              if (j_affect < i_affect) {flag=0;break;}
              else{
                PartnerParticleNew = Sys->mesh_dummy_particles[j_affect];
                PartnerParticleOld = Sys->Particles[i_partner];
                flag=1;break;
              }
            }
          }
          if (flag == 2){PartnerParticleNew = PartnerParticleOld = Sys->Particles[i_partner];}
          if(flag>0){
            if(   deltaRsqr(current_Particle.MidPoint, PartnerParticleNew.MidPoint, Lx, Ly) 
                < MaxParticleCutoffSqr(current_Particle, PartnerParticleNew)
              ){
              if ( CheckOverlap(current_Particle, PartnerParticleNew, Lx, Ly, &tmp) == ALOT){Energies[1]=ALOT;return;}
            }
          } 
        }
      }
    } 
  }
}

double MeshEnergyAndNhatChange(sys *Sys){
  int i0=Sys->i0, j0=Sys->j0;
  double suggestedz = Sys->suggestedz;
  double dxsqr=Sys->dxsqr, dysqr=Sys->dysqr;
  double twodx=2.0*Sys->dx, twody=2.0*Sys->dy;
  
  MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  MeshPoint **OtherMesh    = (Sys->choice == TOP) ? Sys->BotMesh : Sys->TopMesh;
  if( Sys->choice*(OtherMesh[i0][j0].Pos.z - suggestedz) < 0.0 ) return ALOT;
  
  vec *lcCent = &(Sys->lcCent);
  vec *lcLeft = &(Sys->lcLeft), *lcRight = &(Sys->lcRight);
  vec *lcIn   = &(Sys->lcIn),   *lcOut   = &(Sys->lcOut);
  
  vec *nhatRight = &(Sys->nhatRight), *nhatLeft = &(Sys->nhatLeft);
  vec *nhatIn = &(Sys->nhatIn), *nhatOut = &(Sys->nhatOut);
  double nx, ny, norm, Orient=Sys->choice;
  
  int in, out, left, right;
  int in_in, out_out, left_left, right_right;
  
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  lcCent->x = (MeshToModify[left][j0].Pos.z - 2.0*suggestedz + MeshToModify[right][j0].Pos.z)/dxsqr;
  lcCent->y = (MeshToModify[i0][in].Pos.z - 2.0*suggestedz + MeshToModify[i0][out].Pos.z)/dysqr;
  
  MeshPBC(&left_left, &right_right, &out_out, &in_in, Sys, right, j0);
  lcRight->x = (suggestedz - 2.0*MeshToModify[right][j0].Pos.z + MeshToModify[right_right][j0].Pos.z)/dxsqr;
  lcRight->y = MeshToModify[right][j0].lc_y;
  
  nx = (suggestedz - MeshToModify[right_right][j0].Pos.z)/twodx;
  ny = (MeshToModify[right][in_in].Pos.z - MeshToModify[right][out_out].Pos.z)/twody;
  norm = sqrt(nx*nx + ny*ny + 1.0);
  nhatRight->x = Orient*nx/norm;
  nhatRight->y = Orient*ny/norm;
  nhatRight->z = Orient/norm;
  Sys->nNormRight = norm;
    
  MeshPBC(&left_left, &right_right, &out_out, &in_in, Sys, left, j0);
  lcLeft->x = (MeshToModify[left_left][j0].Pos.z - 2.0*MeshToModify[left][j0].Pos.z + suggestedz)/dxsqr;
  lcLeft->y = MeshToModify[left][j0].lc_y;
  
  nx = (MeshToModify[left_left][j0].Pos.z - suggestedz)/twodx;
  ny = (MeshToModify[left][in_in].Pos.z - MeshToModify[left][out_out].Pos.z)/twody;
  norm = sqrt(nx*nx + ny*ny + 1.0);
  nhatLeft->x = Orient*nx/norm;
  nhatLeft->y = Orient*ny/norm;
  nhatLeft->z = Orient/norm;
  Sys->nNormLeft = norm;
    
  MeshPBC(&left_left, &right_right, &out_out, &in_in, Sys, i0, out);
  lcOut->x = MeshToModify[i0][out].lc_x;
  lcOut->y = (suggestedz - 2.0*MeshToModify[i0][out].Pos.z + MeshToModify[i0][out_out].Pos.z)/dysqr;
  
  nx = (MeshToModify[left_left][out].Pos.z - MeshToModify[right_right][out].Pos.z)/twodx;
  ny = (suggestedz - MeshToModify[i0][out_out].Pos.z)/twody;
  norm = sqrt(nx*nx + ny*ny + 1.0);
  nhatOut->x = Orient*nx/norm;
  nhatOut->y = Orient*ny/norm;
  nhatOut->z = Orient/norm;
  Sys->nNormOut = norm;
  
  MeshPBC(&left_left, &right_right, &out_out, &in_in, Sys, i0, in);
  lcIn->x = MeshToModify[i0][in].lc_x;
  lcIn->y = (MeshToModify[i0][in_in].Pos.z - 2.0*MeshToModify[i0][in].Pos.z + suggestedz)/dysqr;
  
  nx = (MeshToModify[left_left][in].Pos.z - MeshToModify[right_right][in].Pos.z)/twodx;
  ny = (MeshToModify[i0][in_in].Pos.z - suggestedz)/twody;
  norm = sqrt(nx*nx + ny*ny + 1.0);
  nhatIn->x = Orient*nx/norm;
  nhatIn->y = Orient*ny/norm;
  nhatIn->z = Orient/norm;
  Sys->nNormIn = norm;
  
  double Eold =  MeshToModify[i0][j0].CurvatureEnergy 
               + MeshToModify[left][j0].CurvatureEnergy + MeshToModify[right][j0].CurvatureEnergy 
               + MeshToModify[i0][in].CurvatureEnergy   + MeshToModify[i0][out].CurvatureEnergy;
  
  double Enew = Sys->MeshEnergyPrefac*
                ( (lcCent->x + lcCent->y)*(lcCent->x + lcCent->y) 
                 +(lcRight->x + lcRight->y)*(lcRight->x + lcRight->y) 
                 +(lcLeft->x + lcLeft->y)*(lcLeft->x + lcLeft->y)
                 +(lcIn->x + lcIn->y)*(lcIn->x + lcIn->y)
                 +(lcOut->x + lcOut->y)*(lcOut->x + lcOut->y)
                );
  return Enew-Eold; 
}

double ComputePMChange(sys *Sys){
  double DeltaEMeshParticle=0.0;
  
  for(int k=0;k<Sys->AffectedParticleNum;k++){
    DeltaEMeshParticle += ComputeDeltaMPEnergy(Sys, &(Sys->mesh_dummy_particles[k]), Sys->Unique_dummy_index[k]);
  }
  
  if(DeltaEMeshParticle<ALOT){
    double Lx = Sys->Lx, Ly = Sys->Ly, dx = Sys->dx, dy = Sys->dy;
    int Nx = Sys->Nx, Ny = Sys->Ny, i0=Sys->i0, j0=Sys->j0;
    
    MeshPoint **MyMesh = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
    vec MyPos;CopyVec(&MyPos, MyMesh[i0][j0].Pos);MyPos.z = Sys->suggestedz;
      
    celllist *CellList = &(Sys->CellList);
    int **MyPopulations = (Sys->choice == TOP) ? CellList->BotPopulations : CellList->TopPopulations;
    int  ***ParticleIds = (Sys->choice == TOP) ? CellList->BotParticleIds : CellList->TopParticleIds;
    
    double left  = Sys->dx*IndexPBC(i0 - 1, Nx);
    double right = Sys->dx*IndexPBC(i0 + 1, Nx);
    double down  = Sys->dy*IndexPBC(j0 - 1, Ny);
    double up    = Sys->dy*IndexPBC(j0 + 1, Ny);
    double diffleft, diffright, diffup, diffdown, z1, z2, z3, z4;
    int inew, jnew, xNeigh, yNeigh;
    
    int MyCx, MyCy, Cx, Cy, p, particleId;
    vec diffvec, v, ParticlePos;
    
    Position_To_Index_PBC(&MyCx, MyPos.x, CellList->dCx, CellList->NumCellsx);
    Position_To_Index_PBC(&MyCy, MyPos.y, CellList->dCy, CellList->NumCellsy);
    int *MyNeighbors = CellList->Neighbors[MyCx][MyCy];
    
    for(int k=0;k<9;k++){
      Cx = MyNeighbors[2*k];Cy = MyNeighbors[2*k+1];
      for(p=0;p<MyPopulations[Cx][Cy];p++){
        particleId = ParticleIds[Cx][Cy][p];
        CopyVec(&ParticlePos, Sys->Particles[particleId].EndPoint);
        
        diffleft  = left  - ParticlePos.x;ScalarNearestImageConv( &diffleft,  Lx );
        diffright = right - ParticlePos.x;ScalarNearestImageConv( &diffright, Lx );
        diffup    = up    - ParticlePos.y;ScalarNearestImageConv( &diffup,    Ly );
        diffdown  = down  - ParticlePos.y;ScalarNearestImageConv( &diffdown,  Ly );
        
        if ( (diffleft < 0.0 && diffright > 0.0) && (diffdown < 0.0 && diffup > 0.0 ) ) {
          Position_To_Index_PBC(&inew, ParticlePos.x, dx, Sys->Nx);
          Position_To_Index_PBC(&jnew, ParticlePos.y, dy, Sys->Ny);
    
          twodvecdiff(&(diffvec), ParticlePos, MyMesh[inew][jnew].Pos);
          NearestImageConv(&(diffvec), Lx, Ly);
          v.x=diffvec.x/Sys->dx;
          v.y=diffvec.y/Sys->dy;
    
          xNeigh = IndexPBC(inew + ((diffvec.x > 0) ? 1 : -1), Nx);
          yNeigh = IndexPBC(jnew + ((diffvec.y > 0) ? 1 : -1), Ny);
    
          z1 = (inew==i0   && jnew == j0)   ? Sys->suggestedz : MyMesh[inew][jnew].Pos.z;
          z2 = (xNeigh==i0 && jnew == j0)   ? Sys->suggestedz : MyMesh[xNeigh][jnew].Pos.z;
          z3 = (inew==i0   && yNeigh == j0) ? Sys->suggestedz : MyMesh[inew][yNeigh].Pos.z;
          z4 = (xNeigh==i0 && yNeigh == j0) ? Sys->suggestedz : MyMesh[xNeigh][yNeigh].Pos.z;
          
          InterpolateProperty( &(diffvec.z), v, z1, z2, z3, z4 );
          diffvec.z-= ParticlePos.z;

          if( Sys->Particles[particleId].MeshId*diffvec.z < 0.0 || diffvec.z*diffvec.z < Sys->Particles[particleId].MeshParticleCutoffsqr){
            Sys->dummy_particle.MeshParticleEnergy = ALOT;
            DeltaEMeshParticle+=ALOT - Sys->Particles[particleId].MeshParticleEnergy;
            Sys->p0 = particleId;
          }
        }
      }
    }
  }
  return DeltaEMeshParticle;
}

void StoreParticleMeshChange(sys *Sys, double *Energies, double DeltaEMeshParticle){
  int p;  
  for(int k=0;k<Sys->AffectedParticleNum;k++){
    p = Sys->Unique_dummy_index[k]; 
    CopyVec(&(Sys->Particles[p].nhat)    , Sys->mesh_dummy_particles[k].nhat);
    CopyVec(&(Sys->Particles[p].Pos)     , Sys->mesh_dummy_particles[k].Pos );
    CopyVec(&(Sys->Particles[p].MidPoint), Sys->mesh_dummy_particles[k].MidPoint );
    CopyVec(&(Sys->Particles[p].EndPoint), Sys->mesh_dummy_particles[k].EndPoint );
    Sys->Particles[p].MeshParticleEnergy = Sys->mesh_dummy_particles[k].MeshParticleEnergy;
    Sys->Particles[p].PhaseFactorEnergy  = Sys->mesh_dummy_particles[k].PhaseFactorEnergy;
  }
  
  if(Sys->p0 !=EMPTY) Sys->Particles[Sys->p0].MeshParticleEnergy = Sys->dummy_particle.MeshParticleEnergy;
  Sys->p0 = EMPTY;
  Sys->TotPhaseFactorEnergy  += Energies[0];
  Sys->TotParticleEnergy     += Energies[1];
  Sys->TotMeshParticleEnergy += DeltaEMeshParticle;
}

void StoreMeshChange(sys *Sys, double DeltaEMesh){
  int i0 = Sys->i0, j0 = Sys->j0;
  MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  int in, out, left, right;
  MeshPBC(&left, &right, &out, &in, Sys, i0, j0);
  
  MeshToModify[i0][j0].Pos.z = Sys->suggestedz;
  MeshToModify[i0][j0].lc_x = Sys->lcCent.x;MeshToModify[i0][j0].lc_y = Sys->lcCent.y;
  MeshToModify[i0][j0].CurvatureEnergy    = Sys->MeshEnergyPrefac*
                                            (Sys->lcCent.x + Sys->lcCent.y)*(Sys->lcCent.x + Sys->lcCent.y);
  
  MeshToModify[right][j0].lc_x = Sys->lcRight.x;
  MeshToModify[right][j0].CurvatureEnergy = Sys->MeshEnergyPrefac*
                                            (Sys->lcRight.x + Sys->lcRight.y)*(Sys->lcRight.x + Sys->lcRight.y);
  MeshToModify[right][j0].nNorm = Sys->nNormRight;
  
                      
  MeshToModify[left][j0].lc_x = Sys->lcLeft.x;
  MeshToModify[left][j0].CurvatureEnergy  = Sys->MeshEnergyPrefac*
                                            (Sys->lcLeft.x + Sys->lcLeft.y)*(Sys->lcLeft.x + Sys->lcLeft.y);
  MeshToModify[left][j0].nNorm = Sys->nNormLeft;
  
  MeshToModify[i0][in].lc_y = Sys->lcIn.y;
  MeshToModify[i0][in].CurvatureEnergy    = Sys->MeshEnergyPrefac*
                                            (Sys->lcIn.x + Sys->lcIn.y)*(Sys->lcIn.x + Sys->lcIn.y);
  MeshToModify[i0][in].nNorm = Sys->nNormIn;
                                            
  MeshToModify[i0][out].lc_y = Sys->lcOut.y;
  MeshToModify[i0][out].CurvatureEnergy   = Sys->MeshEnergyPrefac*
                                            (Sys->lcOut.x + Sys->lcOut.y)*(Sys->lcOut.x + Sys->lcOut.y);
  MeshToModify[i0][out].nNorm = Sys->nNormOut;
                                            
  CopyVec(&(MeshToModify[left][j0].nhat) , Sys->nhatLeft);
  CopyVec(&(MeshToModify[right][j0].nhat), Sys->nhatRight);
  CopyVec(&(MeshToModify[i0][in].nhat)   , Sys->nhatIn);
  CopyVec(&(MeshToModify[i0][out].nhat)  , Sys->nhatOut);
  
  Sys->MeshEnergy+=DeltaEMesh;
}

int PerformMeshMove(sys *Sys){
  int acc=0;
  if(gsl_rng_uniform(Sys->rng) > Sys->meshglobalrate){
    double Energies[2];Sys->MeshShift=0;
    SuggestMeshMove(Sys);
    double DeltaEMesh=MeshEnergyAndNhatChange(Sys);
    if(DeltaEMesh<ALOT){
      PopulateAffectedParticles(Sys);
      double DeltaEMeshParticle = ComputePMChange(Sys);
      if(DeltaEMeshParticle<ALOT){
        MeshParticleChange(Sys, Energies);
        if(Energies[1]<ALOT){
          Energies[1]+=BondMeshMove(Sys);
          double DeltaEBias=MeshBiasChange(Sys);
          double DeltaETot=DeltaEMesh + DeltaEBias + Energies[0] + Energies[1] + DeltaEMeshParticle;
          if ( gsl_rng_uniform(Sys->rng)< exp(-DeltaETot) ){
            StoreMeshChange(Sys, DeltaEMesh);
            UpdateBondLists(Sys, 1);
            StoreParticleMeshChange(Sys, Energies, DeltaEMeshParticle);
            StoreBiasChange(Sys, DeltaEBias, (Sys->choice==BOT ? 0:1));
            Sys->TotEnergy+=DeltaETot;
            acc=1;
          }
        }
      }
    }
  }
  else{acc = GlobalMeshShift(Sys);Sys->MeshShift=1;}
  
  return acc;
}
