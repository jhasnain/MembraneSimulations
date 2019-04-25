#include "mesh.h"
#include "effectiveparticlepotential.h"

int ErrorCheck_MeshEnergies_BiasEnergy(sys *Sys){
  int Err=0;
  if( fabs(Sys->BiasEnergy - PrintBiasPotential(Sys)  )>1e-6 ){printf("\nMeshBias error Stored: %lf Computed: %lf\n", Sys->BiasEnergy, PrintBiasPotential(Sys));Err=2;}
  if( fabs(Sys->MeshEnergy - ReturnVanillaMeshEnergy(Sys)  )>1e-6 ){printf("\nMeshEnergy off Stored: %lf Computed: %lf\n", Sys->MeshEnergy, ReturnVanillaMeshEnergy(Sys));Err=2;}
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
          printf("Error in meshpoint calculation, %s i %d j %d,\n", (m==0 ? "Bot":"Top"), i, j );
          PrintVec(testvec, "Computed");
          PrintVec(TheMesh[i][j].nhat, "Stored\n");Err=1;
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
        Err=1;
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
  InitializeVeff(Sys);
  (*ComputeBiasPotential)(Sys);
}

int GlobalMeshShift(sys *Sys){
  Sys->suggestedz = gsl_ran_flat(Sys->rng, -Sys->MDisp, Sys->MDisp);
  Sys->choice = 2*gsl_rng_uniform_int(Sys->rng, 2) - 1;
  
  int acc=0;
  double DeltaEBias = MeshBiasGlobalChange(Sys);
  double DeltaEVeff = VeffGlobalShiftChange(Sys);
  
  if ( gsl_rng_uniform(Sys->rng)< exp(-(DeltaEBias + DeltaEVeff)) ){
    MeshPoint **MyMesh = (Sys->choice==BOT ? Sys->BotMesh : Sys->TopMesh);
    for(int nx=0; nx<Sys->Nx; nx++){
      for(int ny=0; ny<Sys->Ny; ny++){
        MyMesh[nx][ny].Pos.z+=Sys->suggestedz;
      }
    } 
    StoreBiasChange(Sys, DeltaEBias);
    StoreVeffGlobalShiftChange(Sys, DeltaEVeff);
    Sys->TotEnergy+=DeltaEVeff+DeltaEBias;
    acc=1;  
  }
  
  return acc;
}

void SuggestMeshMove(sys *Sys){
  Sys->i0=gsl_rng_uniform_int(Sys->rng, Sys->Nx);
  Sys->j0=gsl_rng_uniform_int(Sys->rng, Sys->Ny);
  Sys->choice = 2*gsl_rng_uniform_int(Sys->rng, 2) - 1 ;
  
  MeshPoint **MeshToModify = (Sys->choice == TOP) ? Sys->TopMesh : Sys->BotMesh;
  Sys->suggestedz=MeshToModify[Sys->i0][Sys->j0].Pos.z + gsl_ran_flat(Sys->rng, -Sys->MDisp, Sys->MDisp);
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
    Sys->MeshShift=0;
    SuggestMeshMove(Sys);
    double DeltaEMesh=MeshEnergyAndNhatChange(Sys);    
    if(DeltaEMesh <ALOT){
      double DeltaEBias = MeshBiasChange(Sys);
      double DeltaEVeff = VeffChange(Sys);
      double DeltaETot = DeltaEMesh + DeltaEBias + DeltaEVeff;
      
      if ( gsl_rng_uniform(Sys->rng)< exp(-DeltaETot) ){
        StoreMeshChange(Sys, DeltaEMesh);
        StoreBiasChange(Sys, DeltaEBias, (Sys->choice==BOT ? 0:1));
        StoreVeffChange(Sys, DeltaEVeff);
        Sys->TotEnergy+=DeltaETot;
        acc=1;
      }
    }
  }
  
  else{acc = GlobalMeshShift(Sys);Sys->MeshShift=1;}  
  return acc;
}
