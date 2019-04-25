#include "initialize.h"


int CountNumberOfCompletedRuns(char *OutputDir, char *tmp){
  int index=0;
  DIR * dirp;
  struct dirent * entry;
  char RestartDir[1024];
  
  strcpy(RestartDir, OutputDir);
  strcat(RestartDir, "/Restarts/");
//   printf("%s %s\n", RestartDir, OutputDir);
  dirp = opendir(RestartDir);
  if(dirp){
    while ((entry = readdir(dirp)) != NULL) {
      if (entry->d_type == DT_REG && strstr(entry->d_name, tmp) ){ 
        index++;
      }
    }
  }
  else{printf("\nERROR! Could not find %s on this machine\n", RestartDir);exit(1);}
  closedir(dirp);
  return index;
}

void SetOutPutFiles(sys *Sys, struct arguments *arguments ){
//Set filepaths
  char tmp[1024], tmp2[1024];
  if(0 == strcmp(arguments->Tag,"")){tmp[0] = '\0';}
  else{snprintf(tmp, 1024, "%s_", arguments->Tag);}
  
  snprintf(tmp2, 1024, "%sNodes_%d_NumPTypes_%d_kc_%.2lf_A_%.1lf_theta0_%lf", 
                tmp,
                arguments->Nx*arguments->Ny,
                arguments->NumPType,
                arguments->kc,
                arguments->Lx*arguments->Ly,
                arguments->theta0
         );
  
  int index = CountNumberOfCompletedRuns(arguments->OutputDir, tmp2);
  
  snprintf(Sys->Trajfile,    1024, "%s/Trajectories/Traj_%s_Run_%d_.xyz",
          arguments->OutputDir, tmp2, index);
  snprintf(Sys->Datafile,    1024, "%s/DataFiles/Data_%s_Run_%d_.dat",
          arguments->OutputDir, tmp2, index);
  snprintf(Sys->Histfile, 1024, "%s/DataFiles/Hists/Zhist_%s_Run_%d_.hist",
          arguments->OutputDir, tmp2, index);
  snprintf(Sys->Corrfile, 1024, "%s/DataFiles/Hists/Corr_%s_Run_%d_.hist",
          arguments->OutputDir, tmp2, index);
  snprintf(Sys->Restartfile, 1024, "%s/Restarts/Restart_%s_Run_%d_.rest",
          arguments->OutputDir, tmp2, index);
  
  Sys->fileptr_traj    = fopen(Sys->Trajfile, "w");
  Sys->fileptr_data    = fopen(Sys->Datafile, "w");
  
}

void ExtractValuesFromStringList( void *array, char *opt, char *str, int flavornums){
  int i=0;
  if( strcmp(opt, "INT") == 0 ){
    ((int *)array)[0]= atoi(strtok(str , " "));
    for (i=1;i<flavornums;i++){
      ((int *)array)[i]=atoi(strtok(NULL , " "));
    }
  }
  else if( strcmp(opt, "DOUBLE")== 0 ){
   ((double *)array)[0]= atof(strtok(str , " "));
   for (i=1;i<flavornums;i++){
     ((double *)array)[i]=atof(strtok(NULL , " "));
    }
  }
  else{printf("Could not determine what this OPT: %s is\n", opt);}
}

void CreateMeshes(sys *Sys, struct arguments *arguments){
  int i, j, k;
  double x;
  
  Sys->Nx      = arguments->Nx;
  Sys->Ny      = arguments->Ny;
  Sys->MTot    = arguments->Nx*arguments->Ny;
  Sys->twoMTot = 2*Sys->MTot;
    
  Sys->Lx    = arguments->Lx;
  Sys->Ly    = arguments->Ly;
  Sys->Lxon2 = arguments->Lx/2;
  Sys->Lyon2 = arguments->Ly/2;
  Sys->A     = arguments->Lx*arguments->Ly;
  
  Sys->dx    = Sys->Lx/Sys->Nx;
  Sys->dy    = Sys->Ly/Sys->Ny;
  Sys->dxsqr = Sys->dx*Sys->dx;
  Sys->dysqr = Sys->dy*Sys->dy;
  
  Sys->dA = Sys->dx*Sys->dy;
  Sys->kc = arguments->kc;
  Sys->MeshEnergyPrefac = (Sys->kc/2)*Sys->dA;
  
  Sys->TopMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  Sys->BotMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  
  for (i=0;i<Sys->Nx;i++){
    Sys->TopMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->TopMesh[i][j].Pos.x            = x;
      Sys->TopMesh[i][j].Pos.y            = Sys->dy*j;
      Sys->TopMesh[i][j].Pos.z            = arguments->z0/2;
      Sys->TopMesh[i][j].ParticleProbs    = (double **)malloc((arguments->NumPType + 1)*sizeof(double *));
      for(k=0;k<(arguments->NumPType + 1);k++) Sys->TopMesh[i][j].ParticleProbs[k] = (double *)malloc((arguments->NumPType + 1)*sizeof(double));
      Sys->TopMesh[i][j].BondNums         = (double *)malloc((arguments->NumPType)*sizeof(double));
    }
  }
  
  for (i=0;i<Sys->Nx;i++){
    Sys->BotMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->BotMesh[i][j].Pos.x            = x;
      Sys->BotMesh[i][j].Pos.y            = Sys->dy*j;
      Sys->BotMesh[i][j].Pos.z            = -arguments->z0/2;
      Sys->BotMesh[i][j].ParticleProbs    = (double **)malloc((arguments->NumPType + 1)*sizeof(double *));
      for(k=0;k<(arguments->NumPType + 1);k++) Sys->BotMesh[i][j].ParticleProbs[k] = (double *)malloc((arguments->NumPType + 1)*sizeof(double));
      Sys->BotMesh[i][j].BondNums         = (double *)malloc((arguments->NumPType)*sizeof(double));
    }
  }
  
  Sys->MDisp  = arguments->MDisp;
}

void CreateBiasPotential(sys *Sys, struct arguments *arguments){
  if (arguments->BiasStr<1e-6)   Sys->meshglobalrate = 0.0;
  else if (arguments->kc>1000.0) Sys->meshglobalrate = 1.0;
  else Sys->meshglobalrate = 1.0/Sys->twoMTot;
  
  if(strcmp(arguments->BiasType, "Harmonic")==0 /*&& arguments->NumPType==0*/){
    Sys->MeshBias.theta0    = arguments->theta0/2;
    Sys->MeshBias.BiasStr   = arguments->BiasStr;
    ComputeBiasPotential    = ComputeBiasPotentialHarmonic;
    MeshBiasChange          = MeshBiasChangeHarmonic;
    MeshBiasGlobalChange    = MeshBiasGlobalHarmonic;
    StoreBiasChange         = StoreBiasChangeHarmonic;
    PrintBiasPotential      = PrintBiasPotentialHarmonic;
    ErrorCheck_Bias         = ErrorCheckHarmonicBias;
    Sys->MeshBias.ShiftSize = sqrt(2.0/Sys->MeshBias.BiasStr);
    
  }
  
  else if(strcmp(arguments->BiasType, "Harmonic")==0 || strcmp(arguments->BiasType, "HarmonicVeff")){
    Sys->MeshBias.theta0  = arguments->theta0/2;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    Sys->MeshBias.VeffStr = 0.0*arguments->Nx*arguments->Ny;
    ComputeBiasPotential  = ComputeBiasPotentialHarmonicPlusVeff;
    MeshBiasChange        = MeshBiasChangeHarmonicPlusVeff;
    MeshBiasGlobalChange  = MeshBiasGlobalHarmonicPlusVeff;
    StoreBiasChange       = StoreBiasChangeHarmonicPlusVeff;
    PrintBiasPotential    = PrintBiasPotentialHarmonicPlusVeff;
    ErrorCheck_Bias       = ErrorCheckHarmonicPlusVeffBias;
    snprintf(arguments->BiasType, 1024, "HarmonicVeff");
    Sys->MeshBias.ShiftSize = sqrt(2.0/Sys->MeshBias.BiasStr);
  }
  
  else if (strcmp(arguments->BiasType, "Point")==0){
    Sys->MeshBias.theta0  = arguments->z0/2;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialMinDist;
    MeshBiasChange        = MeshBiasChangeMinDist;
    StoreBiasChange       = StoreBiasChangeMinDist; 
    Sys->MeshBias.ShiftSize = 0.4;
  }
  
  else if (strcmp(arguments->BiasType, "Linear")==0){
    Sys->MeshBias.theta0  = Sys->MTot*arguments->z0/2;  
    Sys->MeshBias.BiasStr = arguments->BiasStr/Sys->MTot;
    ComputeBiasPotential  = ComputeBiasPotentialAverage;
    MeshBiasChange        = MeshBiasChangeAverage;
    StoreBiasChange       = StoreBiasChangeAverage;
    PrintBiasPotential    = PrintBiasPotentialAverage;
    Sys->MeshBias.ShiftSize = 0.4;
  }
  
  else{
    printf("Error: Could not parse bias type option!\nImplemented are: Point, Linear, and Harmonic\nYou Supplied: %s\n", arguments->BiasType);
    exit(1);
  }
}

void SetChemPotParams(sys *Sys, struct arguments *arguments){
  char Densities[1024], LPType[1024], BindPType[1024], RanPType[1024];
  int pflavors = Sys->NumPTypes = arguments->NumPType, i;
  double alpha = 1e6*arguments->Nx*arguments->Ny/(arguments->Lx*arguments->Ly);
  
  if(pflavors>0){
    strcpy(Densities, arguments->Densities);
    strcpy(BindPType, arguments->BindPType);
    strcpy(RanPType,  arguments->RanPType);
    strcpy(LPType,    arguments->LPType);
    
    Sys->CurParPop[0]     = (double *)malloc((pflavors+1)*sizeof(double));
    Sys->CurParPop[1]     = (double *)malloc((pflavors+1)*sizeof(double));
    
    Sys->LPType[0]     = (double *)malloc(pflavors*sizeof(double));
    Sys->LPType[1]     = (double *)malloc(pflavors*sizeof(double));
    Sys->Fugacities[0] = (double *)malloc(pflavors*sizeof(double));
    Sys->Fugacities[1] = (double *)malloc(pflavors*sizeof(double)); 
    Sys->BindPType     = (double *)malloc(pflavors*sizeof(double));
    Sys->BondNums      = (double *)malloc(pflavors*sizeof(double));
    Sys->RanPType      = (double *)malloc(pflavors*sizeof(double));
  
    ExtractValuesFromStringList( (void *) Sys->BindPType,"DOUBLE", BindPType, pflavors);
    ExtractValuesFromStringList( (void *) Sys->RanPType, "DOUBLE", RanPType,  pflavors);
  
    double *tmp = (double *)malloc(2*pflavors*sizeof(double));
    ExtractValuesFromStringList( (void *) tmp, "DOUBLE", Densities,  2*pflavors);
    
    double zbar;
    zbar=0.0;
    for(i=0;i<pflavors;i++) zbar+=tmp[i];
    zbar = zbar/(alpha-zbar);
    for(i=0;i<pflavors;i++) Sys->Fugacities[0][i] = (1.0+zbar)*tmp[i]/alpha;
    
    zbar=0.0;
    for(i=pflavors;i<2*pflavors;i++) zbar+=tmp[i];
    zbar = zbar/(alpha-zbar);
    for(i=pflavors;i<2*pflavors;i++) Sys->Fugacities[1][i-pflavors] = (1.0+zbar)*tmp[i]/alpha;

//     printf("zbarbot %lf alphabot %lf\n", zbar, alpha);
//     printf("zbartop %lf alphatop %lf\n", zbar, alpha);
//     printf("z1 %lf z2 %lf\n", Sys->Fugacities[0][0], Sys->Fugacities[1][0]);
//     exit(1);
//     for(i=0;i<pflavors;i++){Sys->Fugacities[0][i] = exp(tmp[i]);}
//     for(i=pflavors;i<2*pflavors;i++){Sys->Fugacities[1][i-pflavors] = exp(tmp[i]);} 
    
    ExtractValuesFromStringList( (void *) tmp, "DOUBLE", LPType,  2*pflavors);
    for(i=0;i<pflavors;i++){Sys->LPType[0][i] = tmp[i];}
    for(i=pflavors;i<2*pflavors;i++){Sys->LPType[1][i-pflavors] = tmp[i];} 
  }
}

void SetMonteCarloParams(sys *Sys, struct arguments *arguments){
  Sys->SweepNum       = arguments->MCSweeps;
  Sys->DoF            = Sys->twoMTot;
  Sys->MCSteps        = (long)arguments->MCSweeps*Sys->DoF;
  
  Sys->AcceptanceRate = Sys->NodeAccept = Sys->ShiftAccept = Sys->ShiftAttempt = 0.0;
  Sys->SweepIndex     = 0;
  
  if(arguments->datnums == 0) Sys->datarate = Sys->SweepNum;
  else{Sys->datarate  = max(1, Sys->SweepNum/arguments->datnums);}
  
  if(arguments->framenums == 0) Sys->framerate = Sys->SweepNum;
  else{Sys->framerate = max(1, Sys->SweepNum/arguments->framenums);}
}

void InitCorrFuncs(sys *Sys, struct arguments *arguments){
  corrfuncs *CorrFuncs=&(Sys->CorrFuncs);
  int compnums  = CorrFuncs->compnums  = Sys->NumPTypes + 1;
  int statenums = CorrFuncs->statenums = compnums*compnums;
  CorrFuncs->Cur_MeanProbs = (double *)malloc(statenums*sizeof(double));
  CorrFuncs->Ave_MeanProbs = (double *)malloc(statenums*sizeof(double));
  for(int i=0;i<statenums;i++) CorrFuncs->Ave_MeanProbs[i] = 0.0;
  
  if(Sys->Hhist.binnums>1){
    if( fabs(Sys->dx-Sys->dy)>1e-6 ){
      printf("ERROR: dx(=%lf) and dy(=%lf) are not the same.\nSupport for variable grid widths not implemented.\nPlx set Hist->binnums to 0.\n", Sys->dx, Sys->dy);
      exit(1);
    }
    
    int Nx = Sys->Nx, Ny = Sys->Ny, int_rcutoffsqr = CorrFuncs->int_rcutoffsqr = 20*20;
    int *tmp = (int *)malloc(Sys->Nx*Sys->Ny*sizeof(int));
    int i, j, dummy, lenrvals=0, totdistnums=0, int_rsqr;
    for(i=0;i<Nx*Ny;i++)tmp[i]=ALOT;
    
    //Fetch all possible distances, and count unique ones
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
        int_rsqr = computeindexdistsqr(i, j, Nx, Ny);
        if(int_rsqr <= int_rcutoffsqr){
          for(dummy=0;dummy<totdistnums;dummy++){
            if(tmp[dummy]==int_rsqr){break;}
          }      
          if(dummy==totdistnums) lenrvals++;
          tmp[totdistnums] = int_rsqr;
          totdistnums++;
        }
      }
    }
    
  //Sort entries in ascending order
    dummy=0;
    while(dummy==0){
      dummy=1;
      for(i=0;i<totdistnums;i++){
        if(tmp[i+1]<tmp[i]){
          j=tmp[i+1];tmp[i+1]=tmp[i];tmp[i] = j;
          dummy=0;
        }
      }
    }
    
    //Allocate and initialize arrays
    
    CorrFuncs->lenrvals           =             lenrvals;
    CorrFuncs->int_rsqrlookup     =     (int  *)malloc(lenrvals*sizeof(int));
    CorrFuncs->rnorms             =  (double  *)malloc(lenrvals*sizeof(double));
    CorrFuncs->rvals              =  (double  *)malloc(lenrvals*sizeof(double));
    CorrFuncs->Cur_Correlations   = (double ***)malloc(lenrvals*sizeof(double**));
    CorrFuncs->Ave_Correlations   = (double ***)malloc(lenrvals*sizeof(double**));
     
    for(dummy=0;dummy<lenrvals;dummy++){
      CorrFuncs->int_rsqrlookup[dummy]   = -1;
      CorrFuncs->rnorms[dummy]           = 0.0;
      CorrFuncs->Cur_Correlations[dummy] = (double **)malloc(statenums*sizeof(double *));
      CorrFuncs->Ave_Correlations[dummy] = (double **)malloc(statenums*sizeof(double *));
      for(i=0;i<statenums;i++){
        CorrFuncs->Cur_Correlations[dummy][i] = (double *)malloc(statenums*sizeof(double));
        CorrFuncs->Ave_Correlations[dummy][i] = (double *)malloc(statenums*sizeof(double));
        for(j=0;j<statenums;j++)CorrFuncs->Ave_Correlations[dummy][i][j] = 0.0;
      }
    }
    
    //Assign each distance to int_rsqrlookup and count multiplicity in rnorms
    
    dummy=0;
    for(i=0;i<totdistnums;i++){
      for(j=0;j<dummy;j++){
        if(CorrFuncs->int_rsqrlookup[j] == tmp[i]){
          CorrFuncs->rnorms[j]++;
          break;
        }
      }
      if(j==dummy){
        CorrFuncs->int_rsqrlookup[dummy] = tmp[i];
        CorrFuncs->rnorms[dummy] = 1.0;
        dummy++;
      }
    }
    
    //Assign physical distance from int_rsqrlookup to rvals for plotting purposes
    for(i=0;i<lenrvals;i++) CorrFuncs->rvals[i] = Sys->dx*sqrt((double)CorrFuncs->int_rsqrlookup[i]);
    
    free(tmp);
  }
}

void CreatedzHists(sys *Sys, struct arguments *arguments){
  double hstd=ReturnHeightStd(Sys);
  Sys->Hhist.binnums=arguments->Histbnums;
  if(Sys->Hhist.binnums>1){
    Sys->Hhist.zmin = max(0.0, Sys->MeshBias.theta0-5.0*hstd);
    Sys->Hhist.zmax = Sys->MeshBias.theta0+5.0*hstd;
    Sys->Hhist.dz   = (Sys->Hhist.zmax - Sys->Hhist.zmin)/max(1, Sys->Hhist.binnums-1);
    Sys->Hhist.histnorm = (Sys->Hhist.dz*Sys->MTot);
    
    Sys->Hhist.Hhist[0] = (double *)malloc(sizeof(double)*Sys->Hhist.binnums);
    Sys->Hhist.Hhist[1] = (double *)malloc(sizeof(double)*Sys->Hhist.binnums);
    
    Sys->Hhist.CompHhist[0] = (double *)malloc(sizeof(double)*Sys->Hhist.binnums);
    Sys->Hhist.CompHhist[1] = (double *)malloc(sizeof(double)*Sys->Hhist.binnums);
    
    for(int i=0;i<Sys->Hhist.binnums;i++){
      Sys->Hhist.CompHhist[0][i] = Sys->Hhist.Hhist[0][i]=Sys->Hhist.zmin + i*Sys->Hhist.dz;
      Sys->Hhist.CompHhist[1][i]=0.0;
    }
  }
}

void InitializeStructs(sys *Sys, struct arguments *arguments){
  Sys->rng=gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(Sys->rng, arguments->RandSeed);

  if ( strcmp(arguments->Inputfile, "None")!=0 ){
    if ((Sys->fileptr_input = fopen(arguments->Inputfile, "r"))){}
    else{printf("ERROR:\n Could not find file: check if,\n%s\nis correct", arguments->Inputfile);exit(1);}
    char tmp[1024];strcpy(tmp, arguments->Inputfile);
    ReadInputCommands(Sys->fileptr_input, arguments);strcpy(arguments->Inputfile, tmp);
  }
  
  CreateMeshes(Sys, arguments);
  SetOutPutFiles(Sys, arguments);
  SetChemPotParams(Sys, arguments);
  CreateBiasPotential(Sys, arguments);
  SetMonteCarloParams(Sys, arguments);
  CreatedzHists(Sys, arguments);
  InitCorrFuncs(Sys, arguments);
  
  if( strcmp(arguments->Inputfile, "None")!=0 ){
    ReadConfiguration(Sys->fileptr_input, Sys);
    fclose(Sys->fileptr_input);
  }
  if(strcmp(arguments->ErrorCheck, "Yes")==0){
    ErrorCheckIntermittent = ErrorCheckConfiguration;
    ErrorCheckEvery = DontErrorCheck;
    printf("ErrorCheck Intermittently\n");
  }
  else if(strcmp(arguments->ErrorCheck, "YES")==0){
    ErrorCheckIntermittent = DontErrorCheck;
    ErrorCheckEvery = ErrorCheckConfiguration;
    printf("ErrorCheck Every step\n");
  }
  else {
    ErrorCheckIntermittent = DontErrorCheck;
    ErrorCheckEvery = DontErrorCheck;
  }
}
