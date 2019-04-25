#include "initialize.h"


void CreateMeshCircles(double xtop, double ytop, double rtop, 
                       double xbot, double ybot, double rbot, 
                       double z0, sys *Sys){
  int i,j;
  MeshPoint **TopMesh = Sys->TopMesh;
  MeshPoint **BotMesh = Sys->BotMesh;
  vec r0top, r0bot, distvec;
  double rtopsqr=rtop*rtop, rbotsqr=rbot*rbot, distsqr;
  
  r0top.x=xtop;r0top.y=ytop;
  r0bot.x=xbot;r0bot.y=ybot;
  
  for (i=0;i<Sys->Nx;i++) {
    for (j=0; j<Sys->Ny; j++){
      twodvecdiff(&distvec, TopMesh[i][j].Pos, r0top);
      NearestImageConv(&distvec, Sys->Lx, Sys->Ly);
      distsqr=twod_NormSqr(distvec);
      if( distsqr < rtopsqr )TopMesh[i][j].Pos.z -= sqrt(rtopsqr - distsqr);
      
      twodvecdiff(&distvec, BotMesh[i][j].Pos, r0bot);
      NearestImageConv(&distvec, Sys->Lx, Sys->Ly);
      distsqr=twod_NormSqr(distvec);
      if( distsqr < rbotsqr )BotMesh[i][j].Pos.z += sqrt(rbotsqr - distsqr);
    }
  }
}

void CreateMeshSine(double AmpTop, double PeriodTop, double kxTop, double kyTop,
                    double AmpBot, double PeriodBot, double kxBot, double kyBot,
                    double z0,  sys *Sys){
  int i,j;
  MeshPoint **TopMesh = Sys->TopMesh;
  MeshPoint **BotMesh = Sys->BotMesh;
  
  vec kbarTop, kbarBot;
  
  kbarTop.x = kxTop*(2.0*M_PI*PeriodTop/Sys->Lx);
  kbarTop.y = kyTop*(2.0*M_PI*PeriodTop/Sys->Ly);
  kbarTop.z = 0.0;
  
  
  kbarBot.x = kxBot*(2.0*M_PI*PeriodBot/Sys->Lx);
  kbarBot.y = kyBot*(2.0*M_PI*PeriodBot/Sys->Ly);
  kbarBot.z = 0.0;
  
  
  for (i=0;i<Sys->Nx;i++) {
    for (j=0; j<Sys->Ny; j++){
      TopMesh[i][j].Pos.z = AmpTop*sin( dotprod(TopMesh[i][j].Pos, kbarTop) ) + z0;
      BotMesh[i][j].Pos.z = AmpBot*sin( dotprod(BotMesh[i][j].Pos, kbarBot) );
    }
  }
}

int CountNumberOfCompletedRuns(char *OutputDir, char *tmp){
  int index=0;
  DIR * dirp;
  struct dirent * entry;
  char RestartDir[1024];
  
  strcpy(RestartDir, OutputDir);
  strcat(RestartDir, "/Restarts/");
  
  dirp = opendir(RestartDir);
  while ((entry = readdir(dirp)) != NULL) {
    if (entry->d_type == DT_REG && strstr(entry->d_name, tmp) ){ 
      index++;
    }
  }
  closedir(dirp);
  return index;
}

void SetOutPutFiles(sys *Sys, struct arguments *arguments ){
//Set filepaths
  char tmp[1024];
  sprintf(tmp, "Nodes_%d_NumPTypes_%d_kc_%.2lf_A_%.1lf_theta0_%lf", 
                arguments->Nx*arguments->Ny,
                arguments->NumPType,
                arguments->kc,
                arguments->Lx*arguments->Ly,
                arguments->theta0
         );
  
  int index = CountNumberOfCompletedRuns(arguments->OutputDir, tmp);
  
  snprintf(Sys->Trajfile,    1024, "%s/Trajectories/Traj_%s_Run_%d_.xyz",
          arguments->OutputDir, tmp, index);
  snprintf(Sys->Datafile,    1024, "%s/DataFiles/Data_%s_Run_%d_.dat",
          arguments->OutputDir, tmp, index);
  snprintf(Sys->Restartfile, 1024, "%s/Restarts/Restart_%s_Run_%d_.rest",
          arguments->OutputDir, tmp, index);
  
  Sys->fileptr_traj    = fopen(Sys->Trajfile, "w");
  Sys->fileptr_data    = fopen(Sys->Datafile, "w");
//   Sys->fileptr_restart = fopen(Sys->Restartfile, "w");
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
  int i, j;
  double x;
  
  Sys->Nx      = arguments->Nx;
  Sys->Ny      = arguments->Ny;
  Sys->MTot    = arguments->Nx*arguments->Ny;
  Sys->twoMTot = 2*Sys->MTot;
  
  
  Sys->Lx   = arguments->Lx;
  Sys->Ly   = arguments->Ly;
  Sys->Lxon2= arguments->Lx/2;
  Sys->Lyon2= arguments->Ly/2;
  Sys->A    = arguments->Lx*arguments->Ly;
  
  Sys->dx    = Sys->Lx/Sys->Nx;
  Sys->dy    = Sys->Ly/Sys->Ny;
  Sys->dxsqr = Sys->dx*Sys->dx;
  Sys->dysqr = Sys->dy*Sys->dy;
  
  Sys->dA = Sys->dx*Sys->dy;
  Sys->kc = arguments->kc;
  Sys->MeshEnergyPrefac=(Sys->kc/2)*Sys->dA;
  
  Sys->TopMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  Sys->BotMesh = (MeshPoint **)malloc(Sys->Nx*sizeof(MeshPoint *));
  
  for (i=0;i<Sys->Nx;i++) {
    Sys->TopMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->TopMesh[i][j].Pos.x=x;
      Sys->TopMesh[i][j].Pos.y=Sys->dy*j;
      Sys->TopMesh[i][j].Pos.z=arguments->z0;
    }
  }
  
  for (i=0;i<Sys->Nx;i++) {
    Sys->BotMesh[i]= (MeshPoint *)malloc(Sys->Ny*sizeof(MeshPoint));
    x=i*Sys->dx;
    for (j=0; j<Sys->Ny; j++){
      Sys->BotMesh[i][j].Pos.x=x;
      Sys->BotMesh[i][j].Pos.y=Sys->dy*j;
      Sys->BotMesh[i][j].Pos.z=0.0;
    }
  }
  
  Sys->MRate  = arguments->MRate;
  Sys->MDisp  = arguments->MDisp;
}

void InitializeParticleProperties(sys *Sys, struct arguments *arguments){
  int p;
  int type=0;
  int partialPSum=Sys->CurParPop[0][type];
  double x = 0.0, y = 0.0;
  
  double delta = Sys->DPType[type];
  double deltax = Sys->dx*ceil(delta / Sys->dx);
  double deltay = Sys->dy*ceil(delta / Sys->dy);
  
  for(p=0;p<Sys->maxpop;p++){
    if(p >= Sys->TotCurPop[0] ){
      Sys->Particles[p].Pos.x          = EMPTY;
      Sys->Particles[p].Pos.y          = EMPTY;
      Sys->Particles[p].Pos.z          = EMPTY;
      Sys->Particles[p].type           = EMPTY;
      Sys->Particles[p].MyIntStr       = EMPTY;
      Sys->Particles[p].MyIntCutoff    = EMPTY;
      Sys->Particles[p].MyIntCutoffsqr = EMPTY;
      Sys->Particles[p].MyLength = -40;
      
    }
    
    else{
      Sys->Particles[p].Pos.x = x;
      Sys->Particles[p].Pos.y = y;
      Sys->Particles[p].Pos.z = 0.0;
      
      Sys->Particles[p].nhat.x = 0.0;
      Sys->Particles[p].nhat.y = 0.0;
      Sys->Particles[p].nhat.z = 1.0;
      
      Sys->Particles[p].type           = type;
      Sys->Particles[p].MyIntStr       = Sys->BindPType[type];
      Sys->Particles[p].MyIntCutoff    = Sys->RanPType[type];
      Sys->Particles[p].MyIntCutoffsqr = Sys->RanPType[type]*Sys->RanPType[type];
      
      Sys->Particles[p].MeshId   = BOT;
      Sys->Particles[p].MyDiam   = Sys->DPType[type];
      Sys->Particles[p].MyRad    = Sys->Particles[p].MyDiam/2;
      Sys->Particles[p].MyRadSqr = (Sys->Particles[p].MyDiam/2)*(Sys->Particles[p].MyDiam/2);
      Sys->Particles[p].MyLength = Sys->LPType[type];
      Sys->Particles[p].MeshParticleCutoffsqr =  Sys->Particles[p].MyRad;
      Sys->Particles[p].MeshParticleCutoffsqr *= Sys->Particles[p].MeshParticleCutoffsqr;
            
      Sys->Particles[p].MidPoint.x = Sys->Particles[p].Pos.x;
      Sys->Particles[p].MidPoint.y = Sys->Particles[p].Pos.y;
      Sys->Particles[p].MidPoint.z = Sys->Particles[p].MyLength/2;
      
      Sys->Particles[p].EndPoint.x = Sys->Particles[p].Pos.x;
      Sys->Particles[p].EndPoint.y = Sys->Particles[p].Pos.y;
      Sys->Particles[p].EndPoint.z = Sys->Particles[p].MyLength;
            
      x+=deltax;
      if(x + deltax/2.0 > Sys->Lx){x=0.0;y+=deltay;}
        
      if(p==partialPSum-1){
        type++;
        partialPSum+=Sys->CurParPop[0][type];
        delta = Sys->DPType[type];
        deltax = Sys->dx*ceil(delta / Sys->dx);
        deltay = Sys->dy*ceil(delta / Sys->dy);
      }
    }
    Sys->Particles[p].MyIntPartner = -1;
  }
  
  type=0;
  partialPSum=Sys->CurParPop[1][type];
  x=0.0;y=0.0;
  delta = Sys->DPType[type];
  deltax = Sys->dx*ceil(delta / Sys->dx);
  deltay = Sys->dy*ceil(delta / Sys->dy);
  
  for(p=Sys->maxpop;p<Sys->twomaxpop;p++){
    if(p - Sys->maxpop >= Sys->TotCurPop[1]){
      Sys->Particles[p].Pos.x = EMPTY;
      Sys->Particles[p].Pos.y = EMPTY;
      Sys->Particles[p].Pos.z = EMPTY;
      Sys->Particles[p].type           = EMPTY;
      Sys->Particles[p].MyIntStr       = EMPTY;
      Sys->Particles[p].MyIntCutoff    = EMPTY;
      Sys->Particles[p].MyIntCutoffsqr = EMPTY;
    }
    
    else{
      Sys->Particles[p].Pos.x= x;
      Sys->Particles[p].Pos.y= y;
      Sys->Particles[p].Pos.z= arguments->z0;
      
      Sys->Particles[p].nhat.x= 0.0;
      Sys->Particles[p].nhat.y= 0.0;
      Sys->Particles[p].nhat.z=-1.0;
      
      Sys->Particles[p].type           = type;
      Sys->Particles[p].MyIntStr       = Sys->BindPType[type];
      Sys->Particles[p].MyIntCutoff    = Sys->RanPType[type];
      Sys->Particles[p].MyIntCutoffsqr = Sys->RanPType[type]*Sys->RanPType[type];
      
      Sys->Particles[p].MeshId   = TOP;
      Sys->Particles[p].MyDiam   = delta;
      Sys->Particles[p].MyRad    = Sys->Particles[p].MyDiam/2;
      Sys->Particles[p].MyRadSqr = (Sys->Particles[p].MyDiam/2)*(Sys->Particles[p].MyDiam/2);
      Sys->Particles[p].MyLength = Sys->LPType[type];
      Sys->Particles[p].MeshParticleCutoffsqr =  Sys->Particles[p].MyRad;
      Sys->Particles[p].MeshParticleCutoffsqr *= Sys->Particles[p].MeshParticleCutoffsqr;
      
      Sys->Particles[p].MidPoint.x= Sys->Particles[p].Pos.x;
      Sys->Particles[p].MidPoint.y= Sys->Particles[p].Pos.y;
      Sys->Particles[p].MidPoint.z=-Sys->Particles[p].MyLength/2;
      
      Sys->Particles[p].EndPoint.x= Sys->Particles[p].Pos.x;
      Sys->Particles[p].EndPoint.y= Sys->Particles[p].Pos.y;
      Sys->Particles[p].EndPoint.z=-Sys->Particles[p].MyLength;
      
      x+=deltax;
      if(x + deltax/2.0 > Sys->Lx){x=0.0;y+=deltay;}
        
      if(p - Sys->maxpop==partialPSum - 1){
        type++;
        partialPSum+=Sys->CurParPop[1][type];
        delta = Sys->DPType[type];
        deltax = Sys->dx*ceil(delta / Sys->dx);
        deltay = Sys->dy*ceil(delta / Sys->dy);
      }
    }
    Sys->Particles[p].MyIntPartner = -1;
  }
}

void DetermineCellNeighbours(int ***List, int i, int j, int Nx, int Ny){
  List[i][j][0] = IndexPBC(i,   Nx);
  List[i][j][1] = IndexPBC(j,   Ny);
  List[i][j][2] = IndexPBC(i-1, Nx);
  List[i][j][3] = IndexPBC(j-1, Ny);
  List[i][j][4] = IndexPBC(i,   Nx);
  List[i][j][5] = IndexPBC(j-1, Ny);
  List[i][j][6] = IndexPBC(i+1, Nx);
  List[i][j][7] = IndexPBC(j-1, Ny);
  List[i][j][8] = IndexPBC(i-1, Nx);
  List[i][j][9] = IndexPBC(j,   Ny);
  List[i][j][10]= IndexPBC(i+1, Nx);
  List[i][j][11]= IndexPBC(j,   Ny);
  List[i][j][12]= IndexPBC(i-1, Nx);
  List[i][j][13]= IndexPBC(j+1, Ny);
  List[i][j][14]= IndexPBC(i,   Nx);
  List[i][j][15]= IndexPBC(j+1, Ny);
  List[i][j][16]= IndexPBC(i+1, Nx);
  List[i][j][17]= IndexPBC(j+1, Ny);
}

void CreateCellLists(sys *Sys){
  int i=1,j,k;
  celllist *CellList = &(Sys->CellList);
  
  double maxlength = -1000.0, minrad = 10000;
  for(k=0;k<Sys->NumPTypes;k++){
    maxlength=max(maxlength, 0.5*Sys->LPType[k]);
    minrad=min(minrad, Sys->DPType[k]/2);
  }
  
  double cutoff = max(maxlength, 3.0*Sys->dx);
  
  while (Sys->Lx/i > cutoff)i++;
  CellList->NumCellsx = max(1, i-1);i=1;
  
  while (Sys->Ly/i > cutoff)i++;
  CellList->NumCellsy = max(i, i-1);
  
  if(fabs(Sys->Lx -Sys->Ly)<1e-5)CellList->NumCellsy = CellList->NumCellsx;
  
  CellList->dCx = Sys->Lx/CellList->NumCellsx;
  CellList->dCy = Sys->Ly/CellList->NumCellsy;
  
  CellList->MaxPopCell = max(10, ceil(1.5*0.9068*CellList->dCx*CellList->dCy/(M_PI*minrad*minrad)));
  
  CellList->Neighbors      = (int ***)malloc(CellList->NumCellsx*sizeof(int **));
  
  CellList->TopPopulations = (int ** )malloc(CellList->NumCellsx*sizeof(int  *));
  CellList->TopParticleIds = (int ***)malloc(CellList->NumCellsx*sizeof(int **));
  
  CellList->BotPopulations = (int ** )malloc(CellList->NumCellsx*sizeof(int  *));
  CellList->BotParticleIds = (int ***)malloc(CellList->NumCellsx*sizeof(int **));
  
  for(i=0;i<CellList->NumCellsx;i++){
    
    CellList->Neighbors[i]   = (int **)malloc(CellList->NumCellsy*sizeof(int *));
    
    CellList->TopPopulations[i] = (int * )malloc(CellList->NumCellsy*sizeof(int ));
    CellList->TopParticleIds[i] = (int **)malloc(CellList->NumCellsy*sizeof(int *));
    
    CellList->BotPopulations[i] = (int * )malloc(CellList->NumCellsy*sizeof(int ));
    CellList->BotParticleIds[i] = (int **)malloc(CellList->NumCellsy*sizeof(int *));
    
    for(j = 0;j<CellList->NumCellsy;j++){
      
      CellList->Neighbors[i][j] = (int *)malloc(18*sizeof(int));
    
      CellList->TopPopulations[i][j] = 0;
      CellList->TopParticleIds[i][j] = (int *)malloc(CellList->MaxPopCell*sizeof(int));
      
      CellList->BotPopulations[i][j] = 0;
      CellList->BotParticleIds[i][j] = (int *)malloc(CellList->MaxPopCell*sizeof(int));
      
      for(k = 0;k<CellList->MaxPopCell;k++){
        CellList->TopParticleIds[i][j][k] = EMPTY;
        CellList->BotParticleIds[i][j][k] = EMPTY;
      }
      
      DetermineCellNeighbours(CellList->Neighbors, i, j, CellList->NumCellsx, CellList->NumCellsy);
    }
  }
}

void CreateParticles(sys *Sys, struct arguments *arguments){
  int i, pflavors = arguments->NumPType;
  char PTypePop[1024], LPType[1024], DPType[1024], BindPType[1024], RanPType[1024];
  double minrad, maxlength;
  
  Sys->NumPTypes = pflavors;
  Sys->vmdDoF = 0;
  if (arguments->NumPType>0){
    Sys->PTypePop  = (int *)   malloc(pflavors*sizeof(int));
    Sys->LPType    = (double *)malloc(pflavors*sizeof(double));
    Sys->DPType    = (double *)malloc(pflavors*sizeof(double));
    Sys->BindPType = (double *)malloc(pflavors*sizeof(double));
    Sys->RanPType  = (double *)malloc(pflavors*sizeof(double));
    Sys->AccessVol = (double *)malloc(pflavors*sizeof(double));
    
    Sys->CurParPop[0] = (int *)malloc(pflavors*sizeof(int));
    Sys->CurParPop[1] = (int *)malloc(pflavors*sizeof(int));
    
    Sys->vmdLength = (int *)   malloc(pflavors*sizeof(int));
    
    strcpy(PTypePop,  arguments->PTypePop);
    strcpy(LPType,    arguments->LPType);
    strcpy(DPType,    arguments->DPType);
    strcpy(BindPType, arguments->BindPType);
    strcpy(RanPType,  arguments->RanPType);
    
    ExtractValuesFromStringList( (void *) Sys->PTypePop, "INT",    PTypePop,  pflavors);
    ExtractValuesFromStringList( (void *) Sys->LPType,   "DOUBLE", LPType,    pflavors);
    ExtractValuesFromStringList( (void *) Sys->DPType,   "DOUBLE", DPType,    pflavors);
    ExtractValuesFromStringList( (void *) Sys->BindPType,"DOUBLE", BindPType, pflavors);
    ExtractValuesFromStringList( (void *) Sys->RanPType, "DOUBLE", RanPType,  pflavors);
  
    Sys->TotCurPop[0] = Sys->TotCurPop[1] = 0;
    
    minrad = maxlength = 1000.0;
    for(i=0;i<pflavors;i++){
      minrad    = min(minrad, Sys->DPType[i]/2);
      maxlength = min(maxlength, Sys->LPType[i]);
    }
    Sys->maxpop = min(MAXPOP, (int) round(0.9068*Sys->A/(M_PI*minrad*minrad)));
    Sys->twomaxpop = 2*Sys->maxpop;

    for(i=0;i<pflavors;i++){
      Sys->vmdLength[i]= (int)round(Sys->LPType[i]/Sys->DPType[i]) + 1;
      Sys->vmdDoF+= 2*Sys->maxpop*Sys->vmdLength[i];
      Sys->CurParPop[0][i]=Sys->PTypePop[i];
      Sys->CurParPop[1][i]=Sys->PTypePop[i];
      Sys->TotCurPop[0] += Sys->PTypePop[i];
      Sys->TotCurPop[1] += Sys->PTypePop[i];
      Sys->AccessVol[i] = Sys->A/(M_PI*(Sys->DPType[i]/2)*(Sys->DPType[i]/2));
    }

    Sys->Particles = (particle *)malloc( Sys->twomaxpop * sizeof(particle) );
    
    InitializeParticleProperties(Sys, arguments);
    
    Sys->PRate=arguments->PRate;
    Sys->PDisp=arguments->PDisp;
  }
  
  else{
    Sys->vmdDoF=0;
    Sys->maxpop = 0;Sys->twomaxpop = 0;
    Sys->PRate=0;Sys->MRate=1.0;Sys->GCRate = 0.0;
  }
}

void CreateBiasPotential(sys *Sys, struct arguments *arguments){
  if (arguments->BiasStr<1e-6) Sys->meshglobalrate = 0.0;
  else Sys->meshglobalrate = 1.0/Sys->twoMTot;
  
  if(strcmp(arguments->BiasType, "Harmonic")==0){
    Sys->MeshBias.theta0  = arguments->theta0;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialHarmonic;
    MeshBiasChange        = MeshBiasChangeHarmonic;
    MeshBiasGlobalChange  = MeshBiasGlobalHarmonic;
    StoreBiasChange       = StoreBiasChangeHarmonic;
    PrintBiasPotential    = PrintBiasPotentialHarmonic;
  }
  
  else if (strcmp(arguments->BiasType, "Point")==0){
    Sys->MeshBias.theta0  = arguments->z0;
    Sys->MeshBias.BiasStr = arguments->BiasStr;
    ComputeBiasPotential  = ComputeBiasPotentialMinDist;
    MeshBiasChange        = MeshBiasChangeMinDist;
    StoreBiasChange       = StoreBiasChangeMinDist;  
  }
  
  else if (strcmp(arguments->BiasType, "Linear")==0){
    Sys->MeshBias.theta0  = Sys->MTot*arguments->z0;  
    Sys->MeshBias.BiasStr = arguments->BiasStr/Sys->MTot;
    ComputeBiasPotential  = ComputeBiasPotentialAverage;
    MeshBiasChange        = MeshBiasChangeAverage;
    StoreBiasChange       = StoreBiasChangeAverage;
    PrintBiasPotential    = PrintBiasPotentialAverage;
  }
  
  else{
    printf("Error: Could not parse bias type option!\nImplemented are: Point, Linear, and Harmonic\nYou Supplied: %s\n", arguments->BiasType);
    exit(1);
  }
}

void SetChemPotParams(sys *Sys, struct arguments *arguments){
  char ChemPots[1024];
  int pflavors = Sys->NumPTypes, i;
  
  Sys->GCRate = arguments->GCRate;
  
  if(pflavors>0){
    Sys->ChemPots[0] = (double *)malloc(pflavors*sizeof(double));
    Sys->ChemPots[1] = (double *)malloc(pflavors*sizeof(double)); 
    strcpy(ChemPots,  arguments->ChemPots);
    
    double *tmp = (double *)malloc(2*pflavors*sizeof(double));
    ExtractValuesFromStringList( (void *) tmp, "DOUBLE", ChemPots,  2*pflavors);
    for(i=0;i<pflavors;i++){Sys->ChemPots[0][i] = tmp[i];}
    for(i=pflavors;i<2*pflavors;i++){Sys->ChemPots[1][i-pflavors] = tmp[i];} 
  }
}

void SetMonteCarloParams(sys *Sys, struct arguments *arguments){
  
  Sys->SweepNum       = arguments->MCSweeps;
  Sys->DoF            = Sys->twoMTot + Sys->TotCurPop[0] + Sys->TotCurPop[1];
  Sys->vmdDoF        += Sys->twoMTot;
  Sys->MCSteps        = (long)arguments->MCSweeps*Sys->DoF;
  
  Sys->AcceptanceRate = Sys->AcceptMesh = Sys->AcceptPart = Sys->AcceptGC = 0.0;
  Sys->MeshAttempt    = Sys->ParticleAttempt = Sys->GCAttempt = 0 ;
  Sys->SweepIndex     = 0;
  
  if(arguments->datnums == 0) Sys->datarate = Sys->SweepNum;
  else{Sys->datarate  = max(1, Sys->SweepNum/arguments->datnums);}
  
  if(arguments->framenums == 0) Sys->framerate = Sys->SweepNum;
  else{Sys->framerate = max(1, Sys->SweepNum/arguments->framenums);}
}

void InitializeStructs(sys *Sys, struct arguments *arguments){
  Sys->rng=gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(Sys->rng, time(NULL));
  
  if ( strcmp(arguments->Inputfile, "None")!=0 ){
    if ((Sys->fileptr_input = fopen(arguments->Inputfile, "r"))){}
    else{printf("ERROR:\n Could not find file: check if,\n%s\nis correct", arguments->Inputfile);exit(1);}
    char tmp[1024];strcpy(tmp, arguments->Inputfile);
    ReadInputCommands(Sys->fileptr_input, arguments);strcpy(arguments->Inputfile, tmp);
  }
  
  CreateMeshes(Sys, arguments);
  CreateParticles(Sys, arguments);
  CreateCellLists(Sys);
  SetOutPutFiles(Sys, arguments);
  CreateBiasPotential(Sys, arguments);
  SetChemPotParams(Sys, arguments);
  SetMonteCarloParams(Sys, arguments);
//   CreateMeshCircles(7.5, 7.5, 10.0, 3.0, 3.0, 0.0, Sys);
//   CreateMeshSine(3.0, 2.0, 1.0, 0.0, 2.0, 4.0, sqrt(2.0)/2.0, sqrt(2.0)/2.0, arguments->z0, Sys);
//   PrintConfiguration(Sys->fileptr_traj, Sys, 1);
//   exit(1);
  if ( strcmp(arguments->Inputfile, "None")!=0 ){
    ReadConfiguration(Sys->fileptr_input, Sys);
    fclose(Sys->fileptr_input);
  }
  
  if(strcmp(arguments->ErrorCheck, "Yes")==0){ErrorCheck = ErrorCheckConfiguration;printf("ErrorChecking");}
  else ErrorCheck = DontErrorCheck;
}

