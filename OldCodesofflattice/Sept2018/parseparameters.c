#include "parseparameters.h"
#include "membranesys.h"
#include "utilities.h"

// 
const char *__restrict__ ArgumentPrint(){
  return "\nMC properties:\n"
  "   MC Sweeps:       %d\n"
  "   Mesh rate:       %lf\n"
  "   Particle rate:   %lf\n"
  "   GC rate:         %lf\n"
  "   Mesh disp:       %lf\n"
  "   Particle disp:   %lf\n"
  "   Bias type    :   %s\n"
  "   Bias strength:   %lf\n"
  "   Bias location:   %lf\n"
  "   InitialHeight:   %lf\n"
  "\nMesh Properties:\n"
  "   Bending Modulus: %lf\n"
  "   Lx:              %lf\n"
  "   Ly:              %lf\n"
  "   Node nums in x:  %d\n"
  "   Node nums in y:  %d\n"
  "\nParticle Propterties:\n"
  "   Ptype num:       %d\n"
  "   Population nums: %s\n"
  "   Lengths:         %s\n"
  "   Diameters:       %s\n"
  "   Binding str:     %s\n"
  "   Binding ranges:  %s\n"
  "   ChemPot vals:    %s\n"
  "\nOutput Options:\n"
  "   Tag:             %s\n"
  "   OutputDir:       %s\n"
  "   Framenums:       %d\n"
  "   Datanums:        %d\n"
  "\nRestart Info:\n"
  "   Inputfilename:   %s\n"
  "   Equil option:    %s\n"
  "\nError Check:        %s\n"
  "\n";
}

//Replace all occurances of %s above with [^n]
const char *__restrict__ ArgumentRead(){
  return "\nMC properties:\n"
  "   MC Sweeps:       %d\n"
  "   Mesh rate:       %lf\n"
  "   Particle rate:   %lf\n"
  "   GC rate:         %lf\n"
  "   Mesh disp:       %lf\n"
  "   Particle disp:   %lf\n"
  "   Bias type    :   %[^\n]\n"
  "   Bias strength:   %lf\n"
  "   Bias location:   %lf\n"
  "   InitialHeight:   %lf\n"
  "\nMesh Properties:\n"
  "   Bending Modulus: %lf\n"
  "   Lx:              %lf\n"
  "   Ly:              %lf\n"
  "   Node nums in x:  %d\n"
  "   Node nums in y:  %d\n"
  "\nParticle Propterties:\n"
  "   Ptype num:       %d\n"
  "   Population nums: %[^\n]"
  "   Lengths:         %[^\n]"
  "   Diameters:       %[^\n]"
  "   Binding str:     %[^\n]"
  "   Binding ranges:  %[^\n]"
  "   ChemPot vals:    %[^\n]"
  "\nOutput Options:\n"
  "   Tag:             %[^\n]"
  "   OutputDir:       %[^\n]"
  "   Framenums:       %d\n"
  "   Datanums:        %d\n"
  "\nRestart Info:\n"
  "   Inputfilename:   %[^\n]"
  "   Equil option:    %[^\n]"
  "\nError Check:        %[^\n]\n"
  "\n";
}

/* Our argp parser. */
static error_t parse_opt (int key, char *arg, struct argp_state *state){
/* Parse each option, one at a time. */
/* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;
  
  switch (key)
    {
    
    case 'a':sscanf(arg, "%lf", &arguments->Lx);break;
    case 'b':sscanf(arg, "%lf", &arguments->Ly);break;
    case 'c':sscanf(arg, "%d",  &arguments->Nx);break;
    case 'd':sscanf(arg, "%d",  &arguments->Ny);break;
    case 'e':sscanf(arg, "%lf", &arguments->z0);break;
    case 'f':sscanf(arg, "%lf", &arguments->MRate);break;
    case 'g':sscanf(arg, "%lf", &arguments->MDisp);break;
    case 'h':sscanf(arg, "%lf", &arguments->kc);break;
    
    case 'i':sscanf(arg, "%d",  &arguments->NumPType);break;
    case 'j':strcpy(arguments->PTypePop,  arg);break;
    case 'k':strcpy(arguments->LPType,    arg);break;
    case 'l':strcpy(arguments->DPType,    arg);break;
    case 'm':strcpy(arguments->BindPType, arg);break;
    case 'n':strcpy(arguments->RanPType,  arg);break;
    case 'o':sscanf(arg, "%lf", &arguments->PRate);break;
    case 'p':sscanf(arg, "%lf", &arguments->PDisp);break;  
    
    case 'q':strcpy(arguments->ChemPots,  arg);break;
    case 'r':sscanf(arg, "%lf", &arguments->GCRate);break;
        
    case 's':strcpy(arguments->OutputDir, arg);break;
    case '3':strcpy(arguments->Tag,  arg);break;
    case 't':sscanf(arg, "%d",  &arguments->framenums);break;
    case 'u':sscanf(arg, "%d",  &arguments->datnums);break;
    
    case 'v':strcpy(arguments->Inputfile, arg);break;
    case 'w':strcpy(arguments->Equil,     arg);break;
    
    case 'x':sscanf(arg, "%d",  &arguments->MCSweeps);break;
    case 'y':sscanf(arg, "%lf", &arguments->BiasStr);break;
    
    case 'z':strcpy(arguments->BiasType, arg);break;
    
    case '1':sscanf(arg, "%lf", &arguments->theta0);break;
    
    case '2':strcpy(arguments->ErrorCheck,  arg);break;
    
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, "", "" };

void SanityCheckArgs(struct arguments *arguments){
  if(arguments->MRate + arguments->PRate + arguments->GCRate!=1.0){
   printf("Sum of MRate, PRate, and GCRate is not unity,\n");
   printf("MRate: %lf PRate: %lf GCRate: %lf\n", arguments->MRate, arguments->PRate, arguments->GCRate);
   exit(1);
//    arguments->PRate=1.0-arguments->MRate;
//    printf("Resetting to:\nMRate: %lf PRate: %lf\n", arguments->MRate, 1.0 - arguments->MRate);
  }
  if(arguments->NumPType==0 ){arguments->MRate=1.0;arguments->PRate=0.0;}
}

void InterpretInputCommands( struct arguments *arguments, int argc, char **argv){
  
//Default values
  
//Mesh properties
  arguments->Lx=300.0;
  arguments->Ly=300.0;
  arguments->Nx=50;
  arguments->Ny=50;
  arguments->kc=20;
  
//Particle properties
  
  arguments->NumPType=0;

//   arguments->NumPType=1;
//   snprintf(arguments->PTypePop,  1024, "1");
//   snprintf(arguments->LPType,    1024, "0.0");
//   snprintf(arguments->DPType,    1024, "10.0");
//   snprintf(arguments->BindPType, 1024, "0.0");
//   snprintf(arguments->RanPType,  1024, "0.0");
//   snprintf(arguments->ChemPots,  1024, "-2.0 -2.0");

//   arguments->NumPType=2;
//   snprintf(arguments->PTypePop,  1024, "0 0");
//   snprintf(arguments->LPType,    1024, "12.0 8.0");
//   snprintf(arguments->DPType,    1024, "4.0 4.0");
//   snprintf(arguments->BindPType, 1024, "0.0 1.0");
//   snprintf(arguments->RanPType,  1024, "0.0 4.0");
//   snprintf(arguments->ChemPots,  1024, "-1.0 -2.0 -1.0 -2.0");

//MC properties
  arguments->MRate=0.98;
  arguments->PRate=0.02;
  arguments->GCRate=0.0;
  
  arguments->MDisp=1.0;
  arguments->PDisp=100.0; 
  
  arguments->z0=50.0;
  arguments->BiasStr=2.0;
  arguments->theta0=50.0;
  snprintf(arguments->BiasType, 1024, "Harmonic");
  arguments->MCSweeps=10000;
  
// Output options  
  arguments->framenums=500;
  arguments->datnums=50000;
  sscanf("/Data/FletcherMembrane/", "%s", arguments->OutputDir);
//   arguments->Tag[0]='\0';
  strncpy(arguments->Tag, "Test", 1024);
  
// Input options
  strncpy(arguments->Inputfile, "None", 1024);
  snprintf(arguments->Equil,  3, "NO ");
  
// ErrorChecking argument
  strncpy(arguments->ErrorCheck, "Naw", 1024);
    
//Parse Arguments
  argp_parse (&argp, argc, argv, 0, 0, arguments);
  SanityCheckArgs(arguments);
}

void PrintInputCommands(FILE *fp, struct arguments *arguments){

  fprintf(fp, ArgumentPrint(),
              arguments->MCSweeps,
              arguments->MRate, arguments->PRate, arguments->GCRate,
              arguments->MDisp, arguments->PDisp,
              arguments->BiasType, arguments->BiasStr, arguments->theta0,
              arguments->z0,
              arguments->kc, arguments->Lx, arguments->Ly,
              arguments->Nx, arguments->Ny,
              arguments->NumPType,
              arguments->PTypePop, arguments->LPType, arguments->DPType, arguments->BindPType, arguments->RanPType,
              arguments->ChemPots,
              arguments->Tag, arguments->OutputDir,
              arguments->framenums, arguments->datnums,
              arguments->Inputfile, arguments->Equil,
              arguments->ErrorCheck
         );
}

void ConfigurationRestart(FILE *fp, sys *Sys){
  int i, j, mytype;
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf\n", 
      Sys->TopMesh[i][j].Pos.x,  Sys->TopMesh[i][j].Pos.y,  Sys->TopMesh[i][j].Pos.z,
      Sys->TopMesh[i][j].nhat.x, Sys->TopMesh[i][j].nhat.y, Sys->TopMesh[i][j].nhat.z
    );
   }
  }
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    fprintf(
      fp, "N %lf %lf %lf %lf %lf %lf\n", 
      Sys->BotMesh[i][j].Pos.x,  Sys->BotMesh[i][j].Pos.y,  Sys->BotMesh[i][j].Pos.z,
      Sys->BotMesh[i][j].nhat.x, Sys->BotMesh[i][j].nhat.y, Sys->BotMesh[i][j].nhat.z
    );
   }
  }
  
  fprintf(fp, "BotConfig:\n");
  for(mytype = 0;mytype<Sys->NumPTypes;mytype++){
    for(j = 0;j<Sys->TotCurPop[0];j++){
      if(Sys->Particles[j].type == mytype ){
        fprintf(
          fp, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
          Sys->Particles[j].type,     Sys->Particles[j].MeshId,
          Sys->Particles[j].Pos.x,    Sys->Particles[j].Pos.y,    Sys->Particles[j].Pos.z,  
          Sys->Particles[j].nhat.x,   Sys->Particles[j].nhat.y,   Sys->Particles[j].nhat.z,
          Sys->Particles[j].MyDiam,   Sys->Particles[j].MyRad,    Sys->Particles[j].MyRadSqr,
          Sys->Particles[j].MyLength, Sys->Particles[j].MyIntStr, Sys->Particles[j].MyIntCutoff,
          Sys->Particles[j].MyIntCutoffsqr, (Sys->BondList.BondableTypes>0 ? Sys->BondList.ParticleBonds[j].currbond:-1)
        );
      }   
    }
  }
  
  fprintf(fp, "TopConfig:\n");
  for(mytype = 0;mytype<Sys->NumPTypes;mytype++){
    for(j = Sys->maxpop;j<Sys->TotCurPop[1] + Sys->maxpop;j++){
      if(Sys->Particles[j].type == mytype ){
        fprintf(
          fp, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
          Sys->Particles[j].type,     Sys->Particles[j].MeshId,
          Sys->Particles[j].Pos.x,    Sys->Particles[j].Pos.y,    Sys->Particles[j].Pos.z,  
          Sys->Particles[j].nhat.x,   Sys->Particles[j].nhat.y,   Sys->Particles[j].nhat.z,
          Sys->Particles[j].MyDiam,   Sys->Particles[j].MyRad,    Sys->Particles[j].MyRadSqr,
          Sys->Particles[j].MyLength, Sys->Particles[j].MyIntStr, Sys->Particles[j].MyIntCutoff,
          Sys->Particles[j].MyIntCutoffsqr, (Sys->BondList.BondableTypes>0 ? Sys->BondList.ParticleBonds[j].currbond:-1)
        );  
      } 
    }
  }
}

void WriteRestart(sys *Sys, struct arguments *arguments){
  Sys->fileptr_restart = fopen(Sys->Restartfile, "w");
  PrintInputCommands(Sys->fileptr_restart, arguments);
  ConfigurationRestart(Sys->fileptr_restart, Sys);
  fclose(Sys->fileptr_restart);
}

void ReadConfiguration(FILE *fp, sys *Sys){
  char line[1024], *dummy;
  int i,j;
  double length;
  
  dummy = fgets(line, 1024, fp);
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    sscanf(
      line, "N %lf %lf %lf %lf %lf %lf\n", 
      &(Sys->TopMesh[i][j].Pos.x),  &(Sys->TopMesh[i][j].Pos.y),  &(Sys->TopMesh[i][j].Pos.z),
      &(Sys->TopMesh[i][j].nhat.x), &(Sys->TopMesh[i][j].nhat.y), &(Sys->TopMesh[i][j].nhat.z)
    );
    dummy = fgets(line, 1024, fp);
   }
  }
  
  for (i=0;i<Sys->Nx;i++){
   for (j=0;j<Sys->Ny;j++){
    sscanf(
      line, "N %lf %lf %lf %lf %lf %lf\n", 
      &(Sys->BotMesh[i][j].Pos.x),  &(Sys->BotMesh[i][j].Pos.y),  &(Sys->BotMesh[i][j].Pos.z),
      &(Sys->BotMesh[i][j].nhat.x), &(Sys->BotMesh[i][j].nhat.y), &(Sys->BotMesh[i][j].nhat.z)
    );
    dummy = fgets(line, 1024, fp);
   }
  }
  
  Sys->TotCurPop[0] = Sys->TotCurPop[1] = 0;
  for(i=0;i<Sys->NumPTypes;i++)Sys->CurParPop[0][i] = Sys->CurParPop[1][i] = 0;
  
  if ( strcmp(line, "BotConfig:\n")!=0 ){printf("Error: Mismatch in inputfile:\n%s\n", line);exit(1);}
  
  dummy = fgets(line, 1024, fp);
  for (i=0;( i<Sys->maxpop && (strcmp(line, "TopConfig:\n") != 0) );i++){
    sscanf(line, 
      "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
      &(Sys->Particles[i].type),     &(Sys->Particles[i].MeshId),
      &(Sys->Particles[i].Pos.x),    &(Sys->Particles[i].Pos.y),    &(Sys->Particles[i].Pos.z),  
      &(Sys->Particles[i].nhat.x),   &(Sys->Particles[i].nhat.y),   &(Sys->Particles[i].nhat.z),
      &(Sys->Particles[i].MyDiam),   &(Sys->Particles[i].MyRad),    &(Sys->Particles[i].MyRadSqr),
      &(Sys->Particles[i].MyLength), &(Sys->Particles[i].MyIntStr), &(Sys->Particles[i].MyIntCutoff),
      &(Sys->Particles[i].MyIntCutoffsqr), &(Sys->BondList.ParticleBonds[i].currbond) 
    );

    length = Sys->Particles[i].MyLength;Sys->Particles[i].MeshParticleCutoffsqr = Sys->Particles[i].MyRadSqr;
      
    Sys->Particles[i].MidPoint.x = Sys->Particles[i].Pos.x + (length/2)*Sys->Particles[i].nhat.x;
    Sys->Particles[i].MidPoint.y = Sys->Particles[i].Pos.y + (length/2)*Sys->Particles[i].nhat.y;
    Sys->Particles[i].MidPoint.z = Sys->Particles[i].Pos.z + (length/2)*Sys->Particles[i].nhat.z;
    PosPBC(&(Sys->Particles[i].MidPoint), Sys->Lx, Sys->Ly);
    
    Sys->Particles[i].EndPoint.x = Sys->Particles[i].Pos.x + length*Sys->Particles[i].nhat.x;
    Sys->Particles[i].EndPoint.y = Sys->Particles[i].Pos.y + length*Sys->Particles[i].nhat.y;
    Sys->Particles[i].EndPoint.z = Sys->Particles[i].Pos.z + length*Sys->Particles[i].nhat.z; 
    PosPBC(&(Sys->Particles[i].EndPoint), Sys->Lx, Sys->Ly);
    
    Sys->Particles[i].MeshParticleEnergy = 0.0;
    Sys->TotCurPop[0]++;
    Sys->CurParPop[0][Sys->Particles[i].type]++;
    dummy = fgets(line, 1024, fp);
  }
  dummy=fgets(line, 1024, fp);
  for (i=Sys->maxpop;(i<Sys->twomaxpop && (dummy != NULL));i++){
    sscanf(line,
      "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", 
      &(Sys->Particles[i].type),     &(Sys->Particles[i].MeshId),
      &(Sys->Particles[i].Pos.x),    &(Sys->Particles[i].Pos.y),    &(Sys->Particles[i].Pos.z),  
      &(Sys->Particles[i].nhat.x),   &(Sys->Particles[i].nhat.y),   &(Sys->Particles[i].nhat.z),
      &(Sys->Particles[i].MyDiam),   &(Sys->Particles[i].MyRad),    &(Sys->Particles[i].MyRadSqr),
      &(Sys->Particles[i].MyLength), &(Sys->Particles[i].MyIntStr), &(Sys->Particles[i].MyIntCutoff),
      &(Sys->Particles[i].MyIntCutoffsqr), &(Sys->BondList.ParticleBonds[i].currbond) 
    );

    length = Sys->Particles[i].MyLength;Sys->Particles[i].MeshParticleCutoffsqr = Sys->Particles[i].MyRadSqr;
      
    Sys->Particles[i].MidPoint.x = Sys->Particles[i].Pos.x + (length/2)*Sys->Particles[i].nhat.x;
    Sys->Particles[i].MidPoint.y = Sys->Particles[i].Pos.y + (length/2)*Sys->Particles[i].nhat.y;
    Sys->Particles[i].MidPoint.z = Sys->Particles[i].Pos.z + (length/2)*Sys->Particles[i].nhat.z;
    PosPBC(&(Sys->Particles[i].MidPoint), Sys->Lx, Sys->Ly);
    
    Sys->Particles[i].EndPoint.x = Sys->Particles[i].Pos.x + length*Sys->Particles[i].nhat.x;
    Sys->Particles[i].EndPoint.y = Sys->Particles[i].Pos.y + length*Sys->Particles[i].nhat.y;
    Sys->Particles[i].EndPoint.z = Sys->Particles[i].Pos.z + length*Sys->Particles[i].nhat.z; 
    PosPBC(&(Sys->Particles[i].EndPoint), Sys->Lx, Sys->Ly);
    
    Sys->Particles[i].MeshParticleEnergy = 0.0;
    Sys->TotCurPop[1]++;
    Sys->CurParPop[1][Sys->Particles[i].type]++;
    dummy = fgets(line, 1024, fp);
  }
  dummy=dummy;
  fclose(fp);
}

void ReadInputCommands(FILE *fp, struct arguments *arguments){
  
  int tmp=fscanf(fp, ArgumentRead(), 
                     &arguments->MCSweeps,
                     &arguments->MRate, &arguments->PRate, &arguments->GCRate,
                     &arguments->MDisp, &arguments->PDisp,
                     &arguments->BiasType, &arguments->BiasStr, &arguments->theta0,
                     &arguments->z0,
                     &arguments->kc, &arguments->Lx, &arguments->Ly,
                     &arguments->Nx, &arguments->Ny,
                     &arguments->NumPType,
                     &arguments->PTypePop, &arguments->LPType, &arguments->DPType, &arguments->BindPType, &arguments->RanPType,
                     &arguments->ChemPots,
                     &arguments->Tag, &arguments->OutputDir,
                     &arguments->framenums, &arguments->datnums,
                     &arguments->Inputfile, &arguments->Equil,
                     &arguments->ErrorCheck
                );
  tmp=tmp;
}
