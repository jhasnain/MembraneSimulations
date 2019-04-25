#ifndef MEMSYS_H
#define MEMSYS_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

#define TOP -1
#define BOT 1
#define EMPTY -100
#define ALOT 100000000.0
#define MAXPOP 1000
#define EPSILON 0.000001

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

typedef struct vec{
  double x, y, z;
} vec;

typedef struct particle{
  vec Pos, EndPoint, MidPoint;
  vec nhat;
  double MeshParticleEnergy, PhaseFactorEnergy;
  
  int type, MeshId;
  double MyDiam, MyRad, MyRadSqr, MyLength;
  double MyIntStr, MyIntCutoff, MyIntCutoffsqr, MeshParticleCutoffsqr;
  
  int MyIntPartner;
} particle;

typedef struct Meshpoint{
  vec Pos;
  vec nhat;
  double nNorm;
  double lc_x, lc_y;
  double CurvatureEnergy;
} MeshPoint;

typedef struct meshbias{
  double BiasStr, theta0;

  //Current state
  double theta;
  int imin, jmin;
  
  //Suggested state
  int idummy, jdummy;
  double thetadummy;
} meshbias;  

typedef struct celllist{
  int MaxPopCell;
  int NumCellsx, NumCellsy;
  double dCx, dCy;
  
  int ***Neighbors;

  int  **BotPopulations;
  int ***BotParticleIds;
  
  int  **TopPopulations;
  int ***TopParticleIds;
  
} celllist;

typedef struct sys{
//MC Ingredients
  long MCSteps;
  long DoF, SweepIndex, SweepNum;
  double PDisp, PRate; 
  double MRate, MDisp;
  double GCRate;
  gsl_rng *rng;
  
//Mesh properties
  int    Nx, Ny, MTot, twoMTot;
  double Lx, Ly, Lxon2, Lyon2, A;
  double kc, MeshEnergyPrefac;
  double dx, dy, dA;
  double dxsqr, dysqr;
  
  MeshPoint **TopMesh;
  MeshPoint **BotMesh;
  
  //Mesh MC move vars
    int i0, j0, choice;
    double suggestedz, meshglobalrate;
    //Local derivatives dummy vars
      vec lcCent, lcRight, lcLeft, lcIn, lcOut;
      vec nhatRight, nhatLeft, nhatIn, nhatOut;
      double nNormRight, nNormLeft, nNormIn, nNormOut;
      
    //Variables to track particle changes due to mesh move
      int AffectedParticleNum, Unique_dummy_index[12];
      particle mesh_dummy_particles[12];
      
  //MeshBiasParams
    meshbias MeshBias;
    
//Particle properties
  int NumPTypes, *PTypePop;
  double *LPType, *DPType, *BindPType, *RanPType, *AccessVol;
  particle *Particles;

  celllist CellList;
  
  //Particle MC move vars
    int p0;
    particle dummy_particle;
    
  //Grand Canonical MC properties
  int TotCurPop[2], *CurParPop[2], maxpop, twomaxpop;
  double *ChemPots[2];
    
  //Grand move vars
  int GC_AddRemove, GC_TopOrBot, GC_PType;
    
//Measureables
  double AcceptanceRate, AcceptMesh, AcceptPart, AcceptGC;
  long MeshAttempt, ParticleAttempt, GCAttempt;
  double MeshEnergy, BiasEnergy, TotMeshParticleEnergy, TotChemoStatEnergy;
  double TotParticleEnergy, TotPhaseFactorEnergy, TotEnergy;
  
//Input/OutPut options
  int datarate, framerate;
  char Trajfile[1024], Datafile[1024], Restartfile[1024];
  FILE *fileptr_traj, *fileptr_data, *fileptr_restart, *fileptr_input;
  int vmdDoF, *vmdLength;
  
} sys;

#endif
