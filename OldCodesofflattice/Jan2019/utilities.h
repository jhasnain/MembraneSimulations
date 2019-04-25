#ifndef UTILITIES_H
#define UTILITIES_H

#include "membranesys.h"
#include <string.h>

void PrintConfiguration(FILE *fp, sys *Sys, int opt);
void ReadConfiguration(FILE *fp, sys *Sys);
void PrintData(FILE *fp, sys *Sys, long currentMCstep);
void PrintDataToScreen(FILE *fp, sys *Sys, long currentMCstep);
void PrintVec(vec target, char *note);
void PrintJustMesh(FILE *fp, FILE *gp, sys *Sys);
void ConfigurationRestart(FILE *fp, sys *Sys);

void twodvecdiff(vec *target, vec v2, vec v1);
void twodvecsum(vec *target, vec v2, vec v1);
double twod_NormSqr(vec v);

void threedvecdiff(vec *target, vec v2, vec v1);
double threed_NormSqr(vec v);
void threedvecsumscalar(vec *target, vec v2, vec v1, double a);
double deltaRsqr(vec v1, vec v2, double Lx, double Ly);
void threed_NormalizeVec(vec *target);
void crossprod(vec *target, vec v1, vec v2);
double dotprod( vec v1, vec v2);
void rescalevec(vec *target, double a);

void Position_To_Index(int *i, double x, double dx);
void Position_To_Index_PBC(int *i, double x, double dx, int maxx);
void Index_To_Position(double *x, int i, double dx);
void Distance_To_IndexRan(int *IR, double R, double dr);

void TranslateParticle(particle *Target, particle *Source, double x, double y, double z);
void InterpolateProperty(double *Target, vec diffvec, double xy, double xpy, double xyp, double xpyp);

void NearestImageConv(vec *myvec, double Lx, double Ly );
void ScalarNearestImageConv(double *val, double L );
int IndexPBC(int index, int Nmax);
void PosPBC(vec *myvec, double Lx, double Ly);
void CopyVec(vec *target, vec source);

void TestCollisionDetection(sys *Sys);

#endif
