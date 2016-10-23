#ifndef _KF_KS_DS_v2_steady_
#define _KF_KS_DS_v2_steady_

#include "KFKSDS-steady.h"

extern "C" void KFKSDS_steady_v2eps (int *dim, double *sy, double *sZ, double *sZtZ, 
  double *sT, double *sH, double *sR, 
  double *V, double *sQ, double *sa0, double *sP0, 
  double *tol, int *maxiter, double *ksconvfactor,
  double *res);

extern "C" void KFKSDS_steady_v2eta (int *dim, double *sy, double *sZ, double *sZtZ, 
  double *sT, double *sH, double *sR, 
  double *V, double *sQ, double *sa0, double *sP0, 
  double *tol, int *maxiter, double *ksconvfactor,
  double *res);

#endif
