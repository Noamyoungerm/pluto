/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Maxwell solver for the Resistive RMHD equations.

  \author  A. Mignone (mignone@to.infn.it)
  \date    Dec 09, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#if PHYSICS == ResRMHD
/* ********************************************************************* */
void CT_MaxwellSolver (const Data *d, const EMF *emf, Grid *grid)
/*!
 *
 *********************************************************************** */
{
  int i, j, k;
  int recE = RECONSTRUCTION;
  int recB = RECONSTRUCTION;
  double lambda;
  DIM_EXPAND(double ***Bx1s = d->Vs[BX1s];  ,
             double ***Bx2s = d->Vs[BX2s];  ,
             double ***Bx3s = d->Vs[BX3s];)
  #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
  DIM_EXPAND(double ***Ex1s = d->Vs[EX1s];  ,
             double ***Ex2s = d->Vs[EX2s];  ,
             double ***Ex3s = d->Vs[EX3s];)
  #endif

  static uint16_t *flag;
  static double *BxL, *ByL, *BzL;
  static double *BxR, *ByR, *BzR;
  static double *ExL, *EyL, *EzL;
  static double *ExR, *EyR, *EzR;

  if (BxL == NULL) {
    BxL = ARRAY_1D(NMAX_POINT, double);
    ByL = ARRAY_1D(NMAX_POINT, double);
    BzL = ARRAY_1D(NMAX_POINT, double);
    BxR = ARRAY_1D(NMAX_POINT, double);
    ByR = ARRAY_1D(NMAX_POINT, double);
    BzR = ARRAY_1D(NMAX_POINT, double);

    ExL = ARRAY_1D(NMAX_POINT, double);
    EyL = ARRAY_1D(NMAX_POINT, double);
    EzL = ARRAY_1D(NMAX_POINT, double);
    ExR = ARRAY_1D(NMAX_POINT, double);
    EyR = ARRAY_1D(NMAX_POINT, double);
    EzR = ARRAY_1D(NMAX_POINT, double);
    flag = ARRAY_1D(NMAX_POINT, uint16_t);
  }

/* --------------------------------------------------------
   X1. Reconstruct along the x1-direction:

     o Ez(j+1/2) and By(j+1/2) --> (i+1/2, j+1/2, k)
     o Bz(j+1/2) and Ey(j+1/2) --> (i+1/2, j+1/2, k)
     
     o Ey(k+1/2) and Bz(k+1/2) --> (i+1/2, j, k+1/2)
     o By(k+1/2) and Ez(k+1/2) --> (i+1/2, j, k+1/2)
   -------------------------------------------------------- */
   
  #if INCLUDE_IDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (j = emf->jbeg; j <= emf->jend; j++){

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){  /* Interface flag */
      flag[i] = d->flag[k][j][i] | d->flag[k][j+1][i];
    }
    #endif
    ArrayReconstruct (emf->ezj, flag, i, j, k, IDIR, EzL, EzR, recE, grid);
    ArrayReconstruct (Bx2s,     flag, i, j, k, IDIR, ByL, ByR, recB, grid);
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Ex3e[k][j][i] = 0.25*(EzL[i] + EzR[i]) + 0.5*(ByR[i] - ByL[i]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Bzj, flag, i, j, k, IDIR, BzL, BzR, recB, grid);
    ArrayReconstruct (Ex2s,     flag, i, j, k, IDIR, EyL, EyR, recE, grid);
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Bx3e[k][j][i] = 0.25*(BzL[i] + BzR[i]) - 0.5*(EyR[i] - EyL[i]);
    }
    #endif
    #endif  /* INCLUDE_JDIR */

    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (i = emf->ibeg; i <= emf->iend+1; i++){  /* Interface flag */
      flag[i] = d->flag[k][j][i] | d->flag[k+1][j][i];
    }
    #endif
    ArrayReconstruct (emf->eyk, flag, i, j, k, IDIR, EyL, EyR, recE, grid); 
    ArrayReconstruct (Bx3s,     flag, i, j, k, IDIR, BzL, BzR, recB, grid);
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Ex2e[k][j][i] = 0.25*(EyL[i] + EyR[i]) - 0.5*(BzR[i] - BzL[i]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Byk, flag, i, j, k, IDIR, ByL, ByR, recB, grid); 
    ArrayReconstruct (Ex3s,     flag, i, j, k, IDIR, EzL, EzR, recE, grid);
    for (i = emf->ibeg; i <= emf->iend; i++){
      emf->Bx2e[k][j][i] = 0.25*(ByL[i] + ByR[i]) + 0.5*(EzR[i] - EzL[i]);
    }
    #endif
    #endif  /* INCLUDE_KDIR */
  }}
  #endif

/* --------------------------------------------------------
   X2. Reconstruct along the x2-direction:

     o Ex(k+1/2) and Bz(k+1/2) --> (i, j+1/2, k+1/2)
     o Bx(k+1/2) and Ez(k+1/2) --> (i, j+1/2, k+1/2)

     o Ez(i+1/2) and Bx(i+1/2) --> (i+1/2, j+1/2, k)
     o Bz(i+1/2) and Ex(i+1/2) --> (i+1/2, j+1/2, k)
   -------------------------------------------------------- */

  #if INCLUDE_JDIR
  for (k = emf->kbeg; k <= emf->kend; k++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){
  
    #if INCLUDE_KDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){  /* Interfce flag */
      flag[j] = d->flag[k][j][i] | d->flag[k+1][j][i];
    }
    #endif
    ArrayReconstruct (emf->exk, flag, i, j, k, JDIR, ExL, ExR, recE, grid); 
    ArrayReconstruct (Bx3s,     flag, i, j, k, JDIR, BzL, BzR, recB, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Ex1e[k][j][i] = 0.25*(ExL[j] + ExR[j]) + 0.5*(BzR[j] - BzL[j]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Bxk, flag, i, j, k, JDIR, BxL, BxR, recB, grid); 
    ArrayReconstruct (Ex3s,     flag, i, j, k, JDIR, EzL, EzR, recE, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Bx1e[k][j][i] = 0.25*(BxL[j] + BxR[j]) - 0.5*(EzR[j] - EzL[j]);
    }
    #endif
    #endif /* INCLUDE_KDIR */
    
    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (j = emf->jbeg; j <= emf->jend+1; j++){  /* Inteface flag */
      flag[j] = d->flag[k][j][i+1] | d->flag[k][j][i];
    }
    #endif
    ArrayReconstruct (emf->ezi, flag, i, j, k, JDIR, EzL, EzR, recE, grid);  
    ArrayReconstruct (Bx1s,     flag, i, j, k, JDIR, BxL, BxR, recB, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Ex3e[k][j][i] += 0.25*(EzL[j] + EzR[j]) - 0.5*(BxR[j] - BxL[j]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Bzi, flag, i, j, k, JDIR, BzL, BzR, recB, grid);  
    ArrayReconstruct (Ex1s,     flag, i, j, k, JDIR, ExL, ExR, recE, grid);
    for (j = emf->jbeg; j <= emf->jend; j++){
      emf->Bx3e[k][j][i] += 0.25*(BzL[j] + BzR[j]) + 0.5*(ExR[j] - ExL[j]);
    }
    #endif
    #endif /* INCLUDE_IDIR */

  }}
  #endif

/* --------------------------------------------------------
   X3. Reconstruct along the x3-direction:

     o Ey(i+1/2) and Bx(i+1/2) --> (i+1/2, j, k+1/2)
     o By(i+1/2) and Ex(i+1/2) --> (i+1/2, j, k+1/2)
     
     o Ex(j+1/2) and By(j+1/2) --> (i, j+1/2, k+1/2)
     o Bx(j+1/2) and Ey(j+1/2) --> (i, j+1/2, k+1/2)
   -------------------------------------------------------- */

  #if INCLUDE_KDIR
  for (j = emf->jbeg; j <= emf->jend; j++){ 
  for (i = emf->ibeg; i <= emf->iend; i++){

    #if INCLUDE_IDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){   /* Interface flag */
      flag[k] = d->flag[k][j][i] | d->flag[k][j][i+1];
    }
    #endif
    ArrayReconstruct (emf->eyi, flag, i, j, k, KDIR, EyL, EyR, recE, grid);
    ArrayReconstruct (Bx1s,     flag, i, j, k, KDIR, BxL, BxR, recB, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      emf->Ex2e[k][j][i] += 0.25*(EyL[k] + EyR[k]) + 0.5*(BxR[k] - BxL[k]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Byi, flag, i, j, k, KDIR, ByL, ByR, recB, grid);
    ArrayReconstruct (Ex1s,     flag, i, j, k, KDIR, ExL, ExR, recE, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      emf->Bx2e[k][j][i] += 0.25*(ByL[k] + ByR[k]) - 0.5*(ExR[k] - ExL[k]);
    }
    #endif
    #endif /* INCLUDE_IDIR */

    #if INCLUDE_JDIR
    #if SHOCK_FLATTENING == MULTID
    for (k = emf->kbeg; k <= emf->kend+1; k++){   /* Interface flag */
      flag[k] = d->flag[k][j][i] | d->flag[k][j+1][i];
    }
    #endif
    ArrayReconstruct (emf->exj, flag, i, j, k, KDIR, ExL, ExR, recE, grid);
    ArrayReconstruct (Bx2s,     flag, i, j, k, KDIR, ByL, ByR, recB, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      emf->Ex1e[k][j][i] += 0.25*(ExL[k] + ExR[k]) - 0.5*(ByR[k] - ByL[k]);
    }
    #if DIVE_CONTROL == CONSTRAINED_TRANSPORT
    ArrayReconstruct (emf->Bxj, flag, i, j, k, KDIR, BxL, BxR, recB, grid);
    ArrayReconstruct (Ex2s,     flag, i, j, k, KDIR, EyL, EyR, recE, grid);  
    for (k = emf->kbeg; k <= emf->kend; k++){ 
      emf->Bx1e[k][j][i] += 0.25*(BxL[k] + BxR[k]) + 0.5*(EyR[k] - EyL[k]);
    }
    #endif
    #endif /* INCLUDE_JDIR */
  }}
  #endif
    
}

#endif /* PHYSICS == ResRMHD */
