/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief GFORCE Riemann solver for the classical MHD equations.

  Solve the Riemann problem for MHD equations using the
  generalized FORCE flux described in Section 3 of
  Toro \& Titarev, JCP (2006) 216, 403.

  \b Reference:
    - "MUSTA Fluxes for systems of conservation laws"
       Toro \& Titarev, JCP (2006) 216, 403

  \authors A. Mignone (mignone@to.infn.it)
  \date    Jul 20, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void GFORCE_Solver (const Sweep *sweep, int beg, int end, 
                    double *cmax, Grid *grid)
/*!
 *
 *********************************************************************** */
{
  int    nv, i;
  static unsigned char *flag;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double dtdx, dxdt, om, cfl;
  double *a2L = stateL->a2, *a2R = stateR->a2;
  double **uL = stateL->u, **uR  = stateR->u;
  double *vL, *vR;
  double **fL = stateL->flux, **fR  = stateR->flux;
  double *pL = stateL->prs,  *pR = stateR->prs;

  static State stateLW;
  static double **fluxLF;

#if DIVB_CONTROL == EIGHT_WAVES
  print ("! GFORCE_Solver(): does not work with Powell's 8-wave\n");
  QUIT_PLUTO(1);
#endif
  
/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (fluxLF == NULL){
    fluxLF = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);

    StateStructAllocate(&stateLW);

  /* ----------------------------------
     Set scalars to zero (scalar flux
     will be computed later)
     ---------------------------------- */
    #if NSCL > 0
    for (i = 0; i < NMAX_POINT; i++){
      NSCL_LOOP(nv) stateLW.u[i][nv] = 0.0;
    }
    #endif
  }

/* --------------------------------------------------------
   1. Set background field and apply GLM pre-conditioner
   -------------------------------------------------------- */

#if BACKGROUND_FIELD == YES
  /* -- Background field computed in stateL->bck = stateR->bck */
  GetBackgroundField (stateL,  beg, end, FACE_CENTER, grid);
  GetBackgroundField (&stateLW, beg, end, FACE_CENTER, grid);
#endif
  
#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif
 
/* --------------------------------------------------------
   2. Compute characteristic velocities
   -------------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  HLL_Speed (stateL, stateR, sweep->SL, sweep->SR, beg, end);
  for (i = beg; i <= end; i++) {
    cmax[i] = MAX(fabs(sweep->SL[i]), fabs(sweep->SR[i]));
  }
      
/* --------------------------------------------------------
   3. Compute fluxes on the original data
   -------------------------------------------------------- */

  Flux (stateL, beg, end);  
  Flux (stateR, beg, end);  

/* --------------------------------------------------------
   4a. Compute LF flux and LW conservative state
   -------------------------------------------------------- */

  for (i = beg; i <= end; i++) {
    dxdt = cmax[i];
    dtdx = 1.0/dxdt;
    NFLX_LOOP(nv){
      fluxLF[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - dxdt*(uR[i][nv] - uL[i][nv]));
      stateLW.u[i][nv] =   0.5*(uL[i][nv] + uR[i][nv]) 
                         - 0.5*dtdx*(fR[i][nv] - fL[i][nv]);
    }
    stateLW.u[i][MXn] -= 0.5*dtdx*(pR[i] - pL[i]);
    #if ENTROPY_SWITCH != NO
    vL = stateL->v[i];
    vR = stateR->v[i];

    double fLe = fL[i][RHO]*vL[ENTR];
    double fRe = fR[i][RHO]*vR[ENTR];
    stateLW.u[i][ENTR] =  0.5*(uL[i][ENTR] + uR[i][ENTR]) 
                        - 0.5*dtdx*(fRe - fLe);
    #endif 
  }

  int status =  ConsToPrim (stateLW.u, stateLW.v, beg, end, flag);
  if (status != 0){
    printLog ("! GFORCE_Solver(): u(LW) not physical\n");
    QUIT_PLUTO(1);
  }
  SoundSpeed2 (&stateLW, beg, end, FACE_CENTER, NULL);
  Flux (&stateLW, beg, end);

/* --------------------------------------------------------
   4b. Take convex average between F(LF) and F(LW)
   -------------------------------------------------------- */

  cfl = RuntimeGet()->cfl;
  for (i = beg; i <= end; i++) {
    #ifndef GFORCE_OMEGA
    om  = 1.0/(1.0 + cfl);
    #else
    om = GFORCE_OMEGA;
    #endif
    #if SHOCK_FLATTENING == MULTID
    if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){
      om = 0.0;
    }
    #endif

    NFLX_LOOP(nv) {
      sweep->flux[i][nv] = om*stateLW.flux[i][nv] + (1.0 - om)*fluxLF[i][nv];
    }
    sweep->press[i] = om*stateLW.prs[i] + (1.0 - om)*0.5*(pL[i] + pR[i]);
  }

/* --------------------------------------------------------
   5. Define point and diffusive fluxes for CT
   -------------------------------------------------------- */
  
#if DIVB_CONTROL == CONSTRAINED_TRANSPORT 
  CT_Flux (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   6. Add CR flux contribution using simplified upwinding.
   -------------------------------------------------------- */

  #if (PARTICLES == PARTICLES_CR) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux (stateL, beg, end);
  Particles_CR_Flux (stateR, beg, end);

  for (i = beg; i <= end; i++){
    for (nv = NFLX; nv--; ) {
      sweep->flux[i][nv] += 0.5*(stateL->fluxCR[i][nv] + stateR->fluxCR[i][nv]);
    }
  }  
  #endif
  
}

