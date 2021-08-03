/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Wrapper function to assign shearing-box boundary conditions.

  SB_Boundary() is a wrapper function used to assign shearing-box 
  boundary conditions on both cell-centered and staggered variables at
  the X1_BEG and X1_END boundaries.
  It is called as a regular boundary condition after periodic boundary
  conditions have already been imposed on the solution arrays.\n
  The function defines, for each type of variable (centered or staggered),
  the corresponding boundary region where shearing conditions must be
  applied using the RBox structure.
  The actual boundary condition is imposed by calling SB_SetBoundaryVar() 
  with the desired array and its box layout.

  \authors A. Mignone (mignone@to.infn.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)
  \date    Aug 21, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double sb_vy;

/* ********************************************************************* */
void SB_Boundary (const Data *d, int side, Grid *grid) 
/*! 
 * Main wrapper function used to assign shearing-box boundary conditions
 * on flow variables.
 *
 * \param d    pointer to the PLUTO Data structure
 * \param side the side of the computational domain (X1_BEG or X1_END) 
 * \param grid pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *
 * \todo Check if sb_vy needs to be global.
 *********************************************************************** */
{
  int    i, j, k, nv;
  double t, Lx;
  RBox   box;

  Lx    = g_domEnd[IDIR] - g_domBeg[IDIR];
  sb_vy = fabs(2.0*SB_A*Lx);

#if INCLUDE_JDIR

/* ----------------------------------------------
   1. Non-axisymmetric case
   ---------------------------------------------- */

  t = g_time;
  #if TIME_STEPPING == RK2
 if (g_intStage == 2) t += g_dt;
  #elif TIME_STEPPING == RK3
 if (g_intStage == 2) t += 0.5*g_dt;
 if (g_intStage == 3) t += g_dt;
  #endif
  #ifdef CTU
  if (g_intStage == 2) t += 0.5*g_dt;
  #endif

/* ----------------------------------------------
   1a. X1 Beg Boundary
   ---------------------------------------------- */

  if (side == X1_BEG){

  /* ---- loop on cell-centered variables ---- */

    box.ibeg = 0; box.iend = IBEG-1;
    box.jbeg = 0; box.jend = NX2_TOT-1;
    box.kbeg = 0; box.kend = NX3_TOT-1;
    for (nv = 0; nv < NVAR; nv++){
      #ifdef STAGGERED_MHD
      DIM_EXPAND(if (nv == BX1) continue;  ,
               if (nv == BX2) continue;  ,
               if (nv == BX3) continue;)
      #endif
      SB_SetBoundaryVar(d->Vc[nv], &box, side, t, grid);
      #ifndef FARGO
      if (nv == VX2) X1_BEG_LOOP(k,j,i) d->Vc[nv][k][j][i] += sb_vy;
      #endif
    }  /* -- end loop on cell-centered variables -- */

    #ifdef STAGGERED_MHD
    box.ibeg =  0; box.iend = IBEG-1;
    box.jbeg = -1; box.jend = NX2_TOT-1;
    box.kbeg =  0; box.kend = NX3_TOT-1;
    SB_SetBoundaryVar(d->Vs[BX2s], &box, side, t, grid);
    #if INCLUDE_KDIR
    box.ibeg =  0; box.iend = IBEG-1;
    box.jbeg =  0; box.jend = NX2_TOT-1;
    box.kbeg = -1; box.kend = NX3_TOT-1;
    SB_SetBoundaryVar(d->Vs[BX3s], &box, side, t, grid);
    #endif /* INCLUDE_KDIR */
    #endif /* STAGGERED_MHD */
  } /* -- END side X1_BEG -- */

/* ----------------------------------------------
   1b. X1 End Boundary
   ---------------------------------------------- */

  if (side == X1_END){

  /* ---- loop on cell-centered variables ---- */

    box.ibeg = IEND+1; box.iend = NX1_TOT-1;
    box.jbeg =      0; box.jend = NX2_TOT-1;
    box.kbeg =      0; box.kend = NX3_TOT-1;
    for (nv = 0; nv < NVAR; nv++){
      #ifdef STAGGERED_MHD
      DIM_EXPAND(if (nv == BX1) continue;  ,
               if (nv == BX2) continue;  ,
               if (nv == BX3) continue;)
      #endif
      SB_SetBoundaryVar(d->Vc[nv], &box, side, t, grid);
      #ifndef FARGO
      if (nv == VX2)  X1_END_LOOP(k,j,i) d->Vc[nv][k][j][i] -= sb_vy;
      #endif
    }  /* -- end loop on cell-centered variables -- */

    #ifdef STAGGERED_MHD
    box.ibeg = IEND+1; box.iend = NX1_TOT-1;
    box.jbeg =     -1; box.jend = NX2_TOT-1;
    box.kbeg =      0; box.kend = NX3_TOT-1;
    SB_SetBoundaryVar(d->Vs[BX2s], &box, side, t, grid);
    #if INCLUDE_KDIR
    box.ibeg = IEND+1; box.iend = NX1_TOT-1;
    box.jbeg =      0; box.jend = NX2_TOT-1;
    box.kbeg =     -1; box.kend = NX3_TOT-1;
    SB_SetBoundaryVar(d->Vs[BX3s], &box, side, t, grid);
    #endif /* INCLUDE_KDIR */
    #endif /* STAGGERED_MHD */
  }

#else

/* ----------------------------------------------
   2. Axisymmetric case
   ---------------------------------------------- */

  #ifndef FARGO
  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i) d->Vc[VX2][k][j][i] += sb_vy;
  }else if (side == X1_END){
    X1_END_LOOP(k,j,i) d->Vc[VX2][k][j][i] -= sb_vy;
  }
  #endif
  
#endif 
}
