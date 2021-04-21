/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Add source terms to the right hand side of relativistic HD/MHD eqns.

  Add source terms to the right hand side of RHD or RMHD equations in
  conservative form.
  These include

  -# Body forces;
  -# Powell's 8-waves source terms;
  

  Care is taken to ensure that gravity components are included even
  when a direction is not active.
  The following table summarizes:
  
  1D (x):
  Sweep| gx | gy | gz |
  -----|----|----|----|
   x   | o  | o  | o  |

  2D (x,y):
  Sweep| gx | gy | gz |
  -----|----|----|----|
   x   | o  |    |    |
   y   |    | o  | o  |

  2D (x,z):
  Sweep| gx | gy | gz |
  -----|----|----|----|
   x   | o  |    |    |
   z   |    | o  | o  |

  3D (x,y,z):
  Sweep| gx | gy | gz |
  -----|----|----|----|
   x   | o  |    |    |
   y   |    | o  |    |
   z   |    |    | o  |
        
  For consistency, the same approach must be used in PrimSource().
  
  \author A. Mignone (mignone@to.infn.it)
  \date   Apr 25, 2020
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifndef iMPHI
 #define iMPHI MX2  /* -- for Cartesian coordinates -- */
#endif

/* *********************************************************************** */
void RightHandSideSource (const Sweep *sweep, timeStep *Dts,
                          int beg, int end, double dt, double *phi_p, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      phi_p  force potential at interfaces
 * \param [in]      grid   pointer to Grid structure
 *
 * \return This function has no return value.
 ************************************************************************* */
{
  int    i, j, k, nv;

  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  double r_1, g[3], fg[3], vg;

  double dtdx, scrh, ct;
  double Sm;
  double *x1   = grid->x[IDIR],  *x2  = grid->x[JDIR],  *x3  = grid->x[KDIR];
  double *x1p  = grid->xr[IDIR], *x2p = grid->xr[JDIR], *x3p = grid->xr[KDIR];
  double *x1m  = grid->xl[IDIR], *x2m = grid->xl[JDIR], *x3m = grid->xl[KDIR];
  double *dx1  = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3 = grid->dx[KDIR];
#if GEOMETRY == SPHERICAL
  double *rt   = grid->rt;
  double *sp   = grid->sp;
  double *sm   = grid->sp-1;
  double *s    = grid->s;
  double *dmu  = grid->dmu;
#endif
  double ***dV = grid->dV;
  double **rhs  = sweep->rhs;
  double **flux = sweep->flux;
  double **vp   = stateL->v;
  double **vm   = stateR->v-1;
  double *p;
  double cl;
  double lor2, vel2, vphi, phi_c;
  double *v;

/* --------------------------
      pointer shortcuts
   -------------------------- */

  p  = sweep->press;
  
  i = g_i;  /* will be redefined during x1-sweep */
  j = g_j;  /* will be redefined during x2-sweep */
  k = g_k;  /* will be redefined during x3-sweep */
  
  if (g_dir == IDIR){

    for (i = beg; i <= end; i++) {
      dtdx = dt/dx1[i];
      v    = stateC->v[i];

    /* --------------------------------------------
       I1. Add geometrical source term
       -------------------------------------------- */

#if GEOMETRY == CARTESIAN

#elif GEOMETRY == CYLINDRICAL

#elif GEOMETRY == POLAR

#elif GEOMETRY == SPHERICAL

#endif

    /* ----------------------------------------------------
       I3. Include body forces all at once during the first
           sweep.
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
      vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
      lor2 = 1.0/(1.0 - vel2);
      BodyForceVector(v, g, x1[i], x2[j], x3[k]);

      vg       = v[VX1]*g[IDIR] + v[VX2]*g[JDIR] + v[VX3]*g[KDIR];
      fg[IDIR] = v[RHO]*lor2*(lor2*v[VX1]*vg + g[IDIR]);
      rhs[i][MX1] += dt*fg[IDIR];
      IF_ENERGY (rhs[i][ENG] += dt*v[VX1]*fg[IDIR];)

      /* ----------------------------------------
         Add non-active dimensions during
         this sweep.
         ---------------------------------------- */

      #if !INCLUDE_JDIR && !INCLUDE_KDIR
      fg[JDIR] = v[RHO]*lor2*(lor2*v[VX2]*vg + g[JDIR]);
      fg[KDIR] = v[RHO]*lor2*(lor2*v[VX3]*vg + g[KDIR]);
      rhs[i][MX2] += dt*fg[JDIR];
      rhs[i][MX3] += dt*fg[KDIR];
      IF_ENERGY (rhs[i][ENG] += dt*v[VX2]*fg[JDIR];)
      IF_ENERGY (rhs[i][ENG] += dt*v[VX3]*fg[KDIR];)
      #endif

      #endif  /* BODY_FORCE & VECTOR */


      #if (BODY_FORCE & POTENTIAL)
      #error Cannot use BodyForcePotential in relativistic flows
      #endif

    }

  } else if (g_dir == JDIR){

    scrh = dt;
#if GEOMETRY == POLAR
    scrh /= x1[i];
    r_1   = 1.0/x1[i];
#elif GEOMETRY == SPHERICAL
    scrh /= rt[i];
    r_1   = 1.0/rt[i];
#endif
    for (j = beg; j <= end; j++) {
      dtdx = scrh/dx2[j];
      v = stateC->v[j];

#if GEOMETRY != SPHERICAL


#elif GEOMETRY == SPHERICAL


#endif  /* GEOMETRY == SPHERICAL */

    /* ----------------------------------------------------
       J3. Include Body force
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
      vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
      lor2 = 1.0/(1.0 - vel2);

      BodyForceVector(v, g, x1[i], x2[j], x3[k]);
      vg       = v[VX1]*g[IDIR] + v[VX2]*g[JDIR] + v[VX3]*g[KDIR];
      fg[JDIR] = v[RHO]*lor2*(lor2*v[VX2]*vg + g[JDIR]);
      rhs[j][MX2] += dt*fg[JDIR];
      IF_ENERGY (rhs[j][ENG] += dt*v[VX2]*fg[JDIR];)

      #if !INCLUDE_KDIR
      fg[KDIR] = v[RHO]*lor2*(lor2*v[VX3]*vg + g[KDIR]);
      rhs[j][MX3] += dt*fg[KDIR];
      IF_ENERGY (rhs[j][ENG] += dt*v[VX3]*fg[KDIR];)
      #endif  /* !INCLUDE_KDIR */
      #endif  /* (BODY_FORCE & VECTOR) */

    }

  }else if (g_dir == KDIR){

    scrh  = dt;
    #if GEOMETRY == SPHERICAL
    scrh *= dx2[j]/(rt[i]*dmu[j]);
    #endif

    for (k = beg; k <= end; k++) {
      dtdx = scrh/dx3[k];
      v    = stateC->v[k];

    /* ----------------------------------------------------
       K3. Include body forces
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
      vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
      lor2 = 1.0/(1.0 - vel2);

      BodyForceVector(v, g, x1[i], x2[j], x3[k]);
      vg       = v[VX1]*g[IDIR] + v[VX2]*g[JDIR] + v[VX3]*g[KDIR];
      fg[KDIR] = v[RHO]*lor2*(lor2*v[VX3]*vg + g[KDIR]);

      rhs[k][MX3] += dt*fg[KDIR];
      IF_ENERGY (rhs[k][ENG] += dt*v[VX3]*fg[KDIR];)
      #if !INCLUDE_JDIR
      fg[JDIR] = v[RHO]*lor2*(lor2*v[VX2]*vg + g[JDIR]);
      rhs[k][MX2] += dt*fg[JDIR];
      IF_ENERGY (rhs[k][ENG]   += dt*v[VX2]*fg[JDIR];)
      #endif
      #endif  /* (BODY_FORCE & VECTOR) */
      
    }
  }

/* --------------------------------------------------
              Powell's source terms
   -------------------------------------------------- */

  #if (PHYSICS == RMHD) && (DIVB_CONTROL == EIGHT_WAVES)
  for (i = beg; i <= end; i++) {
    rhs[i][MX1] += dt*sweep->src[i][MX1];
    rhs[i][MX2] += dt*sweep->src[i][MX2];
    rhs[i][MX3] += dt*sweep->src[i][MX3];

    rhs[i][BX1] += dt*sweep->src[i][BX1];
    rhs[i][BX2] += dt*sweep->src[i][BX2];
    rhs[i][BX3] += dt*sweep->src[i][BX3];
    #if HAVE_ENERGY
    rhs[i][ENG] += dt*sweep->src[i][ENG];
    #endif
  }
  #endif

/* -------------------------------------------------
            Extended GLM source terms
   ------------------------------------------------- */

  #if (defined GLM_MHD) && (GLM_EXTENDED == YES)
   print ("! RightHandSide(): Extended GLM source terms not defined for RMHD\n");
   QUIT_PLUTO(1);
  #endif

}
