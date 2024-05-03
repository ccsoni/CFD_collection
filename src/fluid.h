#pragma once
#include "run_param.h"

struct fluid_1d 
{
  REAL dens, momx, eneg, etrp;
  short high_mach, under_shock;
};

struct fluid_flux 
{
  REAL dens_flux, momx_flux, eneg_flux, etrp_flux;
};

void calc_flux(struct fluid_1d *, struct fluid_flux *, struct run_param*);
REAL muscl_L(REAL, REAL, REAL);
REAL muscl_R(REAL, REAL, REAL);
REAL mp5_R(REAL, REAL, REAL, REAL, REAL);
REAL mp5_L(REAL, REAL, REAL, REAL, REAL);
REAL mp5_bound_R(REAL, REAL, REAL, REAL, REAL, REAL);
REAL mp5_bound_L(REAL, REAL, REAL, REAL, REAL, REAL);
REAL csl5(REAL, REAL, REAL, REAL, REAL, REAL, REAL);
REAL M_p(REAL);
REAL M_m(REAL);
REAL P_p(REAL);
REAL P_m(REAL);

#define GAMMA (5.0/3.0)
#define MUSCL_KAPPA (-2.5)

#define AUSM_alpha (0.1875) /* 3.0/16.0 */
//#define AUSM_alpha (0.0)
#define AUSM_beta  (0.125) /* 1.0/8.0 */
//#define AUSM_beta (0.0)

#define HIGH_MACH_THRESHOLD (0.1)
#define SHOCK_PRESS_THRESHOLD (0.3)

#if defined(__AUSMP__) || defined(__AUSMP_MUSCL__) || defined(__AUSM__MP5) || defined(__HLL__) || defined(__HLLC__) || defined(__HLLC_MUSCL__) || defined(__HLLC_MP5__)
#define __ENTROPY_INTEGRATION_AVAILABLE__
#endif
