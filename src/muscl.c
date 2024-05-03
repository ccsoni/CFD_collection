#include <math.h>
#include "run_param.h"
#include "fluid.h"

REAL muscl_L(REAL u_im1, REAL u_i, REAL u_ip1)
{
  REAL delta_p, delta_m, delta_pp, delta_mm, beta;

  beta = (3.0-MUSCL_KAPPA)/(1.0-MUSCL_KAPPA);

  delta_p = u_ip1-u_i;
  delta_m = u_i-u_im1;

  if(delta_p*beta*delta_m<0.0) {
    delta_pp=0.0;
  }else if(fabs(delta_p) < fabs(beta*delta_m)) {
    delta_pp=delta_p;
  }else{
    delta_pp = beta*delta_m;
  }

  if(delta_m*beta*delta_p<0.0) {
    delta_mm=0.0;
  }else if(fabs(delta_m) < fabs(beta*delta_p)){
    delta_mm=delta_m;
  }else{
    delta_mm=beta*delta_p;
  }

  return (u_i+0.25*(1.0-MUSCL_KAPPA)*delta_mm+0.25*(1.0+MUSCL_KAPPA)*delta_pp);
}

REAL muscl_R(REAL u_im1, REAL u_i, REAL u_ip1)
{
  REAL delta_p, delta_m, delta_pp, delta_mm, beta;

  beta = (3.0-MUSCL_KAPPA)/(1.0-MUSCL_KAPPA);

  delta_p = u_ip1-u_i;
  delta_m = u_i-u_im1;

  if(delta_p*beta*delta_m<0.0) {
    delta_pp=0.0;
  }else if(fabs(delta_p) < fabs(beta*delta_m)) {
    delta_pp=delta_p;
  }else{
    delta_pp = beta*delta_m;
  }

  if(delta_m*beta*delta_p<0.0) {
    delta_mm=0.0;
  }else if(fabs(delta_m) < fabs(beta*delta_p)){
    delta_mm=delta_m;
  }else{
    delta_mm=beta*delta_p;
  }

  return (u_i-0.25*(1.0+MUSCL_KAPPA)*delta_mm-0.25*(1.0-MUSCL_KAPPA)*delta_pp);
}
