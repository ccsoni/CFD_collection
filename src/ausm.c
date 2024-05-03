#include <math.h>
#include "fluid.h"

REAL Mp_beta(REAL M)
{
  return 0.25*SQR(M+1.0) + AUSM_beta*SQR(M*M-1.0);
}

REAL Mm_beta(REAL M)
{
  return -0.25*SQR(M-1.0) - AUSM_beta*SQR(M*M-1.0);
}

REAL M_p(REAL M) 
{
  REAL Mplus;
  if(fabs(M) >= 1.0) {
    Mplus = 0.5*(M+fabs(M));
  }else{
    Mplus = Mp_beta(M);
  }
  
  return Mplus;
}

REAL M_m(REAL M)
{
  REAL Mminus;
  if(fabs(M) >= 1.0) {
    Mminus = 0.5*(M+fabs(M));
  }else{
    Mminus = Mm_beta(M);
  }

  return Mminus;
}

REAL Pp_alpha(REAL M)
{
  return 0.25*SQR(M+1.0)*(2.0-M)+AUSM_alpha*M*SQR(M*M-1.0);
}

REAL Pm_alpha(REAL M)
{
  return 0.25*SQR(M-1.0)*(2.0+M)-AUSM_alpha*M*SQR(M*M-1.0);
}

REAL P_p(REAL M)
{
  REAL Pplus;
  if(fabs(M) >= 1.0) {
    Pplus = 0.5*(1.0+copysign(1.0, M));
  }else{
    Pplus = Pp_alpha(M);
  }
  
  return Pplus;
}

REAL P_m(REAL M)
{
  REAL Pminus;
  if(fabs(M) >= 1.0) {
    Pminus = 0.5*(1.0-copysign(1.0,M));
  }else{
    Pminus = Pm_alpha(M);
  }

  return Pminus;
}
