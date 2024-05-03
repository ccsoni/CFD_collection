#include "run_param.h"

REAL csl5(REAL fuu, REAL fu, REAL f0, REAL fd, REAL fdd, 
	  REAL dx, REAL del)
// dx = |dx|
// del = |a*dt/dx|
{
#if 0
  return (dx*(del*((4.0*(fuu-6.5*fu+23.5*f0)+6.0*(9.0*fd-fdd))
                   +del*(5.0*(15.0*(f0-fd)-(fu-fdd))
                         +del*((10.0*(3.0*fu-4.0*f0+fd)-5.0*(fuu-fdd))
                               +del*(5.0*(-3.0*(f0-fd)+(fu-fdd))
                                     +del*(6.0*f0-4.0*(fu+fd)+(fuu+fdd)))))))/120.0);
#else
  return ((4.0*(fuu-6.5*fu+23.5*f0)+6.0*(9.0*fd-fdd))
               +del*(5.0*(15.0*(f0-fd)-(fu-fdd))
                     +del*((10.0*(3.0*fu-4.0*f0+fd)-5.0*(fuu-fdd))
                           +del*(5.0*(-3.0*(f0-fd)+(fu-fdd))
                                 +del*(6.0*f0-4.0*(fu+fd)+(fuu+fdd))))))/120.0;
#endif
}
