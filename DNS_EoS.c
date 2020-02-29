/* genral EoS funcs */
/* Wolfgang Tichy, 2/2020 */

#include "sgrid.h"
#include "DNSdata.h"


/* global struct with EoS info */
tEoS EoS[1];


/**************************************************************************/
/* some general functions in that hold for any EoS */
/**************************************************************************/

/* return h-1 = epsl + P/rho0 */
double hm1_of_rho0_epsl_P(double rho0, double epsl, double P)
{
  if(rho0==0.) return 0.;
  else         return epsl + P/rho0;
}

/* return h = 1 + epsl + P/rho0 */
double h_of_rho0_epsl_P(double rho0, double epsl, double P)
{
  return hm1_of_rho0_epsl_P(rho0, epsl, P) + 1.;
}
