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


/**************************************************************************/
/* wrappers for some funcs in tab1d_Of_rho0_AtT0.c */
/**************************************************************************/

/* for EoS->vars_from_hm1 */
void tab1d_rho0_P_rhoE_drho0dhm1_from_hm1(double hm1, double *rho0,
                                          double *P, double *rhoE,
                                          double *drho0dhm1)
{
  double epsl;

  tab1d_rho0_epsl_P_drho0dhm1_Of_hm1_AtT0(hm1, rho0, &epsl, P,
                                          drho0dhm1);
  /* total energy density */
  *rhoE = *rho0 + (*rho0) * epsl;
}

/* for EoS->rho0_of_hm1 */
double tab1d_rho0_of_hm1(double hm1)
{
  double rho0, epsl, P, dPdrho0, dPdepsl;

  tab1d_Of_hm1_AtT0(hm1, &rho0, &epsl, &P, &dPdrho0, &dPdepsl);
  return rho0;
}

/* for EoS->hm1_of_P */
double tab1d_hm1_of_P(double P)
{
  double rho0, epsl, dPdrho0, dPdepsl;

  tab1d_Of_P_AtT0(P, &rho0, &epsl, &dPdrho0, &dPdepsl);

  return hm1_of_rho0_epsl_P(rho0, epsl, P);
}

/* for EoS->rho0_rhoE_from_P */
void tab1d_rho0_rhoE_from_P(double P, double *rho0, double *rhoE)
{
  double epsl, dPdrho0, dPdepsl;

  tab1d_Of_P_AtT0(P, rho0, &epsl, &dPdrho0, &dPdepsl);

  /* total energy density */
  *rhoE = *rho0 + (*rho0) * epsl;
}
