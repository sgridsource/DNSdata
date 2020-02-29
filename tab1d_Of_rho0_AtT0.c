/* deal with 1d EoS tables, where we give all EoS vars as funcs of rho0 */
/* Wolfgang Tichy, 2/2020 */

#include "sgrid.h"
#include "DNSdata.h"

#define ln10 2.302585092994046
#define PR 0



/* global struct to hold 1d table of rho at T=0 and Y=0 */
struct TAB1D_rho0 {
  double *Logrho0;  /* primary x-var according to which tab is accessed */
  double *Logepsl;  /* primary y-var that depends on primary x-var */
  double *dLogepsl; /* deriv dy/dx */
  int n;
} tab1d_rho0[1];


/* struct to pass around vars */
typedef struct {
  double rho0;
  double epsl;
  double P;
  double dPdrho0;
  double dPdepsl;
  double P0;
  double hm10;
  double LogP0;
  double Loghm10;
} ttab1d_vars;


/**************************************************************************/
/* functions that deal with 1d EoS tables of rho0 */
/**************************************************************************/

/* alloc table members of tab1d_rho0 */
void EoS_alloc_tab1d_rho0(int n)
{
  if(n)
  {
    tab1d_rho0->Logrho0  = calloc(n, sizeof(double));
    tab1d_rho0->Logepsl  = calloc(n, sizeof(double));
    tab1d_rho0->dLogepsl = calloc(n, sizeof(double));
  }
  tab1d_rho0->n = n;
}

/* free table */
void EoS_free_tab1d_rho0(void)
{
  if(tab1d_rho0->n)
  {
    free(tab1d_rho0->Logrho0);
    free(tab1d_rho0->Logepsl);
    free(tab1d_rho0->dLogepsl);
    tab1d_rho0->Logrho0  = NULL;
    tab1d_rho0->Logepsl  = NULL;
    tab1d_rho0->dLogepsl = NULL;
    tab1d_rho0->n = 0;
  }
}

/* min and max Logrho0 */
void min_max_Logrho0_tab1d(double *Logrho0min, double *Logrho0max)
{
  *Logrho0min = tab1d_rho0->Logrho0[0];
  *Logrho0max = tab1d_rho0->Logrho0[tab1d_rho0->n - 1];
}


/* load a 1d EoS as func of rho0 from file fname,
   keep log10(rho0), log10(eps), d(log(epsl))/d(log(rho0)) in memory
   -----------------------------------------------------------------
   The files must contain the line:
   # nrows = ...
   so that we can know the number of data rows

   Afterwards we need 3 cols:
   rho0  epsl  P
   in units of G=c=1
     col1:  rho0 = baryonic mass density
     col2:  epsl = kinetic energy per baryonic mass
     col3:  P    = pressure
   lines starting with # are comments

   !!! IMPORTANT !!!
   to use routine "tab1d_locate" below tabs have to be ordered in rho0 */
void EoS_tab1d_load_rho0_epsl_P_AtT0(char *fname)
{
  FILE *fp;
  char str[1024], par[1024], val[1024];
  int i, nrows;
  double rho0, epsl, P;

  PRF;
  printf("EoS_tab1d_load_rho0_epsl_P_AtT0:\nread in EoS table from file %s: ",
         fname);

  /* open file */
  fp = fopen(fname, "r");
  if(!fp) errorexits(" could not open file %s", fname);

  /* count data rows */
  nrows = 0;
  while(fgets(str, 1024, fp))
  {
    /* omit comments and pars */
    if(str[0] == '#') continue;
    nrows++;
  }
  printf("nrows=%d\n", nrows);
  EoS_free_tab1d_rho0();
  EoS_alloc_tab1d_rho0(nrows);

  /* read data file line by line */
  rewind(fp);
  i = 0;
  while(fgets(str, 1024, fp))
  {
    /* treat comments and pars */
    if(str[0] == '#')
    {
      if(get_par_from_str(str, par, "=", val, 1023)) printf("%s", str);
      continue;
    }

    /* read data line */
    if(sscanf(str, "%lf %lf %lf", &rho0, &epsl, &P)<3)
      errorexit("problem reading data line");

    /* save log of rho0 and epsl */
    tab1d_rho0->Logrho0[i] = log10(rho0);
    tab1d_rho0->Logepsl[i] = log10(epsl);

    // the following quantity is
    //
    // d Log epsl   rho0 d epsl      P
    // ---------- = ---- ------ = -------
    // d Log rho0   epsl d rho0   epsl rho0
    //
    // the last line uses 1st law of thermo @ T=0
    // d epsl = P/rho0^2 d rho0  <==>  P = rho0^2 d(epsl)/d(rho0)
    //
    tab1d_rho0->dLogepsl[i]= P/(epsl*rho0);

    if(0 && PR)
      printf("%e %e %e  %e %e %e\n", rho0,epsl,P,
                              tab1d_rho0->Logrho0[i], tab1d_rho0->Logepsl[i],
                              tab1d_rho0->dLogepsl[i]);
    i++;
    if(i >= nrows) break;
  }

  /* close file */
  fclose (fp);

  /* check order of table: locate routine assumes monotonic rho */
  for(i=1; i<nrows; i++)
  {
    if((tab1d_rho0->Logrho0[i] - tab1d_rho0->Logrho0[i-1]) < 0.)
      errorexit("EoS table file must be ordered with monotonic rho0");
  }

  /* printf info */
  printf("  table loaded: tab1d_rho0->n=%d:\n", tab1d_rho0->n);
  printf("  table:  Logrho0       Logepsl       dLogepsl\n");
  i = 0;
  printf("%6d:  %+e %+e %+e\n", i,
         tab1d_rho0->Logrho0[i], tab1d_rho0->Logepsl[i],
         tab1d_rho0->dLogepsl[i]);
  printf("     |   ...\n");
  i = tab1d_rho0->n - 1;
  printf("%6d:  %+e %+e %+e\n", i,
         tab1d_rho0->Logrho0[i], tab1d_rho0->Logepsl[i],
         tab1d_rho0->dLogepsl[i]);
}


/*****************************************************************************/
/* interpolation routines with derivatives */

/* Locate the index of the nearest xval element in the array x */
int tab1d_locate(double *x, int Nx, double xval)
{
  int ju,jm,jl;
  int ascnd;

  jl=0;  /* if xval is NAN we return this value */
  ju=Nx;

  if(xval <= x[0])
  {
    if(xval < x[0])
      if(PR) printf("tab1d_locate warning: pt is outside (xval<x[0]).\n");
    return 0;
  }
  else if (xval >= x[Nx-1])
  {
    if(xval > x[Nx-1])
      if(PR) printf("tab1d_locate warning: pt is outside (xval>x[Nx-1]).\n");
    return Nx-1;
  }

  ascnd = (x[Nx-1] >= x[0]);

  while (ju-jl > 1)
  {
    jm = (ju+jl) >> 1;

    if ((xval >= x[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }

  return jl;
}

/* Hermite 2 pts interpolation with derivs */
int interp1d_her4(double *f, double *df, double *x, int Nx,
                  double xv, double *fv_p, double *dfv_p, double *ddfv_p)
{
  /* Given the values in xv, it returns the interpolated values fv and
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x)
  Hermite 2 pts interpolation is used */
  int ret = 0;

  if (Nx < 2) errorexit(" too few points for interpolation");

  int i = tab1d_locate(x,Nx,xv);

  if(i < 0)
  {
    if(PR)
      printf(" too few points on left => interpolation inaccurate! (rho=%e)\n",
             pow(10.,xv));
    i = 0;
    ret = -1;
  }
  if(i > (Nx-2))
  {
    if(PR)
      printf(" too few points on right => interpolation inaccurate! (rho=%e)\n",
             pow(10.,xv));
    i = Nx-2;
    ret = 1;
  }

  double xi    =  x[i];
  double fi    =  f[i];
  double dfi   = df[i];

  double xipo  =  x[i+1];
  double fipo  =  f[i+1];
  double dfipo = df[i+1];

  double w   = (xv-xi)/(xipo-xi);
  double w2  = w*w;
  double w3  = w2*w;
  double h0w = 2.*w3 - 3.*w2 + 1.;
  double h1w = w3 - 2.*w2 + w;

  double mw   = 1. - w;
  double mw2  = mw*mw;
  double mw3  = mw2*mw;
  double h0mw = 2.*mw3 - 3.*mw2 + 1.;
  double h1mw = mw3 - 2.*mw2 + mw;

  double dw = 1./(xipo-xi);

  *fv_p   = fi*h0w + fipo*h0mw + (xipo-xi)*(dfi*h1w - dfipo*h1mw);

  *dfv_p  = 6.*( fi*(w2-w) - fipo*(mw2-mw) )*dw
      + dfi*(3.*w2-4.*w+1.) + dfipo*(3.*mw2-4.*mw+1.);

  *ddfv_p = ( fi*(12.*w-6.) + fipo*(12.*mw-6.) )*dw*dw
      + ( dfi*(6.*w-4.) - dfipo*(6.*mw-4.) )*dw;

  return ret;
}

/* Lagrange 4 pts interpolation with derivs */
int interp1d_lag4(double *f, double *df, double *x, int Nx,
                   double xv, double *fv_p, double *dfv_p, double *ddfv_p)
{
  /* Given the values in xv, it returns the interpolated values fv and
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x)
  Lagrangian 4 pts interpolation is used */
  int ret = 0;

  if (Nx < 4) errorexit(" too few points for interpolation");

  int i = tab1d_locate(x,Nx,xv);

  if(i < 1)
  {
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = 1;
    ret = -1;
  }
  if(i > (Nx-3))
  {
    if (PR) printf(" too few points on the right => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = Nx-3;
    ret = 1;
  }

  double ximo =  x[i-1];
  double xi   =  x[i];
  double xipo =  x[i+1];
  double xipt =  x[i+2];

  double C1   = (f[i] - f[i-1])/(xi - ximo);
  double C2   = (-f[i] + f[i+1])/(-xi + xipo);
  double C3   = (-f[i+1] + f[i+2])/(-xipo + xipt);
  double CC1  = (-C1 + C2)/(-ximo + xipo);
  double CC2  = (-C2 + C3)/(-xi + xipt);
  double CCC1 = (-CC1 + CC2)/(-ximo + xipt);

  *fv_p   = f[i-1] + (-ximo + xv)*(C1 + (-xi + xv)*(CC1 + CCC1*(-xipo + xv)));
  *dfv_p  = C1 - (CC1 - CCC1*(xi + xipo - 2.*xv))*(ximo - xv)
      + (-xi + xv)*(CC1 + CCC1*(-xipo + xv));
  *ddfv_p = 2.*(CC1 - CCC1*(xi + ximo + xipo - 3.*xv));

  return ret;
}

/* Linear interpolation with derivs */
int interp1d_lag1(double *f, double *df, double *x, int Nx,
                   double xv, double *fv_p, double *dfv_p, double *ddfv_p)
{
  /* Given the values in xv, it returns the interpolated values fv and
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x)
  Linear interpolation is used (check and testing) */
  int ret = 0;

  if (Nx < 2) errorexit(" too few points for interpolation");

  int i;
  i = tab1d_locate(x,Nx,xv);

  if(i < 0)
  {
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = 0;
    ret = -1;
  }
  if(i > (Nx-2))
  {
    if (PR) printf(" too few points on the right => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = Nx-2;
    ret = 1;
  }

  double xi   =  x[i];
  double xipo =  x[i+1];

  double fi   =  f[i];
  double fipo =  f[i+1];

  *dfv_p  = (fipo-fi)/(xipo-xi);
  *fv_p   = fi + *dfv_p*(xv - xi);
  *ddfv_p = 0.;

  return ret;
}

/* wrapper for interp to get */
int EoS_interp_tab1d_rho0(double Logrho,
                    double *Logepsl, double *dLogepsl, double *d2Logepsl)
{
  int ret;
  ret=interp1d_her4(tab1d_rho0->Logepsl, tab1d_rho0->dLogepsl,
                    tab1d_rho0->Logrho0, tab1d_rho0->n,
                    Logrho, Logepsl, dLogepsl, d2Logepsl);
  /* interp1d_lag4(tab1d_rho0->Logepsl, tab1d_rho0->dLogepsl,
                   tab1d_rho0->Logrho0, tab1d_rho0->n,
                   Logrho, Logepsl, dLogepsl, d2Logepsl); */
  if(ret)
    ret=interp1d_lag1(tab1d_rho0->Logepsl, tab1d_rho0->dLogepsl,
                      tab1d_rho0->Logrho0, tab1d_rho0->n,
                      Logrho, Logepsl, dLogepsl, d2Logepsl);
  return ret;
}



/*****************************************************************************/
/* interpolate 1d EoS Table (Termodynamically consistent procedure) */
int tab1d_Of_rho0_AtT0(double rho0, double *epsl,
                       double *P, double *dPdrho0, double *dPdepsl)//double *cs2)
{
  double Logepsl, dLogepsl, d2Logepsl; /* tmp Log epsl and derivs */
  double E, DEDrho0, D2EDrho02;        /* tmp epsl and derivs */
  double Logrho0;

  /* take care of vacuum */
  if(rho0<=0.)
  {
    *epsl = *P = *dPdrho0 = *dPdepsl = 0;
    return 1;
  }

  Logrho0 = log10(rho0);
  EoS_interp_tab1d_rho0(Logrho0, &Logepsl, &dLogepsl, &d2Logepsl);

  E         = pow(10., Logepsl);// this is epsl from tab
  DEDrho0   = E*dLogepsl/(rho0);
  D2EDrho02 = (E*d2Logepsl/(ln10) +
              dLogepsl*(DEDrho0*(rho0)-E))/((rho0)*(rho0));

  *P       = (rho0)*(rho0)*DEDrho0;
  *epsl    = E; // set epsl here
  *dPdrho0 = 2.*(rho0)*DEDrho0 + (rho0)*(rho0)*D2EDrho02;
  *dPdepsl = 0.;

  /* can get sound speed^2 cs2 like this: */
  //*cs2 = cs2_of_rho0_epsl_P_dPdrho0_dPdepsl(rho0, E, *P, *dPdrho0, *dPdepsl);
  // BAM had last two args interchanged:
  //*cs2 = eos_cs2_rep(rho0,E, *P,*dPdepsl,*dPdrho0); // BUG in bam???

  if(0&PR) printf("%e %e -> %e %e %e -> %e %e %e\n",
                rho0,*epsl, E,DEDrho0,D2EDrho02, *P,*dPdrho0,*dPdepsl);
  if(!isfinite(rho0) || !isfinite(*epsl)) return 1;
  return 0;
}

/* same as tab1d_Of_rho0_AtT0 but operate with Log = log_10 */
int tab1d_Of_Logrho0_AtT0(double Logrho0, double *Logepsl,
                          double *LogP, double *dLogPdLogrho0)
{
  double epsl, P, dPdrho0, dPdepsl;
  double rho0 = pow(10., Logrho0);
  int ret;

  ret = tab1d_Of_rho0_AtT0(rho0, &epsl, &P, &dPdrho0, &dPdepsl);
  if(epsl>0.)
    *Logepsl = log10(epsl);
  else
    *Logepsl = -300;

  if(P>0.)
  {
    *LogP = log10(P);
    /* dLnP/dLnrho0 = dLnP/dP dP/dLnrho0 = (1/P) dP/drho0 drho0/dLnrho0
                    = (rho0/P) dP/drho0
       Log(x) = Ln(x)/Ln(10)
       dLogP/dLogrho0 = dLnP/dLnrho0 = (rho0/P) dP/drho0 */
    *dLogPdLogrho0 = dPdrho0 * rho0/P;
  }
  else
  {
    *LogP = -300;
    *dLogPdLogrho0 = 1e-6;
  }
  return ret;
}

/*****************************************************************************/
/* Root finders to get EoS vars from P or hm1=h-1 */

/* set P(rho0) - P0 in f and its rho0 deriv in df */
void Pofrho0_minus_P0(double rho0, void *par, double *f, double *df)
{
  double epsl, P, dPdrho0, dPdepsl;
  ttab1d_vars *p = (ttab1d_vars *) par;

  tab1d_Of_rho0_AtT0(rho0, &epsl, &P, &dPdrho0, &dPdepsl);
  *f  = P - p->P0;
  *df = dPdrho0;

  /* return vals also in par */
  //p->rho0 = rho0;
  //p->epsl = epsl;
  //p->P    = P;
  //p->dPdrho0 = dPdrho0;
  //p->dPdepsl = dPdepsl;

  if(PR)
  {
    PRF;printf(": P0=%g\n rho0=%g: P=%g -> *f=%g *df=%g\n",
               p->P0,rho0, P, *f, *df);
  }
}

/* set LogP(Logrho0) - LogP0 in f and its rho0 deriv in df */
void LogPofLogrho0_minus_LogP0(double Logrho0, void *par,
                               double *f, double *df)
{
  ttab1d_vars *p = (ttab1d_vars *) par;
  double Logepsl, LogP, dLogPdLogrho0;

  tab1d_Of_Logrho0_AtT0(Logrho0, &Logepsl, &LogP, &dLogPdLogrho0);
  *f  = LogP - p->LogP0;
  *df = dLogPdLogrho0;

  if(PR)
  {
    PRF;printf(": LogP0=%g\n Logrho0=%g: LogP=%g -> *f=%g *df=%g\n",
               p->LogP0, Logrho0, LogP, *f, *df);
  }
}

/* interpolate 1d EoS Table at P (needs root finder),
   *rho0 is init guess for rho0 */
int tab1d_Of_P_AtT0(double P, double *rho0, double *epsl,
                    double *dPdrho0, double *dPdepsl)
{
  double Logrho0, Logrho0min, Logrho0max, Lrho0min, Lrho0max;
  double Papprox;
  int ret;
  ttab1d_vars p[1];

  /* take care of vacuum */
  if(P<=0.)
  {
    *rho0 = *epsl = *dPdrho0 = *dPdepsl = 0;
    return 1;
  }

  /* set min and max rho0 */
  min_max_Logrho0_tab1d(&Logrho0min, &Logrho0max);
  Lrho0min = Logrho0min - 50.;
  Lrho0max = Logrho0max + 5.;

  /* find Logrho0 */
  if(*rho0 >= 0.)
    Logrho0 = log10(*rho0);
  else
    Logrho0 = 0.5*(Logrho0min+Logrho0max);
  p->LogP0 = log10(P);
  ret = newton1d_brak(&Logrho0, LogPofLogrho0_minus_LogP0,
                      Lrho0min,Lrho0max, p, 100, 1e-10, 1);
  if(PR) printf("newton1d_brak -> ret=%d\n", ret);
  if(ret<0) errorexit("newton1d_brak failed!");
  *rho0 = pow(10., Logrho0);

  return tab1d_Of_rho0_AtT0(*rho0,  epsl, &Papprox, dPdrho0, dPdepsl);

  /* use Pofrho0_minus_P0 to find root */
  //   /* set min and max rho0 */
  //   min_max_Logrho0_tab1d(&rho0min, &rho0max);
  //   rho0min = 0.;
  //   rho0max = pow(10., rho0max) * 1e3;
  //
  //   /* find rho0 */
  //   p->P0 = P;
  //   ret = newton1d_brak(rho0, Pofrho0_minus_P0, rho0min,rho0max, p,
  //                       100, 1e-10, 1);
  //   if(PR) printf("newton1d_brak -> ret=%d\n", ret);
  //   if(ret<0) errorexit("newton1d_brak failed!");
  //   return tab1d_Of_rho0_AtT0(*rho0,  epsl, &Papprox, dPdrho0, dPdepsl);
}


/* set h(rho0) - h0 in f and its rho0 deriv in df */
void hm1ofrho0_minus_hm10(double rho0, void *par, double *f, double *df)
{
  double epsl, P, dPdrho0, dPdepsl, hm1;
  ttab1d_vars *p = (ttab1d_vars *) par;

  tab1d_Of_rho0_AtT0(rho0, &epsl, &P, &dPdrho0, &dPdepsl);
  hm1 = hm1_of_rho0_epsl_P(rho0, epsl, P);

  *f  = hm1 - p->hm10;
  if(rho0>0.)
    *df = dPdrho0/rho0;
  else
    *df = 0;

  if(PR)
  {
    PRF;printf(": hm10=%g\n rho0=%g: hm1=%g -> *f=%g *df=%g\n",
               p->hm10, rho0, hm1, *f, *df);
  }
}

/* set Loghm1(Logrho0) - Loghm10 in f and its rho0 deriv in df */
void Loghm1ofLogrho0_minus_Loghm10(double Logrho0, void *par,
                                   double *f, double *df)
{
  ttab1d_vars *p = (ttab1d_vars *) par;
  double Logepsl, LogP, dLogPdLogrho0;
  double rho0, epsl, P, hm1, Loghm1, dLoghm1dLogrho0;

  tab1d_Of_Logrho0_AtT0(Logrho0, &Logepsl, &LogP, &dLogPdLogrho0);
  rho0 = pow(10., Logrho0);
  epsl = pow(10., Logepsl);
  P    = pow(10., LogP);
  hm1  = hm1_of_rho0_epsl_P(rho0, epsl, P);
  Loghm1 = log10(hm1);

  /* dLnP/dLnrho0 = dLnP/dP dP/dLnrho0 = (1/P) dP/drho0 drho0/dLnrho0
                  = (rho0/P) dP/drho0
     Log(x) = Ln(x)/Ln(10)
     dLogP/dLogrho0 = dLnP/dLnrho0 = (rho0/P) dP/drho0 */
  /* dh = dP/rho0  ==>  dh/drho0 = (dP/drho0)/rho0
     ==> d(h-1)/drho0 = (dP/drho0)/rho0
     dLoghm1/dLogrho0 = (rho0/hm1) dhm1/drho0 = (rho0/hm1) (dP/drho0)/rho0
                      = (1/hm1) (P/rho0) dLogP/dLogrho0 */
  dLoghm1dLogrho0 = dLogPdLogrho0 * P/(hm1*rho0);

  *f  = Loghm1 - p->Loghm10;
  *df = dLoghm1dLogrho0;

  if(PR)
  {
    PRF;printf(": Loghm10=%g\n Logrho0=%g: Loghm1=%g -> *f=%g *df=%g\n",
               p->Loghm10, Logrho0, Loghm1, *f, *df);
  }
}

/* interpolate 1d EoS Table at hm1=h-1 (needs root finder),
   *rho0 is init guess for rho0 */
int tab1d_Of_hm1_AtT0(double hm1, double *rho0, double *epsl,
                      double *P, double *dPdrho0, double *dPdepsl)
{
  double Logrho0, Logrho0min, Logrho0max, Lrho0min, Lrho0max;
  int ret;
  ttab1d_vars p[1];

  /* take care of vacuum */
  if(hm1<=0.)
  {
    *rho0 = *epsl = *P = *dPdrho0 = *dPdepsl = 0;
    return 1;
  }

  /* set min and max rho0 */
  min_max_Logrho0_tab1d(&Logrho0min, &Logrho0max);
  Lrho0min = Logrho0min - 50.;
  Lrho0max = Logrho0max + 5.;

  /* find Logrho0 */
  if(*rho0 >= 0.)
    Logrho0 = log10(*rho0);
  else
    Logrho0 = 0.5*(Logrho0min+Logrho0max);
  p->Loghm10 = log10(hm1);
  ret = newton1d_brak(&Logrho0, Loghm1ofLogrho0_minus_Loghm10,
                      Lrho0min,Lrho0max, p, 100, 1e-10, 1);
  if(PR) printf("newton1d_brak -> ret=%d\n", ret);
  if(ret<0) errorexit("newton1d_brak failed!");
  *rho0 = pow(10., Logrho0);

  return tab1d_Of_rho0_AtT0(*rho0, epsl, P, dPdrho0, dPdepsl);

  //   /* set min and max rho0 */
  //   min_max_Logrho0_tab1d(&rho0min, &rho0max);
  //   rho0min = 0.;
  //   rho0max = pow(10., rho0max) * 1e3;
  //
  //   /* find rho0 */
  //   p->hm10 = hm1;
  //   ret = newton1d_brak(rho0, hm1ofrho0_minus_hm10, rho0min,rho0max, p,
  //                       100, 1e-10, 1);
  //   if(PR) printf("newton1d_brak -> ret=%d\n", ret);
  //   if(ret<0) errorexit("newton1d_brak failed!");
  //   return tab1d_Of_rho0_AtT0(*rho0, epsl, P, dPdrho0, dPdepsl);
}

/* we also need drho0/dh or dh/drho0: *rho0 is init guess for rho0 */
void tab1d_rho0_epsl_P_drho0dhm1_Of_hm1_AtT0(double hm1,
                                             double *rho0, double *epsl,
                                             double *P, double *drho0dhm1)
{
  double dPdrho0, dPdepsl;
  tab1d_Of_hm1_AtT0(hm1, rho0, epsl, P, &dPdrho0, &dPdepsl);

  /* dh = dP/rho0  ==>  dh/drho0 = (dP/drho0)/rho0
     ==> d(h-1)/drho0 = (dP/drho0)/rho0
     ==> drho0/d(h-1) = rho0 / (dP/drho0) */
  if(dPdrho0>0.)
    *drho0dhm1 = (*rho0) / (dPdrho0);
  else
    *drho0dhm1 = 0.;
}
