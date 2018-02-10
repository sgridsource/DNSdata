/* grid_setup.c */
/* Wolfgang Tichy 2008 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow
#define pow2(x)    ((x)*(x))

/* struct types used in root finder newton_linesrch_itsP */
typedef struct T_grid_b_struct {
  tGrid *grid; /* grid */
  int b;      /* box1 */
} t_grid_b_struct;

typedef struct T_grid_box_A_B_icoeffs_innerdom_outerdom_struct {
  tGrid *grid; /* grid */
  tBox *box;   /* box */
  double A;         /* A for q_of_sig_forgiven_ABP */
  double B;         /* B for q_of_sig_forgiven_ABP */
  int    icoeffs;   /* icoeffs for q_of_sig_forgiven_ABP */
  int    innerdom;  /* inner domain for q_of_sig_forgiven_ABP */
  int    outerdom;  /* outer domain for q_of_sig_forgiven_ABP */
} t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct;

/* global vars in this file */
double rf_surf1; /* radius of star1 */
double rf_surf2; /* radius of star2 */
double P_core1;  /* core pressure of star1 */
double P_core2;  /* core pressure of star2 */
/* global vars in this file for root finding */
tBox *DNSdata_q_VectorFunc_box;      /* box for BNdata_q_VectorFunc */
double *DNSdata_q_VectorFunc_coeffs; /* coeffs for BNdata_q_VectorFunc */
double DNSdata_q_VectorFunc_A;       /* B for DNSdata_q_VectorFunc */
double DNSdata_q_VectorFunc_B;     /* phi for DNSdata_q_VectorFunc */
/* global vars in this file for minimization with numrec's powell */
t_grid_b_struct *grid_b_PARS;

/* funs in this file */
void m01_VectorFunc(int n, double *vec, double *fvec);
void m02_VectorFunc(int n, double *vec, double *fvec);
void minimize_dsigma_pm_dB_1_ByAdjusting_sigp_1phi(tGrid *grid, int innerdom);
void set_dsigma_pm_dB_toZero_atB01(tGrid *grid, int innerdom);
void DNSdata_LowPassFilter_with_dsigma_pm_dB_01_EQ0(tGrid *grid, int innerdom);
void DNSdata_LowPassFilter_with_dsigma_pm_dBphi_01_EQ0(tGrid *grid,int innerdom);




/* setup initial boxsizes */
int set_boxsizes(tGrid *grid)
{
  double m1, Phic1, Psic1;
  double m2, Phic2, Psic2;
  double kappa     = Getd("DNSdata_kappa");
  double DNSdata_n = Getd("DNSdata_n");
  double DNSdata_b = Getd("DNSdata_b");
  double m01 = Getd("DNSdata_m01");
  double m02 = Getd("DNSdata_m02");
  double xmin1,xmax1, xmin2,xmax2, xc1, xc2; /* x-positions of stars */
  double DoM, nu; /* distance over total rest mass, and rest mass ratio */
  double DoM3, DoM4, DoM5; /* powers of DoM */
  double xCM, Omega;
  double Fc, qc, Cc, oouzerosqr;
  double vec[2];
  double fvec[2];
  int check, stat;

  printf("set_boxsizes: setting box sizes and coordinates used ...\n");
  prTimeIn_s("WallTime: ");

  /* reset initial DNSdata_m01/2 if needed */
  if(Getv("DNSdata_iterate_m0", "yes"))
  {
    printf(" DNSdata_iterate_m0 = yes : setting:\n");
    /* set DNSdata_desired_m01/2 if needed */
    if(Getd("DNSdata_desired_m01")<0.0) Setd("DNSdata_desired_m01", m01);
    if(Getd("DNSdata_desired_m02")<0.0) Setd("DNSdata_desired_m02", m02);
    if(Getd("DNSdata_desired_kappa")<0.0) Setd("DNSdata_desired_kappa", kappa);
    printf("   DNSdata_desired_m01 = %g\n", Getd("DNSdata_desired_m01"));
    printf("   DNSdata_desired_m02 = %g\n", Getd("DNSdata_desired_m02"));
    printf("   DNSdata_desired_kappa = %g\n", Getd("DNSdata_desired_kappa"));
    /* set initial DNSdata_m01/2 */
    m01 = Getd("DNSdata_init_m01");
    m02 = Getd("DNSdata_init_m02");
    Setd("DNSdata_m01", m01);
    Setd("DNSdata_m02", m02);
    printf("   DNSdata_m01 = DNSdata_init_m01 = m01 = %g\n", m01);
    printf("   DNSdata_m02 = DNSdata_init_m02 = m02 = %g\n", m02);
    printf("   DNSdata_m0change = %g\n", Getd("DNSdata_m0change"));
  }

  /* for TOV we need an EoS so we need to init EoS here */
  if(Getv("DNSdata_EoS_type", "pwp"))
  {     
    if(!Getv("DNSdata_EoS_file", " ")) DNS_pwp_init_file();
    else DNS_pwp_init_parameter();
  }
  else if(Getv("DNSdata_EoS_type", "poly"))  DNS_poly_init();
  else errorexit("unkown DNSdata_EoS_type");

  /* do we set star from mass m01 or from max q called qm1 */
  if(Getd("DNSdata_qm1")<0.0)
  {
    /* many parts of the code only work with m01>0 */
    if(Getd("DNSdata_m01")<=0.0) errorexit("make sure DNSdata_m01 > 0");

    /* find P_core1, s.t. rest mass is m01 */
    printf("find P_core1 of a TOV star, s.t. rest mass is m01=%g\n",
           Getd("DNSdata_m01"));
    vec[1] = 1e-7;   /* initial guess */
    /* do newton_linesrch_its iterations: */
    stat=newton_linesrch_its(vec, 1, &check, m01_VectorFunc, 
   		           Geti("Coordinates_newtMAXITS"),
   		           Getd("Coordinates_newtTOLF") );
    if(check || stat<0)
    {
      printf("WARNING from newton_linesrch_its with MAXITS=%d TOLF=%g:\n",
             Geti("Coordinates_newtMAXITS"), Getd("Coordinates_newtTOLF"));
      printf("  check=%d stat=%d\n", check, stat);
    }
    /* check if we found correct vec or just a local max in m0 */
    m01_VectorFunc(1, vec, fvec);
    if(fabs(fvec[1])>Getd("Coordinates_newtTOLF")*1000)
    {
      int i;
      double v[2], f[2];
      printf("There may be a max in m0 near Pc=%g\n", vec[1]);
      for(i=-15; i<=15; i++)
      { 
        v[1] = vec[1] + i*vec[1]/1000;
        m01_VectorFunc(1, v, f);
      }
    }
    P_core1 = vec[1];
  }
  else  /* get P_core1 from max q */
  {
    double rho0, rhoE, drho0dhm1;
    DNS_polytrope_EoS_of_hm1(Getd("DNSdata_qm1"),
                             &rho0, &P_core1, &rhoE, &drho0dhm1);
  }
  printf("setting: P_core1=%g\n", P_core1);

  /* TOV_init yields m01 for a given P_core1 */
  TOV_init(P_core1, 1, &rf_surf1, &m1, &Phic1, &Psic1, &m01);
  printf(" rf_surf1=%g: m1=%g Phic1=%g Psic1=%g m01=%g\n",
         rf_surf1, m1, Phic1, Psic1, m01);
  if(Getd("DNSdata_qm1")>0.0)
  {
    Setd("DNSdata_m01", m01);
    printf("setting: DNSdata_m01=%g\n", Getd("DNSdata_m01"));
  }

  if(Getd("DNSdata_m02")>0.0 && Getd("DNSdata_qm2")<0.0)
  {
    /* find P_core2, s.t. rest mass is m02 */
    printf("find P_core2 of a TOV star, s.t. rest mass is m02=%g\n",
           Getd("DNSdata_m02"));
    vec[1] = 1e-7;   /* initial guess */
    /* do newton_linesrch_its iterations: */
    stat=newton_linesrch_its(vec, 1, &check, m02_VectorFunc, 
                             Geti("Coordinates_newtMAXITS"),
                             Getd("Coordinates_newtTOLF") );
    if(check || stat<0)
    {
      printf("WARNING from newton_linesrch_its with MAXITS=%d TOLF=%g:\n",
             Geti("Coordinates_newtMAXITS"), Getd("Coordinates_newtTOLF"));
      printf("  check=%d stat=%d\n", check, stat);
    }
    /* check if we found correct vec or just a local max in m0 */
    m02_VectorFunc(1, vec, fvec);
    if(fabs(fvec[1])>Getd("Coordinates_newtTOLF")*1000)
    {
      int i;
      double v[2], f[2];
      printf("There may be a max in m0 near Pc=%g\n", vec[1]);
      for(i=-15; i<=15; i++)
      { 
        v[1] = vec[1] + i*vec[1]/1000;
        m02_VectorFunc(1, v, f);
      }
    }
    P_core2 = vec[1];
    printf("setting: P_core2=%g\n", P_core2);

    /* TOV_init yields m02 for a given P_core2 */
    TOV_init(P_core2, 1, &rf_surf2, &m2, &Phic2, &Psic2, &m02);
    printf(" rf_surf2=%g: m2=%g Phic2=%g Psic2=%g m02=%g\n",
           rf_surf2, m2, Phic2, Psic2, m02);
  }
  else if(Getd("DNSdata_qm2")>0.0) /* get P_core2 from max q */
  {
    double rho0, rhoE, drho0dhm1;
    DNS_polytrope_EoS_of_hm1(Getd("DNSdata_qm2"),
                             &rho0, &P_core2, &rhoE, &drho0dhm1);
    printf("setting: P_core2=%g\n", P_core2);
    /* TOV_init yields m02 for a given P_core2 */
    TOV_init(P_core2, 1, &rf_surf2, &m2, &Phic2, &Psic2, &m02);
    printf(" rf_surf2=%g: m2=%g Phic2=%g Psic2=%g m02=%g\n",
           rf_surf2, m2, Phic2, Psic2, m02);
    Setd("DNSdata_m02", m02);
    printf("setting: DNSdata_m02=%g\n", Getd("DNSdata_m02"));
  }
  else 
  { P_core2=0.0; m2=m02=Phic2=0.0; Psic2=1.0; rf_surf2=rf_surf1; }

  /* Note: sigma1/2 are simply the radii rf_surf1/2 */
  printf(" found: rf_surf1=%g\n", rf_surf1);
  printf(" found: rf_surf2=%g\n", rf_surf2);

  /* find  xmin1,xmax1, xmin2,xmax2, xc1, xc2, and CM for both stars */
  xmax1 = DNSdata_b + rf_surf1;
  xmin1 = DNSdata_b - rf_surf1;
  xmax2 = -DNSdata_b + rf_surf2;
  xmin2 = -DNSdata_b - rf_surf2;
  printf(" => radii of star domains = %g, %g\n", 
         0.5*(xmax1-xmin1), 0.5*(xmax2-xmin2) );
  xc1 = 0.5*(xmax1+xmin1);
  xc2 = 0.5*(xmax2+xmin2);
  printf(" => 1st domain extends from  xin1 = %g  to  xout1 = %g\n", xmin1,xmax1);
  printf(" => 2nd domain extends from  xout2 = %g  to  xin2 = %g\n", xmin2,xmax2);
  printf(" => centers are at: xc1 = %.13g,  xc2 = %.13g\n", xc1, xc2);
  printf(" DNSdata_b = %g\n", Getd("DNSdata_b"));
  DoM = fabs(xc1-xc2)/(m01+m02);
  DoM3 = DoM*DoM*DoM;
  DoM4 = DoM3*DoM;
  DoM5 = DoM4*DoM;
  nu = (m01*m02)/pow(m01+m02, 2.0);

  /* set: DNSdata_x_CM  = DNSdata_x_CM_init 
          DNSdata_Omega = DNSdata_Omega_init
     (unless: DNSdata_x_CM_init  = DNSdata_x_CM
              DNSdata_Omega_init = DNSdata_Omega) */
  if(!Getv("DNSdata_x_CM_init", "DNSdata_x_CM"))
    Sets("DNSdata_x_CM", Gets("DNSdata_x_CM_init"));
  if(!Getv("DNSdata_Omega_init", "DNSdata_Omega"))
    Sets("DNSdata_Omega", Gets("DNSdata_Omega_init"));
  printf(" initializing DNSdata_x_CM and DNSdata_Omega using:\n");
  printf(" DNSdata_x_CM  = %s\n", Gets("DNSdata_x_CM"));
  printf(" DNSdata_Omega = %s\n", Gets("DNSdata_Omega"));

  /* set CM and Omega (taken from PN_ADM_2.m) */
  if(Getv("DNSdata_x_CM", "estimate"))
    xCM = (m01*xc1 + m02*xc2)/(m01+m02);
  else
    xCM = Getd("DNSdata_x_CM");
  if(Getv("DNSdata_Omega", "estimate"))
  {
    Omega = sqrt( 64*DoM3/pow(1 + 2*DoM, 6) +nu/DoM4 +
                 (-5*nu + 8*nu*nu)/(8*DoM5)            )/(m01+m02);
    if(nu<=0.0) Omega=0.0;
  }
  else if(Getv("DNSdata_Omega", "estimate_from_desired_m0"))
  {
    double m1 = Getd("DNSdata_desired_m01");
    double m2 = Getd("DNSdata_desired_m02");
    DoM  = fabs(xc1-xc2)/(m1+m2);
    DoM3 = DoM*DoM*DoM;
    DoM4 = DoM3*DoM;
    DoM5 = DoM4*DoM;
    nu = (m1*m2)/pow(m1+m2, 2.0);
    Omega = sqrt( 64*DoM3/pow(1 + 2*DoM, 6) +nu/DoM4 +
                 (-5*nu + 8*nu*nu)/(8*DoM5)            )/(m01+m02);
    if(nu<=0.0) Omega=0.0;
  }
  else
    Omega = Getd("DNSdata_Omega");
  Setd("DNSdata_x_CM", xCM);
  Setd("DNSdata_Omega", Omega);
  printf(" DNSdata_x_CM and DNSdata_Omega are now set to:\n");
  printf(" DNSdata_x_CM = %g\n", Getd("DNSdata_x_CM"));
  printf(" DNSdata_Omega = %g,  (m01+m02)*DNSdata_Omega = %g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_Omega")*(m01+m02));

  /* set DNSdata_C1 */
  qc = DNS_polytrope_hm1_of_P(P_core1);
  /* oouzerosqr == alpha2 - 
                   Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]), */
  oouzerosqr = exp(Phic1*2.0) - 
               pow(Psic1+m02/(2*fabs(xc1-xc2)), 4) *
               Omega*Omega*(xc1-xCM)*(xc1-xCM);
  if(oouzerosqr<0.0)  oouzerosqr = exp(Phic1*2.0);
  Fc = -sqrt(oouzerosqr);
  /* q == (C1/F - 1.0) */
  Cc = Fc*(qc + 1.0);
  Setd("DNSdata_C1", Cc);
  printf(" DNSdata_C1 = %g\n", Getd("DNSdata_C1"));

  /* set DNSdata_C2 */
  qc = DNS_polytrope_hm1_of_P(P_core2);
  /* oouzerosqr == alpha2 - 
                   Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]), */
  oouzerosqr = exp(Phic2*2.0) - 
               pow(Psic2+m01/(2*fabs(xc1-xc2)), 4) *
               Omega*Omega*(xc2-xCM)*(xc2-xCM);
  if(oouzerosqr<0.0)  oouzerosqr = exp(Phic2*2.0);
  Fc = -sqrt(oouzerosqr);
  /* q == (C2/F - 1.0) */
  Cc = Fc*(qc + 1.0);
  if(m02<=0.0) Cc = 0.0;
  Setd("DNSdata_C2", Cc);
  printf(" DNSdata_C2 = %g\n", Getd("DNSdata_C2"));

  /* set some box pars */
  if(Getv("DNSdata_grid", "36CS_2xyz"))
  {
    double dc = DNSdata_b; 
    /* Sets("nboxes", "38");  // <--cannot be set here, needs to be set earlier */
    printf("using: DNSdata_grid = 36CS_2xyz\n");
    two_spheres_around_two_full_cubes(grid, 0, dc,
                                      0.5*rf_surf1, rf_surf1,
                                      0.5*rf_surf2, rf_surf2,
                                      4.0*dc, 10000.0*dc);
  }

  prTimeIn_s("WallTime: ");
  return 0;
}


/* set attributes in boxes */
int set_DNS_box_attribs(tGrid *grid)
{
  int b;
  //printf("set_DNS_box_attribs:\n");
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    char str[1000];
    double x = 0.0;

    /* left or right? */
    if(box->x_of_X[1]!=NULL)
      x = box->x_of_X[1]((void *)box, -1, 0.5,0.0,0,0);
    if(x>0.0)    box->SIDE = STAR1;
    else (x<0.0) box->SIDE = STAR2;
    else         box->SIDE = ZERO;

    /* default is that there is no star surface */
    box->BOUND = ZERO;

    /* check all coords */
    snprintf(str, 999, "box%d_Coordinates", b);
    if( Getv(str, "Cartesian") )
    {
      double xm = 0.5*(box->bbox[0]+box->bbox[1]);

      box->COORD = CART;
      if(xm>=0.0)  box->SIDE = STAR1;
      else         box->SIDE = STAR2;
      /* if Cartesian it must be inside a star */
      box->MATTR = INSIDE;
    }
    if( Getv(str, "outerCubedSphere") )
    {
      double y, z, r, rs;
      box->COORD = CUBSPH;
      y = box->x_of_X[2]((void *)box, -1, 0.5,0.0,0,0);
      z = box->x_of_X[3]((void *)box, -1, 0.5,0.0,0,0);
      if(box->SIDE == STAR1) { r = x-xmax1;  rs = rf_surf1; }
      else                   { r = x-xmax2;  rs = rf_surf2; }
      r = sqrt(r*r + y*y + z*z)
      if(r<rs)
      {
        box->BOUND = SSURF;   /* star surface */
        box->MATTR = INSIDE;
      }
      else
        box->MATTR = AWAY;
    }
    if( Getv(str, "innerCubedSphere") )
    { 
      double y, z, r, rs;
      box->COORD = CUBSPH;
      x = box->x_of_X[1]((void *)box, -1, 0.0,0.0,0,0);
      y = box->x_of_X[2]((void *)box, -1, 0.0,0.0,0,0);
      z = box->x_of_X[3]((void *)box, -1, 0.0,0.0,0,0);
      if(box->SIDE == STAR1) { r = x-xmax1;  rs = rf_surf1; }
      else                   { r = x-xmax2;  rs = rf_surf2; }
      r = sqrt(r*r + y*y + z*z)
      if(dlesseq(r,rs))
      {
        box->BOUND = SSURF;   /* star surface */
        box->MATTR = TOUCH;
      }
      else
        box->MATTR = AWAY;
    }
    if( Getv(str, "PyramidFrustum") || Getv(str, "CubedShell") )
    { 
      box->COORD = CUBSPH;
      box->MATTR = AWAY;
    }
    if( Getv(str, "stretchedCubedShell") )
    { 
      box->COORD = SCUBSH;
      box->MATTR = AWAY;
    }
  }
  return 0;
}

/* print attributes in DNS boxes */
int pr_DNS_box_attribs(tGrid *grid)
{
  int b;
  printf("pr_DNS_box_attribs:\n");
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    printf("  b=%d:  box  ->SIDE=%d  ->MATTR=%d  ->BOUND=%d  ->COORD=%d\n",
           b, box->SIDE, box->MATTR, box->BOUND, box->COORD);
  }
  return 0;
}


/* funtion to be passed into newton_linesrch_its to find P_core1 from m01 */
void m01_VectorFunc(int n, double *vec, double *fvec)
{
  double m01 = Getd("DNSdata_m01");
  double Pc, m, Phic, Psic, m0;

  Pc = vec[1];
  TOV_init(Pc, 0, &rf_surf1, &m, &Phic, &Psic, &m0);
  printf("   Pc=%+.16e  m0=%.16e\n",Pc,m0);
  fvec[1] = m0-m01;
}
/* funtion to be passed into newton_linesrch_its to find P_core2 from m02 */
void m02_VectorFunc(int n, double *vec, double *fvec)
{
  double m02 = Getd("DNSdata_m02");
  double Pc, m, Phic, Psic, m0;

  Pc = vec[1];
  TOV_init(Pc, 0, &rf_surf2, &m, &Phic, &Psic, &m0);
  printf("   Pc=%+.16e  m0=%.16e\n",Pc,m0);
  fvec[1] = m0-m02;
}



/*****************************/
/* functions to to avoid q<0 */
/*****************************/

/* set q (or other Var) to zero if q<0 or outside stars */
void set_Var_to_Val_if_below_limit_or_outside(tGrid *grid, int vi, 
                                              double Val, double lim)
{
  int b;
  forallboxes(grid, b)
  {
    int i;
    tBox *box = grid->box[b];
    double *Var = box->v[vi];
    forallpoints(box, i)
      if( Var[i]<lim || box->MATTR!=INSIDE )  Var[i] = Val;
  }
}
/* set q (or other Var) to zero at A=0 */
void set_Var_to_Val_atSurface(tGrid *grid, int vi, double Val)
{
  int b;
errorexit("determine boxes to loop over")
  for(b=0; b<=3; b++)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    double *Var = box->v[vi];
    int j,k;

    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)  
      Var[Index(0,j,k)] = Val;
  }
}


/*************************************/
/* functions to adjust star surfaces */
/*************************************/

/* get q at lam by direct computation along a (A,B)=(const1,const2) line */
void DNSdata_q_VectorFunc(int n, double *vec, double *fvec)
{
  tBox *box = DNSdata_q_VectorFunc_box;
  tGrid *grid = box->grid;
  int b = box->b;
  double A = DNSdata_q_VectorFunc_A;
  double B = DNSdata_q_VectorFunc_B;
  double lam = vec[1];  

  /* compute q */
  fvec[1] = DNS_compute_new_centered_q_atXYZ(grid,b, lam,A,B);
}


/* WE NEED to find sigma at A,B such that q(sigma; lam=0or1, A,B) = 0 */
/* q as a func of sig for a given lam=0or1, A, B */
void q_of_sig_forgiven_ABP(int n, double *sigvec, double *qvec, void *p)
{
  double sigp_Bphi = sigvec[1];
  /* get pars */
  t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *pars = 
         (t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct *) p;
  double sigp_1phi = pars->sigp_1phi;
  double B         = pars->B;
  double phi       = pars->phi;
  tGrid *grid      = pars->grid;
  int innerdom     = pars->innerdom;
  int outerdom     = pars->outerdom;
  double AbsCp_Bphi, ArgCp_Bphi, ReCp_Bphi, ImCp_Bphi;
  double AbsCp_1phi, ArgCp_1phi, ReCp_1phi, ImCp_1phi;
  double X,R;
  double Ac,Bc, Acin,Bcin, Acout,Bcout, Acmax, q;
  double vec[3];
  int i, check, stat,statin,statout, dom;
  int guessmode;

  /* check if sigp_Bphi is within range */
  if( (innerdom==0 && sigp_Bphi>=0.0) || (innerdom==3 && sigp_Bphi<=0.0) )
  {
    /* set Cp at B,phi and 1,phi */
    AbsCp_Bphi = sqrt( Abstanh(0.25*sigp_Bphi, 0.25*PI*B) );
    ArgCp_Bphi = 0.5 * Argtanh(0.25*sigp_Bphi, 0.25*PI*B);
    ReCp_Bphi = AbsCp_Bphi * cos(ArgCp_Bphi);
    ImCp_Bphi = AbsCp_Bphi * sin(ArgCp_Bphi);
    AbsCp_1phi = sqrt( Abstanh(0.25*sigp_1phi, 0.25*PI) );
    ArgCp_1phi = 0.5 * Argtanh(0.25*sigp_1phi, 0.25*PI);
    ReCp_1phi = AbsCp_1phi * cos(ArgCp_1phi);
    ImCp_1phi = AbsCp_1phi * sin(ArgCp_1phi);

    /* use Eq. (22), (23) or (24) at A=0 to compute X,R */
    X = ReCp_Bphi - B*ReCp_1phi + B*cos(ArgCp_1phi);
    R = ImCp_Bphi - B*ImCp_1phi + B*sin(ArgCp_1phi);

    /* try 4 ways of finding an initial guess for Ac,Bc */
    for(guessmode=0; guessmode<=3; guessmode++)
    {
      /* set Acin,Bcin, Acout,Bcout, statin,statout to invalid values */
      statin=statout=-1;
      Acin=Acout=Bcin=Bcout=-1.0;
    
      /* find domain and Ac,Bc on current grid, corresponding to X,R */
      dom = innerdom;
      for(i=1; i<=2; i++)
      {
        t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];
        pars->phi = phi;
        pars->box = grid->box[dom];
        pars->X   = X;
        pars->R   = R;
        Acmax  = grid->box[dom]->bbox[1];
    
        /* get initial guess for Ac,Bc in vec[1],vec[2] */
        if(guessmode==0)
        {
          vec[1] = 1e-7; /* initial guess is that Ac,Bc = 0,B */
          vec[2] = B;
        }
        else if(guessmode==1)
        {
          find_nearest_A_B_given_X_R_phi(grid->box[dom], X,R,phi, &Ac,&Bc);
          if(dlesseq(Ac,0.0)) Ac=1e-7;
          if(dequal(Bc,0.0) || dequal(Bc,1.0)) Bc=B;
          if(dequal(B,0.0))   Bc=0.0;
          vec[1] = Ac;
          vec[2] = Bc;
        }
        else if(guessmode==2)
        {
          vec[1] = 0.01;  /* initial guess is that Ac,Bc = 0.01,B */
          vec[2] = B;
        }
        else
        {
          vec[1] = 0.1;   /* initial guess is that Ac,Bc = 0.1,B */
          vec[2] = B;
        }
//printf("q_of_sig_forgiven_ABP: guessmode=%d (A,B)=(%g,%g)\n",
//guessmode, vec[1],vec[2]);
        /* do newton line searches */
        if(dequal(B,0.0))
          stat = newton_linesrch_itsP(vec, 1, &check, DelXR_of_A_forB0_VectorFuncP,
                                      (void *) pars, 1000, 1e-10);
        else
          stat = newton_linesrch_itsP(vec, 2, &check, DelXR_of_AB_VectorFuncP,
                                      (void *) pars, 1000, 1e-10);
        if(check) printf("q_of_sig_forgiven_ABP: check=%d\n", check);  
        Ac = vec[1];
        Bc = vec[2];
    
        /* save vals for later */
        if(dom == innerdom) {  Acin=Ac;  Bcin=Bc;  statin=stat; }
        else                { Acout=Ac; Bcout=Bc; statout=stat; }
//printf("  sigp_Bphi=%g sigp_1phi=%g:\n"
//"       X=%g R=%g: stat=%d dom=%d Ac=%g Bc=%g\n",
//sigp_Bphi,sigp_1phi, X,R, stat,dom, Ac,Bc);
        /* if(stat>=0 && Ac>=0.0 && Ac<=Acmax && Bc>=0.0 && Bc<=1.0) break; */
        if(stat>=0 && dlesseq(0.0,Ac) && dlesseq(Ac,1.0) &&
                      dlesseq(0.0,Bc) && dlesseq(Bc,1.0)   ) break;
        dom = outerdom;
      }
      /* decide which results to use */
      if(dom == outerdom && statin>=0)
      {
        double dA=fabs(Acout)-fabs(Acin);
        double dB=fabs(Bcout)-fabs(Bcin);
        /* switch back to innerdom in some cases */
        if(Ac<0.0)
          if( dA>0.0 && dlesseq(0.0,Bcin) && dlesseq(Bcin,1.0) )
          {dom=innerdom; Ac=0.0; Bc=Bcin; stat=statin;}
        if(Bc<0.0)
          if( dB>0.0 && dlesseq(0.0,Acin) && dlesseq(Acin,1.0) )
          {dom=innerdom; Ac=Acin; Bc=0.0; stat=statin;}
        if(Bc>1.0)
          if(fabs(Bcin)<fabs(Bcout) && dlesseq(0.0,Acin) && dlesseq(Acin,1.0))
          {dom=innerdom; Ac=Acin; Bc=1.0; stat=statin;}
      }
      /* leave guessmode loop if we found a resonable result */
      if(stat>=0 && dgreatereq(Ac,0.0) && dgreatereq(1.0,Ac) &&
                    dgreatereq(Bc,0.0) && dgreatereq(1.0,Bc)   )  break;
    }
    /* check for failure */
    if(stat<0 || dless(Ac,0.0) || dless(1.0,Ac) ||
                 dless(Bc,0.0) || dless(1.0,Bc)   )
    {
      printf("q_of_sig_forgiven_ABP: guessmode=%d\n", guessmode);
      printf("q_of_sig_forgiven_ABP: stat=%d dom=%d Ac=%g Bc=%g\n",
             stat,dom, Ac,Bc);
      printf("statin=%d Acin=%g Bcin=%g  statout=%d Acout=%g Bcout=%g\n",
             statin, Acin,Bcin, statout, Acout,Bcout);
      printf("q_of_sig_forgiven_ABP: X=%g R=%g for: A=0 B=%g phi=%g\n"
             " sigp_Bphi=%g sigp_1phi=%g\n"
             " ReCp_Bphi=%g ImCp_Bphi=%g ReCp_1phi=%g ImCp_1phi=%g\n",
             X,R, B,phi, sigp_Bphi,sigp_1phi,
             ReCp_Bphi,ImCp_Bphi, ReCp_1phi,ImCp_1phi);
      errorexit("q_of_sig_forgiven_ABP: could not find Ac,Bc");
    }
    if(Ac<0.0) Ac=0.0; /* make sure we stay in our box */
    if(Ac>1.0) Ac=1.0;
    if(Bc<0.0) Bc=0.0;
    if(Bc>1.0) Bc=1.0;

    /* if we get point in innerdom with Acmax<Ac<1 we retrun a huge value
       for q, so that newton_linesrch_its thinks it has to backtrack... */
    if(Ac>Acmax && dom==0)      q = 1e6;
    else if(Ac>Acmax && dom==3) q = 1e6;
    else
      /* obtain q at Ac,Bc,phi by direct computation */
      /* grid->box[dom]->v[icoeffs]  contains coeffs of q in box */
      q = DNS_compute_new_centered_q_atXYZ(grid,dom, Ac,Bc,phi);
    qvec[1] = q;
  }
  else /* sigp_Bphi not in range */
  {
    /* make q negative so that we can backtrack */
    if(innerdom==0) q = +(sigp_Bphi-sigp_1phi)*1e6;
    else            q = -(sigp_Bphi-sigp_1phi)*1e6;
    qvec[1] = q;
    printf("q_of_sig_forgiven_ABP: innerdom=%d (B,phi)=(%g,%g)",
          innerdom, B,phi);
    printf(" sigp_1phi=%g sigp_Bphi=%g out of range!\n", sigp_1phi, sigp_Bphi);
    printf(" setting q=%g\n", q);
  }
}
/* q_of_sig_forgiven_ABP wrapper for use with zbrent_itsP */
double q_of_sigp_forgiven_AB_ZP(double sigp, void *p)
{
  double sigvec[2];
  double qvec[2];
  sigvec[1] = sigp;
  q_of_sig_forgiven_ABP(1, sigvec, qvec, p);
//  printf("q_of_sigp_forgiven_Bphi_ZP: sigvec[1]=%g qvec[1]=%g\n",
//         sigvec[1], qvec[1]);
  return qvec[1];
}

/* reset sigma such that the zeros in DNSdata_q are at A=0 */
void reset_Coordinates_CubedSphere_sigma01(tGrid *grid, tGrid *gridnew,
                                          int innerdom,  int outerdom)
{
  int iq = Ind("DNSdata_q");
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int isigma      = Ind("Coordinates_CubedSphere_dsigma01");
  int isigma_dA   = Ind("Coordinates_CubedSphere_dsigma01_dA");
  int isigma_dB   = Ind("Coordinates_CubedSphere_dsigma01_dB");
  double *q_in = grid->box[innerdom]->v[iq];
  double *q_out= grid->box[outerdom]->v[iq];  
  int n1 = grid->box[innerdom]->n1;
  int n2 = grid->box[innerdom]->n2;
  int n3 = grid->box[innerdom]->n3;
  int i,j,k, kk, JK;
  int inz_in;   /* q_in<=0  at i=inz_in (and q_in>0 i=inz_in+1) */
  int inz_out;  /* q_out<=0 at i=inz_out (and q_out>0 i=inz_out-1) */
  int i1, i2, dom; /* zero occurs between index i1 and i2 in domain dom */
  double lam1, lam2;   /* zero occurs between lam=lam1 and lam=lam2 in domain dom */
  double lam0;      /* q=0 at lam=lam0 in domain dom */
  double q1, q2;
  double A,B, x,y,z;
  int itmax = Geti("Coordinates_newtMAXITS");
  double tol = Getd("Coordinates_newtTOLF");
  double vec[2];
  int check, stat;
  int use_last_result_as_new_guess;
  int use_newton_tofind_sigp_Bphi = Getv("DNSdata_sigp_Bphi_FINDER_reset_Coordinates_AnsorgNS_sigma_pm", "newton_linesrch_itsP");
  int b;
  int star;
  if(innerdom==0)      star=STAR1;
  else if(innerdom==3) star=STAR2;

  
  {
    tGrid *grid_p = grid;
    /* how we pick initial guess in the loop below */
    use_last_result_as_new_guess=1;

    /* loop over the remaining j,k i.e. B,phi. 
       NOTE: we assume that n1,n2,n3 are the same in both domains */
    /* (j,k)-loop should be:
    for(j=n2-2; j>0; j--) // we could include j=0 (B=0) here again, so that most sigp_Bphi are found with the same method
    for(k=0; k<n3; k++)  */
    for(JK=n3*(n2-2)-1; JK>=0; JK--) /* JK = n3*(J-1) + K , J=1,...,n2-2, K=n3-1-k, K=0,...,n3-1 */
    {
      /* compute j,k from JK */
      int j = JK/n3 + 1;
      int k = (n3-1) - (JK%n3);
      double vec[2];
      double B,phi, sigp_Bphi, sigp_old;
      double w0, w1, wold;
      int check, stat;
      t_grid_box_XRphi_sigp_1phi_B_icoeffs_innerdom_outerdom_struct pars[1];

      /* set sigp_Bphi = sigp_1phi when we enter loop at j=n2-2, k=0 */
      if(JK==n3*(n2-2)-1) sigp_Bphi = sigp_1phi;

      /* find sigp_Bphi at B,phi such that q(sigp_Bphi; lam=0, B, phi)=0 */
      B   = grid_p->box[dom]->v[iY][Index(0,j,k)];
      phi = grid_p->box[dom]->v[iZ][Index(0,j,k)];
      /* use newton_linesrch_its to find sigp_Bphi */
      pars->sigp_1phi = sigp_1phi;
      pars->B = B;
      pars->phi = phi;
      pars->grid = grid_p;
      pars->innerdom = innerdom;
      pars->outerdom = outerdom;
      if(use_last_result_as_new_guess)
      {
        /* guess for vec[1] is set to sigp_Bphi, the result of previous iteration */
        vec[1] = sigp_Bphi;
      }
      else
      {
        /* old sigp at B,phi and weights */
        sigp_old = grid_p->box[innerdom]->v[isigma][Index(0,j,k)];
        w0 = 1.0-Attenuation01(B*2.0, 2.5, 0.5);
        w1 = Attenuation01(B*2.0-1.0, 2.5, 0.5);
        wold = 1.0 - w0 - w1;
        /* guess for vec[1] is weighted average with weights w0,w1,wold */
        vec[1] = w0*sigp_0phi + w1*sigp_1phi + wold*sigp_old;
      }
//printf("itmax=%d tol=%g vec[1]=%g B=%g phi=%g\n",itmax,tol,vec[1], B,phi);
      if(use_newton_tofind_sigp_Bphi)
      {
        stat=newton_linesrch_itsP(vec, 1, &check, q_of_sig_forgiven_ABP,
                                  (void *) pars, itmax, tol);
      }
      else
      {
        double sigp=vec[1], sigp1, sigp2;

        if( (innerdom==0 && sigp<=0) || (innerdom==3 && sigp>=0) )
          sigp = sigp_1phi;
        if(sigp>0) { sigp1 = sigp*0.99;  sigp2 = sigp*1.01; }
        else       { sigp2 = sigp*0.99;  sigp1 = sigp*1.01; }
        if(zbrac_P(q_of_sigp_forgiven_Bphi_ZP, &sigp1,&sigp2, (void *) pars)<0)
          errorexit("cannot find bracket for q_of_sigp_forgiven_Bphi_ZP");
        if(sigp1*sigp2<=0.0)
        {
          printf("bad bracket: [sigp1,sigp,sigp2]=[%g,%g,%g] -->\n",
                 sigp1,sigp,sigp2);
          if(innerdom==0) { sigp1 = sigp*0.1;  sigp2=sigp*2; }
          else            { sigp2 = sigp*0.1;  sigp1=sigp*2; }
          printf("new bracket: [sigp1,sigp,sigp2]=[%g,%g,%g]\n",
                 sigp1,sigp,sigp2);
        }
        check=0;
        stat=zbrent_itsP(&sigp, q_of_sigp_forgiven_Bphi_ZP,  sigp1,sigp2,
                         (void *) pars, itmax, tol);
        vec[1] = sigp;
      }
      /* If q is nowhere negative newton_linesrch_its may not work. In this
         case we should probably search for the zero in (q - 1e-8). */
//printf("stat=%d\n",stat);
      if(check || stat<0)
        printf("reset_Coordinates_AnsorgNS_sigma_pm: check=%d stat=%d\n",
               check, stat);
      sigp_Bphi = vec[1];

      /* set Coordinates_AnsorgNS_sigma_pm = sigp_Bphi in both domains */
      for(i=0; i<n1; i++)
      {
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] = sigp_Bphi;
      }
//printf("B=%g phi=%g  ", B, phi);
//printf("sigp_Bphi=%g sigp_0phi=%g sigp_1phi=%g\n", sigp_Bphi, sigp_0phi, sigp_1phi);
    } /* end for JK */
  }
  /* make sure that sigma has only one value at B=0 and also at B=1 */
  for(j=0; j<n2; j+=n2-1)
    for(k=1; k<n3; k++)
      for(i=0; i<n1; i++)
        gridnew->box[innerdom]->v[isigma][Index(i,j,k)] =
        gridnew->box[outerdom]->v[isigma][Index(i,j,k)] =
                        gridnew->box[innerdom]->v[isigma][Index(0,j,0)];

  /* compute derivs of sigma in both domains */
  spec_Deriv1(gridnew->box[innerdom], 2, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[innerdom], 3, gridnew->box[innerdom]->v[isigma],
              gridnew->box[innerdom]->v[isigma_dphi]);
  spec_Deriv1(gridnew->box[outerdom], 2, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dB]);
  spec_Deriv1(gridnew->box[outerdom], 3, gridnew->box[outerdom]->v[isigma],
              gridnew->box[outerdom]->v[isigma_dphi]);
}





/******************************************************************/
/* useful functions                                               */
/******************************************************************/

/* compute volume integral of var with index vind in a star (STAR1 or STAR2).
   Here any Psi^6 needs to be already included in the var we integrate. */
double InnerVolumeIntegral(tGrid *grid, int star, int vind)
{
  ...
  return VolInt;
}

/* compute volume integral of var with index vind in domain1 (if b=1)
   or domain2 (if b=2), but it should also work in the other domains.
   Here any Psi^6 needs to be already included in the var we integrate. */
double VolumeIntegral_inDNSgridBox(tGrid *grid, int b, int vind)
{
  double *var;
  double *Integ;
  double *Temp3;
  double *dXdx;
  double *dXdy;
  double *dXdz;
  double *dYdx;
  double *dYdy;
  double *dYdz;
  double *dZdx;
  double *dZdy;
  double *dZdz;
  double VolInt;
  int i;

  var   = grid->box[b]->v[vind];
  dXdx  = grid->box[b]->v[Ind("dXdx")];
  dXdy  = grid->box[b]->v[Ind("dXdx")+1];
  dXdz  = grid->box[b]->v[Ind("dXdx")+2];
  dYdx  = grid->box[b]->v[Ind("dYdx")];
  dYdy  = grid->box[b]->v[Ind("dYdx")+1];
  dYdz  = grid->box[b]->v[Ind("dYdx")+2];
  dZdx  = grid->box[b]->v[Ind("dZdx")];
  dZdy  = grid->box[b]->v[Ind("dZdx")+1];
  dZdz  = grid->box[b]->v[Ind("dZdx")+2];
  Integ = grid->box[b]->v[Ind("DNSdata_temp2")];
  Temp3 = grid->box[b]->v[Ind("DNSdata_temp3")];

  /* set integrand in Integ */
  forallpoints(grid->box[b], i)
  {
    double det = dXdx[i]*dYdy[i]*dZdz[i] + dXdy[i]*dYdz[i]*dZdx[i] +
                 dXdz[i]*dYdx[i]*dZdy[i] - dXdz[i]*dYdy[i]*dZdx[i] -
                 dXdy[i]*dYdx[i]*dZdz[i] - dXdx[i]*dYdz[i]*dZdy[i];
    double jac;

    if(det!=0.0) jac = 1.0/fabs(det);
    /* if det=0 jac should really be infinite, but we hope that the
       integrand goes to zero quickly enough that jac=0 makes difference! */
    else jac = 0.0;

    Integ[i] = var[i] * jac;
  }

  /* integrate */
  VolInt = spec_3dIntegral(grid->box[b], Integ, Temp3);

  return VolInt;
}




/* figure out max A inside stars and adjust boxes4/5 accordingly */
void adjust_box4_5_pars(tGrid *grid)
{
  adjust_box4_5_pars_pm(grid, -1); /* adjust both box4 and 5 */
}




/* set bfaces, set pointlist only if set_fpts=1.  */
int DNSgrid_set_bfaces(tGrid *grid, int set_fpts, int pr)
{
  int b, ob;
  int i,j,k, f, fi, p, dir;
  int Xi = Ind("X");
  int Yi = Xi+1;
  int Zi = Xi+2;
  int xi = Ind("x");
  int yi = xi+1;
  int zi = xi+2;
  int Ai = Ind("DNSdata_A");
  int Bi = Ai+1;
  int phii = Ai+2;

  if(Getv("DNSdata_grid", "4ABphi_2xyz"))
  {
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;

      free_all_bfaces(box);

      if(b==0)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 1; /* face touches box1 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 5; /* face overlaps box5 */
        box->bface[fi]->oXi = xi;
        box->bface[fi]->oYi = yi;
        box->bface[fi]->oZi = zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==3)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 2; /* face touches box2 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 4; /* face overlaps box4 */
        box->bface[fi]->oXi = xi;
        box->bface[fi]->oYi = yi;
        box->bface[fi]->oZi = zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==1)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 0; /* face touches box0 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 2; /* face touches box2 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* infinity edge */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1; /* mark as outer boundary */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) for(k=0;k<n3;k++)
        {
          int ijk = Index(p,0,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==2)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 3; /* face touches box3 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 1; /* face touches box1 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* infinity edge */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1; /* mark as outer boundary */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) for(k=0;k<n3;k++)
        {
          int ijk = Index(p,0,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==4)
      {
        for(f=0; f<6; f++)
        {
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 3; /* face overlaps box3 */
          box->bface[fi]->oXi = Ai;
          box->bface[fi]->oYi = Bi;
          box->bface[fi]->oZi = phii;
          dir = 1+f/2;
          p = ( (n1-1) )*(f%2);
          if(set_fpts) forplaneN(dir, i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else if(b==5)
      {
        for(f=0; f<6; f++)
        {
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 0; /* face overlaps box0 */
          box->bface[fi]->oXi = Ai;
          box->bface[fi]->oYi = Bi;
          box->bface[fi]->oZi = phii;
          dir = 1+f/2;
          p = ( (n1-1) )*(f%2);
          if(set_fpts) forplaneN(dir, i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else errorexit("we must have 0 <= b <= 5");
    }
  }
  if(Getv("DNSdata_grid", "6ABphi_2xyz"))
  {
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      int n1 = box->n1;
      int n2 = box->n2;
      int n3 = box->n3;

      free_all_bfaces(box);

      if(b==0)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 1; /* face touches box1 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 5; /* face overlaps box5 */
        box->bface[fi]->oXi = xi;
        box->bface[fi]->oYi = yi;
        box->bface[fi]->oZi = zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==3)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 2; /* face touches box2 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 4; /* face overlaps box4 */
        box->bface[fi]->oXi = xi;
        box->bface[fi]->oYi = yi;
        box->bface[fi]->oZi = zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
      }
      else if(b==1)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 0; /* face touches box0 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        p = ( (n1-1) )*(f%2);
        /* part of face1 touches box6 */
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 6;
        box->bface[fi]->oXi = Xi;
        box->bface[fi]->oYi = Yi;
        box->bface[fi]->oZi = Zi;
        if(set_fpts)
        {
          double *B = box->v[Ind("Y")];
          double Bmin = grid->box[6]->bbox[2];
          forplane1(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(B[ijk]>Bmin)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
        /* part of face1 is outerbound */
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1;
        if(set_fpts)
        {
          double *B = box->v[Ind("Y")];
          double Bmin = grid->box[6]->bbox[2];
          forplane1(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(B[ijk]<=Bmin)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else if(b==2)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 3; /* face touches box3 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        p = ( (n1-1) )*(f%2);
        /* part of face1 touches box7 */
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 7;
        box->bface[fi]->oXi = Xi;
        box->bface[fi]->oYi = Yi;
        box->bface[fi]->oZi = Zi;
        if(set_fpts)
        {
          double *B = box->v[Ind("Y")];
          double Bmin = grid->box[7]->bbox[2];
          forplane1(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(B[ijk]>Bmin)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
        /* part of face1 is outerbound */
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1;
        if(set_fpts)
        {
          double *B = box->v[Ind("Y")];
          double Bmin = grid->box[7]->bbox[2];
          forplane1(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(B[ijk]<=Bmin)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else if(b==6)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 1; /* face touches box1 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        box->bface[fi]->oXi = Xi;
        box->bface[fi]->oYi = Yi;
        box->bface[fi]->oZi = Zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 7; /* face touches box7 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* outer boundary at face 2 */
        f = 2;
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1; /* mark as outer boundary */
        p = ( (n1-1) )*(f%2);
        if(set_fpts)
        {
          double *A = box->v[Ind("X")];
          double Amax = grid->box[1]->bbox[1];
          forplane2(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(A[ijk]>=Amax)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
        /* interpolation at face 2, if box6 extends into box1 */
        if(box->bbox[0] < grid->box[1]->bbox[1])
        {
          f = 2;
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 1;  /* is inside box1 */
          box->bface[fi]->oXi = Xi;
          box->bface[fi]->oYi = Yi;
          box->bface[fi]->oZi = Zi;
          p = ( (n1-1) )*(f%2);
          if(set_fpts)
          {
            double *A = box->v[Ind("X")];
            double Amax = grid->box[1]->bbox[1];
            forplane2(i,j,k, n1,n2,n3, p)
            {
              int ijk = Index(i,j,k);
              if(A[ijk]<Amax)
                add_point_to_bface_inbox(box, fi, ijk, f);
            }
          }
        }
      }
      else if(b==7)
      {
        /* face0 */
        f = 0;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 2; /* face touches box2 */
        box->bface[fi]->setnormalderiv = 1; /* set normal derivs equal */
        box->bface[fi]->oXi = Xi;
        box->bface[fi]->oYi = Yi;
        box->bface[fi]->oZi = Zi;
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* face1 */
        f = 1;
        fi = add_empty_bface(box, f);
        box->bface[fi]->ob = 6; /* face touches box6 */
        p = ( (n1-1) )*(f%2);
        if(set_fpts) forplane1(i,j,k, n1,n2,n3, p)
        {
          int ijk = Index(i,j,k);
          add_point_to_bface_inbox(box, fi, ijk, f);
        }
        /* outer boundary at face 2 */
        f = 2;
        fi = add_empty_bface(box, f);
        box->bface[fi]->outerbound = 1; /* mark as outer boundary */
        p = ( (n1-1) )*(f%2);
        if(set_fpts)
        {
          double *A = box->v[Ind("X")];
          double Amax = grid->box[2]->bbox[1];
          forplane2(i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            if(A[ijk]>=Amax)
              add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
        /* interpolation at face 2, if box7 extends into box2 */
        if(box->bbox[0] < grid->box[2]->bbox[1])
        {
          f = 2;
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 2;  /* is inside box2 */
          box->bface[fi]->oXi = Xi;
          box->bface[fi]->oYi = Yi;
          box->bface[fi]->oZi = Zi;
          p = ( (n1-1) )*(f%2);
          if(set_fpts)
          {
            double *A = box->v[Ind("X")];
            double Amax = grid->box[2]->bbox[1];
            forplane2(i,j,k, n1,n2,n3, p)
            {
              int ijk = Index(i,j,k);
              if(A[ijk]<Amax)
                add_point_to_bface_inbox(box, fi, ijk, f);
            }
          }
        }
      }
      else if(b==4)
      {
        for(f=0; f<6; f++)
        {
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 3; /* face overlaps box3 */
          box->bface[fi]->oXi = Ai;
          box->bface[fi]->oYi = Bi;
          box->bface[fi]->oZi = phii;
          dir = 1+f/2;
          p = ( (n1-1) )*(f%2);
          if(set_fpts) forplaneN(dir, i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else if(b==5)
      {
        for(f=0; f<6; f++)
        {
          fi = add_empty_bface(box, f);
          box->bface[fi]->ob = 0; /* face overlaps box0 */
          box->bface[fi]->oXi = Ai;
          box->bface[fi]->oYi = Bi;
          box->bface[fi]->oZi = phii;
          dir = 1+f/2;
          p = ( (n1-1) )*(f%2);
          if(set_fpts) forplaneN(dir, i,j,k, n1,n2,n3, p)
          {
            int ijk = Index(i,j,k);
            add_point_to_bface_inbox(box, fi, ijk, f);
          }
        }
      }
      else errorexit("we must have 0 <= b <= 7");
    }
  }

  /* set ofi and bit fields */
  set_ofi_in_all_bfaces(grid);
  set_bits_in_all_bfaces(grid);
  if(pr) forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    printbfaces(box);
    //for(fi=0; fi<box->nbfaces; fi++)
    //  prPointList(box->bface[fi]->fpts);
  }
  return 0;
}



/* return the box index and the coords of the point at the Cartesian x,y,z */
/* initially box, *X,*Y,*Z contain the box and the coords of the
   point on the other grid */
/* NOTE: Interp_Var_From_Grid1_To_Grid2_pm has certain smoothing
         because of DNSgrid_Get_BoxAndCoords_of_xyz (see below). */
int DNSgrid_Get_BoxAndCoords_of_xyz(tGrid *grid1,
                                    double *X1, double *Y1, double *Z1,
                                    tBox *box, int ind, double x, double y, double z)
{
  int b1;

  b1 = b_XYZ_of_xyz(grid1, &X,&Y,&Z, x,y,z);
  return b1;
}


/* Interpolate Var with index vind from grid1 to grid2 
   for domains centered on innerdom */
void Interp_Var_From_Grid1_To_Grid2_pm(tGrid *grid1, tGrid *grid2, int vind,
                                       int star)
{
  int cind = Ind("temp1");
  int Xind = Ind("X");
  int Yind = Ind("Y");
  int Zind = Ind("Z");
  int xind = Ind("x");
  int yind = Ind("y");
  int zind = Ind("z");
  int b,i;

  /* save coeffs of vind on grid1 in cind = Ind("Temp1") */
  forallboxes(grid1, b)
  {
    tBox *box = grid1->box[b];
    /* do nothing if we are on wrong side of grid */
    if(box->SIDE!=star) continue;

    spec_Coeffs(box, box->v[vind], box->v[cind]);
  }

  /* loop over grid2 */
  forallboxes(grid2,b)
  {
    tBox *box = grid2->box[b];
    double *pX = box->v[Xind];
    double *pY = box->v[Yind];
    double *pZ = box->v[Zind];
    double *px = box->v[xind];
    double *py = box->v[yind];
    double *pz = box->v[zind];
    double *pv = box->v[vind];

    /* do nothing if we are on wrong side of grid */
    if(box->SIDE!=star) continue;

    /* here we can use SGRID_LEVEL6_Pragma(omp parallel) */
#undef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
#ifndef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
    SGRID_LEVEL6_Pragma(omp parallel)
    {
      tGrid *grid1_p = make_empty_grid(grid1->nvariables, 0);
      copy_grid(grid1, grid1_p, 0);
#else
    {
      tGrid *grid1_p = grid1;
#endif
      /* start loop */
#ifndef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
      SGRID_LEVEL6_Pragma(omp for) 
#endif
      forallpoints(box,i)
      {
        double X = pX[i];
        double Y = pY[i];
        double Z = pZ[i];
        double x = px[i];
        double y = py[i];
        double z = pz[i];
        double r = sqrt(x*x + y*y + z*z);
        int b1;

        /* get b1, X,Y,Z on grid1_p */
        b1 = DNSgrid_Get_BoxAndCoords_of_xyz(grid1_p, &X,&Y,&Z, box,i, x,y,z);
        if(b1>=0)
        {
          /* get var at point X,Y,Z by interpolation */
          pv[i] = spec_interpolate(grid1_p->box[b1], grid1_p->box[b1]->v[cind], X,Y,Z);
        }
        else /* point px[i],py[i],pz[i] may be beyond outer boundary */
        {
          errorexit("could not find point")
        }

      } /* end forallpoints loop */
#ifdef SERIAL_Interp_Var_From_Grid1_To_Grid2_pm
    }
#else
      /* free local grid copy */
      free_grid(grid1_p);
    }
#endif
  } /* end forallboxes(grid2,b) */
}
/* wrapper for Interpolate_Var_From_Grid1_To_Grid2 with extra dummy argument */
void Interpolate_Var_From_Grid1_To_Grid2_wrapper(tGrid *grid1, tGrid *grid2,
                                                 int vind, int dummy)
{
  Interpolate_Var_From_Grid1_To_Grid2(grid1, grid2, vind);
}






/* compute weighted average of the new q2 on grid2 and the old q1 on grid1 
   at a point X2,Y2,Z2 in grid2 coords */
double DNS_update_q_atXYZ(tGrid *grid2, 
                          int b2, double X2, double Y2, double Z2,
                          double w, tGrid *grid1)
{
  double x,y,z, Xp,Rp;
  int b1;
  double X1,Y1,Z1;
  double q2, q1;
  tBox *box2 = grid2->box[b2];

  /* get q on grid 2 at X2,Y2,Z2 */
  q2 = DNS_compute_new_centered_q_atXYZ(grid2,b2, X2,Y2,Z2);

  /* get q on grid 1 at the same point */
  if(box2->COORD!=CART)
  {
    int domain = get_AnsorgNS_domain(box2);
    xyz_of_AnsorgNS(box2, -1, domain, X2,Y2,Z2, &x,&y,&z, &Xp,&Rp);
  }
  else
    { x=X2;  y=Y2;  z=Z2; }
  X1 = X2;  Y1 = Y2;  Z1 = Z2;
  b1 = DNSgrid_Get_BoxAndCoords_of_xyz(grid1, &X1,&Y1,&Z1, box2,-1,x,y,z);
  q1 = DNS_compute_new_centered_q_atXYZ(grid1,b1, X1,Y1,Z1);

  /* return weighted average */
  return w*q2 + (1.0-w)*q1;
}

/* compute weighted average of the new q2 on grid2 and the old q1 on grid1 
   at a point X2,Y2,Z2 in grid2 coords */
void DNS_update_q(tGrid *grid2, double w, tGrid *grid1)
{
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int iq = Ind("DNSdata_q");
  int iqg= Ind("DNSdata_qg");
  int b2, i;

  forallboxes(grid2,b2)
  {
    tBox *box = grid2->box[b2];
    double *X2 = box->v[iX];
    double *Y2 = box->v[iY];
    double *Z2 = box->v[iZ];
    double *q = box->v[iq];
    double *qg= box->v[iqg];

    forallpoints(box, i)
    {
      q[i] = DNS_update_q_atXYZ(grid2,b2, X2[i],Y2[i],Z2[i], w, grid1);
      qg[i]= q[i];
    }
  }
}


/* set a coord. dependent factor that we can use to multiply eqns that go to
   zero because dX^i/dx^j goes to zero */
int set_DNSdata_CoordFac(tGrid *grid)
{
  int index_dXdx = Ind("dXdx");
  int index_dYdx = Ind("dYdx");
  int index_CoordFac = Ind("DNSdata_CoordFac");
  double CoordFacPower = Getd("DNSdata_CoordFacPower");
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    double *dXdx1 = box->v[index_dXdx + 0];
    double *dXdx2 = box->v[index_dXdx + 1];
    double *dXdx3 = box->v[index_dXdx + 2];
    double *dYdx1 = box->v[index_dYdx + 0];
    double *dYdx2 = box->v[index_dYdx + 1];
    double *dYdx3 = box->v[index_dYdx + 2];
    double *CoordFac = box->v[index_CoordFac];
    int ijk;

    forallpoints(box, ijk)
    {
      if(CoordFacPower != 0. && dXdx1 != NULL && (b==1 || b==2))
      {
        double s11=0.0;
        double s12=0.0;
        double s13=0.0;
        double s21=1.0;
        double s22=1.0;
        double s23=1.0;
        double fac1 = s11*pow2(dXdx1[ijk])+s12*pow2(dXdx2[ijk])+s13*pow2(dXdx3[ijk]);
        double fac2 = s21*pow2(dYdx1[ijk])+s22*pow2(dYdx2[ijk])+s23*pow2(dYdx3[ijk]);

        CoordFac[ijk] = fac1 + fac2;
        CoordFac[ijk] = sqrt(CoordFac[ijk]);
        if(CoordFac[ijk] == 0.) CoordFac[ijk] = 1.0;
        else                    CoordFac[ijk] = pow(CoordFac[ijk], CoordFacPower);
      }
      else
        CoordFac[ijk] = 1.0;
    }
  }
  return 0;
}

/************************************************************************/
/* utilities to load checkpoints with different resolution and pars */
/************************************************************************/
/* for newton_linesrch_itsP: compute error in rest mass */
void m0errOFsigmafac_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int dom;
  double fac = vec[1];
  double m0;
  t_grid_b_struct *pars;

  /* get pars */
  pars = (t_grid_b_struct *) p;
  grid = pars->grid;
  dom  = pars->b;
  /* In here it seems pars->b is the inner box i.e. 0 or 3, it should never
     contain something else. */

  /* scale sigma */
  DNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, fac, dom);
  DNSgrid_init_Coords(grid);
  /* set wB */
  if(dom==0) DNS_set_wB(grid, 1, Getd("DNSdata_xmax1"),0.0,0.0);
  else       DNS_set_wB(grid, 2, Getd("DNSdata_xmax2"),0.0,0.0); 
  /* get new mass */
  m0 =  GetInnerRestMass(grid, dom);
  printf("m0errOFsigmafac_VectorFuncP: fac=%g => m0%d=%g\n",
         fac, (dom>=2)+1, m0);
  fflush(stdout);
  /* set sigma back to what it was */
  DNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, 1.0/fac, dom);
  DNSgrid_init_Coords(grid);
  /* set wB */
  if(dom==0) DNS_set_wB(grid, 1, Getd("DNSdata_xmax1"),0.0,0.0);
  else       DNS_set_wB(grid, 2, Getd("DNSdata_xmax2"),0.0,0.0); 
  /* set fvec */
  if(dom==0)
    fvec[1] = Getd("DNSdata_m01") - m0;
  else
    fvec[1] = Getd("DNSdata_m02") - m0;
  //printf("  => m0_err=%g\n", fvec[1]);
  //fflush(stdout);
}

/* load data from a previous result in a checkpoint file 
   (with e.g. lower res) */
void DNSgrid_load_initial_guess_from_checkpoint(tGrid *grid, char *filename)
{
  tVarList *varlist;
  char *parlist;
  int j, bi;
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");
  char *DNSdata_b_sav = strdup(Gets("DNSdata_b"));
  char *DNSdata_x_CM_sav = strdup(Gets("DNSdata_x_CM"));
  char *DNSdata_xmax1_sav = strdup(Gets("DNSdata_xmax1"));
  char *DNSdata_xmax2_sav = strdup(Gets("DNSdata_xmax2"));

  /* make list of pars we want to read */  
  parlist = "DNSdata_b DNSdata_m01 DNSdata_m02 DNSdata_kappa "
            "DNSdata_Omega DNSdata_x_CM DNSdata_C1 DNSdata_C2 "
            "DNSdata_qmax1 DNSdata_qmax2 DNSdata_xmax1 DNSdata_xmax2 "
            "DNSdata_actual_xmax1 DNSdata_actual_xmax2";

  /* make list of vars we want to read */
  varlist = vlalloc(grid);
  vlpush(varlist, Ind("DNSdata_Psi"));
  vlpush(varlist, Ind("DNSdata_Bx"));
  vlpush(varlist, Ind("DNSdata_alphaP"));
  vlpush(varlist, Ind("DNSdata_Sigma"));
  vlpush(varlist, Ind("DNSdata_q"));
  vlpush(varlist, Ind("Coordinates_AnsorgNS_sigma_pm"));

//printVarList(varlist);

  /* read vars in varlist and pars in parlist*/
  checkpoint_interpolate_Vars_get_Pars(filename, varlist, parlist);

  /* once we have a new Coordinates_AnsorgNS_sigma_pm, we need re-init
     the Coords again. */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");
  init_CoordTransform_And_Derivs(grid);
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  if(strcmp(DNSdata_b_sav, Gets("DNSdata_b"))==0)
  {
    /* set values of A,B,phi in box4/5 */
    set_DNSdata_ABphi(grid);

    /* set wB */
    DNS_set_wB(grid, 1, Getd("DNSdata_actual_xmax1"),0.0,0.0);
    DNS_set_wB(grid, 2, Getd("DNSdata_actual_xmax2"),0.0,0.0); 
  }
  else /* set DNSdata_b=DNSdata_b_sav and adjust other pars */
  {
    double m01 = Getd("DNSdata_m01");
    double m02 = Getd("DNSdata_m02");
    double r, rc, dr, dOm;
    double facvec[2];
    int stat, check;
    t_grid_b_struct pars[1];

    /* rc = (distance from checkpoint) */
    rc = Getd("DNSdata_xmax1")-Getd("DNSdata_xmax2");
    Sets("DNSdata_b", DNSdata_b_sav); /* set b back to saved value */
    /* set some things back to saved values */
    Sets("DNSdata_x_CM", DNSdata_x_CM_sav);
    Sets("DNSdata_xmax1", DNSdata_xmax1_sav);
    Sets("DNSdata_xmax2", DNSdata_xmax2_sav);
    r = Getd("DNSdata_xmax1")-Getd("DNSdata_xmax2");
    dr = r-rc;
    /* adjust Omega according to Kepler's law 
       Om = sqrt(M/r^3)  ==> dOm = -(3/2) Om/r * dr. Here r=2*bc dr=2dx */
    dOm = -1.5*(Getd("DNSdata_Omega")/rc)*dr;
    Setd("DNSdata_Omega", Getd("DNSdata_Omega")+dOm);

    printf("Setting:\n");
    printf(" DNSdata_b = %s\n", Gets("DNSdata_b"));
    printf(" DNSdata_xmax1 = %s\n", Gets("DNSdata_xmax1"));
    printf(" DNSdata_xmax2 = %s\n", Gets("DNSdata_xmax2"));
    printf(" DNSdata_Omega = %s\n", Gets("DNSdata_Omega"));
    printf(" DNSdata_x_CM = %s\n", Gets("DNSdata_x_CM"));
    /* adjust coords such that all is ok, e.g. box5/4 need to cover
       holes in box0/3 */
    DNSgrid_init_Coords(grid);
    /* set wB */
    DNS_set_wB(grid, 1, Getd("DNSdata_xmax1"),0.0,0.0);
    DNS_set_wB(grid, 2, Getd("DNSdata_xmax2"),0.0,0.0); 
    /* Ok now we have complete data, but the var 
       Coordinates_AnsorgNS_sigma_pm should probabaly be changed as well! */

    /* pick Coordinates_AnsorgNS_sigma_pm such that the masses are what we
       want */
    /***********************************************************************************************/
    /* do newton_linesrch_itsP iterations until sigma+ and sigma- are such that masses are correct */
    /***********************************************************************************************/
    printf("Rescaling Coordinates_AnsorgNS_sigma_pm so that masses are correct:\n");
    /* star1 */
    pars->grid = grid;
    pars->b = 0;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    DNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, facvec[1], 0);
    DNSgrid_init_Coords(grid);
    /* set wB */
    DNS_set_wB(grid, 1, Getd("DNSdata_xmax1"),0.0,0.0);

    /* star2 */
    pars->grid = grid;
    pars->b = 3;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m02);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    DNSgrid_scale_Coordinates_AnsorgNS_sigma(grid, facvec[1], 3);
    DNSgrid_init_Coords(grid);
    /* set wB */
    DNS_set_wB(grid, 2, Getd("DNSdata_xmax2"),0.0,0.0); 
  }

  /* check the vars we set for NANs */
  vlpush(varlist, Ind("DNSdata_wBx"));
  printf("checking for NANs in:\n");
  for(j = 0; j < varlist->n; j++) printf("%s ", VarName(varlist->index[j]));
  printf("\n");  
  NumberChecker_CheckIfFinite_VarList(varlist);
  fflush(stdout);

  /* free varlist */
  vlfree(varlist);
  
  /* free saved strings */
  free(DNSdata_b_sav);
  free(DNSdata_x_CM_sav);
  free(DNSdata_xmax1_sav);
  free(DNSdata_xmax2_sav);
}
