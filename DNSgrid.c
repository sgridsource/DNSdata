/* grid_setup.c */
/* Wolfgang Tichy 2008 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow
#define pow2(x)    ((x)*(x))

/* struct types used in root finder newton_linesrch_itsP */
typedef struct T_grid_b_star_A_B_struct {
  tGrid *grid; /* grid */
  int b;       /* box number */
  int b_in;    /* box number of inner box */
  int star;    /* STAR1/2 */
  double A;    /* A coord */
  double B;    /* B coord */
} t_grid_b_star_A_B_struct;

/* global vars in this file */
double rf_surf1; /* radius of star1 */
double rf_surf2; /* radius of star2 */
double P_core1;  /* core pressure of star1 */
double P_core2;  /* core pressure of star2 */

/* funs in this file */
void m0_VectorFuncP(int n, double *vec, double *fvec, void *p);


/* find core pressure in TOV star */
double DNS_find_P_core(double m0, int pr)
{
  double vec[2];
  double fvec[2];
  int check, stat;
  double par[2];
  par[0] = m0;
  par[1] = pr;

  /* guard against unreasonable m0 */
  if(m0<=0.) return 0.;

  /* find P_core, s.t. rest mass is m0 */
  if(pr) printf("find P_core of a TOV star, s.t. rest mass is m0=%g\n", m0);
  vec[1] = 1e-7;   /* initial guess */
  /* do newton_linesrch_its iterations: */
  stat=newton_linesrch_itsP(vec, 1, &check, m0_VectorFuncP,
                            (void *) par,
                            Geti("Coordinates_newtMAXITS"),
                            Getd("Coordinates_newtTOLF") );
  if(check || stat<0)
  {
    printf("WARNING from newton_linesrch_its with MAXITS=%d TOLF=%g:\n",
           Geti("Coordinates_newtMAXITS"), Getd("Coordinates_newtTOLF"));
    printf("  check=%d stat=%d\n", check, stat);
  }

  /* check if we found correct vec or just a local max in m0 */
  m0_VectorFuncP(1, vec, fvec, (void *) par);
  if( (fabs(fvec[1])>Getd("Coordinates_newtTOLF")*1000) && pr )
  {
    int i;
    double v[2], f[2];
    printf("There may be a max in m0 near Pc=%g\n", vec[1]);
    for(i=-15; i<=15; i++)
    { 
      v[1] = vec[1] + i*vec[1]/1000;
      m0_VectorFuncP(1, v, f, (void *) par);
    }
  }
  return vec[1];
}

/* setup initial boxsizes */
int set_DNS_boxsizes(tGrid *grid)
{
  double m1, Phic1, Psic1;
  double m2, Phic2, Psic2;
  //double kappa     = Getd("DNSdata_kappa");
  //double DNSdata_n = Getd("DNSdata_n");
  double DNSdata_b = Getd("DNSdata_b");
  double m01 = Getd("DNSdata_m01");
  double m02 = Getd("DNSdata_m02");
  double xmin1,xmax1, xmin2,xmax2, xc1, xc2; /* x-positions of stars */
  double DoM, nu; /* distance over total rest mass, and rest mass ratio */
  double DoM3, DoM4, DoM5; /* powers of DoM */
  double xCM, Omega;
  double Fc, qc, Cc, oouzerosqr;

  printf("set_DNS_boxsizes: setting box sizes and coordinates used ...\n");
  prTimeIn_s("WallTime: ");

  /* reset initial DNSdata_m01/2 if needed */
  if(Getv("DNSdata_iterate_m0", "yes"))
  {
    printf(" DNSdata_iterate_m0 = yes : setting:\n");
    /* set DNSdata_desired_m01/2 if needed */
    if(Getd("DNSdata_desired_m01")<0.0) Setd("DNSdata_desired_m01", m01);
    if(Getd("DNSdata_desired_m02")<0.0) Setd("DNSdata_desired_m02", m02);
    printf("   DNSdata_desired_m01 = %g\n", Getd("DNSdata_desired_m01"));
    printf("   DNSdata_desired_m02 = %g\n", Getd("DNSdata_desired_m02"));
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
    DNS_pwp_init_parameter();
  }
  else if(Getv("DNSdata_EoS_type", "poly"))  DNS_poly_init();
  else errorexit("unkown DNSdata_EoS_type");

  /* do we set star from mass m01 or from max q called qm1 */
  if(Getd("DNSdata_qm1")<0.0)
  {
    /* many parts of the code only work with m01>0 */
    if(Getd("DNSdata_m01")<=0.0) errorexit("make sure DNSdata_m01 > 0");
    /* find P_core1, s.t. rest mass is m01 */
    P_core1 = DNS_find_P_core(Getd("DNSdata_m01"), 1);
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
    P_core2 = DNS_find_P_core(Getd("DNSdata_m02"), 1);
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

  /* look at DNSdata_grid and maybe setup box resolutions */
  // ...

  prTimeIn_s("WallTime: ");
  return 0;
}


/* setup the cubed sphere boxes */
int DNSdata_setup_boxes(tGrid *grid)
{
  double dc = Getd("DNSdata_b");
  double csize = Getd("DNSdata_InnerCubesSize");
  double ssfac = Getd("DNSdata_OuterShellStart");
  double obfac = Getd("DNSdata_OuterBoundary");
  double xc[4];

  /* enable the vars for sigma01, so they get used in CI set below */
  enablevar(grid, Ind("Coordinates_CubedSphere_sigma01_def"));
  //enablevar(grid, Ind("Coordinates_CubedSphere_dsigma01_dA"));
  //enablevar(grid, Ind("Coordinates_CubedSphere_dsigma01_dB"));

  /* set cubed spheres, for now we disregard the par DNSdata_grid */
  switch(grid->nboxes)
  {
    case 13:
      xc[2] = xc[3] = 0.0;
      xc[1] = dc;
      arrange_1box12CubSph_into_full_cube(grid, 0, xc, 
                                          csize*rf_surf1, rf_surf1, dc);
      break;
    case 26:
      two_full_cubes_touching_at_x0(grid, 0, dc,
                                    csize*rf_surf1, rf_surf1,
                                    csize*rf_surf2, rf_surf2);
      break;
    case 32:
      sphere_around_two_full_cubes_touching_at_x0(grid, 0, dc,
                                                  csize*rf_surf1, rf_surf1,
                                                  csize*rf_surf2, rf_surf2,
                                                  ssfac*dc);
      break;
    case 38:
      printf("using: DNSdata_grid = 36CS_2xyz\n");
      two_spheres_around_two_full_cubes(grid, 0, dc,
                                        csize*rf_surf1, rf_surf1,
                                        csize*rf_surf2, rf_surf2,
                                        ssfac*dc, obfac*dc);
      break;
    default:
      errorexit("nboxes should be 13, 26, 32, or 38");
  }
  return 0;
}

/* set attributes in boxes */
int set_DNS_box_attribs(tGrid *grid)
{
  double DNSdata_b = Getd("DNSdata_b");
  double xmax1 = +DNSdata_b;
  double xmax2 = -DNSdata_b;
  int b;
  //printf("set_DNS_box_attribs:\n");
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    char str[1000];
    double x = 0.0;

    /* left or right? */
    if(box->x_of_X[1]!=NULL)
      x = box->x_of_X[1]((void *)box, -1, 0.5,0.0,0.0);
    if(x>0.0  )    box->SIDE = STAR1;
    else if(x<0.0) box->SIDE = STAR2;
    else           box->SIDE = ZERO;

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
      /* remove sigma01 in this box */
      disable_Coordinates_CubedSphere_sigma01(box);
    }
    if( Getv(str, "CubedSphere") && box->CI->type==outerCubedSphere )
    {
      double y, z, r, rs;
      box->COORD = CUBSPH;
      y = box->x_of_X[2]((void *)box, -1, 0.5,0.0,0.0);
      z = box->x_of_X[3]((void *)box, -1, 0.5,0.0,0.0);
      if(box->SIDE == STAR1) { r = x-xmax1;  rs = rf_surf1; }
      else                   { r = x-xmax2;  rs = rf_surf2; }
      r = sqrt(r*r + y*y + z*z);
      if(r<rs)
      {
        box->BOUND = SSURF;   /* star surface */
        box->MATTR = INSIDE;
      }
      else
      {
        box->MATTR = AWAY;
        /* remove sigma01 in this box */
        disable_Coordinates_CubedSphere_sigma01(box);
      }
    }
    if( Getv(str, "CubedSphere") && box->CI->type==innerCubedSphere )
    { 
      double y, z, r, rs;
      box->COORD = CUBSPH;
      x = box->x_of_X[1]((void *)box, -1, 0.0,0.0,0.0);
      y = box->x_of_X[2]((void *)box, -1, 0.0,0.0,0.0);
      z = box->x_of_X[3]((void *)box, -1, 0.0,0.0,0.0);
      if(box->SIDE == STAR1) { r = x-xmax1;  rs = rf_surf1; }
      else                   { r = x-xmax2;  rs = rf_surf2; }
      r = sqrt(r*r + y*y + z*z);
      if(dlesseq(r,rs))
      {
        box->BOUND = SSURF;   /* star surface */
        box->MATTR = TOUCH;
      }
      else
      {
        box->MATTR = AWAY;
        /* remove sigma01 in this box */
        disable_Coordinates_CubedSphere_sigma01(box);
      }
    }
    if( Getv(str, "CubedSphere") &&
        (box->CI->type==PyramidFrustum || box->CI->type==CubedShell) )
    { 
      box->COORD = CUBSPH;
      box->MATTR = AWAY;
      /* remove sigma01 in this box */
      disable_Coordinates_CubedSphere_sigma01(box);
    }
    if( Getv(str, "stretchedCubedSphere") )
    { 
      box->COORD = SCUBSH;
      box->MATTR = AWAY;
      /* remove sigma01 in this box */
      disable_Coordinates_CubedSphere_sigma01(box);
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


/* funtion to be passed into newton_linesrch_its to find P_core1/2 from m01/2 */
void m0_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double *par = (double *) p;
  double m0A = par[0];
  int    pr  = par[1];
  double Pc, rf_surf, m, Phic, Psic, m0;

  Pc = vec[1];
  TOV_init(Pc, 0, &rf_surf, &m, &Phic, &Psic, &m0);
  if(pr) printf("   Pc=%+.16e  m0=%.16e\n",Pc,m0);
  fvec[1] = m0-m0A;
}


/* get inner and outer edges of both stars */
void DNS_find_xin_xout(tGrid *grid, double *xin1, double *xout1,
                                    double *xin2, double *xout2)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    /* find boxes with matter and star surfaces */
    if(box->MATTR==INSIDE && box->BOUND==SSURF)
    {
      if(box->CI->dom==0)
      {
        double x1 = box->x_of_X[1]((void *) box, -1, 1.0,0.0,0.0);
        if(box->SIDE==STAR1) *xin1  = x1;
        if(box->SIDE==STAR2) *xout2 = x1;
      }
      if(box->CI->dom==1) 
      {
        double x1 = box->x_of_X[1]((void *) box, -1, 1.0,0.0,0.0);
        if(box->SIDE==STAR1) *xout1 = x1;
        if(box->SIDE==STAR2) *xin2  = x1;
      }
    }
  }
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
/* set q (or other Var) to zero at star surface */
void set_Var_to_Val_atSurface(tGrid *grid, int vi, double Val)
{
  int b;

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    double *Var = box->v[vi];
    int i,j,k;

    /* check if we are at surface, and if so set index i there */
    if(box->BOUND!=SSURF) continue;
    if(box->MATTR==INSIDE) i=n1-1;
    else                   i=0;

    for(k=0; k<n3; k++)
    for(j=0; j<n2; j++)  
      Var[Index(i,j,k)] = Val;
  }
}


/*************************************/
/* functions to adjust star surfaces */
/*************************************/

/* WE NEED to find sigma at A,B such that q(sigma; lam=0or1, A,B) = 0 */
/* q as a func of Tlam for a given  A, B
   the transformed lambda:
   Tlam := lam - 1 , if box is b_in
   Tlam := lam,      if box is b    */
double q_of_Tlam_forgiven_AB_ZP(double Tlam, void *p)
{
  t_grid_b_star_A_B_struct *pars = (t_grid_b_star_A_B_struct *) p;
  tGrid *grid = pars->grid;
  double A    = pars->A;
  double B    = pars->B;
  int b;
  double lam, q;

  if(Tlam<0.) { lam = 1. + Tlam;  b = pars->b_in; }
  else        { lam = Tlam;       b = pars->b;    }
  q = DNS_compute_new_centered_q_atXYZ(grid,b, lam, A,B);
  //printf("b=%d lam=%g A=%g B=%g q=%g\n", b, lam, A,B, q);
  return q;
}
/* reset sigma such that the zeros in DNSdata_q are at A=0 */
void reset_Coordinates_CubedSphere_sigma01(tGrid *grid, tGrid *gridnew,
                                           int star)
{
  int iq = Ind("DNSdata_q");
  int iX = Ind("X");
  int iY = Ind("Y");
  int iZ = Ind("Z");
  int itmax = Geti("Coordinates_newtMAXITS");
  double tol = Getd("Coordinates_newtTOLF");
  int outerdom;

  forallboxes(grid, outerdom)
  {
    int i, j, k;
    int n1, n2, n3, innerdom, n1in, n2in, n3in, i0;
    double *q_in, *q_out;
    tBox *box = grid->box[outerdom];
    tBox *boxin, *boxq0, *boxnew, *boxnewin;
    int isigma0, isigma1;
    //int isigma0_dA, isigma0_dB, isigma1_dA, isigma1_dB;
    
    /* do nothing for other star and all boxes that do not touch surface */
    if(box->SIDE  != star)  continue;
    if(box->MATTR != TOUCH) continue;

    n1 = box->n1;
    n2 = box->n2;
    n3 = box->n3;
    innerdom = outerdom - 6; /* works only for my CubSph setup */
    boxin = grid->box[innerdom];
    n1in = boxin->n1;
    n2in = boxin->n2;
    n3in = boxin->n3;
    if(n2!=n2in || n3!=n3in)
      errorexit("reset_Coordinates_CubedSphere_sigma01 needs n2=n2in n3=n3in");

    /* find new boxes also on new grid */
    boxnew   = gridnew->box[outerdom];
    boxnewin = gridnew->box[innerdom];
    isigma0 = boxnew->CI->iFS[0];
    isigma1 = boxnewin->CI->iFS[1];
    //isigma0_dA = boxnew->CI->idSurfdX[0][2];
    //isigma0_dB = boxnew->CI->idSurfdX[0][3];
    //isigma1_dA = boxnewin->CI->idSurfdX[1][2];
    //isigma1_dB = boxnewin->CI->idSurfdX[1][3];

    /* set pointer to q outside and inside star */
    q_out= box->v[iq];  
    q_in = boxin->v[iq];

    /* loop over surface of star touching box */
    forplane1(i0,j,k, n1,n2,n3, 0)
    {
      t_grid_b_star_A_B_struct pars[1];
      int ind = Index(0,j,k);
      int indin = Ind_n1n2(n1in-1,j,k, n1in,n2in);
      double A = box->v[iY][ind];
      double B = box->v[iZ][ind];
      double x, y, z, xc, sig01_AB;
      int dom, stat;
      int inz_in;   /* q_in<=0  at i=inz_in (and q_in>0 i=inz_in+1) */
      int inz_out;  /* q_out<=0 at i=inz_out (and q_out>0 i=inz_out-1) */
      int i1, i2;   /* zero occurs between index i1 and i2 in domain dom */
      double Tlam1, Tlam2; /* zero occurs between Tlam=Tlam1 and Tlam=Tlam2 */
      double Tlam0;        /* q=0 at Tlam=Tlam0, Tlam \in [-1,1] */
      double lam0;         /* q=0 at lam=lam0 in domain dom */

      /* find indices where q_in and q_out switch sign */
      for(i=1; i<n1in; i++) if(q_in[Index(i,j,k)]<=0.0) break;
      inz_in=i;
      for(i=0; i<n1; i++)   if(q_out[Index(i,j,k)]<=0.0) break;
      inz_out=i;

      /* if inz_in<n1in, q has zero in inner domain */
      /* if inz_out<n1,  q is negative in outer domain */
      if(inz_in<n1in)     { i1=inz_in-1;  i2=inz_in;  dom=innerdom; }
      else if(inz_out==0) { i1=i2=0;                  dom=outerdom; }
      else if(inz_out<n1) { i1=inz_out-1; i2=inz_out; dom=outerdom; }
      else
      {
        printf("reset_Coordinates_CubedSphere_sigma01: innerdom=%d  A=%g B=%g  "
               "inz_in=%d inz_out=%d\n", innerdom, A,B, inz_in,inz_out);
        printf("q_in[Index(0,j,k)]=%g\n", q_in[Index(0,j,k)]);
        printf("q_in[Index(n1in-1,j,k)]=%g\n", q_in[Index(n1in-1,j,k)]);
        printf("q_out[Index(0,j,k)]=%g\n", q_out[Index(0,j,k)]);
        printf("C1=%.13g  C2=%.13g\n", Getd("DNSdata_C1"), Getd("DNSdata_C2"));
        quick_Vars_output(box, 
        "Coordinates_CubedSphere_sigma01 DNSdata_q", 666,666);
        errorexit("reset_Coordinates_CubedSphere_sigma01: q>0 everywhere???");
      }

      /* initial bracket for Tlam */
      boxq0 = grid->box[dom];
      Tlam1 = boxq0->v[iX][Index(i1,j,k)];
      Tlam2 = boxq0->v[iX][Index(i2,j,k)];
      if(dom == innerdom)
      {
        Tlam1 = Tlam1 - 1.;
        Tlam2 = Tlam2 - 1.;
      }

      if(i1 == i2)  /* zero is between inner and outer box (see above) */
      {
        Tlam1 = -0.01;
        Tlam2 = +0.01;
      }
      /* set pars */
      pars->grid = grid;
      pars->A = A;
      pars->B = B;
      pars->b = outerdom;
      pars->b_in = innerdom;


      /* use Brent's method to find Tlam0 where q=0 */
      /* zbrac_P may not be needed as Tlam1/2 should already bracket Tlam0 */
      if(zbrac_P(q_of_Tlam_forgiven_AB_ZP, &Tlam1,&Tlam2, (void *) pars)<0)
      {
        printf("reset_Coordinates_CubedSphere_sigma01: outerdom=%d  A=%g B=%g  "
               "inz_in=%d inz_out=%d\n", outerdom, A,B, inz_in,inz_out);
        printf("reset_Coordinates_CubedSphere_sigma01: innerdom=%d  j=%d k=%d "
               "pars->b=%d pars->b_in=%d\n", innerdom,j,k, pars->b,pars->b_in);
        printf("dom=%d i1=%d i2=%d Tlam1=%g Tlam2=%g\n", dom, i1,i2, Tlam1,Tlam2);
        printf("q_in[Index(0,j,k)]=%g\n", q_in[Index(0,j,k)]);
        printf("q_in[Index(n1in-1,j,k)]=%g\n", q_in[Index(n1in-1,j,k)]);
        printf("q_out[Index(0,j,k)]=%g\n", q_out[Index(0,j,k)]);
        printf("new q(Tlam1)=%g, new q(Tlam2)=%g\n",
               q_of_Tlam_forgiven_AB_ZP(Tlam1, (void *) pars),
               q_of_Tlam_forgiven_AB_ZP(Tlam2, (void *) pars));
        printf("new q(-0.1)=%g, new q(0.1)=%g\n",
               q_of_Tlam_forgiven_AB_ZP(-0.1, (void *) pars),
               q_of_Tlam_forgiven_AB_ZP(0.1, (void *) pars));
        printf("C1=%.13g  C2=%.13g\n", Getd("DNSdata_C1"), Getd("DNSdata_C2"));
        quick_Vars_output(box, 
        "Coordinates_CubedSphere_sigma01 DNSdata_q", 666,666);
        errorexit("cannot find bracket for q_of_Tlam_forgiven_AB_ZP");
      }
      stat=zbrent_itsP(&Tlam0, q_of_Tlam_forgiven_AB_ZP,  Tlam1,Tlam2,
                       (void *) pars, itmax, tol);
      if(Tlam0<0.) { lam0 = 1. + Tlam0;  dom = innerdom; }
      else         { lam0 = Tlam0;       dom = outerdom; }

      /* now that we have lam0 we can find x,y,z */
      boxq0 = grid->box[dom];
      x = boxq0->x_of_X[1]((void *)boxq0, -1, lam0, A,B);
      y = boxq0->x_of_X[2]((void *)boxq0, -1, lam0, A,B);
      z = boxq0->x_of_X[3]((void *)boxq0, -1, lam0, A,B);
      xc = boxq0->CI->xc[1];

      /* from x,y,z we now get the new sigma in box outerdom and innerdom */
      sig01_AB = sqrt((x-xc)*(x-xc) + y*y + z*z);
      //printf("lam0=%g x=%g y=%g z=%g sig01_AB=%g\n", lam0, x,y,z, sig01_AB);

      /* set sigma = sig01_AB in both domains */
      boxnew->v[isigma0][ind]     = sig01_AB;
      boxnewin->v[isigma1][indin] = sig01_AB;
    } /* end forplane1 */
  }
  /* ensure that sigma is continuous between boxes on gridnew */
  if(Getv("DNSdata_CubSph_sigma_continuity","yes"))
    DNSgrid_Coordinates_CubSph_sigma_continuity(gridnew, star);

  /* compute sigma and its derives from sigma_def (from CI->iFS) on gridnew */
  DNS_set_sigma01_and_derivs(gridnew, star);
}

/* compute sigma_dA, sigma_dB for all boxes with star surface */
void DNS_set_sigma01_and_derivs(tGrid *grid, int star)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int si;

    /* do nothing for other star and all boxes that do not touch surface */
    if(box->SIDE != star || box->BOUND != SSURF) continue;

    /* find index si sigma in this box */
    if(box->CI->type==innerCubedSphere)      si=0;
    else if(box->CI->type==outerCubedSphere) si=1;
    else errorexit("there should only be outer or inner Cubed Spheres");

    /* set sigma and its derivs */
    init_CubedSphere_from_CI_iFS(box, si);
  }
}



/******************************************************************/
/* useful functions                                               */
/******************************************************************/

/* compute volume integral of var with index vind inside STAR1 or STAR2.
   Here any Psi^6 needs to be already included in the var we integrate. */
double InnerVolumeIntegral(tGrid *grid, int star, int vind)
{
  double VolInt = 0.0;
  int b;
  //pr_DNS_box_attribs(grid);
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if( (box->SIDE == star) && (box->MATTR == INSIDE) )
      VolInt += BoxVolumeIntegral(box, vind);
  }
  return VolInt;
}

/* surface integral over star surface: integrate over lam=0 surface
   of MATTR == TOUCH boxes */
double StarSurfaceIntegral(tGrid *grid, int star, int vind)
{
  double SurfInt = 0.0;
  int b;
  int ig;

  /* do we use physical or flat metric for Surface Integrals? */
  if(Getv("DNSdata_StarSurfaceIntegral_metric","flat")) ig = -1;
  else                                                  ig = Ind("gxx");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    double *Integ = malloc(n1*n2*n3 * sizeof(double));
    double *v     = box->v[vind];

    if( (box->SIDE == star) && (box->MATTR == TOUCH) )
    {
      spec_SurfaceIntegral(box, ig, 1, v, Integ);
      SurfInt += Integ[0]; /* ind=0 is at lam=0 */
    }
    free(Integ);
  }
  return SurfInt;
}

/* set bfaces, set pointlist only if set_fpts=1.  */
int DNSgrid_set_bfaces(tGrid *grid, int set_fpts, int pr)
{
  int b;

  if(Getv("DNSdata_grid","2starcubes"))
  {
    /* set bface structures */
    set_touching_bfaces_of_boxes_with_same_facepoints(grid, 0, grid->nboxes);
    set_all_bfaces_with_ob_minus1_to_outerbound(grid, 0, grid->nboxes);
    /* set ofi and bit fields */
    set_ofi_in_all_bfaces(grid);
    set_bits_in_all_bfaces(grid);
    if(pr) forallboxes(grid, b) printbfaces(grid->box[b]);
  }
  return 0;
}



/* return the box index and the coords of the point at the Cartesian x,y,z */
/* initially box, *X,*Y,*Z contain the box and the coords of the
   point on the other grid */
int DNSgrid_Get_BoxAndCoords_of_xyz(tGrid *grid1,
                                    double *X1, double *Y1, double *Z1,
                                    intList *bl1, double x, double y, double z)
{
  int b1;
  if( bl1==NULL || bl1->n==0 )
    b1 = b_XYZ_of_xyz(grid1, X1,Y1,Z1, x,y,z);
  else
    b1 = b_XYZ_of_xyz_inboxlist(grid1, bl1->e,bl1->n, X1,Y1,Z1, x,y,z);
  return b1;
}


/* Interpolate Var with index vind from grid1 to grid2 */
void Interpolate_Var_From_Grid1_To_Grid2(tGrid *grid1, tGrid *grid2, int vind)
{
  Interp_Var_From_Grid1_To_Grid2_star(grid1, grid2, vind, STAR1);
  Interp_Var_From_Grid1_To_Grid2_star(grid1, grid2, vind, STAR2);
}
/* Interpolate Var with index vind from grid1 to grid2 
   for domains that touch surface of star */
void Interp_Var_From_Grid1_To_Grid2_star(tGrid *grid1, tGrid *grid2, int vind,
                                         int star)
{
  int cind = Ind("temp1");
  int Xind = Ind("X");
  int Yind = Ind("Y");
  int Zind = Ind("Z");
  int xind = Ind("x");
  int yind = Ind("y");
  int zind = Ind("z");
  intList *bl1 = alloc_intList(); /* list that contains boxes to look in */
  int b,i;

  /* save coeffs of vind on grid1 in cind = Ind("temp1") */
  forallboxes(grid1, b)
  {
    tBox *box = grid1->box[b];
    /* do nothing if we are on wrong side of grid or not near star surface */
    if(box->SIDE!=star || box->BOUND!=SSURF) continue;

    spec_Coeffs(box, box->v[vind], box->v[cind]);
    push_intList(bl1, box->b); /* look in box->b */
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

    /* do nothing if we are on wrong side of grid or not near star surface */
    if(box->SIDE!=star || box->BOUND!=SSURF) continue;

    //printCI(box);
    //printf("sigma=%g\n", box->v[37][7]);
    //printf("sigma=%g\n", CubedSphere_sigma(box, 1, 7, -1., -1.));
    //printf("sigma=%g\n", CubedSphere_sigma(box, 1, -1, -1., -1.));

    /* Here we use SGRID_LEVEL6orTOP_Pragma(omp parallel), since both
       DNSgrid_Get_BoxAndCoords_of_xyz and spec_interpolate should be
       thread safe. So we may not want a separate grid copy per thread. */
#undef GRIDperTHREAD_Interp_Var_From_Grid1_To_Grid2_star
    //define GRIDperTHREAD_Interp_Var_From_Grid1_To_Grid2_star 1
    SGRID_LEVEL6orTOP_Pragma(omp parallel)
    {
#ifdef GRIDperTHREAD_Interp_Var_From_Grid1_To_Grid2_star
      tGrid *grid1_p = make_empty_grid(grid1->nvariables, 0);
      copy_grid(grid1, grid1_p, 0);
#else
      tGrid *grid1_p = grid1;
#endif
      /* start loop */
      SGRID_LEVEL6orTOP_Pragma(omp for) 
      forallpoints(box,i)
      {
        double X = pX[i];
        double Y = pY[i];
        double Z = pZ[i];
        double x = px[i];
        double y = py[i];
        double z = pz[i];
        int b1;

        /* get b1, X,Y,Z on grid1_p */
        //printf("b=%d i=%d x=%g y=%g z=%g  guess X=%g Y=%g Z=%g\n",
        //box->b, i ,x,y,z, X,Y,Z);
        b1 = DNSgrid_Get_BoxAndCoords_of_xyz(grid1_p, &X,&Y,&Z, bl1, x,y,z);
        if(b1>=0)
        {
          /* get var at point X,Y,Z by interpolation */
          pv[i] = spec_interpolate(grid1_p->box[b1],
                                   grid1_p->box[b1]->v[cind], X,Y,Z);
          if(!finite(pv[i]))
          {
            printf("box->b=%d x=%g y=%g z=%g  b1=%d X=%.18g Y=%.18g Z=%.18g\n",
                    box->b, x,y,z, b1, X,Y,Z);
            errorexit("!finite(pv[i])");
          }
        }
        else /* point x,y,z may be beyond outer boundary */
        {
          errorexit("could not find point");
        }
        //printf("  b1=%d X=%g Y=%g Z=%g  pv[i]=%g\n", b1, X,Y,Z, pv[i]);
        //printf("  box=%p grid1->box[b1]=%p grid1_p->box[b1=%p\n",
        //       box, grid1->box[b1], grid1_p->box[b1]);
        //if(i>40) exit(88);

      } /* end forallpoints loop */
#ifdef GRIDperTHREAD_Interp_Var_From_Grid1_To_Grid2_star
      /* free local grid copy */
      free_grid(grid1_p);
#endif
    }
  } /* end forallboxes(grid2,b) */
  free_intList(bl1);
}

/* copy variable at i=n1-1 from Box1 to i=0 in Box2 */
/* This can be used to make a var continues across inner and outer domains. */
void copy_Var_Box1ATlam1_to_Box2ATlam0(tGrid *grid, int vind, int b1, int b2)
{
  int i,j,k, ijk1, ijk2; 
  int n1_1 = grid->box[b1]->n1;
  int n2_1 = grid->box[b1]->n2;
  int n3_1 = grid->box[b1]->n3;
  int n1_2 = grid->box[b2]->n1;
  int n2_2 = grid->box[b2]->n2;
  int n3_2 = grid->box[b2]->n3;
  double *v1 = grid->box[b1]->v[vind];
  double *v2 = grid->box[b2]->v[vind];
  if(n2_2!=n2_1 || n3_2!=n3_1)
    errorexit("(n2,n3) has to be equal in both boxes");
  forplane1(i,j,k, n1_2,n2_2,n3_2, 0)
  {
    ijk1 = Ind_n1n2(n1_1-1,j,k, n1_1,n2_1);
    ijk2 = Ind_n1n2(0,j,k, n1_2,n2_2);
    v2[ijk2] = v1[ijk1];
  }
}


/* extrapolate Sigma to outside pof star in a very simple way.
   This will cause jumps in Sigma's derivs! But Sigma should be C0 */
void set_outside_DNSdata_Sigma_by_S0Extrap(tGrid *grid, int star)
{
  int b;
  int iSig = Ind("DNSdata_Sigma");
  int iX   = Ind("X");
  int iY   = Ind("Y");
  int iZ   = Ind("Z");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = grid->box[b]->n1;
    int n2 = grid->box[b]->n2;
    int n3 = grid->box[b]->n3;
    double *X   = box->v[iX];
    double *Y   = box->v[iY];
    double *Z   = box->v[iZ];
    double *Sig = box->v[iSig];
    double *Sig_in;
    int b_in, ijk, ijk_in;
    tBox *box_in;
    int n1_in, n2_in, n3_in;

    /* do nothing if we are away from our star */
    if((box->SIDE != STAR1) || (box->MATTR == AWAY)) continue;

    /* here we could find val of Sigma in center of cube inside star,
       if we need it */
    //if(box->COORD == CART) Sigc=...;

    /* do nothing more if we are not in a box that touches the star */
    if(!(box->MATTR == TOUCH)) continue;

    b_in = b-6; /* works only for my CubSph arrangement */
    box_in = grid->box[b_in];
    n1_in = box_in->n1;
    n2_in = box_in->n2;
    n3_in = box_in->n3;

    if(n2_in!=n2 || n3_in!=n3)
      errorexit("(n2,n3) has to be equal in both boxes");

    Sig_in = box_in->v[iSig];

    forallpoints(box, ijk)
    {
      int k = kOfInd_n1n2(ijk,n1,n2);
      int j = jOfInd_n1n2_k(ijk,n1,n2,k);
      double Si1, Si2, a,b, lam,A,B, r, drdlam, drdlam_in;

      /* our coords in box */
      lam = X[ijk];
      A   = Y[ijk];
      B   = Z[ijk];

      /* get Sigma at both ends of inner box */
      ijk_in = Ind_n1n2(0,j,k, n1_in,n2_in);
      Si1 = Sig_in[ijk_in];
      ijk_in = Ind_n1n2(n1_in-1,j,k, n1_in,n2_in);
      Si2 = Sig_in[ijk_in];

      /* get drdlam, drdlam_in at star surface */
      r_dr_dlam_of_lamAB_CubSph(box, Index(0,j,k), 0.,A,B, &r, &drdlam);
      r_dr_dlam_of_lamAB_CubSph(box_in, ijk_in, 1.,A,B, &r, &drdlam_in);

      /* get slope and intercept */
      a = (Si2 - Si1);        /* dlam here is 1 */
      a *= drdlam/drdlam_in;  /* we want d/dr to be the same in both sides */
      b = Si2;
      Sig[ijk] = a * lam + b;
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
  int star;
  double fac = vec[1];
  double m0;
  t_grid_b_star_A_B_struct *pars;

  /* get pars */
  pars = (t_grid_b_star_A_B_struct *) p;
  grid = pars->grid;
  star = pars->star;

  /* scale sigma */
  DNSgrid_scale_Coordinates_CubSph_sigma(grid, fac, star);
  DNSgrid_init_Coords(grid);
  /* set wB */
  if(star==STAR1) DNS_set_wB(grid, star, Getd("DNSdata_xmax1"),0.0,0.0);
  else            DNS_set_wB(grid, star, Getd("DNSdata_xmax2"),0.0,0.0);
  /* get new mass */
  m0 =  GetInnerRestMass(grid, star);
  printf("m0errOFsigmafac_VectorFuncP: fac=%g => m0%d=%g\n",
         fac, star, m0);
  fflush(stdout);
  /* set sigma back to what it was */
  DNSgrid_scale_Coordinates_CubSph_sigma(grid, 1.0/fac, star);
  DNSgrid_init_Coords(grid);
  /* set wB */
  if(star==STAR1) DNS_set_wB(grid, star, Getd("DNSdata_xmax1"),0.0,0.0);
  else            DNS_set_wB(grid, star, Getd("DNSdata_xmax2"),0.0,0.0);
  /* set fvec */
  if(star==STAR1)
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
  int j;
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
  vlpush(varlist, Ind("Coordinates_CubedSphere_sigma01"));

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
    /* set wB */
    DNS_set_wB(grid, STAR1, Getd("DNSdata_actual_xmax1"),0.0,0.0);
    DNS_set_wB(grid, STAR2, Getd("DNSdata_actual_xmax2"),0.0,0.0);
  }
  else /* set DNSdata_b=DNSdata_b_sav and adjust other pars */
  {
    double m01 = Getd("DNSdata_m01");
    double m02 = Getd("DNSdata_m02");
    double r, rc, dr, dOm;
    double facvec[2];
    int stat, check;
    t_grid_b_star_A_B_struct pars[1];

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
    DNS_set_wB(grid, STAR1, Getd("DNSdata_xmax1"),0.0,0.0);
    DNS_set_wB(grid, STAR2, Getd("DNSdata_xmax2"),0.0,0.0);
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
    pars->star = STAR1;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    DNSgrid_scale_Coordinates_CubSph_sigma(grid, facvec[1], STAR1);
    DNSgrid_init_Coords(grid);
    /* set wB */
    DNS_set_wB(grid, STAR1, Getd("DNSdata_xmax1"),0.0,0.0);

    /* star2 */
    pars->grid = grid;
    pars->star = STAR2;
    facvec[1] = 1.0;
    stat = newton_linesrch_itsP(facvec, 1, &check, m0errOFsigmafac_VectorFuncP,
                               (void *) pars, 100, 0.01*m02);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
    /* use factor to scale sigma */
    DNSgrid_scale_Coordinates_CubSph_sigma(grid, facvec[1], STAR2);
    DNSgrid_init_Coords(grid);
    /* set wB */
    DNS_set_wB(grid, STAR2, Getd("DNSdata_xmax2"),0.0,0.0);
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


/************************************************************************/
/* utilities to manipulate the grid */
/************************************************************************/
/* given a new cubed sphere sigma0/1 and its derivs
   initialize Coordinates. */
void DNSgrid_init_Coords_for_star(tGrid *grid, int star)
{
  int Coordinates_verbose = Getv("Coordinates_verbose", "yes");

  /* avoid too much printf */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "no");

  /* initialize coords on grid, reset x,y,z, dXdx and such */
  init_CoordTransform_And_Derivs(grid);

  /* put back original Coordinates_verbose */
  if(Coordinates_verbose) Sets("Coordinates_verbose", "yes");
}
/* given a new cubed sphere sigma0/1 and its derivs
   initialize Coordinates. */
void DNSgrid_init_Coords(tGrid *grid)
{
  DNSgrid_init_Coords_for_star(grid, -1); /* init both sides of grid */
}


/* scale Coordinates_AnsorgNS_sigma and its deriv on one side
   by a factor fac. */
void DNSgrid_scale_Coordinates_CubSph_sigma(tGrid *grid, double fac, int star)
{
  int b;

  if(star<STAR1 || star>STAR2)
    errorexit("must have STAR1<=star<=STAR2");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    tBox *boxin;
    int innerdom;
    int i,j,k, n1,n2,n3;
    int isigdef0, isigdef1;
    double *sigdef;

    if(box->SIDE!=star || box->MATTR!=TOUCH) continue;

    innerdom = b-6; /* works only for my CubSph arrangement */
    boxin = grid->box[innerdom];
    isigdef0    = box->CI->iFS[0];
    isigdef1    = boxin->CI->iFS[1];

    /* scale sigdef in outer box b */
    sigdef      =  box->v[isigdef0];
    n1=box->n1;
    n2=box->n2;
    n3=box->n3;
    forplane1(i,j,k, n1,n2,n3, 0)
    {
      int ind = Index(i,j,k);
      sigdef[ind] *= fac;
    }

    /* scale sigdef in innerdom */
    sigdef      =  boxin->v[isigdef1];
    n1=boxin->n1;
    n2=boxin->n2;
    n3=boxin->n3;
    forplane1(i,j,k, n1,n2,n3, n1-1)
    {
      int ind = Index(i,j,k);
      sigdef[ind] *= fac;
    }
  }
  /* set sigma and its derivs */
  DNS_set_sigma01_and_derivs(grid, star);
}

/* ensure that sigma is continuous between boxes on gridnew */
/* right now we copy values between bfaces, because averaging does not
   work when a edge or corner point in a bface is not paired with some
   bfaces that contain points from other boxes at the same x,y,z */
void DNSgrid_Coordinates_CubSph_sigma_continuity(tGrid *grid, int star)
{
  int b;

  if(star<STAR1 || star>STAR2)
    errorexit("must have STAR1<=star<=STAR2");

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1,n2,n3;
    int isigma, iosigma;
    double *sigma;
    int fi;

    /* do nothing except for cubed sph. at surface of this star */
    if(box->SIDE!=star || box->BOUND!=SSURF) continue;

    /* find index of sigma in this box */
    if(box->CI->type==innerCubedSphere)
      isigma    = box->CI->iFS[0];
    else if(box->CI->type==outerCubedSphere)
      isigma    = box->CI->iFS[1];
    else
      errorexit("there should only be outer or inner Cubed Spheres");

    /* get sigma pointer */
    sigma    = box->v[isigma];
    n1=box->n1;
    n2=box->n2;
    n3=box->n3;

    forallbfaces(box, fi)
    {
      tBface *bface = box->bface[fi];
      double *osigma;
      int pi, ijk;
      tBface *obface;
      tBox *obox;
      int ob  = bface->ob;
      int ofi = bface->ofi;

      if(ob<0 || ofi<0) continue;
      obox = grid->box[ob];
      obface = obox->bface[ofi];
      if(obface->fpts==NULL)  continue;
      if(bface->same_fpts==0) continue;

      /* find index and poiters of sigma in other box */
      if(obox->CI->type==innerCubedSphere)
        iosigma    = obox->CI->iFS[0];
      else if(obox->CI->type==outerCubedSphere)
        iosigma    = obox->CI->iFS[1];
      else
        continue;
      osigma    = obox->v[iosigma];

      /* go over point pairs */
      forPointList_inbox(bface->fpts, box, pi, ijk)
      {
        int k = kOfInd_n1n2(ijk,n1,n2);
        int j = jOfInd_n1n2_k(ijk,n1,n2,k);
        int i = iOfInd_n1n2_jk(ijk,n1,n2,j,k);

        /* go to next point if it is not in curved inner or outer face */
        if(box->CI->type==outerCubedSphere && i<n1-1) continue;
        if(box->CI->type==innerCubedSphere && i>0)    continue;

        if(j==0 || j==n2-1 || k==0 || k==n3-1)
        {
          int oijk = obface->fpts->point[ob][pi];
          
          /* propagate this sigma to all other boxes */
          osigma[oijk] = sigma[ijk];
        }
      }
    } /* end forallbfaces */
  }
}
