/* DNSdata.c */
/* Wolfgang Tichy 2007 */


#include "sgrid.h"
#include "DNSdata.h"

#define Power pow

/* struct types used in root finder newton_linesrch_itsP */
typedef struct T_grid_bXYZ1_bXYZ2_struct {
  tGrid *grid; /* grid */
  int b1;      /* box1 */
  double X1 ;  /* X pos 1 */
  double Y1 ;  /* Y pos 1 */
  double Z1 ;  /* Z pos 1 */
  int b2;      /* box2 */
  double X2 ;  /* X pos 2 */
  double Y2 ;  /* Y pos 2 */
  double Z2 ;  /* Z pos 2 */
} t_grid_bXYZ1_bXYZ2_struct;

typedef struct T_grid_grid0_m01_m02_struct {
  tGrid *grid;  /* grid where we operate */
  tGrid *grid0; /* grid from which we interpolate vars when domains change */
  double m01;   /* desired m01 */
  double m02;   /* desired m02 */
  int    bi1;   /* location of max1 */
  double Xmax1;
  double Ymax1;
  int    bi2;   /* location of max2 */
  double Xmax2;
  double Ymax2;
  double qm1;   /* desired qm1 = max of q1 */
  double qm2;   /* desired qm2 */
} t_grid_grid0_m01_m02_struct;


/* global vars */
extern double rf_surf1; /* radius of star1 */
extern double rf_surf2; /* radius of star2 */
extern double P_core1;  /* core pressure of star1 */
extern double P_core2;  /* core pressure of star2 */
extern tParameter *pdb;
extern int npdb;
extern int npdbmax;
tGrid *central_q_errors_VectorFunc__grid; /* grid var for central_q_errors_VectorFunc */
tGrid *xmaxs_error_VectorFunc__grid; /* grid var for xmaxs_error_VectorFunc */
double xmaxs_error_VectorFunc__xmax1;  /* xmax1 we currently try to achieve */
double xmaxs_error_VectorFunc__xmax2;  /* xmax2 we currently try to achieve */
tGrid *dqdx_at_Xmax1_2_VectorFunc__grid;  /* grid var for dqdx_at_Xmax1_2_VectorFunc */
int    dqdx_at_Xmax1_2_VectorFunc__bi1; /* boxind of max1 */
int    dqdx_at_Xmax1_2_VectorFunc__bi2; /* boxind of max2 */
double dqdx_at_Xmax1_2_VectorFunc__Xmax1; /* pos. of max1 */
double dqdx_at_Xmax1_2_VectorFunc__Ymax1; /* pos. of max1 */
double dqdx_at_Xmax1_2_VectorFunc__Xmax2; /* pos. of max2 */
double dqdx_at_Xmax1_2_VectorFunc__Ymax2; /* pos. of max2 */
              


/* global var lists */
tVarList *vlu, *vlFu, *vluDerivs;
tVarList *vldu, *vlJdu, *vlduDerivs;


/* functions in this file */
void make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(tGrid *grid,
     tVarList **vlw,  tVarList **vlwDerivs,  tVarList **vlFw, 
     tVarList **vldw, tVarList **vldwDerivs, tVarList **vlJdw, char *Name);
void free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(
     tVarList *vlw,  tVarList *vlwDerivs,  tVarList *vlFw,
     tVarList *vldw, tVarList *vldwDerivs, tVarList *vlJdw);
double GridL2Norm_of_vars_in_string(tGrid *grid, char *str);
double GridL2Norm_of_vars_in_string_withZeroErr_outside(tGrid *grid, char *str);
int DNS_Eqn_Iterator_for_vars_in_string(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr, char *str);
int DNS_ordered_Eqn_Iterator(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
int DNS_Eqn_Iterator(tGrid *grid, int itmax, double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr);
void compute_new_q_and_adjust_domainshapes_InterpFromGrid0(tGrid *grid, 
                                                           tGrid *grid0, 
                                                           int star);
void compute_new_q_and_adjust_domainshapes(tGrid *grid, int star);
void m01_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m02_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m01_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
void m02_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
void qm1_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
void qm2_error_VectorFuncP(int n, double *vec, double *fvec, void *p);
double m01_guesserror_ZP(double C, void *p);
double m02_guesserror_ZP(double C, void *p);
double m01_error_ZP(double C, void *p);
double m02_error_ZP(double C, void *p);
double qm1_error_ZP(double C, void *p);
double qm2_error_ZP(double C, void *p);
int find_Varmax_along_x_axis_in_star(tGrid *grid, int varind, int star,
                                     int *bi, double *X, double *vmax);
int find_qmax_along_x_axis(tGrid *grid, int star,
                           int *bi, double *X, double *qmax);
void central_q_errors_VectorFunc(int n, double *vec, double *fvec);
void estimate_q_errors_VectorFunc(int n, double *vec, double *fvec);
double DNSdata_find_position_of_qmax(tGrid *grid, int star, int *bi,
                                     double *X, double *Y, double *Z);
double DNSdata_find_xyz_of_qmax(tGrid *grid, int star, int *bi,
                                double *x, double *y, double *z);
void set_DNSdata_actual_xyzmax_pars(tGrid *grid);





/* set conformal rotational velocity for one star */
void DNS_set_wB(tGrid *grid, int star, double xc,double yc,double zc)
{
  int corot1 = Getv("DNSdata_rotationstate1","corotation");
  int corot2 = Getv("DNSdata_rotationstate2","corotation");
  int wBfac_Psi6    = Getv("DNSdata_wB_factor","Psi6");
  int wBfac_h       = Getv("DNSdata_wB_factor","h");
  int wBfac_ooalpha = Getv("DNSdata_wB_factor","1/alpha");
  int b;
  double omegax1   = Getd("DNSdata_omegax1");
  double omegay1   = Getd("DNSdata_omegay1");
  double omegaz1   = Getd("DNSdata_omegaz1");
  double omegax2   = Getd("DNSdata_omegax2");
  double omegay2   = Getd("DNSdata_omegay2");
  double omegaz2   = Getd("DNSdata_omegaz2");
  int corot;
  double omegax, omegay, omegaz;

  if(star==STAR1)
  {
    corot=corot1;
    omegax=omegax1;  omegay=omegay1;  omegaz=omegaz1;
  }
  else if(star==STAR2)
  {
    corot=corot2;
    omegax=omegax2;  omegay=omegay2;  omegaz=omegaz2;
  }
  else
    errorexit("DNS_set_wB: star must be 1 or 2.");

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int i;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];
    double *DNSdata_wBx = box->v[Ind("DNSdata_wBx")];
    double *DNSdata_wBy = box->v[Ind("DNSdata_wBx")+1];
    double *DNSdata_wBz = box->v[Ind("DNSdata_wBx")+2];
    double *DNSdata_Psi = box->v[Ind("DNSdata_Psi")];
    double *DNSdata_alphaP = box->v[Ind("DNSdata_alphaP")];
    double *DNSdata_q = box->v[Ind("DNSdata_q")];

    if( box->SIDE!=star ) continue;
    
    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }

      /* set wB */
      if(!corot)
      {
        double Psi1 = DNSdata_Psi[i];
        double Psi4 = Psi1*Psi1*Psi1*Psi1;
        double Psi6 = Psi1*Psi1*Psi4;
        double ooalpha = Psi1/DNSdata_alphaP[i]; /* 1/alpha */
        double h = DNSdata_q[i] + 1.0; /* h = q + 1 */
        // double u0;
        double vx,vy,vz;
        double Att1, wBfac;
        double lam = pX[i];

        if(box->MATTR==AWAY) Att1=0.0; //1.0-Attenuation01((lam-0.1)/0.8, 2.0, 0.5);
        else                 Att1=1.0;

        /* omega cross r-rc */
        vx = ( omegay* (z-zc) - omegaz* (y-yc) )*Att1;
        vy = ( omegaz* (x-xc) - omegax* (z-zc) )*Att1;
        vz = ( omegax* (y-yc) - omegay* (x-xc) )*Att1;
        /* 1/u0^2 ~ alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),*/
        /* u0 = 1.0/sqrt(fabs(exp(2*Phi1) - Psi4*(vx*vx + vy*vy + vz*vz))); */

        /* set wBfac */
        /* wBfac = h*u0 * Psi^6  ???  */
        wBfac = 1.0;
        if(wBfac_Psi6)    wBfac = wBfac * Psi6;
        if(wBfac_h)       wBfac = wBfac * h;
        if(wBfac_ooalpha) wBfac = wBfac * ooalpha; /* u0 ~ 1/alpha */

        DNSdata_wBx[i] = vx * wBfac;
        DNSdata_wBy[i] = vy * wBfac;
        DNSdata_wBz[i] = vz * wBfac;
      }
      else
      {
        DNSdata_wBx[i] = 0.0; 
        DNSdata_wBy[i] = 0.0;
        DNSdata_wBz[i] = 0.0;
      }
    } /* end forallpoints */
  }
}


/* initial shift as in gr-qc/0007028v2*/
void DNSdata_initial_shift(int star, double fac, 
                           double m1, double m2, double Omega, double r12,
                           double rs, double x, double y, double z,
                           double *Bx, double *By, double *Bz)
{
  double F, epsW, epsXi, epsxy, Wy, dWydx,dWydy,dWydz, dXidx,dXidy,dXidz;
  double r = sqrt(x*x + y*y + z*z);
  double nx = x;
  double ny = y;
  double nz = z;
  double xz_on = 0.0;

  /* vector n = x/r */
  if(r>0) { nx=x/r; ny=y/r; nz=z/r; }

  /* function F */
  if(star==1) { F=fac*m1*Omega*r12/(1.0+m1/m2); epsW=+1; epsXi=+1; epsxy=-1; }
  else        { F=fac*m2*Omega*r12/(1.0+m2/m1); epsW=+1; epsXi=+1; epsxy=+1; }

  /* vector W , Wx = Wz = 0.0. And derivs of W and derivs of Xi */
  if(r<rs)
  {
    Wy = (6.0*F/rs)*( 1.0 - r*r/(3*rs*rs) )*epsW;
    dWydx = ( -(4.0*F/(rs*rs))*(r/rs)*nx )*epsW;
    dWydy = ( -(4.0*F/(rs*rs))*(r/rs)*ny )*epsW;
    dWydz = ( -(4.0*F/(rs*rs))*(r/rs)*nz )*epsW;
    dXidx = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * nx )*epsXi;
    dXidy = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * ny
              +((2.0*F)/rs) * (1.0 - 3.0*r*r/(5.0*rs*rs)) )*epsXi;
    dXidz = ( -((12.0*F)/(5.0*rs)) * ((r*y)/(rs*rs)) * nz )*epsXi;
  }
  else
  {
    Wy = (4.0*F/r)*epsW;
    dWydx = ( (-4.0*F/(r*r))*nx )*epsW;    
    dWydy = ( (-4.0*F/(r*r))*ny )*epsW;    
    dWydz = ( (-4.0*F/(r*r))*nz )*epsW;    
    dXidx = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*nx )*epsXi;
    dXidy = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*ny
              +((4.0*F)/5.0) * ((rs*rs)/(r*r*r)) )*epsXi;
    dXidz = ( -((12.0*F)/(5.0*r)) * ((rs*rs)/(r*r)) * ny*nz )*epsXi;
  }

  /* shift */
  *Bx = (          - 0.125*(dXidx + y*dWydx) )*epsxy *xz_on;
  *By = ( 0.875*Wy - 0.125*(dXidy + y*dWydy) )*epsxy;
  *Bz =            - 0.125*(dXidz + y*dWydz) *xz_on;
}

/* initialize DNSdata */
int DNSdata_startup(tGrid *grid)
{
  int rot1 = Getv("DNSdata_rotationstate1","rotation");
  int rot2 = Getv("DNSdata_rotationstate2","rotation");
  int b;
  int TOVav        = Getv("DNSdata_guess", "TOVaverage");
  int TOVprod      = Getv("DNSdata_guess", "TOVproduct");
  int initShift    = Getv("DNSdata_guess", "TaniguchiShift");
  int initFromChkp = Getv("DNSdata_guess", "initialize_from_checkpoint");
  double Omega     = Getd("DNSdata_Omega");
  double omegax1   = Getd("DNSdata_omegax1");
  double omegay1   = Getd("DNSdata_omegay1");
  double omegaz1   = Getd("DNSdata_omegaz1");
  double omegax2   = Getd("DNSdata_omegax2");
  double omegay2   = Getd("DNSdata_omegay2");
  double omegaz2   = Getd("DNSdata_omegaz2");
  double xCM       = Getd("DNSdata_x_CM");
  double DNSdata_b = Getd("DNSdata_b");
  double rs1, m1, Phic1, Psic1, m01;
  double rs2, m2, Phic2, Psic2, m02;
  double xc1 = DNSdata_b;
  double xc2 = -DNSdata_b;
  double xout1 = DNSdata_b + rf_surf1;
  double xin1  = DNSdata_b - rf_surf1;
  double xin2  = -DNSdata_b + rf_surf2;
  double xout2 = -DNSdata_b - rf_surf2;
  double ysh1;

  printf("Initializing DNSdata:\n");
  prTimeIn_s("WallTime: ");

  /* set boundary information: farlimit, falloff, propagation speed */
  VarNameSetBoundaryInfo("DNSdata_Psi",   1, 1, 1.0);
  VarNameSetBoundaryInfo("DNSdata_Bx",    0, 1, 1.0);
  VarNameSetBoundaryInfo("DNSdata_By",    0, 1, 1.0);
  VarNameSetBoundaryInfo("DNSdata_Bz",    0, 1, 1.0);
  VarNameSetBoundaryInfo("DNSdata_alphaP",1, 1, 1.0);
  VarNameSetBoundaryInfo("DNSdata_Sigma", 0, 1, 1.0);

  /* enable all DNSdata vars */
  enablevar(grid, Ind("DNSdata_Psi"));
  enablevar(grid, Ind("DNSdata_Psix"));
  enablevar(grid, Ind("DNSdata_Psixx"));
  enablevar(grid, Ind("DNSdata_Bx"));
  enablevar(grid, Ind("DNSdata_Bxx"));
  enablevar(grid, Ind("DNSdata_Bxxx"));
  enablevar(grid, Ind("DNSdata_alphaP"));
  enablevar(grid, Ind("DNSdata_alphaPx"));
  enablevar(grid, Ind("DNSdata_alphaPxx"));
  enablevar(grid, Ind("DNSdata_Sigma"));
  enablevar(grid, Ind("DNSdata_Sigmax"));
  enablevar(grid, Ind("DNSdata_Sigmaxx"));
  enablevar(grid, Ind("DNSdata_SigmaX"));
  enablevar(grid, Ind("DNSdata_SigmaXX"));
  enablevar(grid, Ind("DNSdata_SigmaXXX"));
  enablevar(grid, Ind("DNSdata_lSigmaX"));
  enablevar(grid, Ind("DNSdata_lSigmaXX"));
  enablevar(grid, Ind("DNSdata_lSigmaXXX"));
  enablevar(grid, Ind("DNSdata_rhobar"));
  enablevar(grid, Ind("DNSdata_wBx"));
  enablevar(grid, Ind("DNSdata_q"));
  enablevar(grid, Ind("DNSdata_qg"));
  enablevar(grid, Ind("DNSdata_wBxx"));
  enablevar(grid, Ind("DNSdata_qx"));
  enablevar(grid, Ind("DNSdata_VRx"));
  enablevar(grid, Ind("DNSdata_temp1"));
  enablevar(grid, Ind("DNSdata_temp2"));
  enablevar(grid, Ind("DNSdata_temp3"));
  enablevar(grid, Ind("DNSdata_temp4"));
  enablevar(grid, Ind("DNSdata_Psiold"));
  enablevar(grid, Ind("DNSdata_Boldx"));
  enablevar(grid, Ind("DNSdata_alphaPold"));
  enablevar(grid, Ind("DNSdata_Sigmaold"));
  enablevar(grid, Ind("DNSdata_qgold"));
  enablevar(grid, Ind("DNSdata_qnocent"));
  //enablevar(grid, Ind("DNSdata_qcorot"));
  //enablevar(grid, Ind("DNSdata_surface_sigma_pm"));

  /* enable some lapse and shift of ADMvars */
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));

  /* set rs, m, Phic, Psic, m0 for both stars */
  TOV_init(P_core1, 1, &rs1, &m1, &Phic1, &Psic1, &m01);
  TOV_init(P_core2, 1, &rs2, &m2, &Phic2, &Psic2, &m02);

  /* set qmax1/2 */
  Setd("DNSdata_qmax1", DNS_polytrope_hm1_of_P(P_core1));
  Setd("DNSdata_qmax2", DNS_polytrope_hm1_of_P(P_core2));
  /* set cart positions of qmax1/2 */
  if(Getd("DNSdata_xmax1")<=0.0) Setd("DNSdata_xmax1", xc1);
  if(Getd("DNSdata_xmax2")>=0.0) Setd("DNSdata_xmax2", xc2);

  /* set the var DNSdata_CoordFacPower */
  enablevar(grid, Ind("DNSdata_CoordFac"));
  set_DNSdata_CoordFac(grid);

  /* set bfaces */
  DNSgrid_set_bfaces(grid, 1, 1);

  /* do nothing if BNSdata_Interpolate_pointsfile exists */
  if(GetsLax("BNSdata_Interpolate_pointsfile")!=0) return 0;

  /* load data from some old checkpoint file */
  if(initFromChkp && GetsLax("outdir_previous_iteration")!=NULL)
  {
    char filename[10000];
    snprintf(filename, 9999, "%s/%s", 
             Gets("outdir_previous_iteration"), "checkpoint.0");
    prdivider(1);
    printf("loading initial guess from %s\n", filename);
    DNSgrid_load_initial_guess_from_checkpoint(grid, filename);
  }
  else if(initFromChkp && strlen(Gets("DNSdata_initfile"))>0)
  {
    char filename[10000];
    snprintf(filename, 9999, "%s", Gets("DNSdata_initfile"));
    prdivider(1);
    printf("loading initial guess from %s\n", filename);
    DNSgrid_load_initial_guess_from_checkpoint(grid, filename);
  }
  else /* use some TOV data */
  {
    /* get yshift1 (for testing) */
    if(m02==0.0) ysh1 = Getd("DNSdata_yshift1");
    else         ysh1 = 0.0;

    /* set initial values in all in boxes */
    forallboxes(grid,b)
    {  
      tBox *box = grid->box[b];
      int i;
      double *pX = box->v[Ind("X")];
      double *pY = box->v[Ind("Y")];
      double *pZ = box->v[Ind("Z")];
      double *px = box->v[Ind("x")];
      double *py = box->v[Ind("y")];
      double *pz = box->v[Ind("z")];
      double *DNSdata_Psi    = box->v[Ind("DNSdata_Psi")];
      double *DNSdata_alphaP = box->v[Ind("DNSdata_alphaP")];
      double *DNSdata_q      = box->v[Ind("DNSdata_q")];
      double *DNSdata_Bx     = box->v[Ind("DNSdata_Bx")];
      double *DNSdata_By     = box->v[Ind("DNSdata_Bx")+1];
      double *DNSdata_Bz     = box->v[Ind("DNSdata_Bx")+2];
      double *DNSdata_Sigma  = box->v[Ind("DNSdata_Sigma")];
      double *DNSdata_wBx    = box->v[Ind("DNSdata_wBx")];
      double *DNSdata_wBy    = box->v[Ind("DNSdata_wBx")+1];
      double *DNSdata_wBz    = box->v[Ind("DNSdata_wBx")+2];
      double r1, m1_r, P1, Phi1, Psi1, q1;
      double r2, m2_r, P2, Phi2, Psi2, q2;
      double Bx1,By1,Bz1;
      double Bx2,By2,Bz2;

      forallpoints(box,i)
      {
        double x = pX[i];
        double y = pY[i];
        double z = pZ[i];

        if(px!=NULL) 
        {
          x = px[i];
          y = py[i];
          z = pz[i];
        }
        /* set Psi, alphaP, q */
        if(TOVav || TOVprod || box->SIDE==STAR1)
        {
          r1 = sqrt((x-xc1)*(x-xc1) + (y-ysh1)*(y-ysh1) + z*z);
          TOV_m_P_Phi_Psi_OF_rf(r1, rs1, m1, P_core1, Phic1, Psic1,
                                &m1_r, &P1, &Phi1, &Psi1);
          q1 = DNS_polytrope_hm1_of_P(P1);
          DNSdata_initial_shift(1, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs1, 
                                -(x-xc1), -(y-ysh1), z,  &Bx1,&By1,&Bz1);
        }
        else { m1_r = P1 = Phi1 = q1 = Bx1=By1=Bz1 = 0.0;  Psi1 = 1.0; }
        if(TOVav || TOVprod || box->SIDE==STAR2)
        {
          r2 = sqrt((x-xc2)*(x-xc2) + y*y + z*z);
          TOV_m_P_Phi_Psi_OF_rf(r2, rs2, m2, P_core2, Phic2, Psic2,
                                &m2_r, &P2, &Phi2, &Psi2);
          q2 = DNS_polytrope_hm1_of_P(P2);
          DNSdata_initial_shift(2, 1.0, m1,m2, Omega, fabs(xc1-xc2), rs2, 
                                x-xc2, y, z,  &Bx2,&By2,&Bz2);
        }
        else { m2_r = P2 = Phi2 = q2 = Bx2=By2=Bz2 = 0.0;  Psi2 = 1.0; }
        
        /* set the data */
        if(TOVprod)
        {
          DNSdata_Psi[i]   = Psi1*Psi2;
          DNSdata_alphaP[i]= exp(Phi1+Phi2)*Psi1*Psi2;
          DNSdata_q[i]     = q1 + q2;
        }
        else /* for TOVav or TOV */
        {
          DNSdata_Psi[i]   = Psi1 + Psi2 - 1.0;
          DNSdata_alphaP[i]= exp(Phi1)*Psi1 + exp(Phi2)*Psi2 - 1.0;
          DNSdata_q[i]     = q1 + q2;
        }
        /* set inertial shift B^i if wanted */
        if(initShift)
        {
          DNSdata_Bx[i] = Bx1 + Bx2;
          DNSdata_By[i] = By1 + By2;
          DNSdata_Bz[i] = Bz1 + Bz2;
        }
        /* set Sigma and wB if needed */
        if( (box->SIDE==STAR1) && (rot1) )
        {
          double Att;
          double lam = pX[i];
          if(box->MATTR==AWAY) Att=0.0; //1.0-Attenuation01((lam-0.1)/0.8, 2.0, 0.5);
          else                 Att=1.0;

          DNSdata_Sigma[i] = Omega*(xc1-xCM) * y * Att;
        }
        if( (box->SIDE==STAR2) && (rot2) )
        {
          double Att;
          double lam = pX[i];
          if(box->MATTR==AWAY) Att=0.0; //1.0-Attenuation01((lam-0.1)/0.8, 2.0, 0.5);
          else                 Att=1.0;

          DNSdata_Sigma[i] = Omega*(xc2-xCM) * y * Att;
        }
      } /* end forallpoints */
    }

    /* set wB in both stars */
    DNS_set_wB(grid, STAR1, xc1,0.0,0.0);
    DNS_set_wB(grid, STAR2, xc2,0.0,0.0);

    /* if NS1 is shifted in y-direc. (for testing) adjust grid on right side */
    if(ysh1 != 0.0 && Getv("DNSdata_adjustdomain01","yes"))
    {
      /* adjust grid so that new q=0 is at A=0 */
      compute_new_q_and_adjust_domainshapes(grid, STAR1);

      /* set q to zero if q<0, and also in region 1 & 2 */
      set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);
    }
  } /* end intialization using TOV data */

  /* print out maxima */
  printf("DNSdata_startup: DNSdata_qmax1 = %g  DNSdata_qmax2 = %g\n"
         "                 DNSdata_xmax1 = %g  DNSdata_xmax2 = %g\n",
         Getd("DNSdata_qmax1"), Getd("DNSdata_qmax2"),
         Getd("DNSdata_xmax1"), Getd("DNSdata_xmax2"));

  /* recompute q from the other fields */
  if(Getv("DNSdata_init_q_fromfields", "yes"))
  {
    /* adjust grid so that new q=0 is at A=0 */
    compute_new_q_and_adjust_domainshapes(grid, STAR1);
    compute_new_q_and_adjust_domainshapes(grid, STAR2);
  
    /* set q to zero if q<0, and also in region 1 & 2 */
    set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);
  }

  /* set qg=q, qgold=qg */
  varcopy(grid, Ind("DNSdata_qg"), Ind("DNSdata_q"));  
  varcopy(grid, Ind("DNSdata_qgold"), Ind("DNSdata_qg"));  

  /* set DNSdata_actual_xyzmax pars */
  set_DNSdata_actual_xyzmax_pars(grid);

/*
for(b=1; b<13; b++)
{
tBox *box = grid->box[b];
int i;
double *pX = box->v[Ind("X")];
double *pY = box->v[Ind("Y")];
double *pZ = box->v[Ind("Z")];
double *px = box->v[Ind("x")];
double *py = box->v[Ind("y")];
double *pz = box->v[Ind("z")];
double *DNSdata_Sigma  = box->v[Ind("DNSdata_Sigma")];
for(i=0; i<6*6*6; i++)
{
double dx = px[i]-20.;
double dy = py[i];
double dz = pz[i];
double r = sqrt(dx*dx + dy*dy + dz*dz);
DNSdata_Sigma[i] = b + 100*b*r;
if(i==0) printf("b=%d r=%g S=%g\n", b, r, DNSdata_Sigma[i]);
}
}
*/

  prTimeIn_s("WallTime: ");
  return 0;
}


/* functions to shift/translate a field/var */
void DNS_translate_var_in_box(tBox *box, int iv, int ivx,
                              double dx, double dy, double dz)
{
  int b = box->b;
  double *var = box->v[iv];
  double *dvdx = box->v[ivx];
  double *dvdy = box->v[ivx+1];
  double *dvdz = box->v[ivx+2];
  int i;

  /* save time if dx=dy=dz=0 */
  if(dequal(dx,0.0) && dequal(dy,0.0) && dequal(dz,0.0))  return;

  /* get gradient of var */
  FirstDerivsOf_S(box, iv, ivx);
  
  /* v_translated(x) = v(x-dx) ~ v(x) - [vq(x)/dx] dx */
  //printf("b=%d: dx=%g dy=%g dz=%g\n", b, dx,dy,dz);
  forallpoints(box, i)
    var[i] = var[i] - (dvdx[i]*dx + dvdy[i]*dy + dvdz[i]*dz);
}
/* translate var on one side of grid specified by star */
void DNS_translate_var_aroundStar(tGrid *grid, int star, int iv, int ivx,
                                  double dx, double dy, double dz)
{
  int b;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    /* continue with next box if on wrong side */
    if(box->SIDE!=star) continue;
    DNS_translate_var_in_box(box, iv, ivx, dx,dy,dz);
  }
}   

/* center Psi, beta^i, alphaP, Sigma around each star */
int DNSdata_center_fields_if_desired(tGrid *grid, int it)
{
  if(    !Getv("DNSdata_center_fields", "no")    &&
     it>=Geti("DNSdata_center_fields_first_at")  && 
         Geti("DNSdata_center_fields_first_at")>=0 )
  {
    int b,i;
    double m01, m02;
    int iPsi = Ind("DNSdata_Psi");
    int iBx  = Ind("DNSdata_Bx");
    int iBy  = Ind("DNSdata_By");
    int iBz  = Ind("DNSdata_Bz");
    int ialP = Ind("DNSdata_alphaP");
    int iSig = Ind("DNSdata_Sigma");
    int itx  = Ind("DNSdata_temp1");
    double dx1,dy1,dz1, dx2,dy2,dz2, x1,y1,z1, x2,y2,z2;
    double fac = Getd("DNSdata_center_fields_fac");
    double maxdx = Getd("DNSdata_center_fields_maxdx");
    double maxdy = Getd("DNSdata_center_fields_maxdy");
    double maxdz = Getd("DNSdata_center_fields_maxdz");
    double xmax1 = Getd("DNSdata_xmax1");
    double xmax2 = Getd("DNSdata_xmax2");
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("DNSdata_center_fields", "center_yz"))  xon=0;

    /* get global max of q in NS1/2 */
    x1 = Getd("DNSdata_actual_xmax1");
    y1 = Getd("DNSdata_actual_ymax1"); 
    z1 = Getd("DNSdata_actual_zmax1");
    x2 = Getd("DNSdata_actual_xmax2");
    y2 = Getd("DNSdata_actual_ymax2");
    z2 = Getd("DNSdata_actual_zmax2");

    printf("Centering fields:  DNSdata_center_fields_fac = %s\n"
           " DNSdata_center_fields = %s\n",
           Gets("DNSdata_center_fields_fac"),
           Gets("DNSdata_center_fields"));

    /* get dx1,dy1,dz1 and center around star1 */
    dx1 = -fac*(x1-xmax1);
    dy1 = -fac*y1;
    dz1 = -fac*z1;
    dx1 = dx1*xon*(fabs(dx1)>maxdx);  /* <-- switch centering on/off */
    dy1 = dy1*yon*(fabs(dy1)>maxdy);
    dz1 = dz1*zon*(fabs(dz1)>maxdz);
    DNS_translate_var_aroundStar(grid,1, iPsi,itx, dx1,dy1,dz1);
    DNS_translate_var_aroundStar(grid,1, iBx ,itx, dx1,dy1,dz1);
    DNS_translate_var_aroundStar(grid,1, iBy ,itx, dx1,dy1,dz1);
    DNS_translate_var_aroundStar(grid,1, iBz ,itx, dx1,dy1,dz1);
    DNS_translate_var_aroundStar(grid,1, ialP,itx, dx1,dy1,dz1);
    DNS_translate_var_aroundStar(grid,1, iSig,itx, dx1,dy1,dz1);

    /* get dx2,dy2,dz2 and center around star2 */
    dx2 = -fac*(x2-xmax2);
    dy2 = -fac*y2;
    dz2 = -fac*z2;
    dx2 = dx2*xon*(fabs(dx2)>maxdx);  /* <-- switch centering on/off */
    dy2 = dy2*yon*(fabs(dy2)>maxdy);
    dz2 = dz2*zon*(fabs(dz2)>maxdz);
    DNS_translate_var_aroundStar(grid,2, iPsi,itx, dx2,dy2,dz2);
    DNS_translate_var_aroundStar(grid,2, iBx ,itx, dx2,dy2,dz2);
    DNS_translate_var_aroundStar(grid,2, iBy ,itx, dx2,dy2,dz2);
    DNS_translate_var_aroundStar(grid,2, iBz ,itx, dx2,dy2,dz2);
    DNS_translate_var_aroundStar(grid,2, ialP,itx, dx2,dy2,dz2);
    DNS_translate_var_aroundStar(grid,2, iSig,itx, dx2,dy2,dz2);

    /* save time if dx1=dy1=dz1=dx2=dy2=dz2=0 */
    if(dequal(dx1,0.0) && dequal(dy1,0.0) && dequal(dz1,0.0) &&
       dequal(dx2,0.0) && dequal(dy2,0.0) && dequal(dz2,0.0)    ) return 0;

    /* recalc q as well? */
    if(Getv("DNSdata_center_fields", "reset_q"))
    {
      if(Getv("DNSdata_center_fields", "adjust_domainshapes"))
      {
        compute_new_q_and_adjust_domainshapes(grid, STAR1);
        compute_new_q_and_adjust_domainshapes(grid, STAR2);
      }
      else
      {
        DNS_compute_new_centered_q(grid, STAR1);
        DNS_compute_new_centered_q(grid, STAR2);
      }

      /* set q to zero if q<0 or in region 1 and 2 */
      set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);
    }
   
    /* set actual positions of maxima again */
    set_DNSdata_actual_xyzmax_pars(grid);

    /* print new masses */
    m01 = GetInnerRestMass(grid, STAR1);
    m02 = GetInnerRestMass(grid, STAR2);
    printf("     => m01=%.19g m02=%.19g\n", m01, m02);
  }
  return 0;
}


/* compute new q on both sides of the grid */
void DNS_compute_new_q(tGrid *grid, int iq)
{
  DNS_compute_new_q_instar(grid, STAR1, iq);
  DNS_compute_new_q_instar(grid, STAR2, iq);
}

/* functions to compute new q from fields, and to also shift them so
   that the max in q is centered where we want it on the x-axis */
void DNS_compute_new_centered_q(tGrid *grid, int star)
{
  int iq = Ind("DNSdata_q");
  int iqg= Ind("DNSdata_qg");

//printf("DNS_compute_new_centered_q: C1=%.13g  C2=%.13g\n",
//Getd("DNSdata_C1"), Getd("DNSdata_C2"));
//quick_Vars_output(grid->box[2], "DNSdata_q",41,41);
  DNS_compute_new_q_instar(grid, star, iq);
//quick_Vars_output(grid->box[2], "DNSdata_q",55,55);

  if(Getv("DNSdata_center_new_q_flag", "yes"))
  {
    int iqx= Ind("DNSdata_qx");
    int b, i;
    int bi1, bi2;
    double dx,dy,dz, xm,ym,zm;
    double fac = Getd("DNSdata_center_new_q_fac");
    double maxdx = Getd("DNSdata_center_new_q_maxdx");
    double maxdy = Getd("DNSdata_center_new_q_maxdy");
    double maxdz = Getd("DNSdata_center_new_q_maxdz");
    double xmax;
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("DNSdata_center_new_q", "center_yz"))  xon=0;

    if(star==STAR1)
    {
      xmax = Getd("DNSdata_xmax1");
      /* get global max of q in NS1/2 */
      xm = Getd("DNSdata_actual_xmax1");
      ym = Getd("DNSdata_actual_ymax1");
      zm = Getd("DNSdata_actual_zmax1");
    }
    else if(star==STAR2)
    {
      xmax = Getd("DNSdata_xmax2");
      /* get global max of q in NS1/2 */
      xm = Getd("DNSdata_actual_xmax2");
      ym = Getd("DNSdata_actual_ymax2");
      zm = Getd("DNSdata_actual_zmax2");
    }
    else
      errorexit("star has to be STAR1 or STAR2");

    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      double *DNSdata_q = box->v[iq];
      double *dqdx = box->v[iqx];
      double *dqdy = box->v[iqx+1];
      double *dqdz = box->v[iqx+2];

      if(box->SIDE != star || box->MATTR == AWAY) continue;

      /* q_centered(x) = q(x+dx) ~ q(x) + [dq(x)/dx] dx */
      dx = xm-xmax; dy = ym; dz = zm;
      dx = dx*xon*(fabs(dx)>maxdx);  /* <-- switch centering on/off */
      dy = dy*yon*(fabs(dy)>maxdy);
      dz = dz*zon*(fabs(dz)>maxdz);

      /* save time if dx=dy=dz=0 */
      if(dequal(dx,0.0) && dequal(dy,0.0) && dequal(dz,0.0))  continue;

      /* get gradient of q */
      FirstDerivsOf_S(box, iq, iqx);

      //printf("b=%d: dx=%g dy=%g dz=%g\n", b, dx,dy,dz);
      forallpoints(box, i)
        DNSdata_q[i] = DNSdata_q[i] + 
                       fac*(dqdx[i]*dx + dqdy[i]*dy + dqdz[i]*dz);
    }
  }
  varcopy(grid, iqg, iq); /* set qg=q */
}
double DNS_compute_new_centered_q_atXYZ(tGrid *grid, int bi,
                                        double X, double Y, double Z)
{
  double q;
  q = DNS_compute_new_q_atXYZ(grid,bi, X,Y,Z);
  if(Getv("DNSdata_center_new_q_flag", "yes"))
  {
    int iq = Ind("DNSdata_temp4");
    int iqx= Ind("DNSdata_qx");
    tBox *box = grid->box[bi];
    double *cx= box->v[Ind("DNSdata_temp1")];
    double *cy= box->v[Ind("DNSdata_temp2")];
    double *cz= box->v[Ind("DNSdata_temp3")];
    double dx,dy,dz, xm,ym,zm;
    double dqdx, dqdy, dqdz;
    double fac = Getd("DNSdata_center_new_q_fac");
    double maxdx = Getd("DNSdata_center_new_q_maxdx");
    double maxdy = Getd("DNSdata_center_new_q_maxdy");
    double maxdz = Getd("DNSdata_center_new_q_maxdz");
    double xmax;
    int xon=1;
    int yon=1;
    int zon=1;
    if(Getv("DNSdata_center_new_q", "center_yz"))  xon=0;

    if(box->SIDE==STAR1)
    {
      /* desired xmax1 */
      xmax = Getd("DNSdata_xmax1");
      /* get global max of q in NS1/2 */
      xm = Getd("DNSdata_actual_xmax1");
      ym = Getd("DNSdata_actual_ymax1"); 
      zm = Getd("DNSdata_actual_zmax1");
    }
    else
    {
      /* desired xmax2 */
      xmax = Getd("DNSdata_xmax2");
      /* get global max of q in NS1/2 */
      xm = Getd("DNSdata_actual_xmax2");
      ym = Getd("DNSdata_actual_ymax2");
      zm = Getd("DNSdata_actual_zmax2");
    }

    /* q_centered(x) = q(x+dx) ~ q(x) + [dq(x)/dx] dx */
    dx = xm-xmax;
    dy = ym;
    dz = zm;
    dx = dx*xon*(fabs(dx)>maxdx);  /* <-- switch centering on/off */
    dy = dy*yon*(fabs(dy)>maxdy);
    dz = dz*zon*(fabs(dz)>maxdz);

    /* save time if dx=dy=dz=0 */
    if(dequal(dx,0.0) && dequal(dy,0.0) && dequal(dz,0.0))
      return q;

    /* get gradient of q and its coeffs, set coeffs of dq in 
       DNSdata_temp1/2/3 */
    DNS_compute_new_q(grid, iq);
    FirstDerivsOf_S(box, iq, iqx);
    spec_Coeffs(box, box->v[iqx], cx);
    spec_Coeffs(box, box->v[iqx+1], cy);
    spec_Coeffs(box, box->v[iqx+2], cz);

    /* find dqdx, dqdy, dqdz by interpolation */
    dqdx = spec_interpolate(box, cx, X,Y,Z);
    dqdy = spec_interpolate(box, cy, X,Y,Z);
    dqdz = spec_interpolate(box, cz, X,Y,Z);

    /* q_centered(x) = q(x+dx) ~ q(x) + [dq(x)/dx] dx */
//printf("bi=%d: dx=%g dy=%g dz=%g\n", bi, dx,dy,dz);
    q = q + fac*(dqdx*dx + dqdy*dy + dqdz*dz);
  }
  return q;
}

/* center q if needed, depending on how pars are set */
int DNSdata_center_q_if_desired(tGrid *grid, int it)
{
  if(    !Getv("DNSdata_center_new_q", "no")    &&
     it>=Geti("DNSdata_center_new_q_first_at")  && 
         Geti("DNSdata_center_new_q_first_at")>=0 )
  {
    int b,i;
    double m01, m02;

    printf("Centering q:  DNSdata_center_new_q_fac = %s\n"
           " DNSdata_center_new_q = %s\n",
           Gets("DNSdata_center_new_q_fac"), Gets("DNSdata_center_new_q"));
    Sets("DNSdata_center_new_q_flag", "yes");  /* activate centering of q */
    if(Getv("DNSdata_center_new_q", "adjust_domainshapes"))
    {
      compute_new_q_and_adjust_domainshapes(grid, STAR1);
      compute_new_q_and_adjust_domainshapes(grid, STAR2);
    }
    else 
    {
      DNS_compute_new_centered_q(grid, STAR1);
      DNS_compute_new_centered_q(grid, STAR2);
    }
    Sets("DNSdata_center_new_q_flag", "no");  /* deactivate centering of q */
    /* set q to zero if q<0 or in region 1 and 2 */
    set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

    /* set actual positions of maxima again ? */
    if(Getv("DNSdata_center_new_q", "set_DNSdata_actual_xyzmax"))
      set_DNSdata_actual_xyzmax_pars(grid);

    /* print new masses */
    m01 = GetInnerRestMass(grid, STAR1);
    m02 = GetInnerRestMass(grid, STAR2);
    printf("     => m01=%.19g m02=%.19g\n", m01, m02);
  }
  return 0;
}


/* find both qmax and reset DNSdata_qmax1/2, DNSdata_xmax1/2 accordingly */
void find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars(tGrid *grid)
{
  int bi1, bi2;
  double xmax1, qmax1, xmax2, qmax2;
  double Xmax1, Xmax2;

  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qmax1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qmax2);
  /* compute qmax1 and qmax2 */
  qmax1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,0.,0.);
  qmax2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,0.,0.);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,0.,0.);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,0.,0.);
  else
    xmax2 = Xmax2;
  /* set qmax1/2 */
  Setd("DNSdata_qmax1", qmax1);
  Setd("DNSdata_qmax2", qmax2);
  /* set cart positions of qmax1/2 */
  Setd("DNSdata_xmax1", xmax1);
  Setd("DNSdata_xmax2", xmax2);
  printf("find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars:\n");
  printf("  resetting: DNSdata_qmax1 = %g  DNSdata_qmax2=%g\n"
         "             DNSdata_xmax1 = %g  DNSdata_xmax2=%g\n",
         Getd("DNSdata_qmax1"), Getd("DNSdata_qmax2"),
         Getd("DNSdata_xmax1"), Getd("DNSdata_xmax2"));
}





/* set the pars DNSdata_desired_VolAvSigma12 to values that min BC violation */
void set_DNSdata_desired_VolAvSigma12_toMinBCerr(tGrid *grid, int index_Sigma)
{
  int b, ijk;
  int AddInnerVolIntToBC=Getv("DNSdata_Sigma_surface_BCs","AddInnerVolIntToBC");
  int InnerVolIntZero  = Getv("DNSdata_Sigma_surface_BCs","InnerVolIntZero");
  int AddInnerSumToBC  = Getv("DNSdata_Sigma_surface_BCs","AddInnerSumToBC");
  int InnerSumZero     = Getv("DNSdata_Sigma_surface_BCs","InnerSumZero");
  int SigmaZeroAtPoint = Getv("DNSdata_Sigma_surface_BCs","ZeroAtPoint");
  double *Sigma;
  double VolAvSigma1=0., VolAvSigma2=0.;

  /* do nothing? */
  if(Getv("DNSdata_set_desired_VolAvSigmas", "no")) return;

  /* set VolAvSigma1/2 */
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int n1 = box->n1;
    int n2 = box->n2;
    int n3 = box->n3;
    int MATTRinside = (box->MATTR== INSIDE);
    int hasSSURF    = (box->BOUND== SSURF);
    //int isXinDom    = (box->CI->dom == box->SIDE - STAR1);
    //int isVolAvBox  = (MATTRinside && hasSSURF && isXinDom);
    int isCube = (box->CI->type == 0);
    int isVolAvBox  = (MATTRinside && isCube);

    if(isVolAvBox)
    {
      double VolAvSigma = 0.;

      Sigma = box->v[index_Sigma];
      /* use VolInt in some cases */
      if (AddInnerVolIntToBC || InnerVolIntZero)
        VolAvSigma = BoxVolumeIntegral(box, index_Sigma);
      else if(AddInnerSumToBC || InnerSumZero)
        forallpoints(box, ijk) VolAvSigma += Sigma[ijk];
      else  /* otherwise value at one point */
      {
        ijk = Index(n1/2, n2/2, n3/2);
        VolAvSigma = Sigma[ijk];
      }
      /* set VolAvSigma1/2 */
      if(box->SIDE==STAR1) VolAvSigma1 = VolAvSigma;
      else                 VolAvSigma2 = VolAvSigma;
    }
  }
  Setd("DNSdata_desired_VolAvSigma1", VolAvSigma1);
  Setd("DNSdata_desired_VolAvSigma2", VolAvSigma2);
  printf(" setting: DNSdata_desired_VolAvSigma1 / 2 = %g / %g\n",
  VolAvSigma1, VolAvSigma2);
}


/* backup grid,pdb to grid_bak,pdb_bak.
   But do it only if DNSdata_domainshape_diff_tol<1e30. 
   Call as:    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak); */
void backup_grid_pdb(tGrid *grid, tParameter *pdb,
                     tGrid *grid_bak, tParameter *pdb_bak)
{
  ///* do nothing if we tolerate large differences */
  //if(Getd("DNSdata_domainshape_diff_tol")>=1e30) return;

  /* make exact copies of grid and pdb */
  copy_grid(grid, grid_bak, 0);
  copy_pdb(pdb, npdb, pdb_bak);
}




/* Adjust C1/2 and thus q by demanding that m01 and m02 stay the same. */
int adjust_C1_C2_q_keep_restmasses(tGrid *grid, int it, double tol)
{
  double Cvec[3];
  double m0errorvec[3];
  double m01, m02, m0_error, dm01, dm02;
  int check, stat, bi, i;
  tGrid *grid_bak;
  tParameter *pdb_bak;
  t_grid_grid0_m01_m02_struct pars[1];

  /* grid and pdb for backups */
  grid_bak = make_empty_grid(grid->nvariables, 0);
  pdb_bak  = make_empty_pdb(npdbmax);

  /* rest masses before adjusting q */
  m01 = GetInnerRestMass(grid, STAR1);
  m02 = GetInnerRestMass(grid, STAR2);
  printf("adjust_C1_C2_q_keep_restmasses: in DNSdata_solve step %d: "
         "WallTime=%gs\n", it, getTimeIn_s());
  printf(" rest mass in inner domains before computing new q:"
         " m01=%g m02=%g\n", m01, m02);

  /* compute error in masses */
  dm01 = m01 - Getd("DNSdata_m01");
  dm02 = m02 - Getd("DNSdata_m02");
  m0_error = dm01*dm01 + dm02*dm02;
  m0_error = sqrt(m0_error)/(Getd("DNSdata_m01")+Getd("DNSdata_m02"));
  printf(" rest mass error = %.4e "
         "(before adjusting C1/2)\n", m0_error);

  /* set desired masses for this iteration */
  pars->m01 = Getd("DNSdata_m01"); // + dm01*0.9;
  pars->m02 = Getd("DNSdata_m02"); // + dm02*0.9;
  printf(" adjusting q,C1,C2 to achieve: m01=%g  m02=%g  tol*0.01=%g\n",
         pars->m01, pars->m02, tol*0.01);

  /* print C1/2 we used before */
  printf(" old: DNSdata_C1=%g DNSdata_C2=%g\n",
         Getd("DNSdata_C1"), Getd("DNSdata_C2"));
  /* printf("     => m01=%g m02=%g\n", m01, m02); */

  /* adjust C1/2 */
  if(!Getv("DNSdata_adjust_C1C2", "no"))
  {
    /* refine guesses for C1,C2? */
    if(Getv("DNSdata_adjust_C1C2", "refineguess"))
    {
      /* choose C1/2 such that rest masses are not too big or too small */
      for(i=0; i<1000; i++)
      {
        DNS_compute_new_centered_q(grid, STAR1);
        DNS_compute_new_centered_q(grid, STAR2);
        m01 = GetInnerRestMass(grid, STAR1);
        m02 = GetInnerRestMass(grid, STAR2);

        check = 0;

        if(m01 > 1.1*Getd("DNSdata_m01"))
          { Setd("DNSdata_C1", 0.999*Getd("DNSdata_C1"));  check=1; }
        else if(m01 < 0.9*Getd("DNSdata_m01"))
          { Setd("DNSdata_C1", 1.002*Getd("DNSdata_C1"));  check=1; }
    
        if(m02 > 1.1*Getd("DNSdata_m02") && Getd("DNSdata_m02")>0)
          { Setd("DNSdata_C2", 0.999*Getd("DNSdata_C2"));  check=1; }
        else if(m02 < 0.9*Getd("DNSdata_m02") && Getd("DNSdata_m02")>0)
          { Setd("DNSdata_C2", 1.002*Getd("DNSdata_C2"));  check=1; }

        if(check==0) break;
      }

      /* refine guess for C1/2 */
      pars->grid = grid;
      Cvec[1] = Getd("DNSdata_C1");
      stat = newton_linesrch_itsP(Cvec, 1, &check, m01_guesserror_VectorFuncP,
                                  (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
      if(check || stat<0)
        printf("  --> C1 guess: check=%d stat=%d\n", check, stat);
      Setd("DNSdata_C1", Cvec[1]);

      Cvec[1] = Getd("DNSdata_C2");
      if(Getd("DNSdata_m02")>0)
        stat = newton_linesrch_itsP(Cvec, 1, &check, m02_guesserror_VectorFuncP,
                                    (void *) pars, 30, max2(m0_error*0.1, tol*0.1));
      if(check || stat<0)
        printf("  --> C2 guess: check=%d stat=%d\n", check, stat);
      Setd("DNSdata_C2", Cvec[1]);

      /* print guess for C1/2 */
      printf(" guess: DNSdata_C1=%g DNSdata_C2=%g\n",
             Getd("DNSdata_C1"), Getd("DNSdata_C2"));
    }
    if(Getv("DNSdata_adjust_C1C2", "zbrent"))
    {
      /***************************************************************/
      /* do zbrent_itsP iterations of Cvec until m01/2 error is zero */
      /***************************************************************/

      /* backup grid,pdb */
      backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
      pars->grid = grid;
      if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
      else  pars->grid0 = grid;

      /* adjust C1 and thus m01 */
      Cvec[1] = Getd("DNSdata_C1"); /* initial guess */
      Cvec[0] = Cvec[1]*1.0001;     /* lower bracket bound (note C<0) */
      Cvec[2] = Cvec[1]*0.9999;     /* upper bracket bound (note C<0) */
      stat = zbrac_P(m01_error_ZP, Cvec, Cvec+2, (void *) pars);
      if(stat<0) errorexit("cannot find bracket for m01_error_ZP");
      stat = zbrent_itsP(Cvec+1, m01_error_ZP, Cvec[0], Cvec[2],
                         (void *) pars, 1000, tol*1e-4);
      if(stat<0) printf("  --> stat=%d\n", stat);
      Setd("DNSdata_C1", Cvec[1]);
      compute_new_q_and_adjust_domainshapes(grid, STAR1);

      /* backup grid,pdb */
      backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
      pars->grid = grid;
      if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
      else  pars->grid0 = grid;

      /* adjust C2 and thus m02 */
      Cvec[1] = Getd("DNSdata_C2"); /* initial guess */
      Cvec[0] = Cvec[1]*1.0001;     /* lower bracket bound (note C<0) */
      Cvec[2] = Cvec[1]*0.9999;     /* upper bracket bound (note C<0) */
      stat = zbrac_P(m02_error_ZP, Cvec, Cvec+2, (void *) pars);
      if(stat<0) errorexit("cannot find bracket for m02_error_ZP");
      stat = zbrent_itsP(Cvec+1, m02_error_ZP, Cvec[0], Cvec[2],
                         (void *) pars, 1000, tol*1e-4);
      if(stat<0) printf("  --> stat=%d\n", stat);
      Setd("DNSdata_C2", Cvec[1]);
      compute_new_q_and_adjust_domainshapes(grid, STAR2);
    }
    else
    {
      /***********************************************************************/
      /* do newton_linesrch_itsP iterations of Cvec until m0errorvec is zero */
      /***********************************************************************/

      /* backup grid,pdb */
      backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
      pars->grid = grid;
      if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
      else  pars->grid0 = grid;

      /* adjust C1 and thus m01 */
      Cvec[1] = Getd("DNSdata_C1");
      stat = newton_linesrch_itsP(Cvec, 1, &check, m01_error_VectorFuncP,
                                  (void *) pars, 1000, tol*0.01);
      if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
      Setd("DNSdata_C1", Cvec[1]);

      /* backup grid,pdb */
      backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
      pars->grid = grid;
      if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
      else  pars->grid0 = grid;

      /* adjust C2 and thus m02 */
      Cvec[1] = Getd("DNSdata_C2");
      if(Getd("DNSdata_m02")>0)
        stat = newton_linesrch_itsP(Cvec, 1, &check, m02_error_VectorFuncP,
                                    (void *) pars, 1000, tol*0.01);
      if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
      Setd("DNSdata_C2", Cvec[1]);
    }

    printf("adjust_C1_C2_q_keep_restmasses:\n");
    printf(" new: DNSdata_C1=%g DNSdata_C2=%g\n",
           Getd("DNSdata_C1"), Getd("DNSdata_C2"));
  }
  else /* keep C1/2, but adjust q */
  {
    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* compute new q in star1 */
    compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, pars->grid0, STAR1);

    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* compute new q in star2 */
    compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, pars->grid0, STAR2);
  
    printf("adjust_C1_C2_q_keep_restmasses: adjusted q, but kept C1/2\n");
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  /* print new masses */
  m01 = GetInnerRestMass(grid, STAR1);
  m02 = GetInnerRestMass(grid, STAR2);
  printf("     => m01=%.19g m02=%.19g\n", m01, m02);

  /* free grid and pdb for backups */
  free_pdb(pdb_bak, npdb);
  free_grid(grid_bak);
        
  return 0;
}

/* Adjust C1/2 and thus q by demanding that qmax1 and qmax2 stay the same. */
int adjust_C1_C2_q_keep_qmax(tGrid *grid, int it, double tol)
{
  double Cvec[3];
  double qmerrorvec[3];
  double qm1, qm2, qm_error, dqm1, dqm2;
  double Xmax1, Xmax2;
  int bi1, bi2;
  int check, stat, bi, i;
  tGrid *grid_bak;
  tParameter *pdb_bak;
  t_grid_grid0_m01_m02_struct pars[1];

  /* grid and pdb for backups */
  grid_bak = make_empty_grid(grid->nvariables, 0);
  pdb_bak  = make_empty_pdb(npdbmax);

  /* max q before adjusting q */
  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qm1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qm2);
  pars->bi1 = bi1;
  pars->Xmax1 = Xmax1;
  pars->Ymax1 = 0.;
  pars->bi2 = bi2;
  pars->Xmax2 = Xmax2;
  pars->Ymax2 = 0.;
  printf("adjust_C1_C2_q_keep_qmax: in DNSdata_solve step %d: "
         "WallTime=%gs\n", it, getTimeIn_s());
  printf(" max q in inner domains before computing new q:"
         " qm1=%g qm2=%g\n", qm1, qm2);
  printf(" in box%d X=%g, box%d X=%g\n", bi1,Xmax1, bi2,Xmax2);

  /* compute error in max q */
  dqm1 = qm1 - Getd("DNSdata_qm1");
  dqm2 = qm2 - Getd("DNSdata_qm2");
  qm_error = dqm1*dqm1 + dqm2*dqm2;
  qm_error = sqrt(qm_error)/(Getd("DNSdata_qm1")+Getd("DNSdata_qm2"));
  printf(" qm error = %.4e "
         "(before adjusting C1/2)\n", qm_error);

  /* set desired qm for this iteration */
  pars->qm1 = Getd("DNSdata_qm1"); // + dqm1*0.9;
  pars->qm2 = Getd("DNSdata_qm2"); // + dqm2*0.9;
  printf(" adjusting q,C1,C2 to achieve: qm1=%g  qm2=%g  tol*0.01=%g\n",
         pars->qm1, pars->qm2, tol*0.01);

  /* print C1/2 we used before */
  printf(" old: DNSdata_C1=%g DNSdata_C2=%g\n",
         Getd("DNSdata_C1"), Getd("DNSdata_C2"));
  /* printf("     => qm1=%g qm2=%g\n", qm1, qm2); */

  /* adjust C1/2 */
  {
    /***********************************************************************/
    /* do newton_linesrch_itsP iterations of Cvec until qmerrorvec is zero */
    /***********************************************************************/
    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* adjust C1 and thus qm1 */
    Cvec[1] = Getd("DNSdata_C1");
    stat = newton_linesrch_itsP(Cvec, 1, &check, qm1_error_VectorFuncP,
                                (void *) pars, 1000, tol*0.01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
    Setd("DNSdata_C1", Cvec[1]);

    /* setup new domain with new q */
    compute_new_q_and_adjust_domainshapes_InterpFromGrid0(pars->grid,
                                                          pars->grid0, STAR1);

    /* backup grid,pdb */
    backup_grid_pdb(grid,pdb, grid_bak,pdb_bak);
    pars->grid = grid;
    if(Getv("DNSdata_m0_error_VectorFuncP_grid0","grid_bak")) pars->grid0=grid_bak;
    else  pars->grid0 = grid;

    /* adjust C2 and thus qm2 */
    Cvec[1] = Getd("DNSdata_C2");
    if(Getd("DNSdata_qm2")>0)
      stat = newton_linesrch_itsP(Cvec, 1, &check, qm2_error_VectorFuncP,
                                  (void *) pars, 1000, tol*0.01);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);  
    Setd("DNSdata_C2", Cvec[1]);

    /* setup new domain with new q */
    if(Getd("DNSdata_qm2")>0)
      compute_new_q_and_adjust_domainshapes_InterpFromGrid0(pars->grid,
                                                            pars->grid0, STAR2);

    printf("adjust_C1_C2_q_keep_qmax:\n");
    printf(" new: DNSdata_C1=%g DNSdata_C2=%g\n",
           Getd("DNSdata_C1"), Getd("DNSdata_C2"));
  }
  
  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  /* print new qm */
  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qm1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qm2);
  printf("     => qm1=%.19g qm2=%.19g\n", qm1, qm2);

  /* free grid and pdb for backups */
  free_pdb(pdb_bak, npdb);
  free_grid(grid_bak);
        
  return 0;
}

int adjust_C1_C2_q_keep_m0_or_qmax(tGrid *grid, int it, double tol)
{
  if(Getd("DNSdata_qm1")<0.0)
    return adjust_C1_C2_q_keep_restmasses(grid, it, tol);
  else
    return adjust_C1_C2_q_keep_qmax(grid, it, tol);  
}










/* for newton_linesrch_its: compute error in xmaxs */
/* if n=1 only DNSdata_Omega is adjusted
   if n=2 both DNSdata_Omega & DNSdata_x_CM are adjusted */
void xmaxs_error_VectorFunc(int n, double *vec, double *fvec)
{
  int bi1, bi2;
  double Xmax1,qmax1, Xmax2,qmax2;
  double xmax1, xmax2;
  tGrid *grid = xmaxs_error_VectorFunc__grid;

//  /* save old q in DNSdata_qgold */
//  varcopy(grid, Ind("DNSdata_qgold"), Ind("DNSdata_qg"));

  /* set DNSdata_Omega & DNSdata_x_CM */
  Setd("DNSdata_Omega", vec[1]);
  if(n>=2) Setd("DNSdata_x_CM",  vec[2]);
  printf("xmaxs_error_VectorFunc: Omega=%g x_CM=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));

  /* compute new q */
  DNS_compute_new_centered_q(grid, STAR1);
  DNS_compute_new_centered_q(grid, STAR2);

  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qmax1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qmax2);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,0.,0.);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,0.,0.);
  else
    xmax2 = Xmax2;
//  /* compute qmax1 and qmax2 */
//  qmax1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,0.,0.);
//  qmax2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,0.,0.);

  printf("xmaxs_error_VectorFunc: Omega=%g x_CM=%g xmax1=%g xmax2=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"),
         xmax1, xmax2);  fflush(stdout);
  prdivider(0);

  /* return errors */
  fvec[1] = xmax1 - xmaxs_error_VectorFunc__xmax1;
  if(n>=2) fvec[2] = xmax2 - xmaxs_error_VectorFunc__xmax2;
}

/* Adjust Omega and x_CM so that point with max q, xmax1/2 stay put */
int adjust_Omega_xCM_keep_xmax(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  double OmxCMvec[3];
  double dxmax_m0[3];
  double dxmax_00[3];
  double dxmax_p0[3];
  double dxmax_0m[3];
  double dxmax_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  int do_lnsrch = Getv("DNSdata_adjust", "always");
  int keepone   = Getv("DNSdata_adjust", "keep_one_xmax");

  /* save old Omega, x_CM */
  Omega = Getd("DNSdata_Omega");
  x_CM  = Getd("DNSdata_x_CM");
  dOmega= Omega * Getd("DNSdata_dOmega_fac");
  dx_CM = Getd("DNSdata_b") * Getd("DNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_xmax: in DNSdata_solve step %d\n"
         "adjust_Omega_xCM_keep_xmax: old Omega = %g  x_CM = %g  tol = %g\n",
         it, Omega, x_CM, tol);

  /* set global vars */
  xmaxs_error_VectorFunc__grid  = grid;
  xmaxs_error_VectorFunc__xmax1 = Getd("DNSdata_xmax1");
  xmaxs_error_VectorFunc__xmax2 = Getd("DNSdata_xmax2");
  printf("adjust_Omega_xCM_keep_xmax: old xmax1 = %g  xmax2 = %g\n",
         xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_m0);
    OmxCMvec[1] = Omega + dOmega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_p0);
    OmxCMvec[1] = Omega;
    /* xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00); */
    printf("adjust_Omega_xCM_keep_xmax: dxmax_m0[1]=%g dxmax_m0[2]=%g\n",
           dxmax_m0[1], dxmax_m0[2]);
    /* printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1]=%g dxmax_00[2]=%g\n",
           dxmax_00[1], dxmax_00[2]); */
    printf("adjust_Omega_xCM_keep_xmax: dxmax_p0[1]=%g dxmax_p0[2]=%g\n",
           dxmax_p0[1], dxmax_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dxmax_m0[1]*dxmax_p0[1]<0.0) && (dxmax_m0[2]*dxmax_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dxmax_m0[1]*dxmax_p0[1]<0.0) || (dxmax_m0[2]*dxmax_p0[2]<0.0)) )
  {
    if(keepone) /* fix one xmax to current value and then set do_lnsrch=1 */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM;
      if(dxmax_m0[1]*dxmax_p0[1]>=0.0)
      {
        xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
        xmaxs_error_VectorFunc__xmax1 += dxmax_00[1];
      }
      if(dxmax_m0[2]*dxmax_p0[2]>=0.0)
      {
        xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
        xmaxs_error_VectorFunc__xmax2 += dxmax_00[2];
      }
      do_lnsrch = 1;
      printf("adjust_Omega_xCM_keep_xmax: "
             "changed: old xmax1 = %g  xmax2 = %g\n",
             xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
    }
    else /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_0m);
      OmxCMvec[2] = x_CM + dx_CM;
      xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_0p);
      OmxCMvec[2] = x_CM;
      /* xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00); */
      printf("adjust_Omega_xCM_keep_xmax: dxmax_0m[1]=%g dxmax_0m[2]=%g\n",
             dxmax_m0[1], dxmax_m0[2]);
      /* printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1]=%g dxmax_00[2]=%g\n",
             dxmax_00[1], dxmax_00[2]); */
      printf("adjust_Omega_xCM_keep_xmax: dxmax_0p[1]=%g dxmax_0p[2]=%g\n",
             dxmax_p0[1], dxmax_p0[2]);
  
      if( (dxmax_m0[1]*dxmax_p0[1]<0.0) && (dxmax_0m[2]*dxmax_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dxmax_m0[2]*dxmax_p0[2]<0.0) && (dxmax_0m[1]*dxmax_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_xmax: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xmax1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_its(OmxCMvec, 2, &check, xmaxs_error_VectorFunc,
                               1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
    printf("adjust_Omega_xCM_keep_xmax: dxmax_00[1] = %g  dxmax_00[2] = %g\n",
           dxmax_00[1], dxmax_00[2]);
  }
  printf("adjust_Omega_xCM_keep_xmax: new Omega = %g  x_CM = %g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));
  prdivider(0);

  /* check where max are and reset DNSdata_xmax1/2 if they are too far off */
  if(Getv("DNSdata_adjust", "reset_xmax_if_problem"))
  {
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    xmaxs_error_VectorFunc__xmax1 = Getd("DNSdata_xmax1");
    xmaxs_error_VectorFunc__xmax2 = Getd("DNSdata_xmax2");
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax_00);
    if(fabs(dxmax_00[1])>5.0*tol) xmaxs_error_VectorFunc__xmax1 += dxmax_00[1];
    if(fabs(dxmax_00[2])>5.0*tol) xmaxs_error_VectorFunc__xmax2 += dxmax_00[2];
    Setd("DNSdata_xmax1", xmaxs_error_VectorFunc__xmax1);
    Setd("DNSdata_xmax2", xmaxs_error_VectorFunc__xmax2);
    /* print current maxs */
    printf("adjust_Omega_xCM_keep_xmax: resetting xmax1 = %g  xmax2 = %g\n",
           xmaxs_error_VectorFunc__xmax1, xmaxs_error_VectorFunc__xmax2);
    prdivider(0);
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_its: compute derivs of q in two domains */
/* if n=1 only DNSdata_Omega is adjusted
   if n=2 both DNSdata_Omega & DNSdata_x_CM are adjusted */
void dqdx_at_Xmax1_2_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = dqdx_at_Xmax1_2_VectorFunc__grid;
  int bi1 = dqdx_at_Xmax1_2_VectorFunc__bi1;
  int bi2 = dqdx_at_Xmax1_2_VectorFunc__bi2;
  double Xmax1 = dqdx_at_Xmax1_2_VectorFunc__Xmax1;
  double Ymax1 = dqdx_at_Xmax1_2_VectorFunc__Ymax1;
  double Xmax2 = dqdx_at_Xmax1_2_VectorFunc__Xmax2;
  double Ymax2 = dqdx_at_Xmax1_2_VectorFunc__Ymax2;
  int b;

//  /* save old q in DNSdata_qgold */
//  varcopy(grid, Ind("DNSdata_qgold"), Ind("DNSdata_qg"));

  /* set DNSdata_Omega & DNSdata_x_CM */
  Setd("DNSdata_Omega", vec[1]);
  if(n>=2) Setd("DNSdata_x_CM",  vec[2]);

  /* compute new q */
  DNS_compute_new_centered_q(grid, STAR1);
  DNS_compute_new_centered_q(grid, STAR2);

  /* get deriv dq of q in box bi1 and bi2 in DNSdata_temp1
     and dq's coeffs c in DNSdata_temp2 */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *q = box->v[Ind("DNSdata_q")];
    double *dq= box->v[Ind("DNSdata_temp1")];
    double *c = box->v[Ind("DNSdata_temp2")];
    spec_Deriv1(box, 1, q, dq);
    spec_Coeffs(box, dq, c);
  }

  /* find dq/dx at the locations Xmax1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("DNSdata_temp2")], Xmax1,Ymax1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("DNSdata_temp2")], Xmax2,Ymax2,0);

  printf("dqdx_at_Xmax1_2_VectorFunc: Omega=%.13g x_CM=%.13g\n"
         "  => dq/dx(xmax1)=%.13g",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"), fvec[1]);

  if(n>=2) printf("  dq/dx(xmax2)=%.13g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that point with max q, xmax1/2 stay put */
int adjust_Omega_xCM_keep_dqdxmax_eq_0(tGrid *grid, int it, double tol)
{
  int check, stat, bi, i;
  int blist[6];
  double OmxCMvec[3];
  double dqdx_m0[3];
  double dqdx_00[3];
  double dqdx_p0[3];
  double dqdx_0m[3];
  double dqdx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double xmax1, xmax2;
  int do_lnsrch = Getv("DNSdata_adjust", "always");
  int pr_curxmax =Getv("DNSdata_adjust", "print_current_xmax");

  /* save old Omega, x_CM */
  Omega = Getd("DNSdata_Omega");
  x_CM  = Getd("DNSdata_x_CM");
  dOmega= Omega * Getd("DNSdata_dOmega_fac");
  dx_CM = Getd("DNSdata_b") * Getd("DNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: in DNSdata_solve step %d\n"
         "adjust_Omega_xCM_keep_dqdxmax_eq_0: old Omega=%g x_CM=%g tol=%g\n",
         it, Omega, x_CM, tol);

  /* save DNSdata_xmax1/2 */
  xmax1 = Getd("DNSdata_xmax1");
  xmax2 = Getd("DNSdata_xmax2");

  /* set global vars */
  dqdx_at_Xmax1_2_VectorFunc__grid  = grid;
  /* for now we assume that the max are in box0/3 at Y=B=0
     or in box4/5 at Y=0 */
  dqdx_at_Xmax1_2_VectorFunc__Ymax1 = 0.0;
errorexit("need other boxes, not 0 and 3");
  blist[0]=0;  blist[1]=5;
  bi=b_X_of_x_forgiven_YZ_inboxlist(grid, blist, 2,
                                    &dqdx_at_Xmax1_2_VectorFunc__Xmax1,
                                    xmax1, 0.0,0.0);
  dqdx_at_Xmax1_2_VectorFunc__bi1 = bi;
  dqdx_at_Xmax1_2_VectorFunc__Ymax2 = 0.0;
  blist[0]=3;  blist[1]=4;
  bi=b_X_of_x_forgiven_YZ_inboxlist(grid, blist, 2,
                                    &dqdx_at_Xmax1_2_VectorFunc__Xmax2,
                                    xmax2, 0.0,0.0);
  dqdx_at_Xmax1_2_VectorFunc__bi2 = bi;
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: xmax1=%g xmax2=%g\n",
         xmax1, xmax2);
  if(dqdx_at_Xmax1_2_VectorFunc__bi1<0 || dqdx_at_Xmax1_2_VectorFunc__bi2<0)
  {
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: failed to find Xmax1 or "
           "Xmax2\n"
           " dqdx_at_Xmax1_2_VectorFunc__bi1=%d "
           " dqdx_at_Xmax1_2_VectorFunc__bi2=%d\n",
           dqdx_at_Xmax1_2_VectorFunc__bi1, dqdx_at_Xmax1_2_VectorFunc__bi2);
    return -1;
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: Xmax1=%g Xmax2=%g\n",
         dqdx_at_Xmax1_2_VectorFunc__Xmax1, dqdx_at_Xmax1_2_VectorFunc__Xmax2);
  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_m0);
    OmxCMvec[1] = Omega + dOmega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_p0);
    OmxCMvec[1] = Omega;
    /* dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00); */
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_m0[1]=%g dqdx_m0[2]=%g\n",
           dqdx_m0[1], dqdx_m0[2]);
    /* printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]); */
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_p0[1]=%g dqdx_p0[2]=%g\n",
           dqdx_p0[1], dqdx_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_m0[2]*dqdx_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dqdx_m0[1]*dqdx_p0[1]<0.0) || (dqdx_m0[2]*dqdx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_0m);
      OmxCMvec[2] = x_CM + dx_CM;
      dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_0p);
      OmxCMvec[2] = x_CM;
      /* dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00); */
      printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_0m[1]=%g dqdx_0m[2]=%g\n",
             dqdx_m0[1], dqdx_m0[2]);
      /* printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
             dqdx_00[1], dqdx_00[2]); */
      printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_0p[1]=%g dqdx_0p[2]=%g\n",
             dqdx_p0[1], dqdx_p0[2]);
  
      if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_0m[2]*dqdx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dqdx_m0[2]*dqdx_p0[2]<0.0) && (dqdx_0m[1]*dqdx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_its iterations until xmax1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_its(OmxCMvec, 2, &check, dqdx_at_Xmax1_2_VectorFunc,
                               1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dqdx_at_Xmax1_2_VectorFunc(2, OmxCMvec, dqdx_00);
    printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]);
  }
  printf("adjust_Omega_xCM_keep_dqdxmax_eq_0: new Omega=%g x_CM=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));
  prdivider(0);

  /* find and print current xmax1/2 */
  if(pr_curxmax)
  {
    double dxmax[3];
    xmaxs_error_VectorFunc__grid  = grid;
    xmaxs_error_VectorFunc__xmax1 = Getd("DNSdata_xmax1");
    xmaxs_error_VectorFunc__xmax2 = Getd("DNSdata_xmax2");
    xmaxs_error_VectorFunc(2, OmxCMvec, dxmax);
  }

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  return 0;
}


/* for newton_linesrch_itsP: compute derivs of Func in two domains */
/* if n=1 only DNSdata_Omega is adjusted
   if n=2 both DNSdata_Omega & DNSdata_x_CM are adjusted */
void dFuncdx_at_Xfm1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int b, bi1, bi2;
  double Xfm1,Yfm1, Xfm2,Yfm2;
  t_grid_bXYZ1_bXYZ2_struct *pars;

  /* get pars */
  pars = (t_grid_bXYZ1_bXYZ2_struct *) p;
  grid = pars->grid;
  bi1 = pars->b1;
  bi2 = pars->b2;
  Xfm1 = pars->X1;
  Yfm1 = pars->Y1;
  Xfm2 = pars->X2;
  Yfm2 = pars->Y2;

  /* set DNSdata_Omega & DNSdata_x_CM */
  Setd("DNSdata_Omega", vec[1]);
  if(n>=2) Setd("DNSdata_x_CM",  vec[2]);

  /* get deriv dFunc of Func in box bi1 and bi2 in DNSdata_temp1
     and dFunc's coeffs c in DNSdata_temp2 */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *Func = box->v[Ind("DNSdata_temp4")];
    double *dFunc= box->v[Ind("DNSdata_temp1")];
    double *c = box->v[Ind("DNSdata_temp2")];
    spec_Deriv1(box, 1, Func, dFunc);
    spec_Coeffs(box, dFunc, c);
  }

  /* find dFunc/dx at the locations Xfm1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("DNSdata_temp2")], Xfm1,Yfm1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("DNSdata_temp2")], Xfm2,Yfm2,0);

  printf("dFuncdx_at_Xfm1_2_VectorFuncP: Omega=%.13g x_CM=%.13g\n"
         "  => dFunc/dX(Xfm1)=%.13g",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"), fvec[1]);

  if(n>=2) printf("  dFunc/dX(Xfm2)=%.13g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that point with max Func, xfm1/2 stay put */
int adjust_Omega_xCM_keep_dFuncdxfm_eq_0(tGrid *grid, int it, double tol)
{
  int check, stat, bi, bi1,bi2, i;
  double OmxCMvec[3];
  double dqdx_m0[3];
  double dqdx_00[3];
  double dqdx_p0[3];
  double dqdx_0m[3];
  double dqdx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double Xfm1,qfm1, Xfm2,qfm2;
  int do_lnsrch = Getv("DNSdata_adjust", "always");
  t_grid_bXYZ1_bXYZ2_struct pars[1];

  /* save old Omega, x_CM */
  Omega = Getd("DNSdata_Omega");
  x_CM  = Getd("DNSdata_x_CM");
  dOmega= Omega * Getd("DNSdata_dOmega_fac");
  dx_CM = Getd("DNSdata_b") * Getd("DNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: in DNSdata_solve step %d\n"
         "  old Omega=%g x_CM=%g tol=%g\n", it, Omega, x_CM, tol);

  /* set Xfm1 and Xfm2, i.e. find max in Func = DNSdata_qcorot */
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_qcorot"), STAR1,
                                   &bi1, &Xfm1, &qfm1);
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_qcorot"), STAR2,
                                   &bi2, &Xfm2, &qfm2);
  /* set global vars */
  pars->grid  = grid;
  pars->b1 = bi1;
  pars->X1 = Xfm1;
  pars->Y1 = 0.0;
  pars->Z1 = 0.0;
  pars->b2 = bi2;
  pars->X2 = Xfm2;
  pars->Y2 = 0.0;
  pars->Z2 = 0.0;

  if(pars->b1<0 || pars->b2<0)
  {
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: failed to find Xfm1 or "
           "Xfm2\n"
           " pars->b1=%d "
           " pars->b2=%d\n",
           pars->b1, pars->b2);
    return -1;
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: Xfm1=%g Xfm2=%g\n",
         pars->X1, pars->X2);
  prdivider(0);
//  /* adjust C1/2 so that masses are correct */
//  adjust_C1_C2_q_keep_m0_or_qmax(grid, it, tol*100.0);
//  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars); */
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_m0[1]=%g dqdx_m0[2]=%g\n",
           dqdx_m0[1], dqdx_m0[2]);
    /* printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]); */
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_p0[1]=%g dqdx_p0[2]=%g\n",
           dqdx_p0[1], dqdx_p0[2]);
    prdivider(0);
  }
  /* check if there is a zero, if so set do_lnsrch=1 */
  if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_m0[2]*dqdx_p0[2]<0.0) )
    do_lnsrch = 1;
  if( (do_lnsrch==0) &&
      ((dqdx_m0[1]*dqdx_p0[1]<0.0) || (dqdx_m0[2]*dqdx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars); */
      printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_0m[1]=%g dqdx_0m[2]=%g\n",
             dqdx_m0[1], dqdx_m0[2]);
      /* printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
             dqdx_00[1], dqdx_00[2]); */
      printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_0p[1]=%g dqdx_0p[2]=%g\n",
             dqdx_p0[1], dqdx_p0[2]);
  
      if( (dqdx_m0[1]*dqdx_p0[1]<0.0) && (dqdx_0m[2]*dqdx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dqdx_m0[2]*dqdx_p0[2]<0.0) && (dqdx_0m[1]*dqdx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xfm1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, dFuncdx_at_Xfm1_2_VectorFuncP,
                               (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dFuncdx_at_Xfm1_2_VectorFuncP(2, OmxCMvec, dqdx_00, (void *) pars);
    printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: dqdx_00[1]=%g dqdx_00[2]=%g\n",
           dqdx_00[1], dqdx_00[2]);
  }
  printf("adjust_Omega_xCM_keep_dFuncdxfm_eq_0: new Omega=%g x_CM=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  return 0;
}

/* for newton_linesrch_itsP: compute derivs of Integrated Euler Eqn 
   in domain 0 and 3, or 4 and 5 */
/* if n=1 only DNSdata_Omega is adjusted
   if n=2 both DNSdata_Omega & DNSdata_x_CM are adjusted */
void dIntegEulerdx_at_X1_2_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  int b, bi1, bi2;
  double X1,Y1, X2,Y2;
  double Om, xcm;
  t_grid_bXYZ1_bXYZ2_struct *pars;

  /* get pars */
  pars = (t_grid_bXYZ1_bXYZ2_struct *) p;
  grid = pars->grid;
  bi1 = pars->b1;
  bi2 = pars->b2;
  X1 = pars->X1;
  Y1 = pars->Y1;
  X2 = pars->X2;
  Y2 = pars->Y2;

  /* set Omega & x_CM */
  Om  = vec[1];
  xcm = vec[2];

  /* compute lnIntegEuler in DNSdata_temp4 and 
     dlnIntegEulerdx/y/z in DNSdata_temp1/2/3 */
  DNS_set_dlnIntegEuler(grid, Ind("DNSdata_temp4"), 
                        Ind("DNSdata_temp1"), Om, xcm);

  /* put coeffs of DNSdata_temp1=dlnIntegEulerdx into DNSdata_temp3=c */
  forallboxes(grid, b) if(b==bi1 || b==bi2)
  {
    tBox *box = grid->box[b];
    double *dlnIE = box->v[Ind("DNSdata_temp1")];
    double *c     = box->v[Ind("DNSdata_temp3")];
    spec_Coeffs(box, dlnIE, c);
  }

  /* find dFunc/dx at the locations X1/2 in NS1/2 and store
     them in fvec[1] and fvec[2] */
  fvec[1] = spec_interpolate(grid->box[bi1],
                       grid->box[bi1]->v[Ind("DNSdata_temp3")], X1,Y1,0);
  if(n>=2) fvec[2] = spec_interpolate(grid->box[bi2], 
                       grid->box[bi2]->v[Ind("DNSdata_temp3")], X2,Y2,0);

  printf("dIntegEulerdx_at_X1_2_VectorFuncP: Om=%.12g xcm=%.12g\n"
         " => dIntEuler/dX(X1)=%.12g", Om, xcm, fvec[1]);

  if(n>=2) printf(" dIntEuler/dX(X2)=%.12g\n", fvec[2]);
  else	   printf("\n");
  fflush(stdout);
  prdivider(0);
}

/* Adjust Omega and x_CM so that force balance is maintained at xmax */
int adjust_Omega_xCM_forcebalance(tGrid *grid, int it, double tol)
{
  int check, stat, bi, bi1,bi2, i;
  double OmxCMvec[3];
  double dIEdxx_m0[3];
  double dIEdxx_00[3];
  double dIEdxx_p0[3];
  double dIEdxx_0m[3];
  double dIEdxx_0p[3];
  double Omega, x_CM;
  double dOmega, dx_CM;
  double Xqm1,qqm1, Xqm2,qqm2;
  int do_lnsrch = Getv("DNSdata_adjust", "always");
  t_grid_bXYZ1_bXYZ2_struct pars[1];

  /* save old Omega, x_CM */
  Omega = Getd("DNSdata_Omega");
  x_CM  = Getd("DNSdata_x_CM");
  dOmega= Omega * Getd("DNSdata_dOmega_fac");
  dx_CM = Getd("DNSdata_b") * Getd("DNSdata_dx_CM_fac");
  prdivider(0);
  printf("adjust_Omega_xCM_forcebalance: in DNSdata_solve step %d\n"
         "  old Omega=%g x_CM=%g tol=%g  WallTime=%gs\n",
         it, Omega, x_CM, tol, getTimeIn_s());

  /* set Xqm1 and Xqm2, i.e. find max in DNSdata_q */
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_q"), STAR1,
                                   &bi1, &Xqm1, &qqm1);
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_q"), STAR2,
                                   &bi2, &Xqm2, &qqm2);
  /* set global vars */
  pars->grid  = grid;
  pars->b1 = bi1;
  pars->X1 = Xqm1;
  pars->Y1 = 0.0;
  pars->Z1 = 0.0;
  pars->b2 = bi2;
  pars->X2 = Xqm2;
  pars->Y2 = 0.0;
  pars->Z2 = 0.0;

  if(pars->b1<0 || pars->b2<0)
  {
    printf("adjust_Omega_xCM_forcebalance: failed to find Xqm1 or "
           "Xqm2\n"
           " pars->b1=%d "
           " pars->b2=%d\n",
           pars->b1, pars->b2);
    return -1;
  }
  printf("adjust_Omega_xCM_forcebalance:\n");
  printf("  box%d (X1,Y1)=(%g,%g) , box%d (X2,Y2)=(%g,%g)\n",
         pars->b1,pars->X1,pars->Y1, pars->b2,pars->X2,pars->Y2);
  prdivider(0);
//  /* adjust C1/2 so that masses are correct */
//  adjust_C1_C2_q_keep_m0_or_qmax(grid, it, tol*100.0);
//  prdivider(0);

  if(do_lnsrch==0)
  {
    /* find deviations for Omega +/- dOmega, x_CM */
    OmxCMvec[2] = x_CM;
    OmxCMvec[1] = Omega - dOmega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_m0, (void *) pars);
    OmxCMvec[1] = Omega + dOmega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_p0, (void *) pars);
    OmxCMvec[1] = Omega;
    /* dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars); */
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_m0[1]=%g dIEdxx_m0[2]=%g\n",
           dIEdxx_m0[1], dIEdxx_m0[2]);
    /* printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
           dIEdxx_00[1], dIEdxx_00[2]); */
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_p0[1]=%g dIEdxx_p0[2]=%g\n",
           dIEdxx_p0[1], dIEdxx_p0[2]);
    prdivider(0);
    /* check if there is a zero, if so set do_lnsrch=1 */
    if( (dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) && (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0) )
      do_lnsrch = 1;
  }
  if( (do_lnsrch==0) &&
      ((dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) || (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0)) )
  {
    /* find deviations for Omega, x_CM +/- dx_CM */
    {
      OmxCMvec[1] = Omega;
      OmxCMvec[2] = x_CM - dx_CM;
      dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_0m, (void *) pars);
      OmxCMvec[2] = x_CM + dx_CM;
      dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_0p, (void *) pars);
      OmxCMvec[2] = x_CM;
      /* dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars); */
      printf("adjust_Omega_xCM_forcebalance: dIEdxx_0m[1]=%g dIEdxx_0m[2]=%g\n",
             dIEdxx_m0[1], dIEdxx_m0[2]);
      /* printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
             dIEdxx_00[1], dIEdxx_00[2]); */
      printf("adjust_Omega_xCM_forcebalance: dIEdxx_0p[1]=%g dIEdxx_0p[2]=%g\n",
             dIEdxx_p0[1], dIEdxx_p0[2]);
  
      if( (dIEdxx_m0[1]*dIEdxx_p0[1]<0.0) && (dIEdxx_0m[2]*dIEdxx_0p[2]<0.0) )
        do_lnsrch = 1;
      if( (dIEdxx_m0[2]*dIEdxx_p0[2]<0.0) && (dIEdxx_0m[1]*dIEdxx_0p[1]<0.0) )
        do_lnsrch = 1;
    }
  }
  printf("adjust_Omega_xCM_forcebalance: do_lnsrch = %d\n", do_lnsrch);
  prdivider(0);
  if(do_lnsrch)
  {
    /**************************************************************************/
    /* do newton_linesrch_itsP iterations until xfm1/2 are where we want them */
    /**************************************************************************/
    OmxCMvec[1] = Omega;
    OmxCMvec[2] = x_CM;
    stat = newton_linesrch_itsP(OmxCMvec, 2, &check, dIntegEulerdx_at_X1_2_VectorFuncP,
                               (void *) pars, 1000, tol*0.5);
    if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  }
  else
  {
    OmxCMvec[1] = Omega;
    dIntegEulerdx_at_X1_2_VectorFuncP(2, OmxCMvec, dIEdxx_00, (void *) pars);
    printf("adjust_Omega_xCM_forcebalance: dIEdxx_00[1]=%g dIEdxx_00[2]=%g\n",
           dIEdxx_00[1], dIEdxx_00[2]);
  }

  /* Nobody has touched the pars DNSdata_Omega and DNSdata_x_CM so far.
     Now set them to what we found */
  Setd("DNSdata_Omega", OmxCMvec[1]);
  Setd("DNSdata_x_CM",  OmxCMvec[2]);
  printf("adjust_Omega_xCM_forcebalance: new Omega=%g x_CM=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  return 0;
}


/* for newton_linesrch_itsP: compute Py_ADM */
void Py_ADM_of_xCM_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  tGrid *grid;
  double xcm;
  double Omega = Getd("DNSdata_Omega");
  t_grid_bXYZ1_bXYZ2_struct *pars;
  int itemp1 = Ind("DNSdata_temp1");
  int itemp2 = Ind("DNSdata_temp2");
  int itemp3 = Ind("DNSdata_temp3");
  double Py_ADM1, Py_ADM2, Py_ADM;

  /* get pars */
  pars = (t_grid_bXYZ1_bXYZ2_struct *) p;
  grid = pars->grid;

  /* set xcm */
  xcm = vec[1];

  /* compute ADM mom. Py_ADM */
  DNS_set_P_ADM_VolInt_integrand_Om_xcm(grid, itemp1,itemp2,itemp3, Omega,xcm);
errorexit("need other boxes, not 0 and 3");
  Py_ADM1 = InnerVolumeIntegral(grid, 0, itemp2);
  Py_ADM2 = InnerVolumeIntegral(grid, 3, itemp2);
  Py_ADM = Py_ADM1 + Py_ADM2;

  printf("Py_ADM_of_xCM_VectorFuncP: xcm=%.12g  Py_ADM=%.12g\n", xcm, Py_ADM);
  fflush(stdout);
  fvec[1] = Py_ADM;
}

/* Adjust Omega and x_CM. First find x_CM s.t. Py_ADM = 0. Then obtain Omega
   from force balance */
int adjust_xCM_Omega_Py0_forcebalance(tGrid *grid, int it, double tol)
{
  int check, stat, bi, bi1,bi2, i;
  double xcmvec[2];
  double OmxCMvec[3];
  double Omega, x_CM, Om1,Om2;
  double Xqm1,qm1, Xqm2,qm2;
  t_grid_bXYZ1_bXYZ2_struct pars[1];

  /* save old Omega, x_CM */
  Omega = Getd("DNSdata_Omega");
  x_CM  = Getd("DNSdata_x_CM");
  prdivider(0);
  printf("adjust_xCM_Omega_Py0_forcebalance: in DNSdata_solve step %d\n"
         "  old Omega=%g x_CM=%g tol=%g  WallTime=%gs\n",
         it, Omega, x_CM, tol, getTimeIn_s());
  prdivider(0);

  /* find x_CM */
  pars->grid  = grid; /* set grid in pars */
  /****************************************************************/
  /* do newton_linesrch_itsP iterations until xcm is what we want */
  /****************************************************************/
  xcmvec[1] = x_CM;
  stat = newton_linesrch_itsP(xcmvec, 1, &check, Py_ADM_of_xCM_VectorFuncP,
                             (void *) pars, 1000, tol*0.5);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);

  /* Nobody has touched the par DNSdata_x_CM so far.
     Now set it to what we found */
  Setd("DNSdata_x_CM",  xcmvec[1]);
  printf("adjust_xCM_Omega_Py0_forcebalance: setting x_CM=%g\n",
         Getd("DNSdata_x_CM"));
  x_CM  = Getd("DNSdata_x_CM");
  prdivider(0);

  /* set Xqm1 and Xqm2, i.e. find max in DNSdata_q */
  bi1=0;  bi2=3;
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_q"), STAR1,
                                   &bi1, &Xqm1, &qm1);
  find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_q"), STAR2,
                                   &bi2, &Xqm2, &qm2);

  if(bi1<0 || bi2<0)
  {
    printf("adjust_xCM_Omega_Py0_forcebalance: failed to find Xqm1 or "
           "Xqm2\n"
           " bi1=%d "
           " bi2=%d\n", bi1, bi2);
    return -1;
  }

  /* set global vars such that we look for Omega such that force balance
     holds at center of star1 at Xqm1,Yqm1 */
  pars->grid  = grid;
  pars->b1 = bi1;
  pars->X1 = Xqm1;
  pars->Y1 = 0.0;
  pars->Z1 = 0.0;
  pars->b2 = bi2;
  pars->X2 = Xqm2;
  pars->Y2 = 0.0;
  pars->Z2 = 0.0;
  printf("adjust_xCM_Omega_Py0_forcebalance:  box%d (Xm,Ym)=(%g,%g)\n",
         pars->b1,pars->X1,pars->Y1);
  /**************************************/
  /* do newton_linesrch_itsP iterations */
  /**************************************/
  OmxCMvec[1] = Omega;
  OmxCMvec[2] = x_CM;
  stat = newton_linesrch_itsP(OmxCMvec, 1, &check, dIntegEulerdx_at_X1_2_VectorFuncP,
                              (void *) pars, 1000, tol*0.5);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  /* save result in Om1 */
  Om1 = OmxCMvec[1];

  /* set global vars such that we look for Omega such that force balance
     holds at center of star2 at Xqm2,Yqm2 */
  pars->grid  = grid;
  pars->b1 = bi2;
  pars->X1 = Xqm2;
  pars->Y1 = 0.0;
  pars->Z1 = 0.0;
  pars->b2 = bi1;
  pars->X2 = Xqm1;
  pars->Y2 = 0.0;
  pars->Z2 = 0.0;
  printf("adjust_xCM_Omega_Py0_forcebalance:  box%d (Xm,Ym)=(%g,%g)\n",
         pars->b1,pars->X1,pars->Y1);
  /**************************************/
  /* do newton_linesrch_itsP iterations */
  /**************************************/
  OmxCMvec[1] = Omega;
  OmxCMvec[2] = x_CM;
  stat = newton_linesrch_itsP(OmxCMvec, 1, &check, dIntegEulerdx_at_X1_2_VectorFuncP,
                              (void *) pars, 1000, tol*0.5);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  /* save result in Om2 */
  Om2 = OmxCMvec[1];

  /* Nobody has touched the par DNSdata_Omega.
     Now set it to what we found */
  Setd("DNSdata_Omega", 0.5*(Om1+Om2)); /* use average of Om1 and Om2 */
  printf("adjust_xCM_Omega_Py0_forcebalance: new Omega=%g x_CM=%g\n",
         Getd("DNSdata_Omega"), Getd("DNSdata_x_CM"));
  prdivider(0);

  /* set q to zero if q<0, and also in region 1 & 2 */
  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  return 0;
}



/* adjust m01 and m02 to let them e.g. grow during iterations */
void adjust_DNSdata_m01_m02(void)
{
  double cm01 = Getd("DNSdata_m01");
  double cm02 = Getd("DNSdata_m02");
  double tm01 = Getd("DNSdata_desired_m01");
  double tm02 = Getd("DNSdata_desired_m02");
  double mCh1 = cm01*Getd("DNSdata_m0change");
  double mCh2 = cm02*Getd("DNSdata_m0change");

  /* change DNSdata_m01 and DNSdata_m02 by mCh1/2 to get closer to
     DNSdata_desired_m01 and DNSdata_desired_m02 */
  if(fabs(tm01-cm01) > mCh1) cm01 += mCh1*signum(tm01-cm01);
  else                       cm01  = tm01;
  if(fabs(tm02-cm02) > mCh2) cm02 += mCh2*signum(tm02-cm02);
  else                       cm02  = tm02;
  Setd("DNSdata_m01", cm01);
  Setd("DNSdata_m02", cm02);
  printf("adjust_DNSdata_m01_m02: DNSdata_m0change=%g mCh1=%g mCh2=%g\n",
         Getd("DNSdata_m0change"), mCh1, mCh2);
  printf(" => DNSdata_m01=%.13g DNSdata_m02=%.13g\n", cm01, cm02);
  printf("    DNSdata_desired_m01=%.13g DNSdata_desired_m02=%.13g\n", tm01, tm02);
}

/* compute weighted average of current and old values,
   and return total error */
double average_current_and_old(double weight, tGrid *grid,
                               tVarList *vlFu, tVarList *vlu,
                               tVarList *vluDerivs, tVarList *vlJdu)
{
  double tm01 = Getd("DNSdata_m01");
  double tm02 = Getd("DNSdata_m02");
  double m01, m02;
  double normresnonlin, L2qdiff, dm01, dm02, m0err, error;
  int b;

  /* if we iterate over the rest masses the true mass goals are different */
  if(Getv("DNSdata_iterate_m0", "yes"))
  {
    tm01 = Getd("DNSdata_desired_m01");
    tm02 = Getd("DNSdata_desired_m02");
  }

  /* reset new values from ell. solve as average between old and new.
     I.e. do: new = weight*new + (1-weight)*old  */
  varadd(grid, Ind("DNSdata_Psi"),   weight,Ind("DNSdata_Psi"),   (1.0-weight),Ind("DNSdata_Psiold"));
  varadd(grid, Ind("DNSdata_alphaP"),weight,Ind("DNSdata_alphaP"),(1.0-weight),Ind("DNSdata_alphaPold"));
  varadd(grid, Ind("DNSdata_Bx"),    weight,Ind("DNSdata_Bx"),    (1.0-weight),Ind("DNSdata_Boldx"));
  varadd(grid, Ind("DNSdata_By"),    weight,Ind("DNSdata_By"),    (1.0-weight),Ind("DNSdata_Boldy"));
  varadd(grid, Ind("DNSdata_Bz"),    weight,Ind("DNSdata_Bz"),    (1.0-weight),Ind("DNSdata_Boldz"));
  varadd(grid, Ind("DNSdata_Sigma"), weight,Ind("DNSdata_Sigma"), (1.0-weight),Ind("DNSdata_Sigmaold"));

  /* compute masses */
  m01 = GetInnerRestMass(grid, STAR1);
  m02 = GetInnerRestMass(grid, STAR2);
  
  /* compute error in masses */
  dm01 = (m01 - tm01)/(tm01+tm02);
  dm02 = (m02 - tm02)/(tm01+tm02);
  m0err = fabs(dm01) + fabs(dm02);
  
  /* evalute residual */
  F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
  normresnonlin = GridL2Norm(vlFu);
  printf("average_current_and_old:  weight=%g\n", weight);
  printf(" => m01=%.19g m02=%.19g\n", m01, m02);

  /* compute error in q. Do it like this: */
  /* compute error in new uncentered q, i.e. the diff between saved 
     DNSdata_qnocent and current uncentered q */
  DNS_compute_new_q(grid, Ind("DNSdata_temp2")); /* set temp2 = new unc. q */
  /* set temp1 = temp2 - qnocent = qnocent_new - qnocent_old */
  varadd(grid, Ind("DNSdata_temp1"),
               1,Ind("DNSdata_temp2"), -1,Ind("DNSdata_qnocent"));
  L2qdiff =0.;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    if(box->MATTR == INSIDE)
      L2qdiff += varBoxL2Norm(box, Ind("DNSdata_temp1"));
  }

  /* compute total error */
  error = normresnonlin + L2qdiff + m0err;
  printf(" => residual=%.4e L2qdiff=%.4e m0err=%.3e => error=%.4e\n",
         normresnonlin, L2qdiff, m0err, error);

  return error;
}

/* get total residual but without contribution in outer boxes */
double normresnonlin_without_DNSdata_Sigma_outside(tGrid *grid)
{
  int b;
  int iSig = Ind("DNSdata_Sigma_Err");

  /* should we switch off all special Sigma BCs when computing vlFu? */
  if(Getv("DNSdata_Sigma_surface_BCs","OutputSurfaceBCres"))
  {
    char *BCsav;
    BCsav = strdup(Gets("DNSdata_Sigma_surface_BCs"));
    Sets("DNSdata_Sigma_surface_BCs", "");
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
    Sets("DNSdata_Sigma_surface_BCs", BCsav);
    free(BCsav);
  }
  else
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);

  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    int i;
    if(box->MATTR==INSIDE) continue; /* do nothing inside stars */
    forallpoints(box, i)  box->v[iSig][i] = 0.0;
  }
  return GridL2Norm(vlFu);
}


/* Solve the Equations */
int DNSdata_solve(tGrid *grid)
{
  int    itmax        = Geti("DNSdata_itmax");
  double tol          = Getd("DNSdata_tol");
  double adjusttol    = max2(tol, Getd("DNSdata_adjust_mintol"));
  double esw          = Getd("DNSdata_esw");
  double esw1         = Getd("DNSdata_esw1");
  int    allow_esw1_it= Geti("DNSdata_allow_esw1_first_at");
  double Sigma_esw    = Getd("DNSdata_Sigma_esw");
  double Sigma_esw1   = Getd("DNSdata_Sigma_esw1");
  int    allow_Sigma_esw1_it = Geti("DNSdata_allow_Sigma_esw1_first_at");
  int    Newton_itmax = itmax;
  double NewtTolFac   = Getd("DNSdata_Newton_tolFac");
  double Newton_tol   = tol*NewtTolFac;
  int    linSolver_itmax  = Geti("DNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("DNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("DNSdata_linSolver_tol");
  double normresnonlin;
  double realnormres = 1e10;
  double realnormres_old;
  double realSigmares, restres;
  int    itsSinceExtraSigma = 0;
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *));
  tVarList *vldummy;
  int it;
  double dOmega = Getd("DNSdata_Omega")*0.1;
  double totalerr1, totalerr;
  char str[1000];

  /* choose linear solver */
  if(Getv("DNSdata_linSolver", "bicgstab"))
    linear_solver=bicgstab;
  else if(Getv("DNSdata_linSolver", "bicgstab_with_fd_UMFPACK_precon"))
    linear_solver=bicgstab_with_fd_UMFPACK_precon;
  else if(Getv("DNSdata_linSolver", "bicgstab_with_BlockJacobi_precon"))
    linear_solver=bicgstab_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "LAPACK"))
    linear_solver=LAPACK_dgesv_wrapper;
  else if(Getv("DNSdata_linSolver", "SOR_Iterator"))
    linear_solver=SOR_Iterator;
  else if(Getv("DNSdata_linSolver", "templates_GMRES"))
    linear_solver=templates_gmres_wrapper;
  else if(Getv("DNSdata_linSolver", "templates_BICGSTAB"))
    linear_solver=templates_bicgstab_wrapper;
  else if(Getv("DNSdata_linSolver", "templates_CGS"))
    linear_solver=templates_cgs_wrapper;
  else if(Getv("DNSdata_linSolver", "templates_SOR"))
    linear_solver=templates_sor_wrapper;
  else if(Getv("DNSdata_linSolver", "UMFPACK"))
    linear_solver=UMFPACK_solve_wrapper;
  else if(Getv("DNSdata_linSolver", "UMFPACK_forSortedVars"))
    linear_solver=UMFPACK_solve_forSortedVars_wrapper;
  else if(Getv("DNSdata_linSolver", "WTsolver"))
    linear_solver=WTsolver;
  else if(Getv("DNSdata_linSolver", "SPQR"))
    linear_solver=SuiteSparseQR_solve_wrapper;
  else if(Getv("DNSdata_linSolver", "templates_GMRES_with_Jacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_Jacobi_precon;
  else if(Getv("DNSdata_linSolver", "templates_GMRES_with_BlockJacobi_precon"))
    linear_solver=templates_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "templates_BICGSTAB_with_BlockJacobi_precon"))
    linear_solver=templates_bicgstab_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "templates_CGS_with_BlockJacobi_precon"))
    linear_solver=templates_cgs_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "templates_QMR_with_BlockJacobi_precon"))
    linear_solver=templates_qmr_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "templates_BICG_with_BlockJacobi_precon"))
    linear_solver=templates_bicg_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "ZIB_GMRES_with_BlockJacobi_precon"))
    linear_solver=ZIB_gmres_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "ZIB_GBIT_with_BlockJacobi_precon"))
    linear_solver=ZIB_gbit_wrapper_with_BlockJacobi_precon;
  else if(Getv("DNSdata_linSolver", "ZIB_PCG_with_BlockJacobi_precon"))
    linear_solver=ZIB_pcg_wrapper_with_BlockJacobi_precon;
  else
    errorexit("DNSdata_solve: unknown DNSdata_linSolver");

  /* allocate varlists */
  vlu  = vlalloc(grid);
  vluDerivs= vlalloc(grid);

  /* add all vars to vlu */
  vlpush(vlu, Ind("DNSdata_Psi"));
  vlpush(vlu, Ind("DNSdata_Bx"));
  vlpush(vlu, Ind("DNSdata_alphaP"));
  vlpush(vlu, Ind("DNSdata_Sigma"));

  /* add derivs to vluDerivs */
  vlpush(vluDerivs, Ind("DNSdata_Psix"));
  vlpush(vluDerivs, Ind("DNSdata_Psixx"));
  vlpush(vluDerivs, Ind("DNSdata_Bxx"));
  vlpush(vluDerivs, Ind("DNSdata_Bxxx"));
  vlpush(vluDerivs, Ind("DNSdata_alphaPx"));
  vlpush(vluDerivs, Ind("DNSdata_alphaPxx"));
  vlpush(vluDerivs, Ind("DNSdata_Sigmax"));
  vlpush(vluDerivs, Ind("DNSdata_Sigmaxx"));

  /* enable vlu, vluDerivs */
  enablevarlist(vlu);
  enablevarlist(vluDerivs); 

  /* now duplicate vlu to get vlFu */  
  vlFu = AddDuplicateEnable(vlu, "_Err");

  /* now duplicate vlFu, vlu and vluDerivs for linarized Eqs. */
  vlJdu      = AddDuplicateEnable(vlFu, "_l");
  vldu       = AddDuplicateEnable(vlu,  "_l");
  vlduDerivs = AddDuplicateEnable(vluDerivs, "_l");

  /* start main iterations */
  prdivider(1);
  printf("DNSdata_solve: starting main iteration loop ...\n");

  /* choose initial Newton_tol */
  F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
  normresnonlin = GridL2Norm(vlFu);
  Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

  /* compute diagnostics like ham and mom */
  DNSdata_verify_solution(grid);

  /* output grid before any iterations are done */
  grid->time  = -itmax;
  write_grid(grid);
  DNSdata_analyze(grid);
  grid->time += 1.0;

  /* main iteration loop, do it until res is small enough */
  for(it=1; it <= itmax; it++)
  {
    int restart;
    char Loop_State[6667];

    printf("DNSdata_solve step %d:\n", it);
    prdivider(1);

    /* save state of all relevant C variables in this loop in 
       DNSdata_Main_Iteration_Loop_State, so a checkpoint saves this */
    snprintf(Loop_State, 6666,
      "%d %d %.15g", 
      it, itsSinceExtraSigma, realnormres);
    Sets("DNSdata_Main_Iteration_Loop_State", Loop_State);

    /* do checkpointing (works only if checkpoint_restart_it<1) */
    restart = checkpoint(grid);

    /* update pars from file, and write new pars */
    if(parameterio_update_pars(grid) == 1 || restart == 1)
    {
      itmax        = Geti("DNSdata_itmax");
      tol          = Getd("DNSdata_tol");
      adjusttol    = max2(tol, Getd("DNSdata_adjust_mintol"));
      esw          = Getd("DNSdata_esw");
      esw1         = Getd("DNSdata_esw1");
      allow_esw1_it= Geti("DNSdata_allow_esw1_first_at");
      Sigma_esw    = Getd("DNSdata_Sigma_esw");
      Sigma_esw1   = Getd("DNSdata_Sigma_esw1");
      allow_Sigma_esw1_it = Geti("DNSdata_allow_Sigma_esw1_first_at");
      Newton_itmax = itmax;
      NewtTolFac   = Getd("DNSdata_Newton_tolFac");
      linSolver_itmax  = Geti("DNSdata_linSolver_itmax");
      linSolver_tolFac = Getd("DNSdata_linSolver_tolFac");
      linSolver_tol    = Getd("DNSdata_linSolver_tol");
      /* reset Newton_tol */
      F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
      normresnonlin = GridL2Norm(vlFu);
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* read some C-vars from DNSdata_Main_Iteration_Loop_State par */
      sscanf(Gets("DNSdata_Main_Iteration_Loop_State"),
        "%d %d %le",
        &it, &itsSinceExtraSigma, &realnormres);
    }
    parameterio_write_current_pars(grid);

    /* If we read a checkpoint that had finished the solve, i.e.
       that has $time = 0, do not solve again! */
    if(grid->time==0.0 && restart == 1) break;

    /* what to do with q at A=0 and q<0 */
    if(Getv("DNSdata_set_negative_q", "zero"))
      set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);
    if(Getv("DNSdata_set_Surface_q", "zero"))
      set_Var_to_Val_atSurface(grid, Ind("DNSdata_q"), 0.0);

    /* center q first, then solve and adjust C1/2, Omega, xCM. */
    if(Getv("DNSdata_center_new_q_timebin", "before_ell_solve"))
      DNSdata_center_q_if_desired(grid, it);

    // Centring of fields is probably not needed! We can center q!:  Yo(42);
    // /* center fields around each star */
    //DNSdata_center_fields_if_desired(grid, it);

    /* save old values before ell. solve */
    varcopy(grid, Ind("DNSdata_Psiold"),    Ind("DNSdata_Psi"));
    varcopy(grid, Ind("DNSdata_alphaPold"), Ind("DNSdata_alphaP"));
    varcopy(grid, Ind("DNSdata_Boldx"),     Ind("DNSdata_Bx"));
    varcopy(grid, Ind("DNSdata_Boldy"),     Ind("DNSdata_By"));
    varcopy(grid, Ind("DNSdata_Boldz"),     Ind("DNSdata_Bz"));
    varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
    varcopy(grid, Ind("DNSdata_qgold"),     Ind("DNSdata_qg"));

    /* set DNSdata_qnocent before ell solve */
    DNS_compute_new_q(grid, Ind("DNSdata_qnocent"));

    /* set wB before we solve */
    DNS_set_wB(grid, STAR1, Getd("DNSdata_actual_xmax1"),0.0,0.0);
    DNS_set_wB(grid, STAR2, Getd("DNSdata_actual_xmax2"),0.0,0.0);

/*
int b;
forallboxes(grid, b)
{
tBox *box = grid->box[b];
int i;
if(b==0 || b==13) continue;
forallpoints(box, i)
{
box->v[37][i] = 6.5 + ((double) rand())*0.5/RAND_MAX;
box->v[38][i] = 6.5 + ((double) rand())*0.5/RAND_MAX;
box->v[38][i] = 6.5 + ((double) rand())*0.5/RAND_MAX;
}
}
quick_Vars_output(grid->box[1], "Coordinates_CubedSphere_sigma01", 66,66);
DNSgrid_Coordinates_CubSph_sigma_continuity(grid, STAR1);
DNSgrid_Coordinates_CubSph_sigma_continuity(grid, STAR2);
quick_Vars_output(grid->box[1], "Coordinates_CubedSphere_sigma01", 77,77);
exit(99);
*/
    /* check if we do another ell. solve for DNSdata_Sigma */
    realnormres_old = realnormres; /* save realnormres */
    realnormres = normresnonlin_without_DNSdata_Sigma_outside(grid);
    itsSinceExtraSigma++;
    prdivider(1);
    printf("DNSdata_solve step %d: itsSinceExtraSigma=%d\n",
           it, itsSinceExtraSigma);
    printf(" realnormres=%g   realnormres_old=%g\n",
           realnormres, realnormres_old);
    if( (realnormres_old <= realnormres*Getd("DNSdata_extraSigmaSolve_fac")) &&
        (itsSinceExtraSigma >= Geti("DNSdata_extraSigmaSolve_every")) )
    {
      if( (!Getv("DNSdata_Sigma_surface_BCs","NoOutsideOnlySolve")) &&
          (!Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside")) )
      {
        /* solve DNSdata_Sigma completely in outer boxes */
        printf("Setting DNSdata_Sigma outside the stars only...\n");
        Sets("DNSdata_KeepInnerSigma", "yes");
        DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol,
             &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
        Sets("DNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        /* reset Sigmaold to take into account new Sigma in outer boxes */
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));

        /* reset Newton_tol, use only other vars */
        normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("DNSdata_CTS_Eqs_Iteration_order"));
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* solve the ell. eqn for Sigma alone */
      DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
             &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
      totalerr1 = average_current_and_old(Sigma_esw, 
                                          grid,vlFu,vlu,vluDerivs, vlJdu);
      if(Sigma_esw<1.0 && it>=allow_Sigma_esw1_it && allow_Sigma_esw1_it>=0)
      {
        /* complete step */
        totalerr = average_current_and_old(Sigma_esw1/Sigma_esw,
                                           grid,vlFu,vlu,vluDerivs,vlJdu);
        /* but go back to Sigma_esw if totalerr is larger */
        if(totalerr>totalerr1)
          totalerr = average_current_and_old(Sigma_esw/Sigma_esw1, 
                                             grid,vlFu,vlu,vluDerivs,vlJdu);
      }

      /* reset Newton_tol, use error norm of all vars in vlu */
      F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
      normresnonlin = GridL2Norm(vlFu);
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* reset Sigmaold so that Sigma does not change when we average later */
      varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));

      /* make sure we do not enter this block in the next iteration */
      itsSinceExtraSigma = 0;
    }

    /* How we solve the coupled ell. eqns */
    if(Getv("DNSdata_EllSolver_method", "allatonce"))
    { 
      /* solve DNSdata_Sigma completely in outer boxes at first iteration */
      if( (grid->time == 1.0-itmax) &&
          (!Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside")) )
      {
        printf("Setting DNSdata_Sigma outside the stars only (using UMFPACK)...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("DNSdata_KeepInnerSigma", "yes");
        DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
        // DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
        //       &normresnonlin, UMFPACK_solve_wrapper, 1, "DNSdata_Sigma");
        Sets("DNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* solve the coupled ell. eqns all together */
      vldummy = vlJdu;
      Newton(F_DNSdata, J_DNSdata, vlu, vlFu, vluDerivs, vldummy,
             Newton_itmax, Newton_tol, &normresnonlin, 1,
             linear_solver, Preconditioner_I, vldu, vlJdu, vlduDerivs, vlu,
             linSolver_itmax, linSolver_tolFac, linSolver_tol);
    }
    else if(Getv("DNSdata_EllSolver_method", "DNS_ordered_Eqn_Iterator"))
    {
      // /* reset Newton_tol, so that we always solve for Sigma */
      // normresnonlin = GridL2Norm_of_vars_in_string(grid, "DNSdata_Sigma");
      // Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* solve completely in outer boxes at first iteration */
      if( (grid->time == 1.0-itmax) &&
          (!Getv("DNSdata_Sigma_surface_BCs","NoOutsideOnlySolve")) &&
          (!Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside")) )
      {
        printf("Setting DNSdata_Sigma outside the stars only...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("DNSdata_KeepInnerSigma", "yes");
        DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
        Sets("DNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
      }
      /* do we want to solve for Sigma? */
      realSigmares =
       GridL2Norm_of_vars_in_string_withZeroErr_outside(grid, "DNSdata_Sigma");
      restres = GridL2Norm_of_vars_in_string(grid,  
                                      Gets("DNSdata_CTS_Eqs_Iteration_order"));
      printf(" realSigmares=%g  restres=%g\n", realSigmares, restres);
      if( realSigmares >= restres * Getd("DNSdata_SigmaSolve_tolFac") )
      {
        /* solve the ell. eqn for Sigma alone */
        DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
        totalerr1 = average_current_and_old(Sigma_esw, 
                                            grid,vlFu,vlu,vluDerivs, vlJdu);
        if(Sigma_esw<1.0 && it>=allow_Sigma_esw1_it && allow_Sigma_esw1_it>=0)
        {
          /* complete step */
          totalerr = average_current_and_old(Sigma_esw1/Sigma_esw,
                                             grid,vlFu,vlu,vluDerivs,vlJdu);
          /* but go back to Sigma_esw if totalerr is larger */
          if(totalerr>totalerr1)
            totalerr = average_current_and_old(Sigma_esw/Sigma_esw1, 
                                               grid,vlFu,vlu,vluDerivs,vlJdu);
        }
        /* reset Sigmaold so that Sigma does not change when we average later */
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
      }

      /* reset Newton_tol */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("DNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* now solve the coupled CTS ell. eqns one after an other */
      DNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("DNSdata_EllSolver_method", "DNS_ordered_Var_Eqn_Iterator"))
    {
      /* solve completely in outer boxes at first iteration */
      if( (grid->time == 1.0-itmax) &&
          (!Getv("DNSdata_Sigma_surface_BCs","NoOutsideOnlySolve")) &&
          (!Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside")) )
      {
        printf("Setting DNSdata_Sigma outside the stars only...\n");
        /* do not touch Sigma in inner boxes, but solve in outer */
        Sets("DNSdata_KeepInnerSigma", "yes");
        DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, Newton_tol, 
               &normresnonlin, linear_solver, 1, "DNSdata_Sigma");
        Sets("DNSdata_KeepInnerSigma", "no");
        totalerr1 = average_current_and_old(1, grid,vlFu,vlu,vluDerivs, vlJdu);
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
        /* reset Newton_tol, use error norm of all vars in vlu */
        F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
        normresnonlin = GridL2Norm(vlFu);
        Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);
        /* reset Sigmaold so that Sigma does not change when we average later */
        varcopy(grid, Ind("DNSdata_Sigmaold"),  Ind("DNSdata_Sigma"));
      }
      /* reset Newton_tol */
      normresnonlin = GridL2Norm_of_vars_in_string(grid, 
                                      Gets("DNSdata_CTS_Eqs_Iteration_order"));
      Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

      /* now solve ell. eqns for vars in DNSdata_CTS_Eqs_Iteration_order
         one after an other */
      DNS_ordered_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                               linear_solver, 1);
    }
    else if(Getv("DNSdata_EllSolver_method", "DNS_Eqn_Iterator"))
    { /* solve the coupled ell. eqns one after an other */
      DNS_Eqn_Iterator(grid, Newton_itmax, Newton_tol, &normresnonlin,
                       linear_solver, 1);
    }
    else
      errorexit("DNSdata_solve: unknown DNSdata_EllSolver_method");

    /* if we iterate rest masses */
    if(Getv("DNSdata_iterate_m0", "yes")) adjust_DNSdata_m01_m02();

    /* reset new values from ell. solve as average between old and new.
       I.e. do: new = esw*new + (1-esw)*old  */
    totalerr1 = average_current_and_old(esw, grid,vlFu,vlu,vluDerivs, vlJdu);
    if(esw<1.0 && it>=allow_esw1_it && allow_esw1_it>=0)
    {
      /* complete step */
      totalerr = average_current_and_old(esw1/esw, grid,vlFu,vlu,vluDerivs,vlJdu);
      /* but go back to esw if totalerr is larger */
      if(totalerr>totalerr1)
        totalerr = average_current_and_old(esw/esw1, grid,vlFu,vlu,vluDerivs,vlJdu);
    }
    else
      totalerr = totalerr1;

    /* compute diagnostics like ham and mom */
    DNSdata_verify_solution(grid);

    /* break if total error is small enough */
    if(totalerr<tol &&
       Getv("DNSdata_break_if_err_below_tol","after_ell_solve")) break;

    /* write after elliptic solve, but before adjusting q */
    grid->time -= 0.5;
    write_grid(grid);
    grid->time += 0.5;

    /* reset DNSdata_qmax1/2, DNSdata_xmax1/2 pars if the iteration 
       # it is contained in the list DNSdata_reset_qmax_xmax_pars_at */
    snprintf(str, 999, "%d", it);
    if( Getv("DNSdata_reset_qmax_xmax_pars_at", str) )
      find_qmaxs_along_x_axis_and_reset_qmaxs_xmaxs_pars(grid);
      
    /* choose how we adjust C1/2, Omega, xCM: */
    if( it>=Geti("DNSdata_adjust_first_at") &&
        Geti("DNSdata_adjust_first_at")>=0 )
    {
      if(Getv("DNSdata_adjust", "findandkeep_xmax")) /* old keep_xmax */
      { /* keep xmax1/2 in place */
        adjust_Omega_xCM_keep_xmax(grid, it, adjusttol);
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
      }
      else if(Getv("DNSdata_adjust", "keep_xmax"))
      { /* keep xmax1/2 in place, by keeping the pos. of dq/dx=0 in place */
        adjust_Omega_xCM_keep_dqdxmax_eq_0(grid, it, adjusttol);
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
      }
      else if(Getv("DNSdata_adjust", "keep_xfm"))
      { /* keep max1/2 of a funtion in place, by keeping the pos. of dFunc/dx=0 in place */
        adjust_Omega_xCM_keep_dFuncdxfm_eq_0(grid, it, adjusttol);
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
      }
      else if(Getv("DNSdata_adjust", "forcebalance"))
      { /* keep max1/2 of a funtion in place, by keeping the pos. of dFunc/dx=0 in place */
        adjust_Omega_xCM_forcebalance(grid, it, adjusttol);
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
      }
      else if(Getv("DNSdata_adjust", "Py0_forcebalance"))
      { /* use Py_ADM = 0 and forcebalance */
        adjust_xCM_Omega_Py0_forcebalance(grid, it, adjusttol);
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
      }
      else /* adjust C1/2, q while keeping restmasses, Omega and xCM */
        adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol);
    }
    else
    { /* adjust C1/2, q while keeping restmasses, Omega and xCM */
      adjust_C1_C2_q_keep_m0_or_qmax(grid, it, adjusttol*100.0); /* *100 because adjust_C1_C2_q_keep_m0_or_qmax multiplies its tol with 0.01 */
    }

    /* compute actual max pos of q and center q if needed and 
       if DNSdata_center_new_q is not "no": */
    /* set actual positions of maxima */
    set_DNSdata_actual_xyzmax_pars(grid);
    if(Getv("DNSdata_center_new_q_timebin", "after_adjusting_Omega_xCM"))
      DNSdata_center_q_if_desired(grid, it);

    /* Set VolAvSigma1/2 so that the next ell. solves all try to 
       achieve a certain Volume Average for DNSdata_Sigma. This also
       results in residuals that do not take into account the arbitrary 
       constant that can be added to DNSdata_Sigma. */
    set_DNSdata_desired_VolAvSigma12_toMinBCerr(grid, Ind("DNSdata_Sigma"));

    /* compute diagnostics like ham and mom */
    DNSdata_verify_solution(grid);

    /* evalute residual and break if it is small enough */
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
    normresnonlin = GridL2Norm(vlFu);
    printf("DNSdata_solve step %d: residual = %e\n", it, normresnonlin);
    totalerr1 = normresnonlin_without_DNSdata_Sigma_outside(grid);
    printf(" with Sigma_Err=0 outside stars: real residual = %e\n", totalerr1);
    prdivider(1);  fflush(stdout);
    if((normresnonlin<tol || totalerr1<tol) && 
       Getv("DNSdata_break_if_err_below_tol","at_iterationend")) break;

    /* set new tol for Newton */
    Newton_tol = max2(normresnonlin*NewtTolFac, tol*NewtTolFac);

    /* write current iteration if we are not done yet and increase counters */
    if(it<itmax)
    {
      write_grid(grid);
      DNS_compute_chi(grid);
      DNSdata_analyze(grid);
    }
    grid->time += 1.0;
  }
  if(it>itmax)
    printf("DNSdata_solve warning: *** Too many steps! ***\n");

  /* do we want to do a final ell. solve for some vars? */
  if( strlen(Gets("DNSdata_FinalEllSolveVars"))>0 )
  {
    double time = grid->time;
    /* write grid once more (at time=-0.1), since we have not written yet */
    grid->time = -0.1;
    write_grid(grid);
    DNSdata_analyze(grid);

    /* ell. solve */
    printf("Final elliptic solve at grid->time = -0.1\n");
    DNS_Eqn_Iterator_for_vars_in_string(grid, Newton_itmax, tol,
                                        &normresnonlin, linear_solver, 1,
                                        Gets("DNSdata_FinalEllSolveVars"));
    /* print residuals */
    DNSdata_verify_solution(grid);
    printf("Final elliptic solve:\n(%s) residual = %e\n",
           Gets("DNSdata_FinalEllSolveVars"), normresnonlin);
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
    normresnonlin = GridL2Norm(vlFu);
    printf("After final elliptic solve: vlu residual = %e\n", normresnonlin);
    totalerr1 = normresnonlin_without_DNSdata_Sigma_outside(grid);
    printf(" with Sigma_Err=0 outside stars: real residual = %e\n", totalerr1);
    prdivider(1);  fflush(stdout);
    grid->time = time; /* restore grid->time */
  }

  /* now we have intial data, set time=0 */
  grid->time = 0.0;

  /* free varlists */     
  VLDisableFree(vldu);
  VLDisableFree(vlduDerivs);
  VLDisableFree(vlJdu);     
  vlfree(vlu);
  vlfree(vluDerivs);
  vlfree(vlFu);
        
  return 0;
}

/* call DNSdata_solve or prepare for interpolation */
int setDNSdata(tGrid *grid)
{
  /* read from checkpoint if BNSdata_Interpolate_pointsfile exists */
  if(GetsLax("BNSdata_Interpolate_pointsfile")!=0)
  {
    /* read checkpoint */
    Sets("checkpoint_restart_it", "0");
    checkpoint(grid); /* works only if checkpoint_restart_it<1 */

    /* now we have read the checkpoint and only want to interpolate,
       so stop any further checkpoints and variable output */
    printf("Switching off further checkpoints and all variable output\n");
    Sets("checkpoint", "no");
    Sets("0doutiter", "-1");
    Sets("1doutiter", "-1");
    Sets("2doutiter", "-1");
    Sets("3doutiter", "-1");
    Sets("0douttime", "-1");
    Sets("1douttime", "-1");
    Sets("2douttime", "-1");
    Sets("3douttime", "-1");
  }
  else /* call solver */
    DNSdata_solve(grid);

  /* enable all ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("rho"));
  enablevar(grid, Ind("jx"));
  enablevar(grid, Ind("Sxx"));

  /* set all ADM vars */
  setADMvars(grid);

  return 0;
}

/* compute absolute error in ANALYSIS */
int DNSdata_verify_solution(tGrid *grid)
{
  /* enable all ADM vars */
  enablevar(grid, Ind("gxx"));
  enablevar(grid, Ind("alpha"));
  enablevar(grid, Ind("betax"));
  enablevar(grid, Ind("Kxx"));
  enablevar(grid, Ind("rho"));
  enablevar(grid, Ind("jx"));
  enablevar(grid, Ind("Sxx"));

  /* set all ADM vars */
  setADMvars(grid);

  /* compute constraint and other things that should be zero */
  computeADMconstraints(grid);

  return 0;
}

/* compute absolute error in ANALYSIS */
int DNSdata_analyze(tGrid *grid)
{
  double DNSdata_b = Getd("DNSdata_b");
  double n     = Getd("DNSdata_n");
  double Gamma = 1.0 + 1.0/n;
  double kappa = Getd("DNSdata_kappa");
  double Omega = Getd("DNSdata_Omega");
  double x_CM  = Getd("DNSdata_x_CM");
  double ecc   = Getd("DNSdata_ecc");
  double rdot  = Getd("DNSdata_rdot");
  double xout1, xin1, xin2, xout2;
  double M_ADM, J_ADM, Jx_ADM,Jy_ADM,Jz_ADM, m01, m02;
  double Jx_ADM1,Jy_ADM1,Jz_ADM1, Jx_ADM2,Jy_ADM2,Jz_ADM2;
  double Px_ADM1,Py_ADM1,Pz_ADM1, Px_ADM2,Py_ADM2,Pz_ADM2;
  double Px_ADM,Py_ADM,Pz_ADM;
  /* double rx_1, rx_2;
     double Sx_ADM1,Sy_ADM1,Sz_ADM1, Sx_ADM2,Sy_ADM2,Sz_ADM2; */
  double Px_1,  Py_1,  Pz_1,  Jx_1, Jy_1, Jz_1, M_1;
  double Rcx_1, Rcy_1, Rcz_1, Sx_1, Sy_1, Sz_1;
  double Px_2,  Py_2,  Pz_2,  Jx_2, Jy_2, Jz_2, M_2;
  double Rcx_2, Rcy_2, Rcz_2, Sx_2, Sy_2, Sz_2;
  int iX = Ind("X");
  int ix = Ind("x");
  int itemp1 = Ind("DNSdata_temp1");
  int itemp2 = Ind("DNSdata_temp2");
  int itemp3 = Ind("DNSdata_temp3");
  double *temp1 = grid->box[1]->v[itemp1];
  double xmax1, qmax1, xmax2, qmax2;
  int bi1, bi2;
  double Xmax1,Ymax1, Xmax2,Ymax2;
  double Zmax1,Zmax2;
  double glob_xmax1, glob_ymax1, glob_zmax1;
  double glob_xmax2, glob_ymax2, glob_zmax2;
  double global_qmax1, global_qmax2;
  double TOV_rf_surf1, TOV_m1, TOV_Phic1, TOV_Psic1, TOV_m01;  /* for TOV */
  double TOV_rf_surf2, TOV_m2, TOV_Phic2, TOV_Psic2, TOV_m02;  /* for TOV */
  double TOV_r_surf1, TOV_r_surf2, TOV_Psis1, TOV_Psis2;
  double chi1 = Getd("DNSdata_mass_shedding1");
  double chi2 = Getd("DNSdata_mass_shedding2");
  FILE *fp;
  char *outdir = Gets("outdir");
  char name[] = "DNSdata_properties.txt";
  char DNS[]  = "DNS";
  char *filename;
  int filenamelen;
  int i, fn;

  /* do nothing if BNSdata_Interpolate_pointsfile exists */
  if(GetsLax("BNSdata_Interpolate_pointsfile")!=0) return 0;
  
  printf("DNSdata_analyze: computing properties of DNS data\n");
  prTimeIn_s("WallTime: ");

  /* get inner and outer edges of both stars */
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
          double xin = box->x_of_X[1]((void *) box, -1, 1.0,0.0,0.0);
          if(box->SIDE==STAR1) xin1 = xin;
          if(box->SIDE==STAR2) xin2 = xin;
        }
        if(box->CI->dom==1) 
        {
          double xout = box->x_of_X[1]((void *) box, -1, 1.0,0.0,0.0);
          if(box->SIDE==STAR1) xout1 = xout;
          if(box->SIDE==STAR2) xout2 = xout;
        }
      }
    }
  }

  /* compute rest masses m01, m02 */
  m01 = GetInnerRestMass(grid, STAR1);
  m02 = GetInnerRestMass(grid, STAR2);
  /* set rest masses if we have fixed qm1/2 */
  if(Getd("DNSdata_qm1")>0.0) Setd("DNSdata_m01", m01);
  if(Getd("DNSdata_qm2")>0.0) Setd("DNSdata_m02", m02);

  /* compute ADM mom. */
  DNS_set_P_ADM_VolInt_integrand(grid, itemp1,itemp2,itemp3);
  Px_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp1);
  Py_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp2);
  Pz_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp3);
  Px_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp1);
  Py_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp2);
  Pz_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp3);
  Px_ADM = Px_ADM1 + Px_ADM2;
  Py_ADM = Py_ADM1 + Py_ADM2;
  Pz_ADM = Pz_ADM1 + Pz_ADM2;
  printf("Px_ADM1=%.9g Py_ADM1=%.9g Pz_ADM1=%.9g "
         "Px_ADM2=%.9g Py_ADM2=%.9g Pz_ADM2=%.9g\n",
         Px_ADM1,Py_ADM1,Pz_ADM1, Px_ADM2,Py_ADM2,Pz_ADM2);

  /* compute ADM ang. mom. */
  DNS_set_J_ADM_VolInt_integrand(grid, itemp1,itemp2,itemp3);
  Jx_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp1);
  Jy_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp2);
  Jz_ADM1 = InnerVolumeIntegral(grid, STAR1, itemp3);
  Jx_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp1);
  Jy_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp2);
  Jz_ADM2 = InnerVolumeIntegral(grid, STAR2, itemp3);
  Jx_ADM = Jx_ADM1 + Jx_ADM2;
  Jy_ADM = Jy_ADM1 + Jy_ADM2;
  Jz_ADM = Jz_ADM1 + Jz_ADM2;
  J_ADM = sqrt(Jx_ADM*Jx_ADM + Jy_ADM*Jy_ADM + Jz_ADM*Jz_ADM);
  printf("Jx_ADM1=%.9g Jy_ADM1=%.9g Jz_ADM1=%.9g "
         "Jx_ADM2=%.9g Jy_ADM2=%.9g Jz_ADM2=%.9g\n",
         Jx_ADM1,Jy_ADM1,Jz_ADM1, Jx_ADM2,Jy_ADM2,Jz_ADM2);

  /* compute ADM mass from Volume int */
  DNS_set_M_ADM_VolInt_integrand(grid, itemp1);
  M_ADM = GridVolumeIntegral(grid, itemp1);

  printf("ADM quantities: M_ADM = %.16g  J_ADM = %.16g\n", M_ADM, J_ADM);

  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qmax1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qmax2);
  /* compute qmax1 and qmax2 */
  qmax1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,0.,0.);
  qmax2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,0.,0.);
  if(grid->box[bi1]->x_of_X[1] != NULL)
    xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,0.0);
  else
    xmax1 = Xmax1;
  if(grid->box[bi2]->x_of_X[1] != NULL)
    xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,0.0);
  else
    xmax2 = Xmax2;
  /* set qmax1/2 */
  Setd("DNSdata_qmax1", qmax1);
  Setd("DNSdata_qmax2", qmax2);
  printf("DNSdata_analyze: DNSdata_qmax1 = %.10g  DNSdata_qmax2=%.10g\n"
         "  at:                    xmax1 = %.10g          xmax2=%.10g\n",
         Getd("DNSdata_qmax1"), Getd("DNSdata_qmax2"), xmax1, xmax2);
  /* set cart positions xmax1/2 of qmax1/2 ? */
  if(Getv("DNSdata_analyze_xmax", "set_DNSdata_xmax"))
  {
    Setd("DNSdata_xmax1", xmax1);
    Setd("DNSdata_xmax2", xmax2);
    printf("  setting:       DNSdata_xmax1 = %.10g  DNSdata_xmax2=%.10g\n",
           Getd("DNSdata_xmax1"), Getd("DNSdata_xmax2"));
  }
  else
  {
    printf("  keeping:       DNSdata_xmax1 = %.10g  DNSdata_xmax2=%.10g\n",
           Getd("DNSdata_xmax1"), Getd("DNSdata_xmax2"));
  }

  /* find global max of q in NS1/2 */
  Ymax1=Ymax2=Zmax1=Zmax2=0.0;
  //global_qmax1 = DNSdata_find_position_of_qmax(grid, STAR1, &bi1, &Xmax1, &Ymax1, &Zmax1);
  //global_qmax2 = DNSdata_find_position_of_qmax(grid, STAR2, &bi2, &Xmax2, &Ymax2, &Zmax2);
  global_qmax1 = DNSdata_find_xyz_of_qmax(grid, STAR1, &bi1, &Xmax1, &Ymax1, &Zmax1);
  global_qmax2 = DNSdata_find_xyz_of_qmax(grid, STAR2, &bi2, &Xmax2, &Ymax2, &Zmax2);

  if(grid->box[bi1]->x_of_X[1] != NULL)
  {
    glob_xmax1 = grid->box[bi1]->x_of_X[1]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
    glob_ymax1 = grid->box[bi1]->x_of_X[2]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
    glob_zmax1 = grid->box[bi1]->x_of_X[3]((void *) grid->box[bi1], -1, Xmax1,Ymax1,Zmax1);
  }
  else
  {
    glob_xmax1 = Xmax1;
    glob_ymax1 = Ymax1;
    glob_zmax1 = Zmax1;
  }
  if(grid->box[bi2]->x_of_X[1] != NULL)
  {
    glob_xmax2 = grid->box[bi2]->x_of_X[1]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
    glob_ymax2 = grid->box[bi2]->x_of_X[2]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
    glob_zmax2 = grid->box[bi2]->x_of_X[3]((void *) grid->box[bi2], -1, Xmax2,Ymax2,Zmax2);
  }
  else
  {
    glob_xmax2 = Xmax2;
    glob_ymax2 = Ymax2;
    glob_zmax2 = Zmax2;
  }
  printf("DNSdata_analyze: global qmax1=%g in box%d\n"
         "                 at (x,y,z)=(%.11g,%.11g,%.11g)\n",
         global_qmax1, bi1, glob_xmax1,glob_ymax1,glob_zmax1);
  printf("DNSdata_analyze: global qmax2=%g in box%d\n"
         "                 at (x,y,z)=(%.11g,%.11g,%.11g)\n",
         global_qmax2, bi2, glob_xmax2,glob_ymax2,glob_zmax2);

  /* compute spin estimates: 
     J_ADM1 = L_ADM1 + S_ADM1,  
     L_ADM1 = r_1 x P_ADM1,  
     S_ADM1 = J_ADM1 - L_ADM1,
     where rx_1 = xmax1-x_CM,  ry_1 = rz_1 = 0  */
  /* Note: this does not work, because we need to use the star's CM and
           and not the point where q is max! */
  /*
  rx_1 = xmax1-x_CM;
  Sx_ADM1 = Jx_ADM1 - 0;
  Sy_ADM1 = Jy_ADM1 - 0; 
  Sz_ADM1 = Jz_ADM1 - rx_1*Py_ADM1;
  rx_2 = xmax2-x_CM;
  Sx_ADM2 = Jx_ADM2 - 0;
  Sy_ADM2 = Jy_ADM2 - 0; 
  Sz_ADM2 = Jz_ADM2 - rx_2*Py_ADM2;
  printf("Sx_ADM1=%.9g Sy_ADM1=%.9g Sz_ADM1=%.9g "
         "Sx_ADM2=%.9g Sy_ADM2=%.9g Sz_ADM2=%.9g\n",
         Sx_ADM1,Sy_ADM1,Sz_ADM1, Sx_ADM2,Sy_ADM2,Sz_ADM2);
  */

  /* compute TOV */
  P_core1 = DNS_find_P_core(m01, 0);
  P_core2 = DNS_find_P_core(m02, 0);
  TOV_init(P_core1, 0,
           &TOV_rf_surf1, &TOV_m1, &TOV_Phic1, &TOV_Psic1, &TOV_m01);
  TOV_init(P_core2, 0,
           &TOV_rf_surf2, &TOV_m2, &TOV_Phic2, &TOV_Psic2, &TOV_m02);
  /* make TOV_rf_surf2 positive to avoid nan in fprintf */
  if(TOV_rf_surf2 <= 0.0) TOV_rf_surf2= 1e-50;
  TOV_Psis1 = 1.0 + TOV_m1/(2.0*TOV_rf_surf1);
  TOV_Psis2 = 1.0 + TOV_m2/(2.0*TOV_rf_surf2);
  TOV_r_surf1 = TOV_rf_surf1 * TOV_Psis1*TOV_Psis1;
  TOV_r_surf2 = TOV_rf_surf2 * TOV_Psis2*TOV_Psis2;
/*
printf("TOV_m1=%g TOV_r_surf1=%g TOV_Psic1=%g\n",
TOV_m1,TOV_r_surf1, TOV_Psic1);
printf("TOV_m1=%g TOV_r_surf1=%g TOV_Psis1=%g\n",
TOV_m1,TOV_r_surf1, TOV_Psis1);
*/

  /* compute star spins from surface integrals */
  setADMvars(grid);
  DNS_set_P_J_SurfInt_integrand(grid, 0, itemp1,itemp2,itemp3); // P integrand
  DNS_StarSurfInt_vector(grid, STAR1, itemp1,itemp2,itemp3, &Px_1,&Py_1,&Pz_1);
  DNS_StarSurfInt_vector(grid, STAR2, itemp1,itemp2,itemp3, &Px_2,&Py_2,&Pz_2);
  //printf("i2=%g\n", grid->box[1]->v[itemp2][0]);
  //printf("gxx=%g\n", grid->box[1]->v[Ind("gxx")][0]);
  //printf("Kxy=%g\n", grid->box[1]->v[Ind("Kxy")][0]);

  DNS_set_P_J_SurfInt_integrand(grid, 1, itemp1,itemp2,itemp3); // J integrand
  DNS_StarSurfInt_vector(grid, STAR1, itemp1,itemp2,itemp3, &Jx_1,&Jy_1,&Jz_1);
  DNS_StarSurfInt_vector(grid, STAR2, itemp1,itemp2,itemp3, &Jx_2,&Jy_2,&Jz_2);

  if(Getv("DNSdata_M_Rc_Integrals","volume"))
  {
    DNS_set_MRc_VolInt_integrand(grid, 0, itemp1,itemp2,itemp3); // M integr.
    M_1 = InnerVolumeIntegral(grid, STAR1, itemp3);
    M_2 = InnerVolumeIntegral(grid, STAR2, itemp3);

    DNS_set_MRc_VolInt_integrand(grid, 1, itemp1,itemp2,itemp3); // Rc integr.
    DNS_InnerVolInt_vector(grid,STAR1,itemp1,itemp2,itemp3, &Rcx_1,&Rcy_1,&Rcz_1);
    DNS_InnerVolInt_vector(grid,STAR2,itemp1,itemp2,itemp3, &Rcx_2,&Rcy_2,&Rcz_2);
  }
  else
  {
    DNS_set_MRc_SurfInt_integrand(grid, 0, itemp1,itemp2,itemp3); // M integr.
    M_1 = StarSurfaceIntegral(grid, STAR1, itemp3);
    M_2 = StarSurfaceIntegral(grid, STAR2, itemp3);

    DNS_set_MRc_SurfInt_integrand(grid, 1, itemp1,itemp2,itemp3); // Rc integr.
    DNS_StarSurfInt_vector(grid,STAR1,itemp1,itemp2,itemp3, &Rcx_1,&Rcy_1,&Rcz_1);
    DNS_StarSurfInt_vector(grid,STAR2,itemp1,itemp2,itemp3, &Rcx_2,&Rcy_2,&Rcz_2);
  }
  Rcx_1 = Rcx_1 / M_1;
  Rcy_1 = Rcy_1 / M_1;
  Rcz_1 = Rcz_1 / M_1;
  Rcx_2 = Rcx_2 / M_2;
  Rcy_2 = Rcy_2 / M_2;
  Rcz_2 = Rcz_2 / M_2;
  DNS_get_Spin( Px_1, Py_1, Pz_1,  Jx_1, Jy_1, Jz_1,
               Rcx_1,Rcy_1,Rcz_1, &Sx_1,&Sy_1,&Sz_1);
  DNS_get_Spin( Px_2, Py_2, Pz_2,  Jx_2, Jy_2, Jz_2,
               Rcx_2,Rcy_2,Rcz_2, &Sx_2,&Sy_2,&Sz_2);
  
  /* write into file */
  filenamelen = strlen(outdir) + strlen(name) + 200;
  filename = cmalloc(filenamelen+1);
  /* 2 filenames: fn=0. DNSdata_properties.txt, fn=1. BNSdata_properties.txt */
  for(fn=0; fn<2; fn++)
  {
    if(fn==1) name[0] = DNS[0]= 'B'; /* BNS ... */
    else      name[0] = DNS[0]= 'D'; /* DNS ... */
    snprintf(filename, filenamelen, "%s/%s", outdir, name);
    fp = fopen(filename, "a");
    if(!fp) errorexits("failed opening %s", filename);

    if(dequal(grid->time, 0.0))
    {
      fprintf(fp, "================================================\n");
      fprintf(fp, "Properties after final iteration are below:\n");
      fprintf(fp, "q_def\t\th-1\n");
      fprintf(fp, "================================================\n\n");
    }
    else fn=2; /* skip file for fn=1 (exits for(fn... loop) */

    fprintf(fp, "%s data properties (time = %g):\n", DNS, grid->time);
    fprintf(fp, "-------------------\n");
    fprintf(fp, "n_list\t\t%s\n", Gets("DNSdata_n"));
    if(EoS_pwp==1)  fprintf(fp, "rho0_list\t%s\n", Gets("DNSdata_pwp_rho0"));
    else            fprintf(fp, "rho0_list\t%s\n", "none");
    fprintf(fp, "kappa\t\t%s\n", Gets("DNSdata_kappa"));
    fprintf(fp, "x_CM\t\t%.19g\n", x_CM);
    fprintf(fp, "Omega\t\t%.19g\n", Omega);
    fprintf(fp, "ecc\t\t%.19g\n", ecc);
    fprintf(fp, "rdot\t\t%.19g\n", rdot);
    fprintf(fp, "m01\t\t%.19g\n", m01);
    fprintf(fp, "m02\t\t%.19g\n", m02);
    fprintf(fp, "J_ADM\t\t%.19g\n", J_ADM);
    fprintf(fp, "M_ADM\t\t%.19g\n", M_ADM);
    /* for(i=b1n1-1; i>=b1n1-5 && i>=0; i--)
       fprintf(fp, "MpADM[%d]\t%.19g\t(A=%.16g  x=%.16g)\n", i, M_ADM_ofA[i],
               grid->box[1]->v[iX][i], grid->box[1]->v[ix][i]); */
    fprintf(fp, "Jx_ADM\t\t%.19g\n", Jx_ADM);
    fprintf(fp, "Jy_ADM\t\t%.19g\n", Jy_ADM);
    fprintf(fp, "Jz_ADM\t\t%.19g\n", Jz_ADM);
    fprintf(fp, "Px_ADM\t\t%.19g\n", Px_ADM);
    fprintf(fp, "Py_ADM\t\t%.19g\n", Py_ADM);
    fprintf(fp, "Pz_ADM\t\t%.19g\n", Pz_ADM);
    fprintf(fp, "\n");
    fprintf(fp, "(m1)_inf\t%.19g\n", TOV_m1);
    fprintf(fp, "(m1/R)_inf\t%.19g\n", TOV_m1/TOV_r_surf1);
    /* fprintf(fp, "(m01/R)_inf\t%.19g\n", m01/TOV_r_surf1); */
    fprintf(fp, "xin1\t\t%+.19g\n", xin1);
    fprintf(fp, "xmax1\t\t%+.19g\n", Getd("DNSdata_actual_xmax1"));
    fprintf(fp, "xout1\t\t%+.19g\n", xout1);
    fprintf(fp, "qmax1\t\t%.19g\n", qmax1);
    fprintf(fp, "chi1\t\t%.19g\n", chi1);
    /* fprintf(fp, "Sx_ADM1\t\t%.19g\n", Sx_ADM1);
       fprintf(fp, "Sy_ADM1\t\t%.19g\n", Sy_ADM1);
       fprintf(fp, "Sz_ADM1\t\t%.19g\n", Sz_ADM1);
       fprintf(fp, "Px_ADM1\t\t%.19g\n", Px_ADM1);
       fprintf(fp, "Py_ADM1\t\t%.19g\n", Py_ADM1);
       fprintf(fp, "Pz_ADM1\t\t%.19g\n", Pz_ADM1); */
    fprintf(fp, "\n");
    fprintf(fp, "M_1\t\t%.19g\n", M_1);
    fprintf(fp, "Px_1\t\t%.19g\n", Px_1);
    fprintf(fp, "Py_1\t\t%.19g\n", Py_1);
    fprintf(fp, "Pz_1\t\t%.19g\n", Pz_1);
    fprintf(fp, "Jx_1\t\t%.19g\n", Jx_1);
    fprintf(fp, "Jy_1\t\t%.19g\n", Jy_1);
    fprintf(fp, "Jz_1\t\t%.19g\n", Jz_1);
    fprintf(fp, "Rcx_1\t\t%.19g\n", Rcx_1);
    fprintf(fp, "Rcy_1\t\t%.19g\n", Rcy_1);
    fprintf(fp, "Rcz_1\t\t%.19g\n", Rcz_1);
    fprintf(fp, "Sx_1\t\t%.19g\n", Sx_1);
    fprintf(fp, "Sy_1\t\t%.19g\n", Sy_1);
    fprintf(fp, "Sz_1\t\t%.19g\n", Sz_1);
    fprintf(fp, "\n");
    fprintf(fp, "(m2)_inf\t%.19g\n", TOV_m2);
    fprintf(fp, "(m2/R)_inf\t%.19g\n", TOV_m2/TOV_r_surf2);
    /* fprintf(fp, "(m02/R)_inf\t%.19g\n", m02/TOV_r_surf2); */
    fprintf(fp, "xin2\t\t%+.19g\n", xin2);
    fprintf(fp, "xmax2\t\t%+.19g\n", Getd("DNSdata_actual_xmax2"));
    fprintf(fp, "xout2\t\t%+.19g\n", xout2);
    fprintf(fp, "qmax2\t\t%.19g\n", qmax2);
    fprintf(fp, "chi2\t\t%.19g\n", chi2);
    /* fprintf(fp, "Sx_ADM2\t\t%.19g\n", Sx_ADM2);
       fprintf(fp, "Sy_ADM2\t\t%.19g\n", Sy_ADM2);
       fprintf(fp, "Sz_ADM2\t\t%.19g\n", Sz_ADM2);
       fprintf(fp, "Px_ADM2\t\t%.19g\n", Px_ADM2);
       fprintf(fp, "Py_ADM2\t\t%.19g\n", Py_ADM2);
       fprintf(fp, "Pz_ADM2\t\t%.19g\n", Pz_ADM2); */
    fprintf(fp, "\n");
    fprintf(fp, "M_2\t\t%.19g\n", M_2);
    fprintf(fp, "Px_2\t\t%.19g\n", Px_2);
    fprintf(fp, "Py_2\t\t%.19g\n", Py_2);
    fprintf(fp, "Pz_2\t\t%.19g\n", Pz_2);
    fprintf(fp, "Jx_2\t\t%.19g\n", Jx_2);
    fprintf(fp, "Jy_2\t\t%.19g\n", Jy_2);
    fprintf(fp, "Jz_2\t\t%.19g\n", Jz_2);
    fprintf(fp, "Rcx_2\t\t%.19g\n", Rcx_2);
    fprintf(fp, "Rcy_2\t\t%.19g\n", Rcy_2);
    fprintf(fp, "Rcz_2\t\t%.19g\n", Rcz_2);
    fprintf(fp, "Sx_2\t\t%.19g\n", Sx_2);
    fprintf(fp, "Sy_2\t\t%.19g\n", Sy_2);
    fprintf(fp, "Sz_2\t\t%.19g\n", Sz_2);
    fprintf(fp, "\n");
    fprintf(fp, "DNSdata_b\t%.19g\n", DNSdata_b);
    fprintf(fp, "\n");
    fclose(fp);
  }
  free(filename);
  prTimeIn_s("WallTime: ");
  return 0;
}


/* evaluate DNSdata eqns for vlu */
void F_DNSdata(tVarList *vlFu, tVarList *vlu,
               tVarList *vluDerivs, tVarList *vlc2)
{
  DNS_CTS(vlFu,vlu,  vlc2,vlc2,vluDerivs, 1);
                   /* ^----^----^--------not used by DNS_CTS if nonlin=1 */
  /* BCs */
  set_DNSdata_BCs(vlFu, vlu, vluDerivs, 1);

  if(Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside"))
  {
    set_Sigma_Omega_r_y_BC(vlFu, vlu, vluDerivs, 1);
    //set_Sigma_0_in1TOUCHATlam1A0B0_BC(vlFu, vlu, vluDerivs, 1);
  }
  else
    set_DNSdata_Sigma_BC(vlFu,vlu,  vlc2,vlc2,vluDerivs, 1);
}

/* evaluate linearized DNSdata eqns */
void J_DNSdata(tVarList *vlJdu, tVarList *vldu,
               tVarList *vlduDerivs, tVarList *vlu)
{
  DNS_CTS(vlJdu,vlu,  vlJdu,vldu,vlduDerivs, 0);
        /* ^--not used by DNS_CTS if nonlin=0 */
  /* BCs */
  set_DNSdata_BCs(vlJdu, vldu, vlduDerivs, 0);

  if(Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside"))
  {
    set_Sigma_Omega_r_y_BC(vlJdu, vldu, vlduDerivs, 0);
    //set_Sigma_0_in1TOUCHATlam1A0B0_BC(vlJdu, vldu, vlduDerivs, 0);
  }
  else
    set_DNSdata_Sigma_BC(vlJdu,vlu,  vlJdu,vldu,vlduDerivs, 0);
}


/* evaluate eqn for a SINGLE one comp. var vlw */
void F_oneComp(tVarList *vlFw, tVarList *vlw,
               tVarList *vlwDerivs, tVarList *vlc2)
{
  /* local varlists with the correct grid, with ind copied from global vl */
  tVarList *lvlFu = vlalloc(vlw->grid);
  tVarList *lvlu  = vlalloc(vlw->grid);
  tVarList *lvlJdu = vlalloc(vlw->grid);
  tVarList *lvldu  = vlalloc(vlw->grid);
  tVarList *lvlduDerivs = vlalloc(vlw->grid);
  vlpushvl(lvlFu, vlFu);
  vlpushvl( lvlu,  vlu);
  vlpushvl(lvlJdu, vlJdu);
  vlpushvl( lvldu,  vldu);
  vlpushvl(lvlduDerivs, vlduDerivs);
 
  /* Note: lvlFu,lvlu contains vlFw,vlw */
  DNS_CTS(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 1);
                    /* ^------^-----^--------not used by DNS_CTS if nonlin=1 */
  /* BCs */
  set_DNSdata_BCs(vlFw, vlw, vlwDerivs, 1);

  if(Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside"))
  {
    set_Sigma_Omega_r_y_BC(vlFw, vlw, vlwDerivs, 1);
    //set_Sigma_0_in1TOUCHATlam1A0B0_BC(vlFw, vlw, vlwDerivs, 1);
  }
  else
    set_DNSdata_Sigma_BC(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 1);
  
  /* free local varlists */
  vlfree(lvlFu);
  vlfree(lvlu);
  vlfree(lvlJdu);
  vlfree(lvldu);
  vlfree(lvlduDerivs);
}

/* evaluate linearized eqn for a SINGLE one comp. var for vldw */
void J_oneComp(tVarList *vlJdw, tVarList *vldw,
               tVarList *vldwDerivs, tVarList *vlw)
{
  /* local varlists with the correct grid, with ind copied from global vl */
  tVarList *lvlFu = vlalloc(vldw->grid);
  tVarList *lvlu  = vlalloc(vldw->grid);
  tVarList *lvlJdu = vlalloc(vldw->grid);
  tVarList *lvldu  = vlalloc(vldw->grid);
  tVarList *lvlduDerivs = vlalloc(vldw->grid);
  vlpushvl(lvlFu, vlFu);
  vlpushvl( lvlu,  vlu);
  vlpushvl(lvlJdu, vlJdu);
  vlpushvl( lvldu,  vldu);
  vlpushvl(lvlduDerivs, vlduDerivs);
  lvldu->vlPars = vldw->vlPars;  /* get pars from vldw into lvldu */

  /* Note: vlJdu,vldu contains vlJdw,vldw */
  DNS_CTS(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 0);
        /* ^--not used by DNS_CTS if nonlin=0 */
  /* BCs */
  set_DNSdata_BCs(vlJdw, vldw, vldwDerivs, 0);

  if(Getv("DNSdata_Sigma_surface_BCs","FakeMatterOutside"))
  {
    set_Sigma_Omega_r_y_BC(vlJdw, vldw, vldwDerivs, 0);
    //set_Sigma_0_in1TOUCHATlam1A0B0_BC(vlJdw, vldw, vldwDerivs, 0);
  }
  else
    set_DNSdata_Sigma_BC(lvlFu,lvlu,  lvlJdu,lvldu,lvlduDerivs, 0);

  /* free local varlists */
  vlfree(lvlFu);
  vlfree(lvlu);
  vlfree(lvlJdu);
  vlfree(lvldu);
  vlfree(lvlduDerivs);
}






/* make var lists that contain VarComp Name, its derivs, its errors,
   and the linearized var and its derivs and errors */
void make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(tGrid *grid,
     tVarList **vlw,  tVarList **vlwDerivs,  tVarList **vlFw, 
     tVarList **vldw, tVarList **vldwDerivs, tVarList **vlJdw, char *Name)
{
  char *str;
  
  str = (char *) calloc(strlen(Name)+20, sizeof(char) );

  /* allocate varlists */
  *vlw       = vlalloc(grid);
  *vlwDerivs = vlalloc(grid);
  *vlFw      = vlalloc(grid);
  *vldw      = vlalloc(grid);
  *vldwDerivs= vlalloc(grid);
  *vlJdw     = vlalloc(grid);

  /* add Name to vlw, ... */
  sprintf(str, "%s", Name);        vlpushone(*vlw,       Ind(str));
  sprintf(str, "%sx", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sy", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sz", Name);       vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxx", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxy", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%sxz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%syy", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%syz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%szz", Name);      vlpushone(*vlwDerivs, Ind(str));
  sprintf(str, "%s_Err", Name);    vlpushone(*vlFw,      Ind(str));
  sprintf(str, "%s_l", Name);      vlpushone(*vldw,       Ind(str));
  sprintf(str, "%sx_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sy_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sz_l", Name);     vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxx_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxy_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%sxz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%syy_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%syz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%szz_l", Name);    vlpushone(*vldwDerivs, Ind(str));
  sprintf(str, "%s_Err_l", Name);  vlpushone(*vlJdw,      Ind(str));
  free(str);
}

/* the var lists vlw, vlwDerivs, vlFw, vldw, vldwDerivs, vlJdw */
void free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(
     tVarList *vlw,  tVarList *vlwDerivs,  tVarList *vlFw,
     tVarList *vldw, tVarList *vldwDerivs, tVarList *vlJdw)
{
  vlfree(vlw);
  vlfree(vlwDerivs);
  vlfree(vlFw);
  vlfree(vldw);
  vlfree(vldwDerivs);
  vlfree(vlJdw);
}          

/* find residual of vars listed in string */
/* this works only it the string contains at most
   "DNSdata_Psi DNSdata_Bx DNSdata_By DNSdata_Bz DNSdata_alphaP DNSdata_Sigma" 
   */
double GridL2Norm_of_vars_in_string(tGrid *grid, char *str)
{
  int pos;
  char *word;
  double norm;
  double sum=0.0;
    
  word = cmalloc(strlen(str) + 10);
  pos=0;
  while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* make new vlw, ... for var in string word */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      
    F_oneComp(vlFw, vlw, vlwDerivs, NULL);
    norm = GridL2Norm(vlFw);
    sum += norm;
    //printf("%s: residual = %g\n", word, norm);
     
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);
  }
  //printf(" => total residual = %g\n", sum);
  free(word);
  return sum;
}

/* find residual of vars listed in string, but set Err to zero outside 
   the stars (e.g. in box1&2) */
/* this works only if the string contains at most
   "DNSdata_Psi DNSdata_Bx DNSdata_By DNSdata_Bz DNSdata_alphaP DNSdata_Sigma" 
   */
double GridL2Norm_of_vars_in_string_withZeroErr_outside(tGrid *grid, char *str)
{
  int pos;
  char *word;
  double norm;
  double sum=0.0;
    
  word = cmalloc(strlen(str) + 10);
  pos=0;
  while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
  {
    int b;
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* make new vlw, ... for var in string word */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      
    F_oneComp(vlFw, vlw, vlwDerivs, NULL);
    /* set Error to zero outside the stars */
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      int iVar  = vlFw->index[0];
      int i;
      if(box->MATTR==INSIDE) continue; /* do nothing inside stars */
      forallpoints(box, i)  box->v[iVar][i] = 0.0;
    }
    norm = GridL2Norm(vlFw);
    sum += norm;
    //printf("%s: residual = %g\n", word, norm);
     
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);
  }
  //printf(" => total residual = %g\n", sum);
  free(word);
  return sum;
}

/* solve some ell. eqns one after an other in the order given in a string */
int DNS_Eqn_Iterator_for_vars_in_string(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr, char *str)
{
  int    Newton_itmax = itmax;    /* Geti("DNSdata_Newton_itmax"); */
  double Newton_tol   = tol*0.1;  /* Getd("DNSdata_Newton_tol"); */
  int    linSolver_itmax  = Geti("DNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("DNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("DNSdata_linSolver_tol");
  double normresnonlin;
  int it, pos;
  int prN=pr;
  char *word;

  /* make sure that rhobar is initialized */
  setADMvars(grid);

  if(pr)
  { 
    printf("DNS_Eqn_Iterator_for_vars_in_string:\n"); 
    printf("%s\n", str);
    printf("  starting iterations, itmax=%d tol=%g\n", itmax, tol);
  }
  for (it = 0; it < itmax; it++)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;
    int Newton_ret, Newton_err=0;

    /* compute error */
    *normres = GridL2Norm_of_vars_in_string(grid, str);

    if(pr)
    {
      prdivider(1);
      printf("DNS_Eqn_Iterator_for_vars_in_string step %d residual = %.4e\n", it, *normres);
      fflush(stdout);
    }
    if (*normres <= tol) break;

    /* set Newton_tol */
    Newton_tol = (*normres)*0.05;

    /* go through DNSdata_Eqn_Iterator_order and solve for all vars in there */
    word = cmalloc(strlen(str) + 10);
    pos=0;
    while( (pos=sscan_word_at_p(str, pos, word)) != EOF)
    {
      int iter;

      /* make new vlw, ... for var in string word */
      make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
               &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw,  word);
      /* call Newton solver for Psi */
      prdivider(1);
      printf("Solving elliptic Eqn for %s:\n", word);
      /* call Newton several times in case we iterate */
      //for(iter=1; iter<=Newton_itmax; iter++)
      //If we solve for Psi with CTSmod, we need to iterate and always update
      //rhobar after each Newton solve. However, this iteration does not seem
      //to converge for m0=2.2 ... So for now we do only 1 iteration, i.e. we
      //don't actually iterate.
      for(iter=1; iter<=1; iter++)
      {
        Newton_ret = 
        Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
               Newton_itmax, Newton_tol, &normresnonlin, prN,
               linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
               linSolver_itmax, linSolver_tolFac, linSolver_tol);
        //if(Getv("DNSdata_CTSmod","yes")) setADMvars(grid); /* sets rhobar */
        //^Don't update rhobar, since this triggers non-conv. it. if m0=2.2

        /* compute current error, after resetting rhobar */
        normresnonlin = GridL2Norm_of_vars_in_string(grid, word);

        /* signal error if Newton_ret>Newton_itmax */
        if(Newton_ret>=Newton_itmax) Newton_err=1;
        if((normresnonlin<=Newton_tol) || Newton_err) break;
      }
      free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                             vldw, vldwDerivs, vlJdw);
    }
    free(word);
    /* quit by setting it=itmax if Newton_err=1 */
    if(Newton_err)
    {
      it=itmax;
      prdivider(1);
      printf("One of the Newton solvers has not converged -> go to step %d\n",
             it+1);
    }
  }
  /* warn if we didn't converge */
  if (it >= itmax)
  {
    *normres = GridL2Norm_of_vars_in_string(grid, str);
    prdivider(1);
    printf("DNS_Eqn_Iterator_for_vars_in_string:\n  *** Too many steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
    printf("DNS_Eqn_Iterator_for_vars_in_string:\n  Residual after %d steps:"
           "  residual = %e\n", it, *normres);
  }
  return it;
}

/* solve the coupled ell. eqns one after an other (in a particular order), 
   and iterate */
int DNS_ordered_Eqn_Iterator(tGrid *grid, int itmax, 
  double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  prdivider(0);
  printf("DNS_ordered_Eqn_Iterator: using DNS_Eqn_Iterator_for_vars_in_string\n");
  return DNS_Eqn_Iterator_for_vars_in_string(grid, itmax, tol, normres,
                 linear_solver, 1, Gets("DNSdata_CTS_Eqs_Iteration_order"));
}

/* solve the coupled ell. eqns one after an other, and iterate */
int DNS_Eqn_Iterator(tGrid *grid, int itmax, double tol, double *normres, 
  int (*linear_solver)(tVarList *x, tVarList *b, 
            tVarList *r, tVarList *c1,tVarList *c2,
	    int itmax, double tol, double *normres,
	    void (*lop)(tVarList *, tVarList *, tVarList *, tVarList *), 
	    void (*precon)(tVarList *, tVarList *, tVarList *, tVarList *)),
  int pr)
{
  int    Newton_itmax = itmax;    /* Geti("DNSdata_Newton_itmax"); */
  double Newton_tol   = tol*0.1;  /* Getd("DNSdata_Newton_tol"); */
  int    linSolver_itmax  = Geti("DNSdata_linSolver_itmax");
  double linSolver_tolFac = Getd("DNSdata_linSolver_tolFac");
  double linSolver_tol    = Getd("DNSdata_linSolver_tol");
  double normresnonlin;
  int it;
  int prN=pr;

  if(pr) printf("DNS_Eqn_Iterator:  starting iterations, itmax=%d tol=%g\n",
                itmax, tol); 
  for (it = 0; it < itmax; it++)
  {
    tVarList *vlw, *vlwDerivs, *vlFw, *vldw, *vldwDerivs, *vlJdw;

    /* compute error vlFu = F(u) */
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
    *normres = GridL2Norm(vlFu);
    if(pr)
    {
      prdivider(1);
      printf("DNS_Eqn_Iterator step %d residual = %.4e\n", it, *normres);
      fflush(stdout);
    }
    if (*normres <= tol) break;

    /* set Newton_tol */
    Newton_tol = (*normres)*0.05;

    /* make new vlw, ... for Psi */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_Psi");
    /* call Newton solver for Psi */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_Psi:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Bx */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_Bx");
    /* call Newton solver for Bx */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_Bx:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for By */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_By");
    /* call Newton solver for By */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_By:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Bz */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_Bz");
    /* call Newton solver for Bz */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_Bz:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for alphaP */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_alphaP");
    /* call Newton solver for alphaP */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_alphaP:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

    /* make new vlw, ... for Sigma */
    make_vl_vlDeriv_vlF_vld_vldDerivs_vlJd_forComponent(grid,
             &vlw,&vlwDerivs,&vlFw,  &vldw,&vldwDerivs,&vlJdw, "DNSdata_Sigma");
    /* call Newton solver for Sigma */
    prdivider(1);
    printf("Solving elliptic Eqn for DNSdata_Sigma:\n");
    Newton(F_oneComp, J_oneComp, vlw, vlFw, vlwDerivs, NULL,
           Newton_itmax, Newton_tol, &normresnonlin, prN,
           linear_solver, Preconditioner_I, vldw, vlJdw, vldwDerivs, vlw,
           linSolver_itmax, linSolver_tolFac, linSolver_tol);
    free_vl_vlDeriv_vlF_vld_vldDerivs_vlJd(vlw, vlwDerivs, vlFw,  
                                           vldw, vldwDerivs, vlJdw);

  }
  /* warn if we didn't converge */
  if (it >= itmax)
  {
    F_DNSdata(vlFu, vlu, vluDerivs, vlJdu);
    *normres = GridL2Norm(vlFu);
    prdivider(1);
    printf("DNS_Eqn_Iterator: *** Too many steps! ");
    if(*normres <= tol) printf("*** \n");
    else		printf("Tolerance goal not reached! *** \n");
    printf("DNS_Eqn_Iterator: Residual after %d steps:"
           "  residual = %e\n", it, *normres);
  }
  return it;
}




/* compute rest mass in STAR1 or STAR2 */
/* rest mass m_0 = \int d^3x \sqrt{g^{(3)}} rho_0 u^0 (-n_0),
   where n_0 = -\alpha, and d^3x \sqrt{g^{(3)}} = d^3x Psi^6  */
double GetInnerRestMass(tGrid *grid, int star)
{
  int iInteg = Ind("DNSdata_temp1");

  /* set DNSdata_temp1 = Integ = rho_0 u^0 alpha * Psi^6 */
  DNS_set_restmassintegrand(grid, iInteg);
  /* Integ is integrand for rest mass. */
  /* get rest masses */
  return InnerVolumeIntegral(grid, star, iInteg);
}



/* compute the new q and adjust the shape of the boundary between domain0/1
   or domain3/2 accordingly. This func modifies grid. It always
   interpolate vars from grid0 which is left untouched. */
void compute_new_q_and_adjust_domainshapes_InterpFromGrid0(tGrid *grid, 
                                                           tGrid *grid0, 
                                                           int star)
{
  tGrid *grid2;
  int iq = Ind("DNSdata_q");
  int interp_qgold = !Getv("DNSdata_new_q", "FromFields");
  int outerdom;

  if(star>STAR2 || star<STAR1)
    errorexit("compute_new_q_and_adjust_domainshapes_InterpFromGrid0: "
              "star is not STAR1 or STAR2");

  /* save q from grid0 on grid */
  copy_gridvar(iq, grid0, grid);

  /* compute new q on grid0 */
  DNS_compute_new_centered_q(grid0, star);

  /* make new grid2, which is an exact copy of grid0 */
  grid2 = make_empty_grid(grid0->nvariables, 0);
  copy_grid(grid0, grid2, 0);

  /* reset sigma such that q=0 at A=0 */
  reset_Coordinates_CubedSphere_sigma01(grid0, grid2, star);

  /* restore q on grid0 from grid */
  copy_gridvar(iq, grid, grid0);

  /* make sure coords on new grid are initialized */
  DNSgrid_init_Coords_for_star(grid2, star);
//quick_Vars_output(grid->box[2], "Coordinates_CubedSphere_sigma01",3,3);
//quick_Vars_output(grid2->box[2], "Coordinates_CubedSphere_sigma01",4,4);

  /* interpolate q (and maybe some other vars) from grid onto new grid2 */
  //  Interp_Var_From_Grid1_To_Grid2_star(grid, grid2, Ind("DNSdata_qg"),star);
  //  Interp_Var_From_Grid1_To_Grid2_star(grid, grid2, Ind("DNSdata_qgold"),star);
  Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_Psi"),star);
  Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_alphaP"),star);
  Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_Bx"),star);
  Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_By"),star);
  Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_Bz"),star);
  if( (star==STAR1 && !Getv("DNSdata_rotationstate1","corotation")) ||
      (star==STAR2 && !Getv("DNSdata_rotationstate2","corotation"))   )
  {
    Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_Sigma"),star);
    Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_wBx"),star);
    Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_wBy"),star);
    Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_wBz"),star);
  }
  if(interp_qgold)
    Interp_Var_From_Grid1_To_Grid2_star(grid0, grid2, Ind("DNSdata_qgold"),star);
//quick_Vars_output(grid->box[2], "DNSdata_Psi",3,3);
//quick_Vars_output(grid2->box[2],"DNSdata_Psi",4,4);

  DNS_compute_new_centered_q(grid2, star);
//quick_Vars_output(grid->box[2], "DNSdata_q",5,5);
//quick_Vars_output(grid2->box[2],"DNSdata_q",6,6);
//quick_Vars_output(grid->box[2], "Coordinates_CubedSphere_dsigma01_dA",7,7);
//quick_Vars_output(grid->box[2], "Coordinates_CubedSphere_dsigma01_dB",7,7);
//quick_Vars_output(grid2->box[2],"Coordinates_CubedSphere_dsigma01_dA",8,8);
//quick_Vars_output(grid2->box[2],"Coordinates_CubedSphere_dsigma01_dB",8,8);

//  /* set q to zero if q<0 or in region 1 and 2 */
//  set_Var_to_Val_if_below_limit_or_outside(grid, Ind("DNSdata_q"), 0.0, 0.0);

  /* copy grid2 back into grid, and free grid2 */
  copy_grid(grid2, grid, 0);
  free_grid(grid2);
}
/* compute the new q and adjust the shape of the boundary between domain0/1
   or domain3/2 accordingly. This is a wrapper that modifies grid */
void compute_new_q_and_adjust_domainshapes(tGrid *grid, int star)
{
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid, star);
}


/* guess error in m01 from inner Volume int., but without adjusting
   surfaces */
void m01_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double m01;
  t_grid_grid0_m01_m02_struct *pars;

  /* get pars */
  pars = (t_grid_grid0_m01_m02_struct *) p;

  Setd("DNSdata_C1", vec[1]);
  DNS_compute_new_centered_q(pars->grid, STAR1);
  m01 = GetInnerRestMass(pars->grid, STAR1);
  fvec[1] = m01 - pars->m01;
}

/* guess error in m02 from inner Volume int., but without adjusting
   surfaces */
void m02_guesserror_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  double m02;
  t_grid_grid0_m01_m02_struct *pars;

  /* get pars */
  pars = (t_grid_grid0_m01_m02_struct *) p;

  Setd("DNSdata_C2", vec[1]);
  DNS_compute_new_centered_q(pars->grid, STAR2);
  m02 = GetInnerRestMass(pars->grid, STAR2);
  fvec[1] = m02 - pars->m02;
}

/* compute difference m01 - DNSdata_m01 */
void m01_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  double m01;

  /* set C1 */
  Setd("DNSdata_C1", vec[1]);

//quick_Vars_output(grid0->box[1], "Coordinates_CubedSphere_sigma01",1,1,0);
//quick_Vars_output(grid->box[1], "Coordinates_CubedSphere_sigma01",2,2,0);

  /* adjust grid so that new q=0 is at A=0 */
  //compute_new_q_and_adjust_domainshapes(grid, STAR1);
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid0, STAR1);

  /*************************************/
  /* compute rest mass error Delta_m01 */
  /*************************************/
  /* get rest mass */
  m01 = GetInnerRestMass(grid, STAR1);

  printf("m01_error_VectorFuncP: C1=%.13g  m01=%.13g\n", vec[1], m01);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);
//quick_Vars_output(grid0->box[1], "Coordinates_CubedSphere_sigma01",11,11,0);
//quick_Vars_output(grid->box[1], "Coordinates_CubedSphere_sigma01",12,12,0);

  fvec[1] = m01 - pars->m01;
}

/* compute difference m02 - DNSdata_m02 */
void m02_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  double m02;

  /* set C2 */
  Setd("DNSdata_C2", vec[1]);

  /* adjust grid so that new q=0 is at A=0 */
  //compute_new_q_and_adjust_domainshapes(grid, STAR2);
  compute_new_q_and_adjust_domainshapes_InterpFromGrid0(grid, grid0, STAR2);

  /*************************************/
  /* compute rest mass error Delta_m01 */
  /*************************************/
  /* get rest mass */
  m02 = GetInnerRestMass(grid, STAR2);

  printf("m02_error_VectorFuncP: C2=%.13g  m02=%.13g\n", vec[1], m02);
  fflush(stdout);
//grid->time=-200;
//write_grid(grid);
  
  fvec[1] = m02 - pars->m02;
}

/* compute difference qm1 - DNSdata_qm1 */
void qm1_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  int bi1;
  double qm1, Xmax1, Ymax1;

  /* set C1 */
  Setd("DNSdata_C1", vec[1]);

  /* get max q locations in NS1 */
  bi1=pars->bi1;
  Xmax1=pars->Xmax1;
  Ymax1=pars->Ymax1;
  /* compute qm1 */
  qm1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,Ymax1,0);

  /******************************/
  /* compute qm error Delta_qm1 */
  /******************************/
  printf("qm1_error_VectorFuncP: C1=%.13g  qm1=%.13g\n", vec[1], qm1);
  fflush(stdout);

  fvec[1] = qm1 - pars->qm1;
}

/* compute difference qm2 - DNSdata_qm2 */
void qm2_error_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_grid_grid0_m01_m02_struct *pars = (t_grid_grid0_m01_m02_struct *) p; /* get pars */
  tGrid *grid = pars->grid;
  tGrid *grid0= pars->grid0;
  int bi2;
  double qm2, Xmax2, Ymax2;

  /* set C2 */
  Setd("DNSdata_C2", vec[1]);

  /* get max q locations in NS2 */
  bi2=pars->bi2;
  Xmax2=pars->Xmax2;
  Ymax2=pars->Ymax2;
  /* compute qm2 */
  qm2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,Ymax2,0);

  /******************************/
  /* compute qm error Delta_qm2 */
  /******************************/
  printf("qm2_error_VectorFuncP: C2=%.13g  qm2=%.13g\n", vec[1], qm2);
  fflush(stdout);

  fvec[1] = qm2 - pars->qm2;
}

/* m01_guesserror_ZP, m02_guesserror_ZP, m01_error_ZP, m02_error_ZP,
   qm1_error_ZP, qm2_error_ZP
   are wrappers around
   m01_guesserror_VectorFuncP, m02_guesserror_VectorFuncP,
   m01_error_VectorFuncP, m02_error_VectorFuncP, qm1_error_VectorFuncP,
   qm2_error_VectorFuncP
   so that we can also use zbrent_itsP */
double m01_guesserror_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  m01_guesserror_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}
double m02_guesserror_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  m02_guesserror_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}
double m01_error_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  m01_error_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}
double m02_error_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  m02_error_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}
double qm1_error_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  qm1_error_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}
double qm2_error_ZP(double C, void *p)
{
  double vec[2], fvec[2];
  vec[1] = C;
  qm2_error_VectorFuncP(1, vec, fvec, p);
  return fvec[1];
}


/* find max q along x-axis, assume that on x-axis: Y=Z=0 */
int find_Varmax_along_x_axis_in_star(tGrid *grid, int varind, int star,
                                      int *bi, double *X, double *vmax)
{
  int stat, b;
  *bi=stat=-1;
  forallboxes(grid, b)
  {
    tBox *box = grid->box[b];
    /* find boxes with matter for star and dom=0 or 1 */
    if(box->MATTR==INSIDE && box->SIDE==star && box->CI->dom<2)
      stat = box_extremum_of_F_in_dir(box, varind, 1, 0.,0., X, vmax);
    if(stat>=0) { *bi=b; break; }
  }
  return stat;
}
int find_qmax_along_x_axis(tGrid *grid, int star,
                           int *bi, double *X, double *qmax)
{
  int ret = find_Varmax_along_x_axis_in_star(grid, Ind("DNSdata_q"), star,
                                             bi, X, qmax);
  if( (*bi<0) && (star==STAR2) ) /* treat case where do not have star2 */
  {
    int b;
    forallboxes(grid, b)
    {
      tBox *box = grid->box[b];
      if( (box->SIDE==STAR2) && (box->COORD==CART) ) break;
    }
    *bi = b;
    *X  = Getd("DNSdata_xmax2");
    *qmax = 0.;
  }
  return ret;
}



/* compute deviation in desired central q value and location */
void central_q_errors_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = central_q_errors_VectorFunc__grid;
  int bi1, bi2;
  double qmax1, Xmax1, xmax1;
  double qmax2, Xmax2, xmax2;

  /* set constants */
  Setd("DNSdata_C1", vec[1]);
  Setd("DNSdata_Omega", vec[2]);
  Setd("DNSdata_C2", vec[3]);
  Setd("DNSdata_x_CM", vec[4]);
  printf("central_q_errors_VectorFunc: C1=%.6g Omega=%.6g C2=%.6g xCM=%.6g\n",
         vec[1], vec[2], vec[3], vec[4]);

  /* compute new q */
  DNS_compute_new_centered_q(grid, STAR1);
  DNS_compute_new_centered_q(grid, STAR2);

  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qmax1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qmax2);

  /* compute qmax1 and qmax2 */
  qmax1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,0.,0.);
  qmax2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,0.,0.);

  /* compute Cartesian xmax1 */
errorexit("need other boxes, not 0 and 3");
  if(bi1==0)
  {
    tBox *box = grid->box[bi1];
    xmax1 = box->x_of_X[1]((void *) box, -1, Xmax1,0.,0.);
  }
  else xmax1 = Xmax1;

  /* compute Cartesian xmax2 */
  if(bi2==3)
  {
    tBox *box = grid->box[bi2];
    xmax2 = box->x_of_X[1]((void *) box, -1, Xmax2,0.,0.);
  }
  else xmax2 = Xmax2;

  /* set fvec */
  fvec[1] = qmax1 - Getd("DNSdata_qmax1");
  fvec[2] = xmax1 - Getd("DNSdata_xmax1");
  fvec[3] = qmax2 - Getd("DNSdata_qmax2");
  fvec[4] = xmax2 - Getd("DNSdata_xmax2");

  /* print results */
  printf(" qmax1=%.6g xmax1=%g X=%.6g  fvec[1]=%g fvec[2]=%g\n",
         qmax1, xmax1, Xmax1, fvec[1], fvec[2]);
  printf(" qmax2=%.6g xmax2=%g X=%.6g  fvec[3]=%g fvec[4]=%g\n",
         qmax2, xmax2, Xmax2, fvec[3], fvec[4]);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);
}

/* compute deviation in q value at xmax */
void estimate_q_errors_VectorFunc(int n, double *vec, double *fvec)
{
  tGrid *grid = central_q_errors_VectorFunc__grid;
  int bi1, bi2;
  double qmax1, Xmax1, xmax1;
  double qmax2, Xmax2, xmax2;

  /* set constants */
  Setd("DNSdata_C1", vec[1]);
  Setd("DNSdata_C2", vec[2]);
  printf("estimate_q_errors_VectorFunc: C1=%.6g C2=%.6g\n", vec[1], vec[2]);

  /* compute new q */
  DNS_compute_new_centered_q(grid, STAR1);
  DNS_compute_new_centered_q(grid, STAR2);

  /* find max q locations xmax1/2 in NS1/2 */
  find_qmax_along_x_axis(grid, STAR1, &bi1, &Xmax1, &qmax1);
  find_qmax_along_x_axis(grid, STAR2, &bi2, &Xmax2, &qmax2);

  /* compute qmax1 and qmax2 */
  qmax1 = DNS_compute_new_centered_q_atXYZ(grid, bi1, Xmax1,0.,0.);
  qmax2 = DNS_compute_new_centered_q_atXYZ(grid, bi2, Xmax2,0.,0.);

  /* compute Cartesian xmax1 */
errorexit("need other boxes, not 0 and 3");
  if(bi1==0)
  {
    tBox *box = grid->box[bi1];
    xmax1 = box->x_of_X[1]((void *) box, -1, Xmax1,0.,0.);
  }
  else xmax1 = Xmax1;

  /* compute Cartesian xmax2 */
  if(bi2==3)
  {
    tBox *box = grid->box[bi2];
    xmax2 = box->x_of_X[1]((void *) box, -1, Xmax2,0.,0.);
  }
  else xmax2 = Xmax2;

  /* set fvec */
  fvec[1] = qmax1 - Getd("DNSdata_qmax1");
  fvec[2] = qmax2 - Getd("DNSdata_qmax2");

  /* print results */
  printf(" qmax1=%.6g xmax1=%g X=%.6g fvec[1]=%g\n",
         qmax1, xmax1, Xmax1, fvec[1]);
  printf(" qmax2=%.6g xmax2=%g X=%.6g fvec[2]=%g\n",
         qmax2, xmax2, Xmax2, fvec[2]);
  fflush(stdout);
//grid->time=-100;
//write_grid(grid);
}


/* find max of Var in star, return Varmax */
double DNSdata_find_position_of_Varmax(tGrid *grid, int vi, int star, int *bi,
                                       double *X, double *Y, double *Z)
{
  double vmax=0;
  int b;
  tBox *box;

  /* first find max on x-axis, to determine box where we search*/
  find_Varmax_along_x_axis_in_star(grid, vi, star, bi, X, &vmax);
  if(*bi<0)
  {
    printf("DNSdata_find_position_of_Varmax: bi<0: "
           "couldn't find max along x-axis\n");
    //errorexit("bi<0: couldn't find max along x-axis");
    return 0.;
  }
  box = grid->box[*bi];

  /* now look anywhere in this box */
  *Y = *Z = 0.;
  box_extremum_of_F(box, vi, X,Y,Z, &vmax);

  return vmax;
}

/* find max q in star, return qmax */
double DNSdata_find_position_of_qmax(tGrid *grid, int star, int *bi,
                                     double *X, double *Y, double *Z)
{
  return DNSdata_find_position_of_Varmax(grid, Ind("DNSdata_q"),
                                         star, bi, X,Y,Z);
}

/* find cart coords of max of q */
double DNSdata_find_xyz_of_qmax(tGrid *grid, int star, int *bi, 
                                double *x, double *y, double *z)
{
  double qmax;

  /* for now we assume max is in Cart. box anyway */
  qmax = DNSdata_find_position_of_qmax(grid, star, bi, x,y,z);
  if(*bi<0)
  {
    find_qmax_along_x_axis(grid, star, bi, x, &qmax);
    *y = *z = 0.;
  }
  if(grid->box[*bi]->COORD != CART)
    errorexit("qmax should be in a Cartesian box");
  return qmax;
}

/* set DNSdata_actual_x/y/z/max1/2 pars */
void set_DNSdata_actual_xyzmax_pars(tGrid *grid)
{
  int bi1, bi2;
  double x1,y1,z1, x2,y2,z2;
  double qmax1, qmax2;
  double xc = Getd("DNSdata_b");
  
  printf("set_DNSdata_actual_xyzmax_pars:  WallTime=%gs\n", getTimeIn_s());
  bi1=0;
  bi2=13;
  x1 = xc;
  x2 =-xc;
  y1 = y2 = z1 = z2 = 0.0;
  qmax1 = DNSdata_find_xyz_of_qmax(grid, STAR1, &bi1, &x1,&y1,&z1);
  qmax2 = DNSdata_find_xyz_of_qmax(grid, STAR2, &bi2, &x2,&y2,&z2);
  Setd("DNSdata_actual_xmax1", x1);
  Setd("DNSdata_actual_ymax1", y1);
  Setd("DNSdata_actual_zmax1", z1);
  Setd("DNSdata_actual_xmax2", x2);
  Setd("DNSdata_actual_ymax2", y2);
  Setd("DNSdata_actual_zmax2", z2);
  printf(" DNSdata_actual_xmax1 = %.13g\n", Getd("DNSdata_actual_xmax1"));
  printf(" DNSdata_actual_ymax1 = %.13g\n", Getd("DNSdata_actual_ymax1"));
  printf(" DNSdata_actual_zmax1 = %.13g\n", Getd("DNSdata_actual_zmax1"));
  printf(" DNSdata_actual_xmax2 = %.13g\n", Getd("DNSdata_actual_xmax2"));
  printf(" DNSdata_actual_ymax2 = %.13g\n", Getd("DNSdata_actual_ymax2"));
  printf(" DNSdata_actual_zmax2 = %.13g\n", Getd("DNSdata_actual_zmax2"));
}


/* set integrand for ADM mass or Center Rc volume integral:
   Rc integrands are in iIntegx, iIntegy, iIntegz if setRc=1, otherwise
   M  integrand is in iIntegz */
void DNS_set_MRc_VolInt_integrand(tGrid *grid, int setRc,
                                  int iIntegx, int iIntegy, int iIntegz)
{
  //int iPsi   = Ind("DNSdata_Psi");
  //int idPsi  = Ind("DNSdata_Psix");
  int iddPsi = Ind("DNSdata_Psixx");
  int ix = Ind("x");
  double xCM = Getd("DNSdata_x_CM");
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    //double *Psi  = box->v[iPsi];
    //double *Psix = box->v[idPsi];
    //double *Psiy = box->v[idPsi+1];
    //double *Psiz = box->v[idPsi+2];
    double *Psixx = box->v[iddPsi];
    double *Psiyy = box->v[iddPsi+3];
    double *Psizz = box->v[iddPsi+5];
    double *x   = box->v[ix];
    double *y   = box->v[ix+1];
    double *z   = box->v[ix+2];
    double *MRx = box->v[iIntegx];
    double *MRy = box->v[iIntegy];
    double *MRz = box->v[iIntegz];
    double LapPsi, oom2PI = -1./(2.*PI);
    double x1,x2,x3;
    //double Psim1,Psim2,Psim4;
    int ijk;

    forallpoints(box, ijk)
    {
      //Psim1 = 1./Psi[ijk];
      //Psim2 = Psim1*Psim1;
      //Psim4 = Psim2*Psim2;
      LapPsi = ( Psixx[ijk] + Psiyy[ijk] + Psizz[ijk] )*oom2PI;
      /* take out Psi^4 so that mass is correct for Schw. in Isotr. coords */
      //LapPsi = LapPsi * Psim4;
      // ^ not a good idea because then SurfInt is diff for M and P
      if(setRc)
      {
        x1 = x[ijk] - xCM;
        x2 = y[ijk];
        x3 = z[ijk];
        MRx[ijk] = x1 * LapPsi;
        MRy[ijk] = x2 * LapPsi;
        MRz[ijk] = x3 * LapPsi;
      }
      else
      {
        MRz[ijk] = LapPsi;
      }
    } /* end forallpoints */
  }
}

/* set Rc by integrating corresponding integrand */
void DNS_InnerVolInt_vector(tGrid *grid, int star,
                            int iIntegx, int iIntegy, int iIntegz,
                            double *Rcx, double *Rcy, double *Rcz)
{
  *Rcx = InnerVolumeIntegral(grid, star, iIntegx);
  *Rcy = InnerVolumeIntegral(grid, star, iIntegy);
  *Rcz = InnerVolumeIntegral(grid, star, iIntegz);
}

/* set integrand for ADM mass or Center Rc surface integral:
   Rc integrands are in iIntegx, iIntegy, iIntegz if setRc=1, otherwise
   M  integrand is in iIntegz */
void DNS_set_MRc_SurfInt_integrand(tGrid *grid, int setRc,
                                   int iIntegx, int iIntegy, int iIntegz)
{
  int iPsi  = Ind("DNSdata_Psi");
  int idPsi = Ind("DNSdata_Psix");
  int ix = Ind("x");
  double xCM = Getd("DNSdata_x_CM");
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    double *Psi  = box->v[iPsi];
    double *Psix = box->v[idPsi];
    double *Psiy = box->v[idPsi+1];
    double *Psiz = box->v[idPsi+2];
    double *x   = box->v[ix];
    double *y   = box->v[ix+1];
    double *z   = box->v[ix+2];
    double *MRx = box->v[iIntegx];
    double *MRy = box->v[iIntegy];
    double *MRz = box->v[iIntegz];
    double n[4];
    double ndPsi, oom2PI = -1./(2.*PI);
    double x1,x2,x3;
    //double Psim1,Psim2,Psim4;
    int ijk;

    forallpoints(box, ijk)
    {
      //Psim1 = 1./Psi[ijk];
      //Psim2 = Psim1*Psim1;
      //Psim4 = Psim2*Psim2;
      boxface_normal_at_ijk(box, 1, ijk, n); /* normal is in n[i] */
      ndPsi = ( Psix[ijk]*n[1] + Psiy[ijk]*n[2] + Psiz[ijk]*n[3] )*oom2PI;
      /* take out Psi^4 so that mass is correct for Schw. in Isotr. coords */
      //ndPsi = ndPsi * Psim4;
      // ^ not a good idea because then SurfInt is diff for M and P
      if(setRc)
      {
        x1 = x[ijk] - xCM;
        x2 = y[ijk];
        x3 = z[ijk];
        MRx[ijk] = x1 * ndPsi;
        MRy[ijk] = x2 * ndPsi;
        MRz[ijk] = x3 * ndPsi;
      }
      else
      {
        MRz[ijk] = ndPsi;
      }
    } /* end forallpoints */
  }
}

/* set integrand for P or J surface integral, set J if setJ=1 */
void DNS_set_P_J_SurfInt_integrand(tGrid *grid, int setJ,
                                   int iIntegx, int iIntegy, int iIntegz)
{
  int iK = Ind("Kxx");
  int ix = Ind("x");
  int iPsi = Ind("DNSdata_Psi");
  double xCM = Getd("DNSdata_x_CM");
  int SurfElemFlat = Getv("DNSdata_StarSurfaceIntegral_metric","flat");
  int b;

  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    double *K11 = box->v[iK];
    double *K12 = box->v[iK+1];
    double *K13 = box->v[iK+2];
    double *K22 = box->v[iK+3];
    double *K23 = box->v[iK+4];
    double *K33 = box->v[iK+5];
    double *x   = box->v[ix];
    double *y   = box->v[ix+1];
    double *z   = box->v[ix+2];
    double *Psi = box->v[iPsi];
    double *PJx = box->v[iIntegx];
    double *PJy = box->v[iIntegy];
    double *PJz = box->v[iIntegz];
    double nf[4];
    double Kn1,Kn2,Kn3, oo8PI=1.0/(8.*PI);
    double Psi2, fac, x1,x2,x3, Knphi1,Knphi2,Knphi3;
    int ijk;

    forallpoints(box, ijk)
    {
      Psi2 = Psi[ijk] * Psi[ijk];
      boxface_normal_at_ijk(box, 1, ijk, nf); /* normal is in nf[i] */
      /* NOTE: nf[i] is normlized wrt. flat metric \delta_ij
         We need the k in  n^i = k nf[i], 
         s.t. 1 = n^i n^j g_ij = k^2 nf[i] nf[j] Psi^4 \delta_ij = k^2 Psi^4.
         Thus k=Psi^{-2}. ==> n^i = Psi^{-2} nf[i]
         ALSO: if SurfElemFlat=1 we need to muliply by another Psi^4
               to get the physical surface element! */
      if(SurfElemFlat) fac = Psi2;
      else             fac = 1./Psi2;
      Kn1 = fac*( K11[ijk]*nf[1] + K12[ijk]*nf[2] + K13[ijk]*nf[3] )*oo8PI;
      Kn2 = fac*( K12[ijk]*nf[1] + K22[ijk]*nf[2] + K23[ijk]*nf[3] )*oo8PI;
      Kn3 = fac*( K13[ijk]*nf[1] + K23[ijk]*nf[2] + K33[ijk]*nf[3] )*oo8PI;
      if(setJ)
      {
        x1 = x[ijk] - xCM;
        x2 = y[ijk];
        x3 = z[ijk];
        Knphi1 = x2 * Kn3 - x3 * Kn2;
        Knphi2 = x3 * Kn1 - x1 * Kn3;
        Knphi3 = x1 * Kn2 - x2 * Kn1;
        PJx[ijk] = Knphi1;
        PJy[ijk] = Knphi2;
        PJz[ijk] = Knphi3;
      }
      else
      {
        PJx[ijk] = Kn1;
        PJy[ijk] = Kn2;
        PJz[ijk] = Kn3;
      }
    } /* end forallpoints */
  }
}

/* compute P1/2, J1/2, Rc1/2, S1/2 on star surface */
/* set P or J by inegrating corresponding integrand */
void DNS_StarSurfInt_vector(tGrid *grid, int star,
                            int iIntegx, int iIntegy, int iIntegz,
                            double *PJx, double *PJy, double *PJz)
{
  *PJx = StarSurfaceIntegral(grid, star, iIntegx);
  *PJy = StarSurfaceIntegral(grid, star, iIntegy);
  *PJz = StarSurfaceIntegral(grid, star, iIntegz);
}

/* get center of star Rc and star spin S */
void DNS_get_Spin(double Px, double Py, double Pz,
                  double Jx, double Jy, double Jz,
                  double Rcx, double Rcy, double Rcz,
                  double  *Sx, double  *Sy, double  *Sz)
{
  /* S = J - Rc /times P */
  *Sx = Jx - (Rcy * Pz - Rcz * Py);
  *Sy = Jy - (Rcz * Px - Rcx * Pz);
  *Sz = Jz - (Rcx * Py - Rcy * Px);
}
